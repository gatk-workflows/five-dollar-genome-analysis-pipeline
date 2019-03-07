version 1.0
## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for germline variant discovery of human whole-genome or exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
task HaplotypeCaller_GATK35_GVCF {
  input {
    File input_bam
    File interval_list
    String gvcf_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Int preemptible_tries
    Int hc_scatter
  }

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Int disk_size = ceil(((size(input_bam, "GB") + 30) / hc_scatter) + ref_size) + 20

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  # We use interval_padding 500 below to make sure that the HaplotypeCaller has context on both sides around
  # the interval because the assembly uses them.
  #
  # Using PrintReads is a temporary solution until we update HaploypeCaller to use GATK4. Once that is done,
  # HaplotypeCaller can stream the required intervals directly from the cloud.
  command {
    /usr/gitc/gatk4/gatk-launch --javaOptions "-Xms2g" \
      PrintReads \
      -I ~{input_bam} \
      --interval_padding 500 \
      -L ~{interval_list} \
      -O local.sharded.bam \
    && \
    java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m \
      -jar /usr/gitc/GATK35.jar \
      -T HaplotypeCaller \
      -R ~{ref_fasta} \
      -o ~{gvcf_basename}.vcf.gz \
      -I local.sharded.bam \
      -L ~{interval_list} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -variant_index_parameter 128000 \
      -variant_index_type LINEAR \
      -contamination ~{default=0 contamination} \
      --read_filter OverclippedRead
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    preemptible: preemptible_tries
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_gvcf = "~{gvcf_basename}.vcf.gz"
    File output_gvcf_index = "~{gvcf_basename}.vcf.gz.tbi"
  }
}

task HaplotypeCaller_GATK4_VCF {
  input {
    File input_bam
    File interval_list
    String vcf_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Float? contamination
    Boolean make_gvcf
    Int preemptible_tries
    Int hc_scatter
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
  }

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Int disk_size = ceil(((size(input_bam, "GB") + 30) / hc_scatter) + ref_size) + 20

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command <<<
    set -e
    gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -L ~{interval_list} \
      -O ~{vcf_basename}.g.vcf.gz \
      -contamination ~{default=0 contamination} \
      -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
      -new-qual \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      ~{true="-ERC GVCF" false="" make_gvcf}
  >>>
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "6.5 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_vcf = "~{vcf_basename}.g.vcf.gz"
    File output_vcf_index = "~{vcf_basename}.g.vcf.gz.tbi"
  }
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
task MergeVCFs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name
    Int preemptible_tries
  }

  Int disk_size = ceil(size(input_vcfs, "GB") * 2.5) + 10

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command {
    java -Xms2000m -jar /usr/gitc/picard.jar \
      MergeVcfs \
      INPUT=~{sep=' INPUT=' input_vcfs} \
      OUTPUT=~{output_vcf_name}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    preemptible: preemptible_tries
    memory: "3 GB"
    disks: "local-disk ~{disk_size} HDD"
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task HardFilterVcf {
  input {
    File input_vcf
    File input_vcf_index
    String vcf_basename
    File interval_list
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
  }

  Int disk_size = ceil(2 * size(input_vcf, "GB")) + 20
  String output_vcf_name = vcf_basename + ".filtered.vcf.gz"

  command {
     gatk --java-options "-Xms3000m" \
      VariantFiltration \
      -V ~{input_vcf} \
      -L ~{interval_list} \
      --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0" \
      --filter-name "HardFiltered" \
      -O ~{output_vcf_name}
  }
  output {
      File output_vcf = "~{output_vcf_name}"
      File output_vcf_index = "~{output_vcf_name}.tbi"
    }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
  }
}

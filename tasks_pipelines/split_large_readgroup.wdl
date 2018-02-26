## Copyright Broad Institute, 2018
##
## This WDL pipeline implements a split of large readgroups for human whole-genome and exome sequencing data.
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

import "./tasks_pipelines/alignment.wdl" as Alignment
import "./tasks_pipelines/bam_processing.wdl" as Processing
import "./tasks_pipelines/utilities.wdl" as Utils

workflow split_large_readgroup {
  File input_bam

  String bwa_commandline
  String bwa_version
  String output_bam_basename
  File ref_fasta
  File ref_fasta_index
  File ref_dict

  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative".
  File ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  Int additional_disk
  Int compression_level
  Int preemptible_tries
  Int reads_per_file = 48000000

  Float bwa_ref_size
  Float disk_multiplier

  Float unmapped_bam_size

  call Alignment.SamSplitter as SamSplitter {
    input :
      input_bam = input_bam,
      n_reads = reads_per_file,
      # Since the output bams are less compressed than the input bam we need a disk multiplier
      # that's larger than 2.
      disk_size = ceil(disk_multiplier * unmapped_bam_size + additional_disk),
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }

  scatter(unmapped_bam in SamSplitter.split_bams) {
    Float current_unmapped_bam_size = size(unmapped_bam, "GB")
    String current_name = basename(unmapped_bam, ".bam")

    call Alignment.SamToFastqAndBwaMemAndMba as SamToFastqAndBwaMemAndMba {
      input:
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = current_name,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_alt = ref_alt,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        bwa_version = bwa_version,
        # The merged bam can be bigger than only the aligned bam,
        # so account for the output size by multiplying the input size by 2.75.
        disk_size = current_unmapped_bam_size + bwa_ref_size + (disk_multiplier * current_unmapped_bam_size) + additional_disk,
        compression_level = compression_level,
        preemptible_tries = preemptible_tries
    }

    Float current_mapped_size = size(SamToFastqAndBwaMemAndMba.output_bam, "GB")
  }

  call Utils.SumFloats as SumSplitAlignedSizes {
    input:
      sizes = current_mapped_size,
      preemptible_tries = preemptible_tries
  }

  call Processing.GatherBamFiles as GatherMonolithicBamFile {
    input:
      input_bams = SamToFastqAndBwaMemAndMba.output_bam,
      disk_size = ceil((2 * SumSplitAlignedSizes.total_size) + additional_disk),
      output_bam_basename = output_bam_basename,
      preemptible_tries = preemptible_tries,
      compression_level = compression_level
  }
  output {
    File aligned_bam = GatherMonolithicBamFile.output_bam
  }
}

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data processing according to the GATK Best Practices (June 2016)
## for human whole-genome and exome sequencing data.
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
import "./tasks_pipelines/split_large_readgroup.wdl" as SplitRG
import "./tasks_pipelines/qc.wdl" as QC
import "./tasks_pipelines/bam_processing.wdl" as Processing
import "./tasks_pipelines/utilities.wdl" as Utils

# WORKFLOW DEFINITION
workflow to_bam_workflow {

  File contamination_sites_ud
  File contamination_sites_bed
  File contamination_sites_mu
  File? fingerprint_genotypes_file
  File? fingerprint_genotypes_index
  File? haplotype_database_file
  File wgs_coverage_interval_list

  String sample_name
  String base_file_name
  Array[File] flowcell_unmapped_bams
  String unmapped_bam_suffix

  Int read_length = 250

  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File ref_alt
  File ref_bwt
  File ref_sa
  File ref_amb
  File ref_ann
  File ref_pac

  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  Int preemptible_tries
  Int agg_preemptible_tries

  Float cutoff_for_large_rg_in_gb = 20.0

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB
  Int additional_disk = select_first([increase_disk_size, 20])
  # Sometimes the output is larger than the input, or a task can spill to disk. In these cases we need to account for the
  # input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
  Float bwa_disk_multiplier = 2.5
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data
  # so it needs more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a
  # larger multiplier
  Float sort_sam_disk_multiplier = 3.25

  # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom
  Float md_disk_multiplier = 2.25

  String bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"

  String recalibrated_bam_basename = base_file_name + ".aligned.duplicates_marked.recalibrated"

  Int compression_level = 2

  # Get the version of BWA to include in the PG record in the header of the BAM produced
  # by MergeBamAlignment.
  call Alignment.GetBwaVersion

  # Get the size of the standard reference files as well as the additional reference files needed for BWA
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
  Float bwa_ref_size = ref_size + size(ref_alt, "GB") + size(ref_amb, "GB") + size(ref_ann, "GB") + size(ref_bwt, "GB") + size(ref_pac, "GB") + size(ref_sa, "GB")
  Float dbsnp_size = size(dbSNP_vcf, "GB")

  # Align flowcell-level unmapped input bams in parallel
  scatter (unmapped_bam in flowcell_unmapped_bams) {

    Float unmapped_bam_size = size(unmapped_bam, "GB")

    String sub_strip_path = "gs://.*/"
    String sub_strip_unmapped = unmapped_bam_suffix + "$"
    String sub_sub = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "")

    # QC the unmapped BAM
    call QC.CollectQualityYieldMetrics as CollectQualityYieldMetrics {
      input:
        input_bam = unmapped_bam,
        metrics_filename = sub_sub + ".unmapped.quality_yield_metrics",
        disk_size = unmapped_bam_size + additional_disk,
        preemptible_tries = preemptible_tries
    }

    if (unmapped_bam_size > cutoff_for_large_rg_in_gb) {
      # Split bam into multiple smaller bams,
      # map reads to reference and recombine into one bam
      call SplitRG.split_large_readgroup as SplitRG {
        input:
          input_bam = unmapped_bam,
          bwa_commandline = bwa_commandline,
          bwa_version = GetBwaVersion.version,
          output_bam_basename = sub_sub + ".aligned.unsorted",
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          ref_alt = ref_alt,
          ref_amb = ref_amb,
          ref_ann = ref_ann, 
          ref_bwt = ref_bwt,
          ref_pac = ref_pac,
          ref_sa = ref_sa,
          additional_disk = additional_disk,
          compression_level = compression_level,
          preemptible_tries = preemptible_tries,
          bwa_ref_size = bwa_ref_size,
          disk_multiplier = bwa_disk_multiplier,
          unmapped_bam_size = unmapped_bam_size
      }
    }

    if (unmapped_bam_size <= cutoff_for_large_rg_in_gb) {
      # Map reads to reference
      call Alignment.SamToFastqAndBwaMemAndMba as SamToFastqAndBwaMemAndMba {
        input:
          input_bam = unmapped_bam,
          bwa_commandline = bwa_commandline,
          output_bam_basename = sub_sub + ".aligned.unsorted",
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          ref_alt = ref_alt,
          ref_bwt = ref_bwt,
          ref_amb = ref_amb,
          ref_ann = ref_ann,
          ref_pac = ref_pac,
          ref_sa = ref_sa,
          bwa_version = GetBwaVersion.version,
          # The merged bam can be bigger than only the aligned bam,
          # so account for the output size by multiplying the input size by 2.75.
          disk_size = unmapped_bam_size + bwa_ref_size + (bwa_disk_multiplier * unmapped_bam_size) + additional_disk,
          compression_level = compression_level,
          preemptible_tries = preemptible_tries
      }
    }

    File output_aligned_bam = select_first([SamToFastqAndBwaMemAndMba.output_bam, SplitRG.aligned_bam])

    Float mapped_bam_size = size(output_aligned_bam, "GB")

    # QC the aligned but unsorted readgroup BAM
    # no reference as the input here is unsorted, providing a reference would cause an error
    call QC.CollectUnsortedReadgroupBamQualityMetrics as CollectUnsortedReadgroupBamQualityMetrics {
      input:
        input_bam = output_aligned_bam,
        output_bam_prefix = sub_sub + ".readgroup",
        disk_size = mapped_bam_size + additional_disk,
        preemptible_tries = preemptible_tries
    }
  }

  # Sum the read group bam sizes to approximate the aggregated bam size
  call Utils.SumFloats as SumFloats {
    input:
      sizes = mapped_bam_size,
      preemptible_tries = preemptible_tries
  }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call Processing.MarkDuplicates as MarkDuplicates {
    input:
      input_bams = output_aligned_bam,
      output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs
      # and the merged output.
      disk_size = (md_disk_multiplier * SumFloats.total_size) + additional_disk,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  Float agg_bam_size = size(MarkDuplicates.output_bam, "GB")

  # Sort aggregated+deduped BAM file and fix tags
  call Processing.SortSam as SortSampleBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
      # This task spills to disk so we need space for the input bam, the output bam, and any spillage.
      disk_size = (sort_sam_disk_multiplier * agg_bam_size) + additional_disk,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  if (defined(haplotype_database_file)) {
    # Check identity of fingerprints across readgroups
    call QC.CrossCheckFingerprints as CrossCheckFingerprints {
      input:
        input_bams = SortSampleBam.output_bam,
        input_bam_indexes = SortSampleBam.output_bam_index,
        haplotype_database_file = haplotype_database_file,
        metrics_filename = base_file_name + ".crosscheck",
        disk_size = agg_bam_size + additional_disk,
        preemptible_tries = agg_preemptible_tries
    }
  }

  # Create list of sequences for scatter-gather parallelization
  call Utils.CreateSequenceGroupingTSV as CreateSequenceGroupingTSV {
    input:
      ref_dict = ref_dict,
      preemptible_tries = preemptible_tries
  }

  # Estimate level of cross-sample contamination
  call Processing.CheckContamination as CheckContamination {
    input:
      input_bam = SortSampleBam.output_bam,
      input_bam_index = SortSampleBam.output_bam_index,
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_prefix = base_file_name + ".preBqsr",
      disk_size = agg_bam_size + ref_size + additional_disk,
      preemptible_tries = agg_preemptible_tries,
      contamination_underestimation_factor = 0.75
  }

  # We need disk to localize the sharded input and output due to the scatter for BQSR.
  # If we take the number we are scattering by and reduce by 3 we will have enough disk space
  # to account for the fact that the data is not split evenly.
  Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping)
  Int potential_bqsr_divisor = num_of_bqsr_scatters - 10
  Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call Processing.BaseRecalibrator as BaseRecalibrator {
      input:
        input_bam = SortSampleBam.output_bam,
        recalibration_report_filename = base_file_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = known_indels_sites_VCFs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        # We need disk to localize the sharded bam due to the scatter.
        disk_size = (agg_bam_size / bqsr_divisor) + ref_size + dbsnp_size + additional_disk,
        preemptible_tries = agg_preemptible_tries
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  # The reports are always the same size
  call Processing.GatherBqsrReports as GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = base_file_name + ".recal_data.csv",
      disk_size = additional_disk,
      preemptible_tries = preemptible_tries
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
    # Apply the recalibration model by interval
    call Processing.ApplyBQSR as ApplyBQSR {
      input:
        input_bam = SortSampleBam.output_bam,
        output_bam_basename = recalibrated_bam_basename,
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        # We need disk to localize the sharded bam and the sharded output due to the scatter.
        disk_size = ((agg_bam_size * 3) / bqsr_divisor) + ref_size + additional_disk,
        compression_level = compression_level,
        preemptible_tries = agg_preemptible_tries
    }
  }

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call Processing.GatherBamFiles as GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = base_file_name,
      # Multiply the input bam size by two to account for the input and output
      disk_size = (2 * agg_bam_size) + additional_disk,
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries
  }

  #BQSR bins the qualities which makes a significantly smaller bam
  Float binned_qual_bam_size = size(GatherBamFiles.output_bam, "GB")

  # QC the final BAM (consolidated after scattered BQSR)
  call QC.CollectReadgroupBamQualityMetrics as CollectReadgroupBamQualityMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      output_bam_prefix = base_file_name + ".readgroup",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      disk_size = binned_qual_bam_size + ref_size + additional_disk,
      preemptible_tries = agg_preemptible_tries
  }

  # QC the final BAM some more (no such thing as too much QC)
  call QC.CollectAggregationMetrics as CollectAggregationMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      output_bam_prefix = base_file_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      disk_size = binned_qual_bam_size + ref_size + additional_disk,
      preemptible_tries = agg_preemptible_tries
  }

  if (defined(haplotype_database_file) && defined(fingerprint_genotypes_file)) {
    # Check the sample BAM fingerprint against the sample array
    call QC.CheckFingerprint as CheckFingerprint {
      input:
        input_bam = GatherBamFiles.output_bam,
        input_bam_index = GatherBamFiles.output_bam_index,
        haplotype_database_file = haplotype_database_file,
        genotypes = fingerprint_genotypes_file,
        genotypes_index = fingerprint_genotypes_index,
        output_basename = base_file_name,
        sample = sample_name,
        disk_size = binned_qual_bam_size + additional_disk,
        preemptible_tries = agg_preemptible_tries
    }
  }

  # QC the sample WGS metrics (stringent thresholds)
  call QC.CollectWgsMetrics as CollectWgsMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      metrics_filename = base_file_name + ".wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      read_length = read_length,
      disk_size = binned_qual_bam_size + ref_size + additional_disk,
      preemptible_tries = agg_preemptible_tries
  }

  # QC the sample raw WGS metrics (common thresholds)
  call QC.CollectRawWgsMetrics as CollectRawWgsMetrics {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      metrics_filename = base_file_name + ".raw_wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      read_length = read_length,
      disk_size = binned_qual_bam_size + ref_size + additional_disk,
      preemptible_tries = agg_preemptible_tries
  }

  # Generate a checksum per readgroup in the final BAM
  call QC.CalculateReadGroupChecksum as CalculateReadGroupChecksum {
    input:
      input_bam = GatherBamFiles.output_bam,
      input_bam_index = GatherBamFiles.output_bam_index,
      read_group_md5_filename = recalibrated_bam_basename + ".bam.read_group_md5",
      disk_size = binned_qual_bam_size + additional_disk,
      preemptible_tries = agg_preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = CollectQualityYieldMetrics.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = CollectUnsortedReadgroupBamQualityMetrics.insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = CollectUnsortedReadgroupBamQualityMetrics.insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_metrics

    File read_group_alignment_summary_metrics = CollectReadgroupBamQualityMetrics.alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = CollectReadgroupBamQualityMetrics.gc_bias_detail_metrics
    File read_group_gc_bias_pdf = CollectReadgroupBamQualityMetrics.gc_bias_pdf
    File read_group_gc_bias_summary_metrics = CollectReadgroupBamQualityMetrics.gc_bias_summary_metrics

    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.cross_check_fingerprints_metrics

    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination

    File calculate_read_group_checksum_md5 = CalculateReadGroupChecksum.md5_file

    File agg_alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics
    File agg_bait_bias_detail_metrics = CollectAggregationMetrics.bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = CollectAggregationMetrics.bait_bias_summary_metrics
    File agg_gc_bias_summary_metrics = CollectAggregationMetrics.gc_bias_summary_metrics
    File agg_gc_bias_detail_metrics = CollectAggregationMetrics.gc_bias_detail_metrics
    File agg_gc_bias_pdf = CollectAggregationMetrics.gc_bias_pdf
    File agg_insert_size_histogram_pdf = CollectAggregationMetrics.insert_size_histogram_pdf
    File agg_insert_size_metrics = CollectAggregationMetrics.insert_size_metrics
    File agg_pre_adapter_detail_metrics = CollectAggregationMetrics.pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = CollectAggregationMetrics.pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = CollectAggregationMetrics.quality_distribution_pdf
    File agg_quality_distribution_metrics = CollectAggregationMetrics.quality_distribution_metrics

    File? fingerprint_summary_metrics = CheckFingerprint.summary_metrics
    File? fingerprint_detail_metrics = CheckFingerprint.detail_metrics

    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    File output_bqsr_reports = GatherBqsrReports.output_bqsr_report

    File raw_wgs_metrics = CollectRawWgsMetrics.metrics
    File wgs_metrics = CollectWgsMetrics.metrics

    File output_bam = GatherBamFiles.output_bam
    File? output_bam_index = GatherBamFiles.output_bam_index
  }
}

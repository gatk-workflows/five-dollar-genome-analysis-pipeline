## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "./tasks_pipelines/unmapped_bam_to_aligned_bam.wdl" as ToBam
import "./tasks_pipelines/germline_variant_discovery.wdl" as Calling
import "./tasks_pipelines/qc.wdl" as QC
import "./tasks_pipelines/utilities.wdl" as Utils

# WORKFLOW DEFINITION
workflow germline_single_sample_workflow {

  File contamination_sites_ud
  File contamination_sites_bed
  File contamination_sites_mu
  File? fingerprint_genotypes_file
  File? fingerprint_genotypes_index
  File? haplotype_database_file
  File wgs_evaluation_interval_list
  File wgs_coverage_interval_list

  String sample_name
  String base_file_name
  String final_gvcf_base_name
  Array[File] flowcell_unmapped_bams
  String unmapped_bam_suffix

  File wgs_calling_interval_list
  Int haplotype_scatter_count
  Int break_bands_at_multiples_of
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

  call ToBam.to_bam_workflow {
    input:
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      haplotype_database_file = haplotype_database_file,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      sample_name = sample_name,
      base_file_name = base_file_name,
      flowcell_unmapped_bams = flowcell_unmapped_bams,
      unmapped_bam_suffix = unmapped_bam_suffix,
      read_length = read_length,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      ref_alt = ref_alt,
      ref_bwt = ref_bwt,
      ref_sa = ref_sa,
      ref_amb = ref_amb,
      ref_ann = ref_ann,
      ref_pac = ref_pac,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      known_indels_sites_VCFs = known_indels_sites_VCFs,
      known_indels_sites_indices = known_indels_sites_indices,
      preemptible_tries = preemptible_tries,
      agg_preemptible_tries = agg_preemptible_tries
  }

  # ValidateSamFile runs out of memory in mate validation on crazy edge case data, so we want to skip the mate validation
  # in those cases.  These values set the thresholds for what is considered outside the normal realm of "reasonable" data.
  Float max_duplication_in_reasonable_sample = 0.30
  Float max_chimerism_in_reasonable_sample = 0.15

  # Convert the final merged recalibrated BAM file to CRAM format
  call Utils.ConvertToCram as ConvertToCram {
    input:
      input_bam = to_bam_workflow.output_bam,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_basename = base_file_name,
      preemptible_tries = agg_preemptible_tries
  }

  # Check whether the data has massively high duplication or chimerism rates
  call QC.CheckPreValidation as CheckPreValidation {
    input:
      duplication_metrics = to_bam_workflow.duplicate_metrics,
      chimerism_metrics = to_bam_workflow.agg_alignment_summary_metrics,
      max_duplication_in_reasonable_sample = max_duplication_in_reasonable_sample,
      max_chimerism_in_reasonable_sample = max_chimerism_in_reasonable_sample,
      preemptible_tries = agg_preemptible_tries
 }

  # Validate the CRAM file
  call QC.ValidateSamFile as ValidateCram {
    input:
      input_bam = ConvertToCram.output_cram,
      input_bam_index = ConvertToCram.output_cram_index,
      report_filename = base_file_name + ".cram.validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ignore = ["MISSING_TAG_NM"],
      max_output = 1000000000,
      is_outlier_data = CheckPreValidation.is_outlier_data,
      preemptible_tries = agg_preemptible_tries
  }

  # Break the calling interval_list into sub-intervals
  # Perform variant calling on the sub-intervals, and then gather the results
  call Utils.ScatterIntervalList as ScatterIntervalList {
    input:
      interval_list = wgs_calling_interval_list,
      scatter_count = haplotype_scatter_count,
      break_bands_at_multiples_of = break_bands_at_multiples_of
  }

  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by 20 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
  Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  # Call variants in parallel over WGS calling intervals
  scatter (index in range(ScatterIntervalList.interval_count)) {
    # Generate GVCF by interval
    call Calling.HaplotypeCaller_GATK35_GVCF as HaplotypeCaller {
      input:
        contamination = to_bam_workflow.contamination,
        input_bam = to_bam_workflow.output_bam,
        interval_list = ScatterIntervalList.out[index],
        gvcf_basename = base_file_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        hc_scatter = hc_divisor,
        preemptible_tries = agg_preemptible_tries
     }
  }

  # Combine by-interval GVCFs into a single sample GVCF file
  call Calling.MergeVCFs as MergeVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_gvcf,
      input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
      output_vcf_name = final_gvcf_base_name + ".g.vcf.gz",
      preemptible_tries = agg_preemptible_tries
  }

  Float gvcf_size = size(MergeVCFs.output_vcf, "GB")

  # Validate the GVCF output of HaplotypeCaller
  call QC.ValidateGVCF as ValidateGVCF {
    input:
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      wgs_calling_interval_list = wgs_calling_interval_list,
      preemptible_tries = agg_preemptible_tries
  }

  # QC the GVCF
  call QC.CollectGvcfCallingMetrics as CollectGvcfCallingMetrics {
    input:
      input_vcf = MergeVCFs.output_vcf,
      input_vcf_index = MergeVCFs.output_vcf_index,
      metrics_basename = final_gvcf_base_name,
      dbSNP_vcf = dbSNP_vcf,
      dbSNP_vcf_index = dbSNP_vcf_index,
      ref_dict = ref_dict,
      wgs_evaluation_interval_list = wgs_evaluation_interval_list,
      preemptible_tries = agg_preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] quality_yield_metrics = to_bam_workflow.quality_yield_metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = to_bam_workflow.unsorted_read_group_base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = to_bam_workflow.unsorted_read_group_base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = to_bam_workflow.unsorted_read_group_insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = to_bam_workflow.unsorted_read_group_insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = to_bam_workflow.unsorted_read_group_quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = to_bam_workflow.unsorted_read_group_quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = to_bam_workflow.unsorted_read_group_quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = to_bam_workflow.unsorted_read_group_quality_distribution_metrics

    File read_group_alignment_summary_metrics = to_bam_workflow.read_group_alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = to_bam_workflow.read_group_gc_bias_detail_metrics
    File read_group_gc_bias_pdf = to_bam_workflow.read_group_gc_bias_pdf
    File read_group_gc_bias_summary_metrics = to_bam_workflow.read_group_gc_bias_summary_metrics

    File? cross_check_fingerprints_metrics = to_bam_workflow.cross_check_fingerprints_metrics

    File selfSM = to_bam_workflow.selfSM
    Float contamination = to_bam_workflow.contamination

    File calculate_read_group_checksum_md5 = to_bam_workflow.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = to_bam_workflow.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = to_bam_workflow.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = to_bam_workflow.agg_bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = to_bam_workflow.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = to_bam_workflow.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = to_bam_workflow.agg_gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = to_bam_workflow.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = to_bam_workflow.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = to_bam_workflow.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = to_bam_workflow.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = to_bam_workflow.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = to_bam_workflow.agg_quality_distribution_metrics

    File? fingerprint_summary_metrics = to_bam_workflow.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = to_bam_workflow.fingerprint_detail_metrics

    File wgs_metrics = to_bam_workflow.wgs_metrics
    File raw_wgs_metrics = to_bam_workflow.raw_wgs_metrics

    File duplicate_metrics = to_bam_workflow.duplicate_metrics
    File output_bqsr_reports = to_bam_workflow.output_bqsr_reports

    File gvcf_summary_metrics = CollectGvcfCallingMetrics.summary_metrics
    File gvcf_detail_metrics = CollectGvcfCallingMetrics.detail_metrics

    File output_cram = ConvertToCram.output_cram
    File output_cram_index = ConvertToCram.output_cram_index
    File output_cram_md5 = ConvertToCram.output_cram_md5

    File validate_cram_file_report = ValidateCram.report

    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_index = MergeVCFs.output_vcf_index
  }
}

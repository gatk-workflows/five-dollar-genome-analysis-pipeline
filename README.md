# five-dollar-genome-analysis-pipeline
Workflows used for germline short variant discovery in WGS data
### germline_single_sample_workflow :
This WDL pipeline implements data pre-processing and initial variant calling (GVCF
generation) according to the GATK Best Practices (June 2016) for germline SNP and
Indel discovery in human whole-genome sequencing data.

Note: For those users interested in running this wdl on FireCloud (FC), the FC
version has been provided as fc_germline_single_sample_workflow.wdl. Please visit the 
FC featured methods and workspaces for more GATK Best Practices pipelines.

#### Requirements/expectations
- Human whole-genome paired-end sequencing data in unmapped BAM (uBAM) format
- One or more read groups, one per uBAM file, all belonging to a single sample (SM)
- Input uBAM files must additionally comply with the following requirements:
- - filenames all have the same suffix (we use ".unmapped.bam")
- - files must pass validation by ValidateSamFile
- - reads are provided in query-sorted order
- - all reads must have an RG tag
- Reference genome must be Hg38 with ALT contigs
#### Outputs 
- Cram, cram index, and cram md5 
- GVCF and its gvcf index 
- BQSR Report
- Several Summary Metrics 

### Software version requirements :
Cromwell version support 
- Successfully tested on v30.2
- Does not work on versions < v23 due to output syntax

#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2020

This software is a computer program whose purpose is to
analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms
of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful,
but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards
their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge
of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/

/*
========================================================================================
                         @git_repo_name@
========================================================================================
 @git_repo_name@ analysis Pipeline.
 #### Homepage / Documentation
 @git_url@
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// File with text to display when a developement version is used
devMessageFile = file("$projectDir/assets/devMessage.txt")

def helpMessage() {
  
  if ("${workflow.manifest.version}" =~ /dev/ ){
     log.info devMessageFile.text
  }

  log.info """
    ==============================================
    EUCANCAN GOLDEN DATASET BENCHMARKING PIPELINE 
    ==============================================
	    
    Usage:
    Run the pipeline with default parameters:
	    nextflow run main.nf -profile docker

    Run with user parameters:

 	    nextflow run main.nf -profile docker --input {variant.calling.file} --public_ref_dir {validation.reference.file} --participant_id {tool.name} --goldstandard_dir {gold.standards.dir} --assess_dir {benchmark.data.dir} --results_dir {output.dir}

    Mandatory arguments:
      --input [file]        Variant calling file(s)
      --community_id [id]   Name or OEB permanent ID for the benchmarking community
      --public_ref_dir      Directory with list of cancer genes used to validate the predictions
      --participant_id  		Name of the tool used for prediction
      --goldstandard_dir 		Dir that contains metrics reference datasets
      --challenges_ids  		List of types of cancer selected by the user, separated by spaces
      --assess_dir			Dir where the data for the benchmark are stored

    Other options:
      --validation_result		The output directory where the results from validation step will be saved
      --assessment_results	The output directory where the results from the computed metrics step will be saved
      --outdir				The output directory where the consolidation of the benchmark will be saved
      --statsdir				The output directory with nextflow statistics
      --data_model_export_dir	The output dir where json file with benchmarking data model contents will be saved
      --otherdir					The output directory where custom results will be saved (no directory inside)

    Flags:
      --help			Display this message
	    
  @git_repo_name@ version: ${workflow.manifest.version}
  ======================================================================

  Usage:
  nextflow run main.nf --reads '*_R{1,2}.fastq.gz' --genome 'hg19' -profile conda
  nextflow run main.nf --samplePlan samplePlan --genome 'hg19' -profile conda

  Mandatory arguments:
    --reads [file]                Path to input data (must be surrounded with quotes)
    --samplePlan [file]           Path to sample plan input file (cannot be used with --reads)
    --genome [str]                Name of genome reference
    -profile [str]                Configuration profile to use. test / conda / multiconda / path / multipath / singularity / docker / cluster (see below)
  
  Inputs:
    --design [file]               Path to design file for extended analysis  
    --singleEnd [bool]            Specifies that the input is single-end reads

  Skip options: All are false by default
    --skipSoftVersion [bool]      Do not report software version
    --skipMultiQC [bool]          Skips MultiQC

  Other options:
    --outDir [file]               The output directory where the results will be saved
    -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
 
  =======================================================
  Available Profiles

    -profile test                Set up the test dataset
    -profile conda               Build a single conda for with all tools used by the different processes before running the pipeline
    -profile multiconda          Build a new conda environment for each tools used by the different processes before running the pipeline
    -profile path                Use the path defined in the configuration for all tools
    -profile multipath           Use the paths defined in the configuration for each tool
    -profile docker              Use the Docker containers for each process
    -profile singularity         Use the singularity images for each process
    -profile cluster             Run the workflow on the cluster, instead of locally

  """.stripIndent()
}

/**********************************
 * SET UP CONFIGURATION VARIABLES *
 **********************************/

// Show help message
if (params.help){
  helpMessage()
  exit 0
}

// Configurable reference genomes

// ADD HERE ANY ANNOTATION



params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false

if ( params.fasta ){
Channel.fromPath(params.fasta)
  .ifEmpty { exit 1, "Reference annotation not found: ${params.fasta}" }
  .set { fastaCh }
}else{
  fastaCh = Channel.empty()
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  customRunName = workflow.runName
}

// Stage config files
multiqcConfigCh = Channel.fromPath(params.multiqcConfig)
outputDocsCh = file("$projectDir/docs/output.md", checkIfExists: true)
outputDocsImagesCh = file("$projectDir/docs/images/", checkIfExists: true)

/************
 * CHANNELS *
 ************/
// input files

input_file = file(params.input)
ref_dir = Channel.fromPath( params.public_ref_dir, type: 'dir' )
tool_name = params.participant_id.replaceAll("\\s","_")
gold_standards_dir = Channel.fromPath(params.goldstandard_dir, type: 'dir' )
benchmark_data = Channel.fromPath(params.assess_dir, type: 'dir' )
community_id = params.community_id

// output
validation_out = file(params.validation_result)
assessment_out = file(params.assessment_results)
aggregation_dir = file(params.outdir)
data_model_export_dir = file(params.data_model_export_dir)
other_dir = file(params.otherdir)

/*******************
 * Header log info *
 *******************/

if ("${workflow.manifest.version}" =~ /dev/ ){
  log.info devMessageFile.text
}

log.info """\
=======================================================
@git_repo_name@ workflow version: ${workflow.manifest.version}
======================================================="""

summary = [
    'Max Memory': params.maxMemory,
    'Max CPUs': params.maxCpus,
    'Max Time': params.maxTime,
    'Container Engine': workflow.containerEngine,
    'Current home': "$USER",
    'Current path': "$PWD",
    'Working dir': workflow.workDir,
    'Output dir': params.outDir,
    'Config Profile': workflow.profile
].findAll{ it.value != null }

log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

// ADD YOUR NEXTFLOW SUBWORFLOWS HERE

// Workflows
// QC : check design and factqc
include { validation } from './nf-modules/local/subworkflow/validation'
include { metrics } from './nf-modules/local/subworkflow/metrics'
include { consolidation } from './nf-modules/local/subworkflow/metrics'

// Processes
include { outputDocumentation } from './nf-modules/local/process/outputDocumentation'

workflow {
    main:

     validation(
       designCheckCh,
       samplePlanCh,
       rawReadsFastqcCh
     )

     oneToFiveCh = Channel.of(1..5)
     metrics(
       oneToFiveCh
     )

     /*********************
      * Software versions *
      *********************/
     consolidation(
       metrics.out.version.first().ifEmpty([])
     )

     /********************
      * Workflow summary *
      ********************/
     workflowSummaryMqc(summary)

     /****************
      * Sub-routines *
      ****************/
     outputDocumentation(
       outputDocsCh,
       outputDocsImagesCh
     )
}

workflow.onComplete {


  // pipelineReport.html
  def reportFields = [:]
  reportFields['pipeline'] = workflow.manifest.name
  reportFields['version'] = workflow.manifest.version
  reportFields['runName'] = customRunName ?: workflow.runName
  reportFields['success'] = workflow.success
  reportFields['dateComplete'] = workflow.complete
  reportFields['duration'] = workflow.duration
  reportFields['exitStatus'] = workflow.exitStatus
  reportFields['errorMessage'] = (workflow.errorMessage ?: 'None')
  reportFields['errorReport'] = (workflow.errorReport ?: 'None')
  reportFields['commandLine'] = workflow.commandLine
  reportFields['projectDir'] = workflow.projectDir
  reportFields['summary'] = summary
  reportFields['summary']['Date Started'] = workflow.start
  reportFields['summary']['Date Completed'] = workflow.complete
  reportFields['summary']['Pipeline script file path'] = workflow.scriptFile
  reportFields['summary']['Pipeline script hash ID'] = workflow.scriptId
  if(workflow.repository) reportFields['summary']['Pipeline repository Git URL'] = workflow.repository
  if(workflow.commitId) reportFields['summary']['Pipeline repository Git Commit'] = workflow.commitId
  if(workflow.revision) reportFields['summary']['Pipeline Git branch/tag'] = workflow.revision

  // Render the TXT template
  def engine = new groovy.text.GStringTemplateEngine()
  def tf = new File("$projectDir/assets/workflowOnCompleteTemplate.txt")
  def txtTemplate = engine.createTemplate(tf).make(reportFields)
  def reportTxt = txtTemplate.toString()

  // Render the HTML template
  def hf = new File("$projectDir/assets/workflowOnCompleteTemplate.html")
  def htmlTemplate = engine.createTemplate(hf).make(reportFields)
  def reportHtml = htmlTemplate.toString()

  // Write summary HTML to a file
  def outputSummaryDir = new File( "${params.summaryDir}/" )
  if( !outputSummaryDir.exists() ) {
    outputSummaryDir.mkdirs()
  }
  def outputHtmlFile = new File( outputSummaryDir, "pipelineReport.html" )
  outputHtmlFile.withWriter { w -> w << reportHtml }
  def outputTxtFile = new File( outputSummaryDir, "pipelineReport.txt" )
  outputTxtFile.withWriter { w -> w << reportTxt }

  // workflowOnComplete file
  File woc = new File("${params.outDir}/workflowOnComplete.txt")
  Map endSummary = [:]
  endSummary['Completed on'] = workflow.complete
  endSummary['Duration']     = workflow.duration
  endSummary['Success']      = workflow.success
  endSummary['exit status']  = workflow.exitStatus
  endSummary['Error report'] = workflow.errorReport ?: '-'
  String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
  println endWfSummary
  String execInfo = "Execution summary\n${endWfSummary}\n"
  woc.write(execInfo)

  // final logs
  if(workflow.success){
    log.info "Pipeline Complete"
  }else{
    log.info "FAILED: $workflow.runName"
  }
}

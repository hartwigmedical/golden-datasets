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

// File with text to display when a developement version is used
devMessageFile = file("$projectDir/assets/devMessage.txt")

def helpMessage() {
  if ("${workflow.manifest.version}" =~ /dev/ ){
     log.info devMessageFile.text
  }

  log.info """
  @git_repo_name@ version: ${workflow.manifest.version}
  ======================================================================

    Usage:
    nextflow run main.nf -t, --truth truth_snv_file.vcf
        -s, --snv snv.vcf
        -i, --indel indel.vcf
        -m, --snvindel indels_and_vcfs.vcf
        -e, --target_bed bed file containing regions to assess
        -v, --sv sv.vcf
        -u, --truth_sv truth_sv.vcf
        -f, --fasta ref_fasta.fa
        -d, --outdir /OUTPUT_DIR/PATH
        -o, --outname output file name
        -n, --sname vcf sample name
        -a, --truth_sv_sname sample name for truth sv if different
        -b, --truth_snv_sname sample name for truth snv if different
        -p, --no_pass keep the pass variants and other variants in test file
        -p, --no_pass_truth keep the pass variants and other variants in truth file
        -c, --cpu number of threads
        -k, --keep (to keep intermediates files)


  Inputs:
    --snv [file]
    --indel [file]               Path to design file for extended analysis
    --sv [file]                  Specifies that the input is single-end reads
    --truth [file]
    --truthSv [file]
    --snvIndel [file]
    --fasta [file]
    --snvSname []
    --svSname
    --truthSnvSname
    --truthSvSname
    --targetBed [file]

  Main Options:
    --noPass [bool]
    --noPassTruth [bool]
    --keep [bool]


  Other options:
    --outDir [file]               The output directory where the results will be saved
    -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    --cpu

  ======================================================================
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


// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
customRunName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  customRunName = workflow.runName
}

/************
 * CHANNELS *
 ************/
// TODO: make sample plan work
samplePlanPathCh = params.samplePlan ? Channel.fromPath(params.samplePlan): Channel.empty()
samplePlanCh = samplePlanPathCh
     .splitCsv(header:true)
     .map { row ->
        [
            type: row.TYPE,
            sampleFile: file(row.PATH),
            sampleName: row.SNAME,
            truthFile: file(row.TRUTH),
            truthSampleName: row.TRUTH_SNAME
        ]
    }.dump(tag : "samplePlanCh")

// Branching samplePlan according to their compression state
samplePlanCh.branch{ row ->
    rawSamplesCh: row.sampleFile=~ /.*vcf$/ || row.truthFile=~ /.*vcf$/
    compressSamplesCh: row.sampleFile=~ /.*gz$/ && row.truthFile=~ /.*gz$/
}.set { samplePlanForks }

(rawSamplesCh, compressSamplesCh) = [
    samplePlanForks.rawSamplesCh.dump(tag: "rawSamplesCh"),
    samplePlanForks.compressSamplesCh.dump(tag: "compressSamplesCh")
]

snvCh = params.snv ? Channel.value(file(params.snv)) : Channel.empty()
indelCh = params.indel ? Channel.value(file(params.indel)) : Channel.empty()
snvIndelCh = params.snvIndel ? Channel.value(file(params.snvIndel)) : Channel.empty()
svCh = params.sv ? Channel.value(file(params.sv)) : Channel.empty()
truthCh = params.truth ? Channel.value(file(params.truth)) : Channel.empty()
truthSvCh = params.truthSv ? Channel.value(file(params.truthSv)) : Channel.empty()

snvSnameCh = params.snvSname ? Channel.value(params.snvSname) : Channel.empty()
indelSnameCh = params.indelSname ? Channel.value(params.indelSname) : Channel.empty()
snvIndelSnameCh = params.snvIndelSname ? Channel.value(params.snvIndelSname) : Channel.empty()
svSnameCh = params.svSname ? Channel.value(params.svSname) : Channel.empty()
truthSnvSnameCh = params.truthSnvSname ? Channel.value(params.truthSnvSname) : Channel.empty()
truthSvSnameCh = params.truthSvSname ? Channel.value(params.truthSvSname) : Channel.empty()

fastaCh = params.fasta ? Channel.value(file(params.fasta)) : Channel.empty()
targetBedCh = params.targetBed ? Channel.value(file(params.targetBed)) : "null"

outNameCh = params.outName ? Channel.value(params.outName) : Channel.empty()
noPassCh = params.noPassCh ? : false
noPassTruthCh = params.noPassTruthCh ? : false
noPassTruthSvCh = params.noPassTruthSvCh ? : false
//keepCh = params.keep ? : false

outNameCh.dump(tag: "outNameCh")
/*******************
 * Header log info *
 *******************/

if ("${workflow.manifest.version}" =~ /dev/ ){
  log.info devMessageFile.text
}

log.info """\
=======================================================
@git_repo_name@ workflow version: ${workflow.manifest.version}
=======================================================\
""".stripIndent()

summary = [
  'Pipeline Release': workflow.revision ?: null,
  'Run Name': customRunName,
  'Max Resources': "${params.maxMemory} memory, ${params.maxCpus} cpus, ${params.maxTime} time per job",
  'Container': workflow.containerEngine && workflow.container ? "${workflow.containerEngine} - ${workflow.container}" : null,
  'Fasta': params.fasta ?: null,
  'Target BED': params.targetBed ?: null,
  //'SV ':
  //'SNV':
  'Output dir': params.outDir,
  'Launch dir': workflow.launchDir,
  'Working dir': workflow.workDir,
  'Script dir': workflow.projectDir,
  'User': workflow.userName,
  'Config Profile': workflow.profile,
  'Config Description': params.configProfileDescription ?: null,
  'Config Contact': params.configProfileContact ?: null,
  'Config URL': params.configProfileUrl ?: null,
].findAll{ it.value != null }

log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "======================================================="

/************
 * PROCESS *
 ************/

 /********************
  * SNV Benchmarking *
  ********************/

// TODO: branch sample plan for bgzip

// process compressVcf {
//     label "onlyLinux"
//     label "lowCpu"
//     label "lowMem"
//
//     input:  $snvindel $truth
//     output: $snvindel $truth
//     script:
//
//     """
//     gzip $OUTPUT_DIR/snv_indel_temp.vcf
//     snvindel=$OUTPUT_DIR/snv_indel_temp.vcf.gz
//
//     gzip $OUTPUT_DIR/truth_temp.vcf
//     truth=$OUTPUT_DIR/truth_temp.vcf.gz
//     """
// }
//
// process mergeSnvIndel {
//     input: $snv $indel
//     output: snvindel=$OUTPUT_DIR/"snv_indel_temp.vcf.gz"
//     script:
//     """
//     bgzip -@ $CPU -c $snv > $snv".gz"
//     bgzip -@ $CPU -c $indel > $indel".gz"
//     snv=$snv".gz"
//     indel=$indel".gz"
//
//     bcftools index --threads $CPU -f -o $snv".csi" $snv
//     bcftools index --threads $CPU -f -o $indel".csi" $indel
//
//     bcftools concat -a $snv $indel -O z -o $OUTPUT_DIR/"snv_indel_temp.vcf.gz" --threads $CPU
//     snvindel=$OUTPUT_DIR/"snv_indel_temp.vcf.gz"
//     """
//
// }
//
process chrHandling {
    input:
    file(snvIndel) from snvIndelCh
    file(truth) from truthCh

    output:
    file("snv_indel_temp.vcf.gz") into snvIndelChrCh
    file("truth_temp.vcf.gz") into truthChrCh

    script:
    """
    zcat ${snvIndel} | awk '{if(\$0 !~ /^#/) print "chr"\$0; else print \$0}' | awk '{gsub(/contig=\\<ID=/,"contig=<ID=chr"); print}' | awk '{gsub(/chrchr/,"chr"); print}' > snv_indel_temp.vcf
    gzip snv_indel_temp.vcf

    zcat ${truth} | awk '{if(\$0 !~ /^#/) print "chr"\$0; else print \$0}'| awk '{gsub(/contig=\\<ID=/,"contig=<ID=chr"); print}' | awk '{gsub(/chrchr/,"chr"); print}' > truth_temp.vcf
    gzip truth_temp.vcf
    """
}


process splitMultiSample {
    label "onlyLinux"
    label "medCpu"
    label "lowMem"

    input:
    file(snvIndel) from snvIndelChrCh
    file(truth) from truthChrCh
    val(sname) from snvIndelSnameCh
    val(truthSnvSname) from truthSnvSnameCh

    output:
    file("snv_indel_temp.sample.vcf.gz") into snvIndelSampleCh
    file("truth_temp.sample.vcf.gz") into truthSampleCh

    script:
    """
    echo ${sname}
    bcftools view -c1 --threads ${task.cpus} -O z -s ${sname} -o snv_indel_temp.sample.vcf.gz ${snvIndel}

    bcftools view -c1 --threads ${task.cpus} -O z -s ${truthSnvSname} -o truth_temp.sample.vcf.gz ${truth}
    """
}

process filterPASS {
    input:
    file(snvIndel) from snvIndelSampleCh
    file(truth) from truthSampleCh
    file(truthSv) from truthSvCh

    output:
    file("snv_indel.pass.vcf") into snvIndelPassCh
    file("truth_temp.pass.vcf") into truthSamplePassCh
    file("truth_temp_sv.pass.vcf") into truthSvPassCh

    script:
    """
    zcat ${snvIndel} | grep "PASS\\|#" > snv_indel.pass.vcf

    zcat ${truth} | grep "PASS\\|#" > truth_temp.pass.vcf

    zcat ${truthSv} | grep "PASS\\|#" > truth_temp_sv.pass.vcf
    """
}

process PrepareAndNormalize {
    label "onlyLinux"
    label "medCpu"
    label "lowMem"

    input:
    file(snvIndel) from snvIndelPassCh
    file(truth) from truthSamplePassCh
    file(fasta) from fastaCh

    output:
    file("snv_indel.pass.sort.prep.norm.vcf.gz") into snvIndelNormCh
    file("truth_temp.sort.prep.norm.vcf.gz") into truthSampleNormCh

    script:
    """
    bcftools sort ${snvIndel} -o snv_indel.pass.sort.vcf.gz -O z
    bcftools sort ${truth} -o truth_temp.pass.sort.vcf.gz -O z

    snvindel="snv_indel.pass.sort.vcf.gz"
    truth="truth_temp.pass.sort.vcf.gz"

    bcftools index -f -o \$snvindel".csi" \$snvindel --threads ${task.cpus}
    bcftools index -f -o \$truth".csi" \$truth --threads ${task.cpus}

    export HGREF=${fasta}

    multimerge \$snvindel -r \$HGREF -o snv_indel.pass.sort.prep.vcf.gz --calls-only=1
    multimerge \$truth -r \$HGREF -o truth_temp.sort.prep.vcf.gz --calls-only=1

    snvindel="snv_indel.pass.sort.prep.vcf.gz"
    truth="truth_temp.sort.prep.vcf.gz"

    pre.py \$snvindel snv_indel.pass.sort.prep.norm.vcf.gz -L --decompose --somatic --threads ${task.cpus}

    pre.py \$truth truth_temp.sort.prep.norm.vcf.gz -L --decompose --somatic --threads ${task.cpus}
    """
}

process benchmarkSNV {
    //TODO: add CPU label
    input:
    file(snvIndel) from snvIndelNormCh
    file(truth) from truthSampleNormCh
    file(target) from targetBedCh
    file(fasta) from fastaCh
    val(outName) from outNameCh

    output:
    file("*.stats.csv") into sompySnvCh

    script:
    if (params.targetBed)
    """
    export HGREF=${fasta}

    som.py ${truth} ${snvIndel} -o ${outName} -N -R ${target}
    """
    else
    """
    export HGREF=${fasta}

    som.py ${truth} ${snvIndel} -o ${outName} -N
    """
}

/********************
 * SV Benchmarking *
 ********************/
 process ingestSV {
     input:
     file(sv) from svCh
     val(svSname) from svSnameCh

     output:
     file("sv_dataframe.csv") into svDf

     script:
     svSnameArg = params.svSname ? " -samplename ${svSname} ": ''
     """
     ingest.py ${sv} ${svSnameArg} -outputfile "sv_dataframe.csv" -filter
     """
 }

 process ingestSvTruth {
     input:
     file(truth) from truthSvCh
     val(truthSvSname) from truthSvSnameCh

     output:
     file("truth_sv_dataframe.csv") into truthSvDf

     script:
     truthSvSnameArg = params.truthSvSname ? " -samplename ${truthSvSname} ": ''
     """
     ingest.py ${truth} ${truthSvSnameArg} -outputfile "truth_sv_dataframe.csv" -filter
     """
 }

 process compareTruthSv {
     input:
     file(svDf) from svDf
     val(truthSvDf) from truthSvDf

     output:
     file("SV_benchmark_results.csv") into truthSvDf

     script:
     truthSvSnameArg = params.truthSvSname ? " -samplename ${truthSvSname} ": ''
     """
     compare_node_to_truth.py ${svDf} ${truthSvDf} -metrics SV_benchmark_results.csv
     """
 }

/**********
 * FastQC *
 **********/

// process fastqc {
//   label 'fastqc'
//   label 'lowMem'
//   label 'lowCpu'
//
//   tag "${prefix}"
//   publishDir "${params.outDir}/fastqc", mode: 'copy'
//
//   input:
//   set val(prefix), file(reads) from rawReadsFastqcCh
//
//   output:
//   file "*_fastqc.{zip,html}" into fastqcResultsCh
//   file "v_fastqc.txt" into fastqcVersionCh
//
//   script:
//   """
//   fastqc -q $reads
//   fastqc --version > v_fastqc.txt
//   """
// }


workflow.onComplete {

  // final logs
  if(workflow.success){
    log.info "Pipeline Complete"
  }else{
    log.info "FAILED: $workflow.runName"
  }
}

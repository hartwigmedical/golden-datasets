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
    nextflow run main.nf \
        -profile multiconda \
        --snvIndel myvariants.vcf.gz \
        --truth truth.snv_indel.vcf.gz \
        --snvIndelSname mySname \
        --truthSnvSname myTruthSname \
        --sv SV.vcf.gz \
        --svSname svSname \
        --truthSv truth.sv.vcf.gz \
        --truthSvSname truthSvSname \
        --fasta hg19.fa \
        --targetBed targets.bed \
        --outName test \
        --condaCacheDir condaEnvsFolder \
        -resume

  Inputs:
    --snv [file]                Snv file if snv & indel are splitted (.vcf)
    --indel [file]              Indel file if snv & indel are splitted (.vcf)
    --snvIndel [file]           Snv & Indel file (.vcf)
    --sv [file]                 Structural variant file (.vcf or .tsv)
    --truth [file]              Truth file for SNV & Indel Calling (.vcf)
    --truthSv [file]            Truth file for SV Calling (.vcf or .tsv)
    --fasta [file]              Reference genome file (.fa)
    --snvSname [str]            Sample Name found in the snvIndel file
    --svSname [str]             Sample Name found in the SV vcf
    --truthSnvSname [str]       Sample Name found in the snvIndel truth file
    --truthSvSname [str]        Sample Name found in the SV truth file
    --targetBed [file]          Target file if WES / Capture (.bed)

  Main Options:
    --noPass [bool]             Keep all variants from test files even those without PASS tag
    --noPassTruth [bool]        Keep all variants from truth files even those without PASS tag


  Other options:
    --outDir [file]               The output directory where the results will be saved
    -name [str]                   Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    --cpu [str]                 Number of CPU available for analysis

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
// samplePlanPathCh = params.samplePlan ? Channel.fromPath(params.samplePlan): Channel.empty()
// samplePlanCh = samplePlanPathCh
//      .splitCsv(header:true)
//      .map { row ->
//         [
//             type: row.TYPE,
//             sampleFile: file(row.PATH),
//             sampleName: row.SNAME,
//             truthFile: file(row.TRUTH),
//             truthSampleName: row.TRUTH_SNAME
//         ]
//     }.dump(tag : "samplePlanCh")
//
// // Branching samplePlan according to their compression state
// samplePlanCh.branch{ row ->
//     rawSamplesCh: row.sampleFile=~ /.*vcf$/ || row.truthFile=~ /.*vcf$/
//     compressSamplesCh: row.sampleFile=~ /.*gz$/ && row.truthFile=~ /.*gz$/
// }.set { samplePlanForks }
//
// (rawSamplesCh, compressSamplesCh) = [
//     samplePlanForks.rawSamplesCh.dump(tag: "rawSamplesCh"),
//     samplePlanForks.compressSamplesCh.dump(tag: "compressSamplesCh")
// ]

snvCh = params.snv ? Channel.of(file(params.snv)) : Channel.empty()
indelCh = params.indel ? Channel.of(file(params.indel)) : Channel.empty()
snvIndelCh = params.snvIndel ? Channel.of(file(params.snvIndel)) : Channel.empty()
svCh = params.sv ? Channel.of(file(params.sv)) : Channel.empty()
truthCh = params.truth ? Channel.of(file(params.truth)) : Channel.empty()
truthSvCh = params.truthSv ? Channel.of(file(params.truthSv)) : Channel.empty()

snvSnameCh = params.snvSname ? Channel.value(params.snvSname) : Channel.empty()
indelSnameCh = params.indelSname ? Channel.value(params.indelSname) : Channel.empty()
snvIndelSnameCh = params.snvIndelSname ? Channel.value(params.snvIndelSname) : Channel.empty()
svSnameCh = params.svSname ? Channel.value(params.svSname) : Channel.empty()
truthSnvSnameCh = params.truthSnvSname ? Channel.value(params.truthSnvSname) : Channel.empty()
truthSvSnameCh = params.truthSvSname ? Channel.value(params.truthSvSname) : Channel.empty()

fastaCh = params.fasta ? Channel.value(file(params.fasta)) : Channel.empty()
targetBedCh = params.targetBed ? Channel.value(file(params.targetBed)) : "null"

outNameCh = params.outName ? Channel.value(params.outName) : Channel.empty()
//keepCh = params.keep ? : false

//outNameCh.dump(tag: "outNameCh")
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
annotateInputCh = snvIndelCh.combine(["snvIndel"]).mix(snvCh.combine(["snv"])).mix(indelCh.combine(["indel"])).mix(truthCh.combine(["truth"])).dump(tag:"testCh")

//(annotateSnvAndIndelCh, annotateSnvIndelCh) = annotateInputCh.into(2)
annotateSnvAndIndelCh = annotateInputCh

// annotateSnameInputCh = snvIndelSnameCh.combine(["snvIndel"]).mix(snvSnameCh.combine(["snv"])).mix(indelSnameCh.combine(["indel"])).mix(truthSnvSnameCh.combine(["truth"])).dump(tag:"testSnameCh")

// (annotateSnvAndIndelSnameCh, annotateSnvIndelSnameCh) = annotateSnameInputCh.into(2)

renameChrCh = Channel.fromPath("assets/rename_chr.txt")

//(renameChrSnvAndIndelCh, renameChrSnvIndelCh) = renameChrCh.into(2)
renameChrSnvAndIndelCh = renameChrCh

// process annotate {
//     label "snv"
//
//     input:
//     tuple file(vcf), val(annotateType) from annotateInputCh
//     file (chrs) from renameChrCh
//
//     output:
//     file("snvIndel.annotate.vcf.gz") into snvIndelAnnotateCh
//     file("snv.annotate.vcf.gz") into snvAnnotateCh
//     file("indel.annotate.vcf.gz") into indelAnnotateCh
//     file("truth.annotate.vcf.gz") into truthAnnotateCh
//
//     script:
//     """
//     bcftools annotate -x INFO,^FORMAT/GT -o ${annotateType}.annotate.vcf.gz -Oz --rename-chrs rename_chr.txt ${vcf}
//     """
// }

//annotateSnvAndIndelCh = annotateSnvAndIndelCh.dump(tag: "annotateSnvAndIndelCh")

process checkSnvAndIndel {
    label "snv"

    input:
    tuple file(vcf), val(annotateType), file(chrs) from annotateSnvAndIndelCh.combine(renameChrSnvAndIndelCh)

    output:
    tuple val(annotateType), file("${annotateType}.annotate.vcf.gz") into annotateCh

    script:
    """
    bcftools annotate -x INFO,^FORMAT/GT -o ${annotateType}.annotate.vcf.gz -Oz --rename-chrs rename_chr.txt ${vcf}
    """
}

//annotateCh = annotateCh.dump(tag: "annotateCh")

annotateCh.branch{ row ->
    snvIndelCh: row[0]=~/snvIndel/
        return row[1]
    snvCh: row[0]=~/snv/
        return row[1]
    indelCh: row[0]=~/indel/
        return row[1]
    truthCh: row[0]=~/truth/
        return row[1]
    otherCh: true
}.set { annotateForks }

(snvCheckedCh, indelCheckedCh, snvIndelCheckedCh, truthCheckedCh, otherCheckedCh) = [annotateForks.snvCh,annotateForks.indelCh,annotateForks.snvIndelCh,annotateForks.truthCh,annotateForks.otherCh]

process concat {
    label "snv"

    input:
    file(snv) from snvCheckedCh
    file(indel) from indelCheckedCh

    output:
    file("snv_indel.vcf.gz") into snvIndelMergedCh

    script:
    """
    bcftools index ${snv}
    bcftools index ${indel}
    bcftools concat -a ${snv} ${indel} -O z -o snv_indel.vcf.gz --threads ${task.cpus}
    """
}

// process checkSnvIndel {
//     label "snv"
//
//     input:
//     tuple file(vcf), val(annotateType) from annotateSnvIndelCh
//     file (chrs) from renameChrSnvIndelCh
//
//     output:
//     file("snvIndel.annotate.vcf.gz") into snvIndelAnnotateCh
//     file("truth.annotate.vcf.gz") into truthAnnotateCh2
//
//     when:
//     params.snvIndel
//
//     script:
//     """
//     bcftools annotate -x INFO,^FORMAT/GT -o ${annotateType}.annotate.vcf.gz -Oz --rename-chrs rename_chr.txt ${vcf}
//     """
// }
// process chrHandling {
//     label "snv"
//
//     input:
//     file(snvIndel) from snvIndelCh
//     file(truth) from truthCh
//
//     output:
//     file("snv_indel_temp.vcf.gz") into snvIndelChrCh
//     file("truth_temp.vcf.gz") into truthChrCh
//
//     script:
//     """
//     zcat ${snvIndel} | awk '{if(\$0 !~ /^#/) print "chr"\$0; else print \$0}' | awk '{gsub(/contig=\\<ID=/,"contig=<ID=chr"); print}' | awk '{gsub(/chrchr/,"chr"); print}' > snv_indel_temp.vcf
//     gzip snv_indel_temp.vcf
//
//     zcat ${truth} | awk '{if(\$0 !~ /^#/) print "chr"\$0; else print \$0}'| awk '{gsub(/contig=\\<ID=/,"contig=<ID=chr"); print}' | awk '{gsub(/chrchr/,"chr"); print}' > truth_temp.vcf
//     gzip truth_temp.vcf
//     """
// }

// TODO: code error message if params.SNV &  params.Indel & params.SnvIndel

//snvIndelMergedCh = snvIndelMergedCh.dump(tag: "snvIndelMergedCh")
//snvIndelCheckedNewCh = snvIndelCheckedCh.ifEmpty([]).dump(tag: "snvIndelCheckedCh")
//annotateCh = annotateCh.dump(tag: "annotateCh")

//(snvIndelCheckedCh,snvIndelCheckedTestCh) = snvIndelCheckedCh.into(2)
//(snvIndelMergedCh,snvIndelMergedTestCh) = snvIndelMergedCh.into(2)

//snvIndelMergedTestCh.concat(snvIndelCheckedTestCh).dump(tag: "mixCh")


process splitMultiSample {
    label "medCpu"
    label "lowMem"
    label "snv"

    input:
    file(snvIndel) from snvIndelMergedCh.concat(snvIndelCheckedCh)
    file(truth) from truthCheckedCh
    val(sname) from snvIndelSnameCh
    val(truthSnvSname) from truthSnvSnameCh

    output:
    file("snv_indel.sample.vcf.gz") into snvIndelSampleCh
    file("truth.sample.vcf.gz") into truthSampleCh

    script:
    """
    echo ${sname}
    bcftools view -c1 --threads ${task.cpus} -O z -s ${sname} -o snv_indel.sample.vcf.gz ${snvIndel}

    bcftools view -c1 --threads ${task.cpus} -O z -s ${truthSnvSname} -o truth.sample.vcf.gz ${truth}
    """
}

process filterPASS {
    label "snv"

    input:
    file snvIndel name "snv_indel_temp.sample.vcf.gz" from snvIndelSampleCh
    file truth name "truth_temp.sample.vcf.gz" from truthSampleCh
    file truthSv name "truth_sv.vcf.gz" from truthSvCh

    output:
    file(params.noPass ? "snv_indel_temp.sample.vcf.gz" : "snv_indel.pass.vcf") into snvIndelPassCh
    file(params.noPassTruth ? "truth_temp.sample.vcf.gz":"truth_temp.pass.vcf") into truthSamplePassCh
    file(params.noPassTruthSv ? "truth_sv.vcf.gz" :"truth_temp_sv.pass.vcf") into truthSvPassCh

    script:
    passSnvCmd = params.noPass ? '': "zcat ${snvIndel} | grep 'PASS\\|#' > snv_indel.pass.vcf"
    passTruthCmd = params.noPassTruth ? '': "zcat ${truth} | grep 'PASS\\|#' > truth_temp.pass.vcf"
    passTruthSvCmd = params.noPassTruthSv ? '':"zcat ${truthSv} | grep 'PASS\\|#' > truth_temp_sv.pass.vcf"
    """
    ${passSnvCmd}
    ${passTruthCmd}
    ${passTruthSvCmd}
    """
}


process PrepareAndNormalize {
    label "medCpu"
    label "lowMem"
    label "snv"

    input:
    file(snvIndel) from snvIndelPassCh.dump(tag: "prepNormsnv")
    file(truth) from truthSamplePassCh.dump(tag: "prepNormtruth")
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
    label "snv"

    publishDir "${params.outDir}/SNV", mode: params.publishDirMode

    //TODO: add CPU label
    input:
    file(snvIndel) from snvIndelNormCh
    file(truth) from truthSampleNormCh
    file(target) from targetBedCh
    file(fasta) from fastaCh
    val(outName) from outNameCh

    output:
    file("*.stats.csv") into sompySnvCh,plotSnvCh

    script:
    targetArg = params.targetBed ? " -R ${target} ": ''
    noPassArg = params.noPass ? " -P ": ''
    """
    export HGREF=${fasta}

    som.py ${truth} ${snvIndel} -o ${outName} -N  ${noPassArg}  ${targetArg}
    """
}

/********************
 * SV Benchmarking *
 ********************/
 process ingestSV {
     label "sv"

     publishDir "${params.outDir}/SV", mode: params.publishDirMode

     input:
     file(sv) from svCh
     val(svSname) from svSnameCh

     output:
     file("sv_dataframe.csv") into svDf,plotSvTestCh

     script:
     svSnameArg = params.svSname ? " -samplename ${svSname} ": ''
     """
     ingest.py ${sv} ${svSnameArg} -outputfile "sv_dataframe.csv" -filter
     """
 }

 process ingestSvTruth {
     label "sv"

     publishDir "${params.outDir}/SV", mode: params.publishDirMode

     input:
     file(truth) from truthSvPassCh
     val(truthSvSname) from truthSvSnameCh

     output:
     file("truth_sv_dataframe.csv") into truthSvDf,plotSvTruthCh

     script:
     truthSvSnameArg = params.truthSvSname ? " -samplename ${truthSvSname} ": ''
     """
     ingest.py ${truth} ${truthSvSnameArg} -outputfile "truth_sv_dataframe.csv" -filter
     """
 }

 process benchmarkSV {
     label "sv"

     publishDir "${params.outDir}/SV", mode: params.publishDirMode

     input:
     file(svDf) from svDf
     file(truthSvDf) from truthSvDf

     output:
     file("SV_benchmark_results.csv") into svResultsCh,plotSvCh

     script:
     truthSvSnameArg = params.truthSvSname ? " -samplename ${params.truthSvSname} ": ''
     """
     compare_node_to_truth.py ${svDf} ${truthSvDf} -metrics SV_benchmark_results.csv
     """
 }

 process plotSNV {
     label "plots"

     publishDir "${params.outDir}/SNV", mode: params.publishDirMode

     input:
     file(statsSNV) from plotSnvCh
     val(outName) from outNameCh

     output:
     file("*") into snvPlotsCh

     script:
     """
     plots_benchmarking_snv.R -b ${statsSNV} -o ${outName}
     """
 }

 process plotSV {
     label "plots"

     publishDir "${params.outDir}/SV", mode: params.publishDirMode

     input:
     file(statsSV) from plotSvCh
     file(truthSvDf) from plotSvTruthCh
     val(outName) from outNameCh

     output:
     file("*") into svPlotsCh

     script:
     """
     plots_benchmarking_SV.R -b ${statsSV} -t ${truthSvDf} -o ${outName}
     """
 }

workflow.onComplete {

  // final logs
  if(workflow.success){
    log.info "Pipeline Complete"
  }else{
    log.info "FAILED: $workflow.runName"
  }
}

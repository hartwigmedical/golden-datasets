/*
Copyright Institut Curie 2020-2021

This software is a computer program whose purpose is to
provide a demo with geniac.
You can use, modify and/ or redistribute the software under the terms
of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful,
but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards
their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge
of the license and that you accept its terms.

*/


/***********
 * MultiQC *
 ***********/

process multiqc {
  label 'multiqc'
  label 'minCpu'
  label 'minMem'
  publishDir "${params.outDir}/MultiQC", mode: 'copy'

  when:
  !params.skipMultiQC

  input:
  val customRunName
  path splan 
  path multiqcConfig 
  path ('fastqc/*') 
  path metadata 
  path ('software_versions/*') 
  path ('workflow_summary/*') 

  output:
  path splan
  path "*report.html", emit: multiqc_report
  path "*_data"

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_report" : "--filename report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  designOpts= params.design ? "-d ${params.design}" : ""
  modulesList = "-m custom_content -m fastqc"
  """
  apMqcHeader.py --splan ${splan} --name "${workflow.manifest.name}" --version "${workflow.manifest.version}" ${metadataOpts} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c multiqc-config-header.yaml -c $multiqcConfig $modulesList
  """
}


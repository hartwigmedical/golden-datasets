# EUCANCan Nextflow Benchmarking pipeline

This is a pipeline to benchmark Single Nucleotide Variants (SNV) and Structural Variants (SV) and create a visualization of the results.

## Pipeline summary:

**SNV:**
1. Check files validity:
2. Concat SNV and Indel calls if they are in two different vcf files
3. Split multisample VCF in order to consider only the somatic sample
4. Extract only PASS variants
5. normalize test and truth files
6. Benchmark test file using som.py tool
7. Visualization of the results with R scriptss

**SV:**

1. transform SV vcf or tsv truth and test files into a well defined table used for benchmarkinkg
2. compare test and truth tables and apply tier rules
3. Visualization with R scripts

## Quick help:
```
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
--keep [bool]


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
```

## Notes:
- you can run the SNV part or the SV paft or both at the same time
- if SNV and Indel files are separated, be sure that the sname inside the vcf are the same


## Defining the '-profile'
By default (whithout any profile), Nextflow will excute the pipeline locally, expecting that all tools are available from your `PATH` variable.
In addition, we set up a few profiles that should allow you i/ to use containers instead of local installation, ii/ to run the pipeline on a cluster instead of on a local architecture.
The description of each profile is available on the help message (see above).

Here are a few examples of how to set the profile option. See the [full documentation](docs/profiles) for details.

```
## Run the pipeline locally, using a global environment where all tools are installed (build by conda for instance)
-profile path --globalPath INSTALLATION_PATH

## Run the pipeline on the cluster, using the Singularity containers
-profile cluster,singularity --singularityPath SINGULARITY_PATH

## Run the pipeline on the cluster, building a new conda environment
-profile cluster,conda --condaCacheDir CONDA_CACHE

```

### Full Documentation

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Profiles](docs/profiles.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

#### Credits

This pipeline has been written by the WP2 -EUCANCan team (F. Allain, T. Gutman, D. VanBeek, À.Ferriz) and was funded by the European Union’s Horizon 2020 research and innovation programme and the Canadian Institutes of Health Research under the grant agreement No 825835 in the framework of the [European-Canadian Cancer Network](https://eucancan.com/)

#### Contacts

For any question, bug or suggestion, please send an issue or contact the WP2.

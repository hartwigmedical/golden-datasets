# Configuration variables for the Docker Nexflow wrapper
# 
# Copy this file as `config.sh` to a directory on your laptop and adjust as necessary
# Expose it to the Docker container at `/config` by including in your `docker run` command:
#   docker run -v /path/to/config_dir_on_your_laptop:/config IMAGE_ID ...
#
# MAke sure all the locations refer to the path they will be mapped to on the Docker container.
# eg if they are in the same directory as the edited config.sh make sure the path is `/config/`.
#
# Additional volume mounts can be provided to `docker run` if needed.
PIPELINE_DIR=/home/daphnevanbeek/golden-datasets/pipeline
CONDA=/home/daphnevanbeek/conda-cache-nextflow
FASTA=/home/daphnevanbeek/data/daphne-eucancan-benchmark-data/hg19_curie.fa
SNVINDEL=/home/daphnevanbeek/data/daphne-eucancan-benchmark-data/COLO_RESULTS/original/Hartwig/5-22/purple/COLO829T.purple.somatic.vcf
TRUTH=/home/daphnevanbeek/data/daphne-eucancan-benchmark-data/COLO_SNV_referencedata/EGAF00001181170/ListforNatureReports.IndelsandSNVs.final.Suppl1.snpEff.validated.vcf
SNVINDELSNAME=COLO829T
TRUTHSNVSNAME=COLO_829_Illumina
SV=/home/daphnevanbeek/data/daphne-eucancan-benchmark-data/COLO_RESULTS/original/Hartwig/5-22/purple/COLO829T.purple.sv.vcf.gz
SVSNAME=COLO829T
TRUTHSV=/home/daphnevanbeek/data/daphne-eucancan-benchmark-data/COLO_RESULTS/truth_file/COLO829_truth.tsv
#TRUTHSVSNAME=XX
OUTNAME=COLO829_hartwig_522_nextflow
OUTDIR=/home/daphnevanbeek/test_nf_bench


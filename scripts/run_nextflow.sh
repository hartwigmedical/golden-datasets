# General parameters
PIPELINE_DIR=/pipeline/location
CONDA=/Place/to/store/conda/images
FASTA=/Reference/Genome/fasta
OUTNAME=
OUTDIR=

# SNV related parameters
SNVINDEL=/PATH/SNV_INDEL_MIXED_FILE.vcf
#OR
# (add --snv --indel params in command line and remove --snvindel)
SNV=/PATH/SNV_FILE.vcf
INDEL=/PATH/INDEL_FILE.vcf

TRUTH=/TRUTH_FILE.vcf

SNVINDELSNAME=/PATH/SNV_INDEL_TUMOR_SAMPLE_NAME
TRUTHSNVSNAME=/PATH/TRUTH_SAMPLE_NAME

# SV related parameters
SV=/PATH/NV_FILE #(tsv,vcf)
SVSNAME=/PATH/SV_SAMPLE_NAME #(not required if file is tsv)

TRUTHSV=/TRUTH_SV_FILE #(vcf,tsv)
TRUTHSVSNAME=XX #(not required if file is tsv, add --truthSvSname in command line)

# Command line
nextflow run $PIPELINE_DIR/main.nf \
        -profile multiconda \
        --snvIndel $SNVINDEL \
        --truth $TRUTH \
        --snvIndelSname $SNVINDELSNAME \
        --truthSnvSname $TRUTHSNVSNAME \
        --sv $SV \
        --truthSv $TRUTHSV \
        --svSname $SVSNAME \
        --fasta $FASTA \
        --outName $OUTNAME \
        --outDir $OUTDIR \
        -w $OUTDIR/$OUTNAME/work \
        --condaCacheDir $CONDA \
        -resume

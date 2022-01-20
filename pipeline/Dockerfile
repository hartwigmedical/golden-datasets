# Set the base image to
FROM google/cloud-sdk:latest
ADD Miniconda3-latest-Linux-x86_64.sh Miniconda3-latest-Linux-x86_64.sh
RUN yes yes | bash Miniconda3-latest-Linux-x86_64.sh
RUN curl -s https://get.nextflow.io | bash

ENTRYPOINT ./nextflow run $PIPELINE_DIR/main.nf \
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
        --condaCacheDir $CONDA
#!/bin/bash --login
# The --login ensures the bash configuration is loaded,
set -euo pipefail
#EUCANCAN SNV & INDEL vcf handling
KEEP=false
NOPASS=false
NOPASS_TRUTH=false
CPU=1
while [ $# -gt 0 ] ; do
  case $1 in
      -t | --truth) truth="$2";;
      -s | --snv) snv="$2";;
      -i | --indel) indel="$2";;
      -m | --snvindel) snvindel="$2";;
      -e | --target_bed) target="$2";;
      -v | --sv) sv="$2";;
      -u | --truth_sv) truth_sv="$2";;
      -f | --fasta) FASTA="$2";;
      -d | --outdir) OUTPUT_DIR="$2";;
      -o | --outname) OUT_NAME="$2";;
      -n | --sname) SAMPLE_NAME="$2";;
      -a | --truth_sv_sname) SV_SAMPLE_NAME="$2";;
      -b | --truth_snv_sname) SNV_SAMPLE_NAME="$2";;
      -p | --no_pass) NOPASS=true;;
      -q | --no_pass_truth) NOPASS_TRUTH=true;;
      -c | --cpu) CPU="$2";;
      -k | --keep) KEEP=true;;
      -h | --help)  echo "Usage:"
          echo "bash ingest_snv.sh -t, --truth truth_snv_file.vcf"
          echo "                   -s, --snv snv.vcf"
          echo "                   -i, --indel indel.vcf"
          echo "                   -m, --snvindel indels_and_vcfs.vcf"
          echo "                   -e, --target_bed bed file containing regions to assess"
          echo "                   -v, --sv sv.vcf"
          echo "                   -u, --truth_sv truth_sv.vcf"
          echo "                   -f, --fasta ref_fasta.fa"
          echo "                   -d, --outdir /OUTPUT_DIR/PATH"
          echo "                   -o, --outname output file name"
          echo "                   -n, --sname vcf sample name"
          echo "                   -a, --truth_sv_sname sample name for truth sv if different"
          echo "                   -b, --truth_snv_sname sample name for truth snv if different"
          echo "                   -p, --no_pass keep the pass variants and other variants in test file"
          echo "                   -p, --no_pass_truth keep the pass variants and other variants in truth file"
          echo "                   -c, --cpu number of threads"
          echo "                   -k, --keep (to keep intermediates files)"
          exit

  esac
  shift
done

echo " "
echo -e "[General Information]:\n"
echo -e "[General parameters]:"
echo "output path:" $OUTPUT_DIR
echo "output file Name:" $OUT_NAME
echo "keep intermediate files ?:" $KEEP
echo "Number of threads:" $CPU
echo "keep pass & other variants in all test file?" $NOPASS
echo "keep pass & other variants in all truth file?" $NOPASS_TRUTH

echo -e "\n[SNV Benchmarking]:"
echo "Truth File:" $truth
echo "SNV vcf file:" $snv
echo "INDEL vcf file:" $indel
echo "SNV + INDEL vcf file:" $snvindel
echo "vcf sample name:" $SAMPLE_NAME
echo "truth snv sample name:" $SNV_SAMPLE_NAME
echo "target bed:" $target
echo "Reference fasta file:" $FASTA

echo -e "\n[SV Benchmarking]:"
echo "SV vcf file:" $sv
echo "TRUTH SV file:" $truth_sv
echo "truth sv sample name:" $SV_SAMPLE_NAME
echo " "

# Create output dir:

mkdir -p $OUTPUT_DIR/$OUT_NAME
OUTPUT_DIR=$OUTPUT_DIR/$OUT_NAME

# Source correct Conda env
source ~/miniconda3/etc/profile.d/conda.sh

# Warnings:

if [[ -n "$truth" && ( -n "$snv" && -n "$indel" || -n "$snvindel") ]];then

    if [[ -n "$snv" && -n "$indel" && -n "$snvindel" ]];then
        echo "[WARNING] too many arguments! Specify SNV and INDEL or SNVINDEL."
        exit
    elif [[ (! -n "$snv" || ! -n "$indel") && ! -n "$snvindel" ]];then
        echo "[WARNING] argument SNV or INDEL missing."
        exit
    elif [[ (-n "$snv" || -n "$indel") &&  -n "$snvindel" ]];then
        echo "[WARNING] too many arguments! Specify SNV and INDEL or SNVINDEL."
        exit
    fi


    # Load conda env:
    #conda env create -n eucancan -f golden-datasets/scripts/environment_snv.yml


    conda activate eucancan
    #source activate eucancan
    # If SNV and INDEL in two files

    if [[ ! -z "$snv" && ! -z "$indel" ]]; then
        echo "[INFO] snv and indel not empty"
        if [[ $snv == *.vcf ]]; then
            bgzip -@ $CPU -c $snv > $snv".gz"
            snv=$snv".gz"
        fi
        if [[ $indel == *.vcf ]]; then
            bgzip -@ $CPU -c $indel > $indel".gz"

            indel=$indel".gz"
        fi

        echo -e "[Running Information]: indexing vcf files\n"
        bcftools index --threads $CPU -f -o $snv".csi" $snv
        bcftools index --threads $CPU -f -o $indel".csi" $indel

        echo -e "[Running Information]: concatenate vcf files\n"
        bcftools concat -a $snv $indel -O z -o $OUTPUT_DIR/"snv_indel_temp.vcf.gz" --threads $CPU
        snvindel=$OUTPUT_DIR/"snv_indel_temp.vcf.gz"
    fi


    # If SNV_INDEL:

    ## if snvindel vcf
    if [[ $snvindel == *.vcf ]]; then
        cat $snvindel > $OUTPUT_DIR/snv_indel_temp.vcf
        gzip $OUTPUT_DIR/snv_indel_temp.vcf
        snvindel=$OUTPUT_DIR/snv_indel_temp.vcf.gz
    fi

    ## if truth vcf
    if [[ $truth == *.vcf ]]; then
        cat $truth > $OUTPUT_DIR/truth_temp.vcf
        gzip $OUTPUT_DIR/truth_temp.vcf
        truth=$OUTPUT_DIR/truth_temp.vcf.gz
    fi

    # Check if vcf is mono sample:

    echo -e "[Running Information]: checking if vcf is multisample\n"

    if [[ `bcftools query -l $snvindel | wc -l` -gt 1  && ! -z "$SAMPLE_NAME" ]]; then
        echo $SAMPLE_NAME
        bcftools view -c1 -O z -s $SAMPLE_NAME -o $OUTPUT_DIR/snv_indel_temp.sample.vcf.gz $snvindel --threads $CPU
        snvindel=$OUTPUT_DIR/snv_indel_temp.sample.vcf.gz
    else [[ `bcftools query -l $snvindel | wc -l` -gt 1 ]];
        echo "[ERROR]" $snvindel "is a multisample"
        echo "[ERROR] sample name must be specified in -n parameter"
        exit
    fi

    if [[ `bcftools query -l $truth | wc -l` -gt 1 ]]; then
        echo "[ERROR]" $truth "is a multisample"
        echo "[ERROR] sample name must be specified in -b parameter"
        exit
    elif [[ `bcftools query -l $truth | wc -l` -gt 1  && ! -z "$SNV_SAMPLE_NAME" ]]; then
        echo $SNV_SAMPLE_NAME
        bcftools view -c1 -O z -s $SNV_SAMPLE_NAME -o $OUTPUT_DIR/truth_temp.sample.vcf.gz $truth --threads $CPU
        truth=$OUTPUT_DIR/truth_temp.sample.vcf.gz
    fi

    ## Replace the 'chr' with '' in the VCFs

    echo -e "[Running Information]: replacing "" by "chr"\n"

    zcat $snvindel | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | awk '{gsub(/contig=\<ID=/,"contig=<ID=chr"); print}' | awk '{gsub(/chrchr/,"chr"); print}' > $OUTPUT_DIR/snv_indel_temp.vcf
    zcat $truth | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' | awk '{gsub(/contig=\<ID=/,"contig=<ID=chr"); print}' | awk '{gsub(/chrchr/,"chr"); print}' > $OUTPUT_DIR/truth_temp.vcf

    snvindel=$OUTPUT_DIR/"snv_indel_temp.vcf"
    truth=$OUTPUT_DIR/"truth_temp.vcf"

    ## Filtering test file PASS variants:
    if [[ "$NOPASS" = false && ! -z "$snvindel" ]]; then
        echo -e "[Running Information]: keeping only PASS variants in snv test file\n"

        grep "PASS\|#" $snvindel > $OUTPUT_DIR/"snv_indel.pass.vcf"

        snvindel=$OUTPUT_DIR/"snv_indel.pass.vcf"
    fi

    ## Filtering truth file PASS variants:
    if [[ "$NOPASS_TRUTH" = false && ! -z "$truth" ]]; then
        echo -e "[Running Information]: keeping only PASS variants in snv truth file\n"

        grep "PASS\|#" $truth > $OUTPUT_DIR/"truth_temp.pass.vcf"

        truth=$OUTPUT_DIR/"truth_temp.pass.vcf"
    fi

    if [[ "$NOPASS_TRUTH" = false && ! -z "$truth_sv" && $truth_sv == *.vcf ]]; then
        echo -e "[Running Information]: keeping only PASS variants in sv truth file\n"

        grep "PASS\|#" $truth_sv> $OUTPUT_DIR/"truth_temp_sv.pass.vcf"
        truth_sv=$OUTPUT_DIR/"truth_temp_sv.pass.vcf"
    elif [[ "$NOPASS_TRUTH" = false && ! -z "$truth_sv" && $truth_sv == *.vcf.gz ]];then
        echo -e "[Running Information]: keeping only PASS variants in sv truth file\n"

        zcat $truth_sv | grep "PASS\|#"> $OUTPUT_DIR/"truth_temp_sv.pass.vcf"
        truth_sv=$OUTPUT_DIR/"truth_temp_sv.pass.vcf"
    fi

    # Sorting vcf files
    echo -e "[Running Information]: sorting vcf files\n"

    bcftools sort $snvindel -o $OUTPUT_DIR/"snv_indel.pass.sort.vcf.gz" -O z
    bcftools sort $truth -o $OUTPUT_DIR/"truth_temp.pass.sort.vcf.gz" -O z

    snvindel=$OUTPUT_DIR/"snv_indel.pass.sort.vcf.gz"
    truth=$OUTPUT_DIR/"truth_temp.pass.sort.vcf.gz"

    # indexing sorted vcf files
    echo -e "[Running Information]: indexing vcf files\n"

    bcftools index -f -o $snvindel".csi" $snvindel --threads $CPU
    bcftools index -f -o $truth".csi" $truth --threads $CPU

    # Preparing for normalization
    echo -e "[Running Information]: preparing for normalizing\n"
    echo -e "[Running Information]: multimerge...]"

    export HGREF=$FASTA

    multimerge $snvindel -r $HGREF -o $OUTPUT_DIR/"snv_indel.pass.sort.prep.vcf.gz" --calls-only=1
    multimerge $truth -r $HGREF -o $OUTPUT_DIR/truth_temp.sort.prep.vcf.gz --calls-only=1

    snvindel=$OUTPUT_DIR/"snv_indel.pass.sort.prep.vcf.gz"
    truth=$OUTPUT_DIR/"truth_temp.sort.prep.vcf.gz"

    # Normalizing vcfs:
    echo -e "[Running Information]: normalizing test file\n"

    pre.py $snvindel $OUTPUT_DIR/"snv_indel.pass.sort.prep.norm.vcf.gz" -L --decompose --somatic --threads $CPU

    echo -e "\n[Running Information]: normalizing truth file\n"

    pre.py $truth $OUTPUT_DIR/truth_temp.sort.prep.norm.vcf.gz -L --decompose --somatic --threads $CPU

    snvindel=$OUTPUT_DIR/"snv_indel.pass.sort.prep.norm.vcf.gz"
    truth=$OUTPUT_DIR/"truth_temp.sort.prep.norm.vcf.gz"

    # Running SNV ingestion script:
    echo -e "[Running Information]: Running ingestion_snv.py script \n"
    # Define directory of git repo with scripts
    DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

    python $DIR/ingest_snv.py -samplename "SAMPLE" -o $OUTPUT_DIR/"snv_indel.pass.sort.prep.norm" $snvindel

    python $DIR/ingest_snv.py -samplename "SAMPLE" -o $OUTPUT_DIR/"truth_temp.sort.prep.norm" $truth

    # Check if we can correctly ingest the SV Truth file
    # No separate sample name for SV file given: this means that it is single sample OR we should reuse samplename OR user should have provided -a argument.
    if [[ `bcftools query -l $sv | wc -l` -gt 1  && -z "$SV_SAMPLE_NAME" ]]; then
      echo "[INFO]" $snvindel "is a multisample"
      echo "[INFO] We will try to use the samplename given by the -n toggle. If this does not work, please provide the correct samplename for the SV file by using the -a toggle."
      SV_SAMPLE_NAME=$SAMPLE_NAME
    fi

    snvindel=$OUTPUT_DIR/"snv_indel.pass.sort.prep.norm.filtered.vcf"
    truth=$OUTPUT_DIR/"truth_temp.sort.prep.norm.filtered.vcf"

    # Running som.py:
    #conda activate eucancan_sv
    #source activate eucancan
    echo -e "[Running Information]: Running som.py evaluation script\n"
    if [[ "$NOPASS" = true ]]; then
        if [[ ! -z "$target" ]];then
            echo "NOPASS + Target"
            som.py $truth $snvindel -o $OUTPUT_DIR/$OUT_NAME --verbose -N -P -R $target
        else
            echo "NOPASS"
            som.py $truth $snvindel -o $OUTPUT_DIR/$OUT_NAME --verbose -N -P
        fi
    else
        if [[ ! -z "$target" ]];then
            echo "Target"
            som.py $truth $snvindel -o $OUTPUT_DIR/$OUT_NAME --verbose -N -R $target
        else
            echo "Basic"
            som.py $truth $snvindel -o $OUTPUT_DIR/$OUT_NAME --verbose -N
        fi
    fi

    echo -e "[Running Information]: SNV becnhmarking ended successfully\n"
    conda deactivate

fi

if [[ -n "$sv" && "$truth_sv" ]];then

    # Running SV part script:

    #conda env create -n eucancan_sv -f ~/golden-datasets/scripts/environment_sv.yml

    conda activate eucancan_sv
    #source activate eucancan_sv

    echo -e "[Running Information]: Running SV ingest.py script \n"
    sv_dataframe=$OUTPUT_DIR/"sv_dataframe.csv"
    truth_sv_dataframe=$OUTPUT_DIR/"truth_sv_dataframe.csv"

    ## Filtering PASS variants:
    if [[ "$NOPASS" = false ]]; then
        echo -e "[Running Information]: keeping only SV PASS variants\n"
        python $DIR/ingest.py $sv -samplename $SAMPLE_NAME -outputfile $sv_dataframe -filter
        if [[ -z "$SV_SAMPLE_NAME" ]]; then
          python $DIR/ingest.py $truth_sv -outputfile $truth_sv_dataframe -filter
        else
          python $DIR/ingest.py $truth_sv -samplename $SV_SAMPLE_NAME -outputfile $truth_sv_dataframe -filter
        fi
    else
        python $DIR/ingest.py $sv -samplename $SAMPLE_NAME -outputfile $sv_dataframe
    fi

    echo -e "[Running Information]: Running compare_node_to_truth.py script\n"
    metrics=$OUTPUT_DIR/"SV_benchmark_results.csv"

    python $DIR/compare_node_to_truth.py $sv_dataframe $truth_sv_dataframe -metrics $metrics

    conda deactivate
    echo -e "[Running Information]: SV becnhmarking ended successfully\n"
fi

# Ploting with R
## SNVs

conda activate benchmarking_plots
#source activate benchmarking_plots

if [[ -n "$truth" && ( -n "$snv" && -n "$indel" || -n "$snvindel") ]];then
    echo -e "[Running Information]: Running plots_benchmarking_snv.R script\n"
    Rscript $DIR/plots_benchmarking_snv.R -b $OUTPUT_DIR/$OUT_NAME".stats.csv" -o $OUTPUT_DIR/$OUT_NAME
fi

if [[ -n "$sv" && "$truth_sv" ]];then
    echo -e "[Running Information]: Running plots_benchmarking_SV.R script\n"
    Rscript $DIR/plots_benchmarking_SV.R -b $OUTPUT_DIR/"SV_benchmark_results.csv" -t $truth_sv -o $OUTPUT_DIR/$OUT_NAME
fi

conda deactivate

# Cleaning
echo -e "[Running Information]: Cleaning files\n"

rm $OUTPUT_DIR/*{tbi,csi,json}

if [[ "$KEEP" = false ]];then
    rm $OUTPUT_DIR/*{vcf,vcf.gz}
fi

echo -e "[Running Information]: The end !\n"

# Example commands
# SCRIPT_DIR=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/golden-datasets-fork-tom/golden-datasets/scripts
# OUT_DIR=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/test_ingestions_sh
# TRUTH=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/truth_dream/truth.snvs.synthetic.challenge.set1.chr.vcf.gz
# SNV=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/colo829/curie_colo829_snps.sample.vcf.gz
# INDEL=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/colo829/curie_colo829_indels.sample.vcf.gz
# SNV_INDEL=/bioinfo/users/tgutman/Documents/Tom/EUCANCan/Benchmark/colo829/head_filtered_Mutect2_ERR2752450_vs_ERR2752449_snpEff.ann.vcf.gz
#
# ## with mixed indels & snvs
# bash $SCRIPT_DIR/temp_ingest_snv.sh --truth $TRUTH --snvindel $SNV_INDEL --outdir $OUT_DIR --outname test_multisample --sname COLO829T  --keep -f /data/annotations/pipelines/Human/hg19/genome/hg19.fa
#
# # with split snvs & indels
# bash $SCRIPT_DIR/temp_ingest_snv.sh --truth $TRUTH --snv $SNV --indel $INDEL --outdir $OUT_DIR --outname test_multisample --keep -f /data/annotations/pipelines/Human/hg19/genome/hg19.fa

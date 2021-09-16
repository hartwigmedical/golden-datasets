# Scripts

Scripts in this folder have two main purpose:
- parsing and benchmarking snv
- parsing and benchmarking sv

## List of scripts:
### SV related:
- `compare_node_to_truth.py`
- `ingest.py`
- `parse_sv_tsv.py`
- `parse_sv_vcf.py`
- `reclassification.py`

### SNV related:

- `ingest_snv.py`
- `parse_snv.py`
- `ingest_snv.sh`

## Tool summary

<center>
<img src="../docs/flowchart.png" alt="flowchart" width="700px"/>
</center>


### SV Benchmark

First of all, clone the devel branch of the repo:

```shell
git clone -b devel https://github.com/EUCANCan/golden-datasets.git
```

**Requirements**

This tool is implemented in `python3` and require several libraries such as `pandas`,`pyVCF` and `argsParse`.
A `conda` environment is provided. Therefore you need to have `conda` installed in your machine.
To use it, simply go to the `scripts` folder of the repo and execute:

```shell
conda env create --name ingestion python=3.6 -f environment_sv.yml
conda activate ingestion
```

Add any new requirement for your code to [environment.yml](https://github.com/EUCANCan/golden-datasets/blob/devel/scripts/environment.yml) (if it is compatible)

#### Usage

**1) Ingestion**


Creates dataframes from VCF/TSV files. This step is required both for test files and truth files.
Flag `-samplename` is required for VCF files and must contain the Tumor sample name

```shell
python ingest.py TEST_FILE -outputfile DATAFRAME_TEST_FILE
python ingest.py TRUTH_FILE -outputfile DATAFRAME_TRUTH_FILE
```


**2) Benchmarking**

Computes metrics for SV calls. Optional: `-metrics` flag to print result to a file

```shell
python compare_node_to_truth.py DATAFRAME_TEST_FILE DATAFRAME_TRUTH_FILE -metrics OUTPUT_METRICS_FILE
```


### SNV Benchmark

**requirements**

This tool is implemented in `bash` and `python3` and require several associated tools and libraries such as `bcftools`,`hap.py` and `pyVCF`.
A `conda` environment is provided.
To use it, simply use:

```shell
conda create -n snv_bench -f golden-datasets/scripts/environment_snv.yml
```
or
```
conda env create -n snv_bench -f golden-datasets/scripts/environment_snv.yml
```
And then:
```
conda activate snv_bench
```

**implementation**

To benchmark snv calls from different centers, multipe steps are required:

- perparing and normalizing test and truth vcfs with `bcftools`, `multimerge` and `pre.py`
- extracting indels > 50bp so they can be assessed as SV with `ingest_snv.py`
- benchmarking each snv file to the truth file with `som.py`

**Usage**


```
bash wrapper.sh  -h

Usage:
bash wrapper.sh -t, --truth truth_file.vcf
                   -s, --snv snv.vcf
                   -i, --indel indel.vcf
                   -m, --snvindel indels_and_vcfs.vcf
                   -v, --sv sv.vcf
                   -u, --truth_sv truth_sv.vcf
                   -f, --fasta ref_fasta.fa
                   -d, --outdir /OUTPUT_DIR/PATH
                   -o, --outname output file name
                   -n, --sname vcf sample name
                   -p, --pass keep only the pass variants
                   -c, --cpu number of threads
                   -k, --keep (to keep intermediates files)
```
Note: Truth file and snv/indel file should have different names

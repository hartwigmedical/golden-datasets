Available companion scripts for pipeline nextflow:

**SNV related:**
- ingest_snv.py : ingestion script to remove indels > 50bp (not used anymore)
- parse_snv.py : dependency of ingest_snv.py

**SV related:**
- compare_node_to_truth.py : comparison between SV test table & truth table
- ingest.py :  ingestion script to create table from test and truth vcf or tsv
- parse_sv_tsv.py : dependency of ingest.py
- parse_sv_vcf.py : dependency of ingest.py

**Visualization related:**

- plots_benchmarking_snv.R : create plots from SNV benchmarking results
- plots_benchmarking_SV.R : create plots from SV benchmarking results

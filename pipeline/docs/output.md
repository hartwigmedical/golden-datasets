# Outputs

Assume that the pipeline has been installed in `${HOME}/tmp/myPipeline/install`. When you run the pipeline, it will create inside the `${HOME}/tmp/myPipeline/install/results` folder  which contains the following subfolders:

* SNV
* SV

## SNV folder:

In this folder you will find results from the SNV benchmark

- `${outName}barplotPrecisionRecallF1.png` and `${outName}barplotTPFPFN.png` are the visualization of the results
- `${outName}snvTable.csv` : is the table summary of the benchmark
- `${outName}.stats.csv` : is the full table summary of the benchmark

## SV folder:

In this folder you will find results from the SV benchmark

- `${outName}barplot_SV_F1_tier3_all.png` `${outName}barplot_SV_general_tier3_all.png` `${outName}barplot_SV_TP_tier3_all.png` `${outName}barplot_SV_truth_types.png` : figures of the results of SV benchmarking
- `SV_benchmark_results.csv` : table with all the SV results with the different tiers
- `sv_dataframe.csv`: standardized table created by ingest.py and used by compare_node_to_truth.py
- `${outName}svTable.csv` : table of results created by the visualization script
- `truth_sv_dataframe.csv` : standardized table created by ingest.py and used by compare_node_to_truth.py

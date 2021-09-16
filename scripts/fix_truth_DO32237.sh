# Solving issue with vcf from ICGC /TCGA
## Goal create a valide vcf with a good header.

# Go to destination dir
cd /home/daphne/data/daphne-eucancan-benchmark-data/ICGC_results/truth/
# Unzip files

gunzip -c 5e68f783-790f-4aed-b7c5-56c34538fb09.pilot50.20160609.somatic.sv.vcf.gz > 5e68f783-790f-4aed-b7c5-56c34538fb09.pilot50.20160609.somatic.sv.vcf

# Add format field:
awk '{if($0 !~ /^#/) print $0"\tGT:AD\t0/1:19,0"; else print $0}' 5e68f783-790f-4aed-b7c5-56c34538fb09.pilot50.20160609.somatic.sv.vcf > temp_sv.vcf

#bien vérifier dans atom si les espaces et tabulations sont bien respectées !

# reaheader
/opt/tools/bcftools/1.9/bcftools reheader -h /home/daphne/golden-datasets/scripts/header_sv.txt -o sv_fixed.vcf -s PILOT50 temp_sv.vcf

# bgzip
bgzip -c sv_fixed.vcf > sv_fixed.vcf.gz

# index:
bcftools index sv_fixed.vcf.gz
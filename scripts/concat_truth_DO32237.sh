# Solving issue with vcf from ICGC /TCGA
## Goal create a valide vcf with a good header.

# Unzip files

gunzip -c /data/tmp/tgutman/VCF_archive/original_files/DO32237T_vs_DO32237N_Mutect2_filtered_pass_norm.vcf.gz > DO32237T_vs_DO32237N_Mutect2_filtered_pass_norm.vcf
gunzip -c /data/tmp/tgutman/VCF_archive/original_files/5e68f783-790f-4aed-b7c5-56c34538fb09.pilot50.20160609.somatic.indel.vcf.gz >5e68f783-790f-4aed-b7c5-56c34538fb09.pilot50.20160609.somatic.indel.vcf
gunzip -c /data/tmp/tgutman/VCF_archive/original_files/5e68f783-790f-4aed-b7c5-56c34538fb09.pilot50.20160609.somatic.snv_mnv.vcf.gz > 5e68f783-790f-4aed-b7c5-56c34538fb09.pilot50.20160609.somatic.snv_mnv.vcf

# Add format field:
awk '{if($0 !~ /^#/) print $0"\tGT:AD\t0/1:19,0"; else print $0}' 5e68f783-790f-4aed-b7c5-56c34538fb09.pilot50.20160609.somatic.indel.vcf > temp_indel.vcf

awk '{if($0 !~ /^#/) print $0"\tGT:AD\t0/1:19,0"; else print $0}' 5e68f783-790f-4aed-b7c5-56c34538fb09.pilot50.20160609.somatic.snv_mnv.vcf > temp_snv.vcf

#bien vérifier dans atom si les espaces et tabulations sont bien respectées !

# reaheader
bcftools reheader -h /data/tmp/tgutman/VCF_archive/original_files/header.txt -o snv.vcf -s PILOT50 temp_snv.vcf

bcftools reheader -h /data/tmp/tgutman/VCF_archive/original_files/header.txt -o indel.vcf -s PILOT50 temp_indel.vcf

# bgzip

bgzip -c snv.vcf > snv.vcf.gz
bgzip -c indel.vcf > indel.vcf.gz

# index:
bcftools index snv.vcf.gz
bcftools index indel.vcf.gz

# concat both files
bcftools concat -a snv.vcf.gz indel.vcf.gz > test.vcf.gz

bcftools sort -o pilot50.snv_indel.vcf.gz -Oz test.vcf.gz

# For the SV part:

gunzip 5e68f783-790f-4aed-b7c5-56c34538fb09.pilot50.20160609.somatic.sv.vcf.gz

# Add format field:
awk '{if($0 !~ /^#/) print $0"\tGT:AD\t0/1:19,0"; else print $0}' 5e68f783-790f-4aed-b7c5-56c34538fb09.pilot50.20160609.somatic.sv.vcf > temp_sv.vcf

# reaheader
## use the header_sv.txt
bcftools reheader -h /bioinfo/users/tgutman/Documents/Tom/EUCANCan/golden-datasets/scripts/header_sv.txt -o sv.vcf -s PILOT50 temp_sv.vcf

bgzip -c sv.vcf> sv.vcf.gz

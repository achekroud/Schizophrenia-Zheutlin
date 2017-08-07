
# generate scores in cnp, swedish twins for each cog var
# 6 cog vars - residual performance (no age, sex)

cogvars=( CVLT CRT DS VOC VR1 VR2 )

for i in "${cogvars[@]}"
do
plink --bfile 77SNPs-cnp --score 77SNP-$i.txt --out 77SNP-cnp-$i
plink --bfile /data/swe_gwas/ABZ/imputed_data/OMNI/call_genos/omni_genos_merged --score 77SNP-$i.txt --out 77SNP-swe-$i
done

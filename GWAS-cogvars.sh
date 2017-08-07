# GWAS to determine beta-weights for each cog var
# 6 cog vars - residual performance (no age, sex)

cogvars=( CVLT_TOTCOR CRT_TIME1 DS_TOTALRAW VOC_TOTALRAW VR1IR_TOTALRAW VR2DR_TOTALRAW )

for i in "${cogvars[@]}"
do
plink --bfile 77SNPs-cnp --linear --pheno $i-resid-cnp.txt --out 77SNPassoc-$i
done

for chr in $(seq 1 33); do
    plink --bfile ./hwe_filtered/snp50_c --chr $chr --make-bed  --nonfounders --cow --out ./chrsplit/snp50_c.chr$chr ;
done

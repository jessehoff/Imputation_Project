for chr in $(seq 1 33); do
    plink --bfile ./merged_files/merged --chr $chr --make-bed  --nonfounders --cow --out ./chrsplit/merged.chr$chr ;
done

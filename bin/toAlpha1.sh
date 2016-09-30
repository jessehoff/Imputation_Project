for chr in $(seq 1 32); do
	plink --bfile ./chrsplit/snp50_c.chr$chr  --recode 01 --nonfounders --cow --out ./toAlpha/peds/snp50_c.chr$chr --output-missing-genotype 3;
done

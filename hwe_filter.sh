#!/bin/bash

# This command removes alleles with HWE P value under specified number

plink --bfile ./individual_filtered/snp50_a --cow --nonfounders --hwe 0.01 --make-bed --out ./hwe_filtered/snp50_a

plink --bfile ./individual_filtered/snp50_b --cow --nonfounders --hwe 0.01 --make-bed --out ./hwe_filtered/snp50_b

plink --bfile ./individual_filtered/snp50_c --cow --nonfounders --hwe 0.01 --make-bed --out ./hwe_filtered/snp50_c

plink --bfile ./individual_filtered/ggpf250 --cow --nonfounders --hwe 0.01 --make-bed --out ./hwe_filtered/ggpf250

plink --bfile ./individual_filtered/hd --cow --nonfounders --hwe 0.01 --make-bed --out ./hwe_filtered/hd


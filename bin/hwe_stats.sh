#!/bin/bash

plink --bfile ./individual_filtered/snp50_a --cow --nonfounders --hardy --out ./hwe_stats/snp50_a

plink --bfile ./individual_filtered/snp50_b --cow --nonfounders --hardy --out ./hwe_stats/snp50_b

plink --bfile ./individual_filtered/snp50_c --cow --nonfounders --hardy --out ./hwe_stats/snp50_c

plink --bfile ./individual_filtered/ggpf250 --cow --nonfounders --hardy --out ./hwe_stats/ggpf250

plink --bfile ./individual_filtered/hd --cow --nonfounders --hardy --out ./hwe_stats/hd


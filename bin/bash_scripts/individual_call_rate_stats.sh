#!/bin/bash

plink --bfile ./allele_filtered/snp50_a --cow --missing --out ./individual_stats/snp50_a

plink --bfile ./allele_filtered/snp50_b --cow --missing --out ./individual_stats/snp50_b

plink --bfile ./allele_filtered/snp50_c --cow --missing --out ./individual_stats/snp50_c

plink --bfile ./allele_filtered/ggpf250 --cow --missing --out ./individual_stats/ggpf250

plink --bfile ./allele_filtered/hd --cow --missing --out ./individual_stats/hd


#!/bin/bash

plink --bfile ./allele_filtered/snp50_a --cow --mind 0.05 --make-bed --out ./individual_filtered/snp50_a

plink --bfile ./allele_filtered/snp50_b --cow --mind 0.05 --make-bed --out ./individual_filtered/snp50_b

plink --bfile ./allele_filtered/snp50_c --cow --mind 0.05 --make-bed --out ./individual_filtered/snp50_c

plink --bfile ./allele_filtered/ggpf250 --cow --mind 0.05 --make-bed --out ./individual_filtered/ggpf250

plink --bfile ./allele_filtered/hd --cow --mind 0.05 --make-bed --out ./individual_filtered/hd


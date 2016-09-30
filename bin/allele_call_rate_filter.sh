#!/bin/bash

plink --file ./raw_genotypes/58336.160906.100.test_snp50_A --map ./maps/9913_SNP50.map --cow --nonfounders --not-chr 0 --geno 0.05 --make-bed --out ./allele_filtered/snp50_a

plink --file ./raw_genotypes/58336.160906.100.test_snp50_B --map ./maps/9913_SNP50.map --cow --nonfounders --not-chr 0 --geno 0.05 --make-bed --out ./allele_filtered/snp50_b

plink --file ./raw_genotypes/58336.160906.100.test_snp50_C --map ./maps/9913_SNP50.map --cow --nonfounders --not-chr 0 --geno 0.05 --make-bed --out ./allele_filtered/snp50_c

plink --file ./raw_genotypes/227234.160906.100.test_ggpf250_A --map ./maps/9913_GGPF250.map --cow --nonfounders --not-chr 0 --geno 0.05 --make-bed --out ./allele_filtered/ggpf250

plink --file ./raw_genotypes/777962.160906.100.test_hd_A --map ./maps/9913_HD.map --cow --nonfounders --not-chr 0 --geno 0.05 --make-bed --out ./allele_filtered/hd

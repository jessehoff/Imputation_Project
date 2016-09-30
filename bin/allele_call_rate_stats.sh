#!/bin/bash

plink --file ./raw_genotypes/58336.160906.100.test_snp50_A --map ./maps/9913_SNP50.map --cow --nonfounders --freq --out ./allele_stats/snp50_a

plink --file ./raw_genotypes/58336.160906.100.test_snp50_B --map ./maps/9913_SNP50.map --cow --nonfounders --freq --out ./allele_stats/snp50_b

plink --file ./raw_genotypes/58336.160906.100.test_snp50_C --map ./maps/9913_SNP50.map --cow --nonfounders --freq --out ./allele_stats/snp50_c

plink --file ./raw_genotypes/227234.160906.100.test_ggpf250_A --map ./maps/9913_GGPF250.map --cow --nonfounders --freq --out ./allele_stats/ggpf250

plink --file ./raw_genotypes/777962.160906.100.test_hd_A --map ./maps/9913_HD.map --cow --nonfounders --freq --out ./allele_stats/hd


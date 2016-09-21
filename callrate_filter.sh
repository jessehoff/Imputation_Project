#!/bin/bash

plink --file 58336.160906.100.test_snp50_A --map 9913_SNP50.map --cow --geno 0.05 --mind 0.05 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/allele_individual_filtered/SNP_50_A

plink --file 58336.160906.100.test_snp50_B --map 9913_SNP50.map --cow --geno 0.06 --mind 0.06 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/allele_individual_filtered/SNP_50_B

plink --file 58336.160906.100.test_snp50_C --map 9913_SNP50.map --cow --geno 0.06 --mind 0.06 --make-bed  --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/allele_individual_filtered/SNP_50_C

plink --file 227234.160906.100.test_ggpf250_A --map 9913_GGPF250.map --cow --geno 0.06 --mind 0.15 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/allele_individual_filtered/GGPF250

plink --file 777962.160906.100.test_hd_A --map 9913_HD.map --cow --geno 0.06 --mind 0.05 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/allele_individual_filtered/HD


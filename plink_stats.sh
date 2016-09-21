#!/bin/bash

plink --file 58336.160906.100.test_snp50_A --map 9913_SNP50.map --cow --nonfounders --freq --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/SNP_50_A

plink --file 58336.160906.100.test_snp50_B --map 9913_SNP50.map --cow --nonfounders --freq --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/SNP_50_B

plink --file 58336.160906.100.test_snp50_C --map 9913_SNP50.map --cow --nonfounders --freq --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/SNP_50_C

plink --file 227234.160906.100.test_ggpf250_A --map 9913_GGPF250.map --cow --nonfounders --freq --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/GGPF250

plink --file 777962.160906.100.test_hd_A --map 9913_HD.map --cow --nonfounders --freq --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/HD



plink --file 58336.160906.100.test_snp50_A --map 9913_SNP50.map --cow --missing --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Individual_Call_Rates/SNP_50_A_missing  

plink --file 58336.160906.100.test_snp50_B --map 9913_SNP50.map --cow --missing --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Individual_Call_Rates/SNP_50_B_missing  

plink --file 58336.160906.100.test_snp50_C --map 9913_SNP50.map --cow --missing --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Individual_Call_Rates/SNP_50_C_missing 

plink --file 227234.160906.100.test_ggpf250_A --map 9913_GGPF250.map --cow --missing --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Individual_Call_Rates/GGPF250_missing

plink --file 777962.160906.100.test_hd_A --map 9913_HD.map --cow --missing --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Individual_Call_Rates/HD_missing

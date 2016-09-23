#!/bin/bash
#this filters multiple files first for missing genotypes, and then uses resulting genotypes to filter out individuals with low call rates. 
plink --file 58336.160906.100.test_snp50_A --map 9913_SNP50.map --cow --geno 0.05 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/SNP_50_A_allele

plink --file 58336.160906.100.test_snp50_B --map 9913_SNP50.map --cow --geno 0.06 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/SNP_50_B_allele

plink --file 58336.160906.100.test_snp50_C --map 9913_SNP50.map --cow --geno 0.06 --make-bed  --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/SNP_50_C_allele

plink --file 227234.160906.100.test_ggpf250_A --map 9913_GGPF250.map --cow --geno 0.06 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/GGPF250_allele

plink --file 777962.160906.100.test_hd_A --map 9913_HD.map --cow --geno 0.06 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/HD_allele



plink --bfile /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/SNP_50_A_allele --cow --mind 0.05 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Individual_Call_Rates/SNP_50_A_allele_individual

plink --bfile /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/SNP_50_B_allele --cow --mind 0.06 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Individual_Call_Rates/SNP_50_B_allele_individual

plink --bfile /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/SNP_50_C_allele --cow --mind 0.06 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Individual_Call_Rates/SNP_50_C_allele_individual

plink --bfile /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/GGPF250_allele --cow --mind 0.85 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Individual_Call_Rates/GGPF250_allele_individual

plink --bfile /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Allele_Frequencies/HD_allele --cow --mind 0.05 --make-bed --out /CIFS/MUG01_N/taylorjerr/JLH/160906_imputation_test/Individual_Call_Rates/HD_allele_individual


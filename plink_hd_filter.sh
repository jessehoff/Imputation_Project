#!/bin/bash 

plink --file 777962.160906.100.test_hd_A --map 9913_HD.map --cow --geno 0.06 --mind 0.05 --make-bed --out HD_cleaned


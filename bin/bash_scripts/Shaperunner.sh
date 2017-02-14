shapeit --input-bed ./chrsplit/merged.chr1.bed ./chrsplit/merged.chr1.bim ./chrsplit/merged.chr1.fam \
        --duohmm \
	--output-max ./shapetest/merged.phased.haps ./shapetest2/merged.phased.sample

        #--input-map ./maps/9913_SNP50.map no map needed, automatically produced with default recombination rates. 

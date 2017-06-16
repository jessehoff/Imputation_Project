DATA =['f250', 'ggpld', 'hd', 'snp50'] #new file names -- these are files that have had ref-alt conversions

rule targ:
	input:
		hd_or_f250  = {'snp50':"correct_sex/777962.170519.1970.HD",'f250':"correct_sex/227234.170519.1970.GGPF250", 'ggpld':"correct_sex/777962.170519.1970.HD", 'hd':"correct_sex/777962.170519.1970.HD"}
def bedchoice(WC):
	return hd_or_f250[WC.assay]

pull_or_not = {'snp50':"--extract dataprepper/snp50_snps.txt",'f250':"", 'ggpld':"--extract dataprepper/ggpld_snps.txt", 'hd':""}
def snpset(WC):
	return pull_or_not[WC.assay]

def runchoice(WC):
	run_dict = {'1':'merged_chrsplit','6':'merged_chrsplit', '2':'assay_chrsplit/','7':'assay_chrsplit/'}
	r = WC.run
	chrom = WC.chr
	if r == '1':
		location = run_dict[r] + '/run' + r+'/hol_testset.merge.chr' + chrom +'.bed'
	if r =='2':
		location = run_dict[r] + WC.sample + '.list1.chr' + chrom + '.bed'
	if r =='6':
		location = run_dict[r] + '/run' + r+'/hol_testset.merge.chr' + chrom +'.bed'
	if r =='7':
		location = run_dict[r] + WC.sample + '.list2.chr' + chrom + '.bed'
	return location
def listchoice(WC):
	r = WC.run
	chrom = WC.chr
	if r =='2':
		location =  'assay_chrsplit/' + WC.sample +'.list1.chr' + chrom
	if r =='7':
		location =  'assay_chrsplit/'+ WC.sample + '.list2.chr' + chrom
	if r =='12':
		location =  'assay_chrsplit/' + WC.sample +'.list3.chr' + chrom
	return location


rule downsample:
	params:
		idslist = "--keep dataprepper/{assay}_ids.list{list}.txt",
		extract = snpset,
		bfile = bedchoice,
		oprefix = "downsample/{assay}.list{list}"
	benchmark:
		"benchmarks/downsample/{assay}.list{list}.txt"
	log:
		"logs/downsample/{assay}.list{list}.log"
	output:
		bed = "downsample/{assay}.list{list}.bed"
	shell:
		"(plink --bfile {params.bfile}  --real-ref-alleles {params.idslist} {params.extract}  --make-bed  --cow --out {params.oprefix})> {log}"

rule assay_chrsplit:
	input:
		bed = expand("downsample/{assay}.list{{list}}.bed", assay = DATA),
	params:
		inprefix = "downsample/{sample}.list{list}",
		oprefix = "assay_chrsplit/{sample}.list{list}.chr{chr}",
		chr = "{chr}"
	benchmark:
		"benchmarks/assay_chrsplit/{sample}.chr{chr}.txt"
	log:
		"logs/eagle_split_chromosomes/{sample}.chr{chr}.log"
	output:
		bed = "assay_chrsplit/{sample}.list{list}.chr{chr}.bed",
		bim = "assay_chrsplit/{sample}.list{list}.chr{chr}.bim",
		fam = "assay_chrsplit/{sample}.list{list}.chr{chr}.fam",
		log = "assay_chrsplit/{sample}.list{list}.chr{chr}.log"
	shell:
		"(plink --bfile {params.inprefix}  --real-ref-alleles --chr {params.chr} --make-bed  --nonfounders --cow --out {params.oprefix})> {log}"

rule eagle_phased_assays:
	input:
		bed = runchoice
	params:
		inprefix = listchoice,
		oprefix = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased"
	benchmark:
		"benchmarks/eagle_phased_assays/run{run}/{sample}.chr{chr}.benchmark.txt"
	threads: 8
	log:
		"logs/eagle_phased_assays/run{run}/{sample}.chr{chr}.log"
	output:
		sample = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.sample",
		haps = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.haps.gz"
	shell:
		"(eagle --bfile={params.inprefix} --geneticMapFile=USE_BIM --maxMissingPerSnp 1 --maxMissingPerIndiv 1 --numThreads 8 --outPrefix {params.oprefix}) > {log}"

rule decompress_single_chrom:
		input:
			gzhaps = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.haps.gz"
		log:
			"logs/decompress/run{run}/{sample}.chr{chr}.log"
		benchmark:
			"benchmarks/decompress/run{run}/{sample}.chr{chr}.benchmark.txt"
		output:
			haps = temp("eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.haps")
		shell:
			"(gunzip -c {input.gzhaps} > {output.haps}) > {log}"

rule hap_leg:
	input:
		haps = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.haps",
		sample = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased.sample"
	params:
		inprefix = "eagle_phased_assays/run{run}/{sample}.chr{chr}.phased",
		oprefix = "impute_input/run{run}/{sample}.chr{chr}.phased"
	log:
		"logs/hap_leg/run{run}/{sample}.chr{chr}.phased.log"
	benchmark:
		"benchmarks/hap_leg/run{run}/{sample}.chr{chr}.phased.benchmark.txt"
	output:
		hap = "impute_input/run{run}/{sample}.chr{chr}.phased.haplotypes",
		leg = "impute_input/run{run}/{sample}.chr{chr}.phased.legend",
		log = "impute_input/run{run}/logs/{sample}.chr{chr}.phased.log"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-log {output.log} --output-ref {params.oprefix}) > {log}"

def findmergelist(wc):
	listdict = {'6':'2','1':'1','2':'1','3':'1'}
	run = wc.run
	list = listdict[run]
	locate = 'assay_chrsplit/hol_testset.list' + list+'.chr'+wc.chr+'.txt'
	#print(locate)
	return locate

rule make_merge_list: #makes a merge list of raw genotypes, doesn't need to have a run because its just going to be a list related thing.
	input:
		filelist = expand("assay_chrsplit/{assay}.list{{list}}.chr{{chr}}.bed", assay= DATA )
	log:
		"logs/{assay}.chr{chr}.txt"
	output:
		"assay_chrsplit/hol_testset.list{list}.chr{chr}.txt"
	shell:
		"python ./bin/merge_file_maker.py {input.filelist} {output}"

rule merged_chrsplit:
	input:
		findmergelist
	params:
		outprefix="merged_chrsplit/run{run}/hol_testset.merge.chr{chr}"
	benchmark:
		"benchmarks/merged_chrsplit/run{run}/hol_testset.merg.chr{chr}.txt"
	log:
		"logs/merged_chrsplit/run{run}/hol_testset.merge.chr{chr}.log"
	output:
		bed="merged_chrsplit/run{run}/hol_testset.merge.chr{chr}.bed",
		animalset = "merged_chrsplit/run{run}/hol_testset.merge.chr{chr}.fam"
	shell:
		"(plink --merge-list {input} --nonfounders --real-ref-alleles --cow --make-bed --out {params.outprefix}) > {log}"

rule eagle_merged:
	input:
		bed = runchoice
	params:
		bed="merged_chrsplit/run{run}/hol_testset.merge.chr{chr}",
		out="eagle_merged/run{run}/hol_testset.merge.chr{chr}"
	threads: 16
	benchmark:
		"benchmarks/eagle_merged/run{run}/hol_testset.merge.chr{chr}.benchmark.txt"
	log:
		"logs/eagle_merged/run{run}/hol_testset.merge.chr{chr}.log"
	output:
		sample = "eagle_merged/run{run}/hol_testset.merge.chr{chr}.sample",
		haps = "eagle_merged/run{run}/hol_testset.merge.chr{chr}.haps.gz"
	shell:
		"(eagle --bfile={params.bed}  --geneticMapFile=USE_BIM --maxMissingPerSnp .99  --maxMissingPerIndiv .99 --numThreads 16 --outPrefix {params.out})> {log} "

rule decompress:
		input:
			gzhaps = "eagle_merged/run{run}/{sample}.chr{chr}.haps.gz"
		log:
			"logs/decompress/run{run}/{sample}.chr{chr}.log"
		benchmark:
			"benchmarks/decompress/run{run}/{sample}.chr{chr}.benchmark.txt"
		output:
			haps = temp("eagle_merged/run{run}/{sample}.chr{chr}.haps")
		shell:
			"(gunzip -c {input.gzhaps} > {output.haps}) > {log}"

rule eagle_merged_vcf:
	input:
		haps="eagle_merged/run{run}/hol_testset.merge.chr{chr}.haps"
	log:
			"logs/eagle_merged_vcf/logs/run{run}/hol_testset.merge.{chr}.log"
	benchmark:
			"benchmarks/eagle_merged_vcf/run{run}/hol_testset.merge.{chr}.benchmark.txt"
	params:
		haps="eagle_merged/run{run}/hol_testset.merge.chr{chr}"
	output:
		vcf="eagle_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf",
		log="eagle_merged_vcf/logs/run{run}/hol_testset.merge.chr{chr}.log"
	shell:
		"(shapeit -convert --input-haps {params.haps} --output-log {output.log} --output-vcf {output.vcf}) > {log}"

rule bgzip_vcf:
	input:
		vcf="eagle_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf",
	log:
		"logs/bgzip_vcf/run{run}/hol_testset.merge.chr{chr}"
	benchmark:
		"benchmarks/bgzip_vcf/run{run}/hol_testset.merge.chr{chr}"
	output:
		vcfgz="eagle_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf.gz",
		index="eagle_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf.gz.tbi"
	shell:
		"(bgzip {input.vcf}; tabix {output.vcfgz}) > {log};"

def pickidsforlist(wc):
	listdict = {'6':'2','1':'1','2':'1','3':'1'}
	run = wc.run
	list = listdict[run]
	return list

rule make_vcf_extract_lists: # this rule doesn't include the list, but that is determined by the run
	input:
		animalset = "merged_chrsplit/run{run}/hol_testset.merge.chr25.fam"
	params:
		list = pickidsforlist
	benchmark:
		"benchmarks/make_vcf_extract_lists/run{run}/{assay}.benchmark.txt"
	log:
		"logs/make_vcf_extract_lists/run{run}/{assay}.benchmark.txt"
	output:
		keep_ids = "merged_chrsplit/run{run}/phased_{assay}.keepvcf",
		keep_maps = "merged_chrsplit/run{run}/phased_{assay}.vcfregion"
	shell:
		"python ./bin/vcfextraction_for_joint_phase.py {input.animalset} {params.list}"



rule vcf_per_assay: #filter the vcfs on a per assay basis
	input:
		vcfgz="eagle_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf.gz",
		index="eagle_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf.gz.tbi",
		keep_ids = "merged_chrsplit/run{run}/phased_{assay}.keepvcf",
		keep_maps = "merged_chrsplit/run{run}/phased_{assay}.vcfregion"
	benchmark:
		"benchmarks/vcf_per_assay/run{run}/{assay}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_per_assay/run{run}/{assay}.chr{chr}.log"
	output:
		vcf = "vcf_per_assay/run{run}/{assay}.chr{chr}.vcf",
	shell:
		"(bcftools view {input.vcfgz} -R {input.keep_maps}  -S {input.keep_ids} -o {output.vcf}) > {log}"

rule vcf_to_hap:
	input:
		vcf = "vcf_per_assay/run{run}/{assay}.chr{chr}.vcf",
	params:
		chr = "{chr}",
		oprefix ="vcf_to_hap/run{run}/{assay}.chr{chr}"
	benchmark:
		"benchmarks/vcf_to_hap/run{run}/{assay}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_hap/run{run}/{assay}.chr{chr}.log"
	output:
		legend = "vcf_to_hap/run{run}/{assay}.chr{chr}.phased.legend",
		haps = "vcf_to_hap/run{run}/{assay}.chr{chr}.phased.haplotypes",
		sample = "vcf_to_hap/run{run}/{assay}.chr{chr}.phased.samples"
	shell:
		"(bcftools convert {input.vcf} --haplegendsample {output.haps},{output.legend},{output.sample}) > {log}" #Updated these to use the naming conventions that the shapeit tool outputs the run2 hap,leg,samples, so naming is consistent in Impute2.

rule vcf_to_haps: #doesn't approrpriately name "haps" haps
	input:
		vcf = "vcf_per_assay/run{run}/{assay}.chr{chr}.vcf",
	params:
		chr = "{chr}",
		oprefix ="vcf_to_haps/run{run}/{assay}.chr{chr}.phased"
	benchmark:
		"benchmarks/vcf_to_haps/run{run}/{assay}.chr{chr}.benchmark.txt"
	log:
		"logs/vcf_to_haps/run{run}/{assay}.chr{chr}.log"
	output:
		hap = "vcf_to_haps/run{run}/{assay}.chr{chr}.phased.haps",
		sample = "vcf_to_haps/run{run}/{assay}.chr{chr}.phased.sample"
	shell:
		"(bcftools convert {input.vcf} --hapsample {output.hap},{output.sample} ) > {log}"

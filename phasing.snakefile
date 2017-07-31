DATA =['f250', 'ggpld', 'hd', 'snp50'] #new file names -- these are files that have had ref-alt conversions
hd_or_f250  = {'snp50':"correct_sex/777962.170519.1970.HD",'f250':"correct_sex/227234.170519.1970.GGPF250", 'ggpld':"correct_sex/777962.170519.1970.HD", 'hd':"correct_sex/777962.170519.1970.HD"}

rule targ:
	input:
		#hd_or_f250  = {'snp50':"correct_sex/777962.170519.1970.HD",'f250':"correct_sex/227234.170519.1970.GGPF250", 'ggpld':"correct_sex/777962.170519.1970.HD", 'hd':"correct_sex/777962.170519.1970.HD"}
		phasehd = expand("eagle_phased_assays/run21/hd.chr{chr}.phased.sample", chr = range(1,30))

def bedchoice(WC):
	return hd_or_f250[WC.assay]

pull_or_not = {'snp50':"--extract dataprepper/snp50_snps.txt",'f250':"", 'ggpld':"--extract dataprepper/ggpld_snps.txt", 'hd':""}
def snpset(WC):
	return pull_or_not[WC.assay]

def runchoice(WC):
	run_dict = {'1':'merged_chrsplit','22':'merged_chrsplit','6':'merged_chrsplit', '2':'assay_chrsplit/','13':'assay_chrsplit/','7':'assay_chrsplit/', '9':'assay_chrsplit','12':'merged_chrsplit', '14':'assay_chrsplit','21':'assay_chrsplit/'}
	r = WC.run
	chrom = WC.chr
	if r == '1':
		location = run_dict[r] + '/run' + r+'/hol_testset.merge.chr' + chrom +'.bed'
	if r ==('2'):# or ('5'):#the sample by sample phasing files identified here need to have their sample referenced in the name, and the combine phasing samples do not have a rule of phasing.
		location = run_dict[r] + WC.sample + '.list1.chr' + chrom + '.bed'
	if r =='6':
		location = run_dict[r] + '/run' + r+'/hol_testset.merge.chr' + chrom +'.bed'
	if r ==('7'):
		location = run_dict[r] + WC.sample + '.list2.chr' + chrom + '.bed'
	if r ==('13'):
		location = run_dict[r] + WC.sample + '.list2.chr' + chrom + '.bed'
	if r ==('9'):
		location = run_dict[r] + WC.sample + '.list2.chr' + chrom + '.bed'
	if r ==('12'):
		location = run_dict[r] + '/run' + r + '/hol_testset.merge.chr' + chrom +'.bed'
	if r ==('14'):
		location = run_dict[r] + WC.sample + '.list3.chr' + chrom + '.bed'
	if r ==('21'):
		if ((WC.sample == 'snp50') or (WC.sample == 'ggpld')):
			location = run_dict[r] +  WC.sample + '.list1.chr' + chrom + '.bed'
		if ((WC.sample == 'f250') or (WC.sample == 'hd')):
			location = run_dict[r] +  WC.sample + '.list4.chr' + chrom + '.bed'
	if r == '22':
		location = run_dict[r] + '/run' + r+'/hol_testset.merge.chr' + chrom +'.bed'
	return location

def listchoice(WC):
	r = WC.run
	chrom = WC.chr
	if r ==('2'):# or ('5'):
		location =  'assay_chrsplit/' + WC.sample +'.list1.chr' + chrom
	if r ==('7') :
		location =  'assay_chrsplit/'+ WC.sample + '.list2.chr' + chrom
	if r ==('10'):
		location =  'assay_chrsplit/'+ WC.sample + '.list2.chr' + chrom
	if r =='13': #not sure this is the right one?
		location =  'assay_chrsplit/' + WC.sample +'.list3.chr' + chrom
	if r ==('21'):
		if ((WC.sample == 'snp50') or (WC.sample == 'ggpld')):
			location =  'assay_chrsplit/' + WC.sample +'.list1.chr' + chrom
		if ((WC.sample == 'f250') or (WC.sample == 'hd')):
			location =  'assay_chrsplit/' + WC.sample +'.list4.chr' + chrom
	return location

rule supplement_extract:
	input:
		bed = "correct_sex/{sample}.bed"
	params:
		inprefix = "correct_sex/{sample}",
		oprefix = "supplement_extract/{sample}"
	benchmark:
		"benchmarks/supplement_extract/{sample}.benchmark.txt"
	log:
		"logs/supplement_extract/{sample}.log"
	output:
		bed = "supplement_extract/{sample}.bed"
	shell:
		"(plink --bfile {params.inprefix} --cow --real-ref-alleles --remove raw_genotypes/hol.1970.txt --make-bed --out {params.oprefix})>{log}"

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
		"(plink --bfile {params.bfile} --real-ref-alleles {params.idslist} {params.extract} --make-bed  --cow --out {params.oprefix})> {log}"

files = ['correct_sex/227234.170619.1255.129_A', 'correct_sex/227234.170619.12703.100_A', 'correct_sex/227234.170619.1351.101_A', 'correct_sex/227234.170619.1667.550_A', 'correct_sex/227234.170619.172.103_A', 'correct_sex/227234.170619.219.102_A', 'correct_sex/227234.170619.442.112_A', 'correct_sex/227234.170619.500.104_A', 'correct_sex/227234.170619.74.124_A', 'supplement_extract/227234.170619.1994.200_A','correct_sex/777962.170619.11.550_A', 'correct_sex/777962.170619.136.124_A', 'correct_sex/777962.170619.1681.100_A', 'correct_sex/777962.170619.213.129_A', 'correct_sex/777962.170619.241.102_A', 'correct_sex/777962.170619.26.103_A', 'correct_sex/777962.170619.315.112_A', 'correct_sex/777962.170619.40.129_B', 'correct_sex/777962.170619.408.550_B', 'correct_sex/777962.170619.41.101_B', 'correct_sex/777962.170619.417.100_B', 'correct_sex/777962.170619.477.104_A', 'correct_sex/777962.170619.552.101_A', 'correct_sex/777962.170619.99.103_B', 'supplement_extract/777962.170619.2779.200_A', 'supplement_extract/777962.170619.411.200_B']

mergelists = {'hd':"hd_assays.txt", "f250":"f250_assays.txt"}
def listpicker(WC):
	assay = WC.sample
	return mergelists[assay]

oldnewlist = {'4':'1', '5':'2', '6':'3','1':'1'}
def oldnew(WC):
	name = "downsample/" + WC.sample + ".list" + oldnewlist[WC.list]
	return name

def oldnewbed(WC):
	name = "downsample/" + WC.sample + ".list" + oldnewlist[WC.list] + '.bed'
	return name

rule across_breed_assay_combine: #makes a single file for each assay, combining animals across breed, and the holstein test data. 
	input:
		sup = expand("{sample}.bed", sample = files), #all files, just to trigger generation of other rules. 
		down = oldnewbed
	params:
		list = listpicker,
		inprefix = oldnew,
		oprefix = "across_breed_assay_combine/{sample}.list{list}"
	benchmark:
		"benchmarks/across_breed_assay_combine/{sample}.benchmark.txt"
	log:
		"logs/across_breed_assay_combine/{sample}.log"
	output:
		combined = "across_breed_assay_combine/{sample}.list{list}.bed"
	shell:
		"plink --bfile {params.inprefix} --cow --real-ref-alleles --merge-list {params.list} --make-bed --out {params.oprefix}"

def downsample_orcombine(WC):
	if WC.sample == 'snp50':
		bed  = 'downsample/snp50.list' + WC.list + '.bed'
	if WC.sample == 'ggpld':
		bed  = 'downsample/ggpld.list' + WC.list + '.bed'
	if WC.sample == 'hd':
		bed  = 'across_breed_assay_combine/hd.list' + WC.list + '.bed'
	if WC.sample == 'f250':
		bed  = 'across_breed_assay_combine/f250.list' + WC.list + '.bed'
	if WC.sample == 'bigholsnp50':
		bed = 'bigholsnp50/bigholsnp50.list' + WC.list + '.bed'
	##plink --merge-list correct_sex/bighol.list --nonfounders -real-ref-alleles --cow --make-bed --out across_breed_assay_combine/bigholsnp50.list5
	#to make the bighol combined 
	return bed

rule assay_chrsplit:
	input:
		downsample_orcombine
	params:
		inprefix = "across_breed_assay_combine/{sample}.list{list}",
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
	priority: 100
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
		hap = temp("impute_input/run{run}/{sample}.chr{chr}.phased.haplotypes"),
		leg = temp("impute_input/run{run}/{sample}.chr{chr}.phased.legend"),
		log = "impute_input/run{run}/logs/{sample}.chr{chr}.phased.log"
	shell:
		"(shapeit -convert --input-haps {params.inprefix} --output-log {output.log} --output-ref {params.oprefix}) > {log}"

def findmergelist(wc):
	listdict = {'6':'2','1':'1','2':'1','3':'1','12': '3','21':'4','22':'5'}
	run = wc.run
	list = listdict[run]
	locate = 'assay_chrsplit/hol_testset.list' + list+'.chr'+wc.chr+'.txt'
	#print(locate)
	return locate

#for the data in the bigphase going to include stuff from list5 and list1
#basically its going to be list 1 + the 50k holsteins (list 5). List 1 data taken directly from the original directoy

DATA =['f250', 'ggpld', 'hd', 'snp50'] 		

rule make_merge_list: #makes a merge list of raw genotypes, doesn't need to have a run because its just going to be a list related thing.
	input:
		filelist = "assay_chrsplit/bigholsnp50.list{list}.chr{chr}.bed" ,
		filelist2 = expand("assay_chrsplit/{assay}.list1.chr{{chr}}.bed", assay= DATA )
	log:
		"logs/make_merge_list/list{list}.txt"
	output:
		"assay_chrsplit/hol_testset.list{list}.chr{chr}.txt"
	shell:
		"python ./bin/merge_file_maker.py  {input.filelist2} {input.filelist}  {output}"

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
	threads: 10
	priority: 30
	benchmark:
		"benchmarks/eagle_merged/run{run}/hol_testset.merge.chr{chr}.benchmark.txt"
	log:
		"logs/eagle_merged/run{run}/hol_testset.merge.chr{chr}.log"
	output:
		sample = "eagle_merged/run{run}/hol_testset.merge.chr{chr}.sample",
		haps = "eagle_merged/run{run}/hol_testset.merge.chr{chr}.haps.gz"
	shell:
		"(eagle --bfile={params.bed}  --geneticMapFile=USE_BIM --maxMissingPerSnp .99  --maxMissingPerIndiv .99 --numThreads 10 --outPrefix {params.out})> {log} "

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
		vcf=temp("eagle_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf"),
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
		vcfgz=temp("eagle_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf.gz"),
		index=temp("eagle_merged_vcf/run{run}/hol_testset.merge.chr{chr}.phased.vcf.gz.tbi")
	shell:
		"(bgzip {input.vcf}; tabix {output.vcfgz}) > {log};"

def pickidsforlist(wc):
	listdict = {'6':'2','1':'1','2':'1','3':'1','12':'3'}
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
		vcf = temp("vcf_per_assay/run{run}/{assay}.chr{chr}.vcf"),
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
		legend = temp("vcf_to_hap/run{run}/{assay}.chr{chr}.phased.legend"),
		haps = temp("vcf_to_hap/run{run}/{assay}.chr{chr}.phased.haplotypes"),
		sample = temp("vcf_to_hap/run{run}/{assay}.chr{chr}.phased.samples")
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
		hap = temp("vcf_to_haps/run{run}/{assay}.chr{chr}.phased.haps"),
		sample = temp("vcf_to_haps/run{run}/{assay}.chr{chr}.phased.sample")
	shell:
		"(bcftools convert {input.vcf} --hapsample {output.hap},{output.sample} ) > {log}"

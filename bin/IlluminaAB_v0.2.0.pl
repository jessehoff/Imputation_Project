#!/usr/bin/perl
# See: http://www.illumina.com/documents/products/technotes/technote_topbot.pdf
# This program takes an input file of SNPs and outputs the bases (A/C/G/T) 
# corresponding to the Illumina A/B format as well as the ref/alt alleles
# and IUPAC code for the SNP.
# Robert Schnabel <schnabelr@missouri.edu.
# University of Missouri
# v0.1.0 11/09/2016

# v0.2.0 02/21/2017
# - Added $RefShort to the input commands so that we can print out
#   the short name of the assembly to the output file.  This output
#   file gets uploaded to markers_[taxon_id]_alleles table in db.
#   This allows us to have mappings for multiple assemblies for each variant.

use autodie;
use Getopt::Long;
use strict;
no strict 'refs';

my %IUPAC;					# Hash for IUPAC codes
my %Comp;					# Hash for complement of bases
my %Chroms;					# Hash for chromosomes
my @Chroms;					# Array of chromosomes present in the input file

# REQUIRED OPTIONS
my $Ref = ();				# Prefix for reference genome
my $Infile = ();			# Input file name
my $KeepTMP = 0;			# 1=keep temp input files, 0=delete temp files. Used for debugging.
my $RefShort = ();			# Short abbreviation for reference from /REF_GENOMES/GENOMES.README
my $help;

# VARIABLES
my $RefIN;					# Reference file
my $fs;						# Field delimiter in input file
my $Rh;						# Reference file header line
my $warning;				# Iterator for the number of warnings found while processing data
my $chr;					# Integer chromosome number from input file
my $LocusName;				# Variant name from input file
my $VChr;					# Integer chromosome from input file
my $Pos;					# Variant position from input file
my $Variant;				# Variant [A/C] or A/C from input file
my $NumBases;				# Number of bases in the reference for each chromosome
my $NumSnp;					# Number of variants on a chromosome
my $TB;						# TOP/BOT
my $A;						# Illumina "A" allele
my $B;						# Illumina "B" allele
my $Step;					# The "walking step" distance for ambiguous variants needed to find an unambiguous position
my $Allele1;				# Allele1 of variant
my $Allele2;				# Allele2 of variant
my $LeftSeq;				# Left sequence flanking variant
my $RightSeq;				# Right sequence flanking variant
my $AltAllele;				# Alternate allele
my $iupac;					# IUPAC code for variant
my @REF;					# Chrom sequence as an array with no delimiter between elements
my @N;						# Locus names
my @C;						# Chromosomes
my @P;						# positions
my @V;						# Variant alleles "$Allele1,$Allele2,$Variant"

# Reference genome fasta files should have the same prefix and the X,Y,MT,XY chromosomes
# should have the leters replaced with the last autosome +1, +2, +3, +4
my $options = GetOptions (
	'help|?'					=> \$help,
	'ref=s'						=> \$Ref,
	'input=s'					=> \$Infile,
	'ref_abbv=s'				=> \$RefShort,
	'keep_tmp=s'				=> \$KeepTMP
);

if ($help) {
print <<HELP;
See: http://www.illumina.com/documents/products/technotes/technote_topbot.pdf
This program takes an input file of SNPs and outputs the bases (A/C/G/T) 
corresponding to the Illumina A/B format as well as the ref/alt alleles
and IUPAC code for the SNP.
Robert Schnabel <schnabelr[at]missouri.edu.
University of Missouri
v1.0 11/09/2016
*******************
REQUIRED OPTIONS
--ref		Reference genome name without the 'NN.fa'.
			This does not work with ChrU where all of
			the Unk contigs are in a multisequence fasta file.
			Reference file names should end with NN.fa,
			where NN is the chromosome number and is padded
			with a 0 (zero) when chrom < 10.
			Example file name: bt_ref_Bos_taurus_UMD_3.1_chr01.fa
			Example flag: --ref bt_ref_Bos_taurus_UMD_3.1_chr
			
--input		Input file name. Input file should be in subdirectory
			named 'input' and an 'output' directory shoudl be present.
			This will be used as the prefix
			for all output file names. The variant input file
			is tab or comma delimited and has the format:
			locus_name, chr, chr_pos, variant
			Where variant is [Allele1/Allele2]
			OR Allele1/Allele2 without the brackets. But the
			'/' must be present.  INDELS not supported.

--ref_abbv	Short abbreviation for reference from /REF_GENOMES/GENOMES.README
OPTIONAL
--keep_tmp	1=keep temp input files, 0=delete temp files. Used for debugging.					
*******************
HELP
exit;
}

##########
# Check to make sure required options are set
if (!$Ref or !$Infile or !$RefShort) {
	print "You must enter --ref --input --ref_abbv\n";
	print "run with --help to see options.";
	exit;
}

@IUPAC{('A/A','A/C','A/G','A/T','C/A','C/C','C/G','C/T','G/A','G/C','G/G','G/T','T/A','T/C','T/G','T/T')} = ('A','M','R','W','M','C','S','Y','R','S','G','K','W','Y','K','T');
%Comp = ( A => 'T', C => 'G', G => 'C', T => 'A');

#######################################
# Read in variant file to determine which chromosomes are present
# and make sure there are no issues with the file.
# We then re-read the input file and parse the variants into temporary
# files so that we can process them chromosome by chromosome.
# There are obviously more elegant ways to do this but this was quick and easy.
# locus_name, chr, chr_pos, variant
# Where variant is [Allele1/Allele2] or Allele1/Allele2 without the brackets.

open VIN, "<input/${Infile}";
print "\nReading input file: $Infile\n";
# Determine what the field seperator is in the input file
my $h = (<VIN>);
if ($h =~ m/,/) {
	$fs = ",";
	print "Field seperator: \"$fs\"\n";
}
elsif ($h =~ m/\t/) {
	$fs = "\t";
	print "Field seperator: \"$fs\"\n";
}
else {
	print "\nUnidentified field seperator\n";
	print "The field seperator should be comma or tab.\n";
	exit;
}

while (<VIN>) {
	chomp $_;
	($LocusName,$VChr,$Pos,$Variant) = split(/$fs/,$_);
	$Chroms{$VChr} = $Pos;
}
close VIN;

foreach (keys %Chroms) {
	push @Chroms, $_;
}

@Chroms = sort @Chroms;

foreach (@Chroms) {
	print "Creating TEMP file for Chrom:$_\n";
	open "OUT${_}", ">input/${Infile}.${_}.TMP";
}

open VIN, "<input/$Infile";
open ERR, ">output/${Infile}.ERR";
print "\nReading input file: $Infile\n";
$h = <VIN>;														# Header row from variant input file
while (<VIN>) {
	chomp $_;
	($LocusName,$VChr,$Pos,$Variant) = split(/$fs/,$_);
	$Variant =~ s/\[|\]//g;										# Remove variant brackets[] if present
	$Allele1 = substr $Variant, 0, 1;							# Grab first allele
	$Allele2 = substr $Variant, 2, 1;							# Grab second allele
	# Check to make sure the variant is valid. If not valid then write to ERR file
	if ( not exists $Comp{$Allele1} or not exists $Comp{$Allele2} ) {
		print "INVALID genotype $Allele1/$Allele2 at $VChr:$Pos\n";
		print ERR "INVALID genotype $Allele1/$Allele2 at $VChr:$Pos\n";
	}
	else {
		print {"OUT${VChr}"} "$LocusName,$VChr,$Pos,$Variant\n";
	}
}

foreach (@Chroms) {
	print "Closing TEMP file for Chrom:$_\n";
	close "OUT${_}";
}

$" = "";	# Set array delimiter to none so that when we print the top genomic seq there are no spaces between characters

#my $chr = $StartChr;
#while ($chr <= $LastChr) {

# One output file for all variants and print header row
open VOUT, ">output/${Infile}.OUT"; 
print VOUT "locus_name,assembly,chr,pos,left_seq,snp,right_seq,step,ref_base,alt_allele,a_allele,b_allele,top_bot,iupac\n";

CHROM:
foreach (@Chroms) {
	$chr = $_;
	$warning = 0;
	# Check to see if the TMP input file has zero length and move to next Chrom if no variants
	# This could happen if there are a small number of input variants and they got removed
	# during the initial screening. By moving to the next chromosome we save the expensive 
	# step of reading the reference file.
	my $FileSize = `du "input/${Infile}.${chr}.TMP" | cut -f1`;
	if ( $FileSize == 0 ) {
		print "Chromosome $chr had no variants after filtering so skipping\n";
		system ("rm input/${Infile}.${chr}.TMP");
		next CHROM;
	}
	
	#######################################
	# Read in reference file and push sequence into array
	# Arrayes are indexed starting at zero so the genome coordinates minus 1
	# will be equal to the array index position for that location in the genome
	if ($chr < 10) { $RefIN = "${Ref}0${chr}.fa"; }
	else { $RefIN = "${Ref}${chr}.fa"; }
	open RIN, "<$RefIN";

	print "\nReading input file: $RefIN\n";
	$Rh = <RIN>;	# Header line from fasta file which will be the sequence identifier
	while (<RIN>) {
		chomp $_;
		push (@REF, split("", $_));
		# Check to make sure that there is only one sequence identifier
		if ($_ =~ m/>/) {
			print "Additional sequence identifiers identified...\n";
			print "LINE $. : $_\n";
		}
	}
	close RIN;

	$NumBases = scalar @REF;
	print "Chr:$chr\tNumBases:$NumBases\tNumLines:$.\n";

	#######################################
	# Read in variant file
	# locus_name, chr, chr_pos, variant
	# Where variant is [Allele1/Allele2] or Allele1/Allele2 without the brackets.

	open VIN, "<input/${Infile}.${chr}.TMP";
	print "\nReading input file: ${Infile}.${chr}.TMP\n";
	while (<VIN>) {
		chomp $_;
		($LocusName,$VChr,$Pos,$Variant) = split(/,/,$_);				# TMP input file is comma delimited
		$Variant =~ s/\[|\]//g;											# Remove variant brackets[] if present
		$Allele1 = substr $Variant, 0, 1;								# Grab first allele
		$Allele2 = substr $Variant, 2, 1;								# Grab second allele
		# Check to make sure the variant is valid
		# If neither of the SNP alleles match the reference then complement them
		if ((@REF[$Pos - 1] ne $Allele1) and (@REF[$Pos - 1] ne $Allele2)) {
			$Allele1 = $Comp{$Allele1};									# Complement of first allele
			$Allele2 = $Comp{$Allele2};									# Complement of second allele
		}
		# After complementing the alleles, make sure at least one of them match the reference.
		# If not, print to ERR file and throw a warning. This check will be triggered if there are
		# indels present in the input file.
		if ((@REF[$Pos - 1] ne $Allele1) and (@REF[$Pos - 1] ne $Allele2)) {
			print "WARNING: $LocusName,Chr:$VChr,Pos:$Pos,$Variant,ref:@REF[$Pos - 1],A1:$Allele1,A2:$Allele2\n";
			$warning++;
		}
		# Check to make sure the reference base is A/C/G/T
		elsif ( not exists $Comp{@REF[$Pos - 1]} ) {
			print "INVALID REFERENCE BASE \"@REF[$Pos - 1]\" AT $VChr:$Pos\n";
			print ERR "INVALID REFERENCE BASE \"@REF[$Pos - 1]\" AT $VChr:$Pos\n";
			$warning++;
		}
		else {
			push (@N,$LocusName);						# locus names
			push (@C,$VChr);							# chromosomes
			push (@P,$Pos);								# positions
			push (@V,"$Allele1,$Allele2,$Variant");		# variant alleles
		}
		# Exit if there are too many discrepencies between the reference base indicated in
		# the variant file and the reference base from the reference sequence or if there
		# are too many reference positions that are invalid bases.
		if ($warning >= 100) {
			print "\n\nERROR: Exiting program due to too many reference base discrepencies.\n";
			print "See error file. Are there indels in the input file?  These are not supported.\n";
			exit;
		}
	}
	close VIN;
	$NumSnp = scalar @N;

	print "\nChromosome: $chr\n";
	print "Number of SNP: $NumSnp\n\n";

	##############################################################################
	# Extract the flanking sequence for all of the valid variants
	my $i = 0;
	foreach (@P) {
		$Pos = $_;
		$TB = ();
		$A = ();
		$B = ();
		$Step = 0;
		($Allele1,$Allele2,$Variant) = split(/,/,@V[$i]);
		# The Alleles need to be in alpha order T/A --> A/T, G/C --> C/G, T/C --> C/T etc.
		if ($Allele1 gt $Allele2) {
			my $Temp1 = $Allele1;
			$Allele1 = $Allele2;
			$Allele2 = $Temp1;
		}
		# Unambiguous cases
		if ($Allele1 eq "A" and ($Allele2 eq "C" or $Allele2 eq "G")) {		# A/(C or G)
			# Ref base = A 
			$TB = "TOP";
			$A = "A";
			$B = "$Allele2";
		}
		elsif ( $Allele2 eq "T" and ($Allele1 eq "C" or $Allele1 eq "G") ) {	# (C or G)/T
			$TB = "BOT";
			$A = "T";
			$B = "$Allele1";
		}
		# Ambiguous cases
		else {
			$Step = 0;
			STEP:
			while ($TB == 0) {
				$Step++;
				my $left = @REF[($Pos - 1 - $Step)];
				my $right = @REF[($Pos - 1 + $Step)];
				if ( $left =~ m/[AT]/ ) {				# Left is A/T 
					if ( $right =~ m/[AT]/ ) { 			# Left is A/T and right is A/T so move to next base
						next STEP;
					}
					$TB = "TOP";						# Left is A/T, Right is not so TOP designation
					$A = "$Allele1";
					$B = "$Allele2";
					last STEP;
				}
				elsif ( $right =~ m/[AT]/ ) {			# Left NOT A/T and right IS A/T so BOT designation
					$TB = "BOT";
					$A = "$Allele2";
					$B = "$Allele1";
					last STEP;
				}
				else { next STEP; }						# Left NOT A/T and right NOT A/T so move to next base
			}
		}

		# Assign alternate allele
		if ($A eq @REF[($Pos - 1)] and $B ne @REF[($Pos - 1)]) {
			$AltAllele = $B;
		}
		else {
			$AltAllele = $A;
		}
		$iupac = $IUPAC{"$A\/$B"};
		
		# We only need to pull flanking sequence for ambiguous alleles
		$LeftSeq = ();
		$RightSeq = ();
		if ( $Step > 0 ) {
			$LeftSeq = "@REF[($Pos - 1 - $Step)..($Pos - 2)]";
			$RightSeq = "@REF[($Pos)..($Pos - 1 + $Step)]";
		}

#		print VOUT "locus_name,assembly,chr,pos,left_seq,snp,right_seq,step,ref_base,alt_allele,a_allele,b_allele,top_bot\n";

		print VOUT "@N[$i],";									# Locus_Name
		print VOUT "$RefShort,"									# Assembly
		print VOUT "@C[$i],";									# Chromosome
		print VOUT "$Pos,";										# Coordinate
		print VOUT "$LeftSeq,";									# Left flanking sequence
		print VOUT "${Variant},";								# Input SNP
		print VOUT "$RightSeq,";								# Right flanking sequence
		print VOUT "$Step,";									# Number of bases needed to find unambiguous
		print VOUT "@REF[($Pos - 1)],";							# Reference base
		print VOUT "$AltAllele,";								# Alt allele
		print VOUT "$A,";										# A allele
		print VOUT "$B,";										# B allele
		print VOUT "$TB,";										# TOP/BOT designation
		print VOUT "$iupac\n";									# IUPAC designation
		$i++;
	}
	@REF = ();
	@N = (); 
	@C = (); 
	@P = (); 
	@V = (); 
	
	# Delete the TMP input files
	if ( $KeepTMP == 0 ) { system ("rm input/${Infile}.${chr}.TMP"); }
}
close VOUT;
exit;

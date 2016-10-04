= BOX 1: Input Files ==========================================================
PedigreeFile				,pedtest.txt
GenotypeFile				,snp50_c.chr2.Alpha
TrueGenotypeFile			,None
= BOX 2: Sex Chromosome ========================================================
SexChrom				,No
= BOX 3: SNPs ==================================================================
NumberSnp				,2526
MultipleHDPanels			,0
NumberSnpxChip                          ,0,0
HDAnimalsThreshold			,2.0
= BOX 4: Internal Editing =======================================================
InternalEdit				,No
EditingParameters			,0.0,0.0,0.0,AllSnpOut
= BOX 5: Phasing ================================================================
NumberPhasingRuns			,2
CoreAndTailLengths			,200,400
CoreLengths				,100,300
PedigreeFreePhasing			,No
GenotypeError				,0.0
NumberOfProcessorsAvailable		,4
LargeDatasets                           ,No,200,1
= BOX 6: Imputation =========================================================
InternalIterations			,5
ConservativeHaplotypeLibraryUse		,No
WellPhasedThreshold			,99.0
= BOX 7: Hidden Markov Model ================================================
HMMOption				,No
TemplateHaplotypes			,200
BurnInRounds				,5
Rounds					,20
ParallelProcessors			,1
Seed					,-123456789
ThresholdForMissingAlleles		,50.0
ThresholdImputed			,90.0
= BOX 8: Running options ====================================================
PreprocessDataOnly			,No
PhasingOnly				,No
UserDefinedAlphaPhaseAnimalsFile	,None
PrePhasedFile				,None
BypassGeneProb				,No
RestartOption				,1

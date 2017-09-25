
###############################################################################
# Purpose:  SnakeMake File TEMPLATE to check for Coverage and QC on the Bam files
# Authors: Nick Brazeau  & Dayne Filer   & Christian Parobek                                               
#Given: Realn.Bam
#Return: coverage, qc       
#Remarks: As before, you will need to change your paths 
#      	  There is a clash issue with snakemake using python 3 and multiqc python 2. Need to run from command line      
#		  Note, the heatmap coverage plot and the bed map plot are not very generalizable and need to be updated. Don't use...depreciated and not up to date at this time   
###############################################################################


####### Working Directory and Project Specifics ############
WRKDIR = '/proj/ideel/meshnick/users/NickB/Projects/Rhea_Plasmepsin_Rx/'
readWD = '/proj/ideel/meshnick/users/NickB/Projects/Rhea_Plasmepsin_Rx/'
#SAMPLES, = glob_wildcards(WRKDIR + 'symlinks/{sample}_R1.fastq.gz')
FQS, = glob_wildcards(WRKDIR + 'symlinks/{fq}.fastq.gz')
SAMPLES, = glob_wildcards(WRKDIR + 'aln/{sample}.realn.bam')

####### Turn on for Pv ##########
#REF = '/proj/meshnick/Genomes/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta'
#GFF = '/proj/meshnick/Genomes/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1.gff'

##### Turn on for Pf PlasmoDB Version 13.0 ##########
#REF = '/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta'
#GFF = '/proj/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7.gff'

##### Turn on for Pf Sanger Version 13.0 BETA ##########
REF = '/proj/ideel/meshnick/Genomes/Pf3D7_Sanger2017/Pfalciparum.genome.fasta'
GFF = '/proj/ideel/meshnick/Genomes/Pf3D7_Sanger2017/Pfalciparum.gff3'
# Note this version is compatible with the Pf3k Pipeline 

######## Tools to Call #########
PICARD = '/nas/longleaf/home/nfb/.linuxbrew/Cellar/picard-tools/2.9.0/share/java/picard.jar'
GATK = '/nas02/home/n/f/nfb/.linuxbrew/Cellar/gatk/3.6/share/java/GenomeAnalysisTK.jar'
FLASH = '/nas/longleaf/home/nfb/.linuxbrew/Cellar/flash/1.2.11/bin/flash'
TMPDIR = '/pine/scr/n/f/nfb/PicardandGATKscratch'


##########################################################################################

############################
####### Target #############
############################
#Rule all checks to see if file is the same and follows directions up to specified point
rule all:
# SUMMARYSTATS
#	input: 'symlinks/FASTQC/Log.txt'
	input: expand('SumSTATsandQC/ValidateBAM/{sample}.validatebams.txt', sample = SAMPLES)
#	input: expand('SumSTATsandQC/AlignSummary/{sample}.AlignSummary.Metrics', 'SumSTATsandQC/FlagStats/{sample}.samtools.flagstats', 'SumSTATsandQC/coverage/data/{sample}.cov', 'SumSTATsandQC/coverage/d_bedtools/{sample}.dbedtools.cov', sample = SAMPLES)
#	input: expand('SumSTATsandQC/AlignSummary/{sample}.AlignSummary.Metrics', sample = SAMPLES)
#	input: expand('SumSTATsandQC/FlagStats/{sample}.samtools.flagstats', sample = SAMPLES)
#	input: expand('SumSTATsandQC/coverage/data/{sample}.cov', sample = SAMPLES)
#	input: expand('SumSTATsandQC/coverage/d_bedtools/{sample}.dbedtools.cov', sample = SAMPLES)
#	input: 'SumSTATsandQC/coverage/{params.prefix}cov_heatmap.pdf'
#	input: 'SumSTATsandQC/coverage/{params.prefix}_individual_bedgraph_plot.pdf'
#	input: 'multiqc.report'

###############################################################################


#####################################
########   MultiQC Summary   ########
#####################################
rule multiqc: 
	input: alignmet='SumSTATsandQC/AlignSummary/', flagstats='SumSTATsandQC/FlagStats/', fastqc='symlinks/FASTQC/'
	output: 'multiqc.report'
	shell: 'multiqc {input.alignment} {input.flagstats} {input.fastqc}'
#This is going to clash with python versions. Need to figure out later

######################################
########   Quality Control   #########
#####################################

rule RunQC:
	output: 'QCCompleted.tab.txt'
	shell: 'bash QCdiscovery.sh ; echo "QC Performed and completed" > {output}'



######################################
######## Summary Statistics #########
#####################################
# rule plot_bedgraph:
# 	input:	expand('SumSTATsandQC/coverage/d_bedtools/{sample}.dbedtools.cov', sample = SAMPLES)
# 	params:	idir = 'SumSTATsandQC/coverage/d_bedtools/', 
# 		prefix = '1877_Plasmepsin_Samples_', 
# 		odir = 'SumSTATsandQC/coverage/',
# 		jackpot = '50',
# 		window = '1000',
# 		yaxis = 'names/1877Sample.list',
# 		title= '1877 Samples',
# 	output:	'SumSTATsandQC/coverage/{params.prefix}_individual_bedgraph_plot.pdf', 
# 		'SumSTATsandQC/coverage/{params.prefix}_stacked_bedgraph.pdf'
# 	shell:	'Coverage_Analyzer_Bedgrapher -I {params.idir} -O {params.odir} -F {params.prefix} -J {params.jackpot} -W {params.window} -Y {params.yaxis} -T {params.title} '

rule plot_coverage:
	input:	expand('SumSTATsandQC/coverage/data/{sample}.cov', sample = SAMPLES)
	params:	idir = 'SumSTATsandQC/coverage/data/', 
		prefix = '1877_Plasmepsin_', 
		odir = 'SumSTATsandQC/coverage/'
	output:	'SumSTATsandQC/coverage/{params.prefix}cov_plot.pdf', 
		'SumSTATsandQC/coverage/{params.prefix}cov_heatmap.pdf'
	shell:	'covPlotter -I {params.idir} -O {params.odir} -F {params.prefix}'

# rule d_bedtools_cov_forbedgraph:
# 	input: 'aln/{sample}.recal.realn.bam'
# 	output: 'SumSTATsandQC/coverage/d_bedtools/{sample}.dbedtools.cov'
# 	shell: 'bedtools genomecov -d -g {GFF} -ibam {input} > {output}'


rule calculate_cov_forheatmap:
	input:	'aln/{sample}.recal.realn.bam'
	output:	'SumSTATsandQC/coverage/data/{sample}.cov'
	shell:	'bedtools genomecov \
		-ibam {input} -g {GFF} -max 500 | grep "genome\|PFC10_API\|M76\|Pf3" \
		> {output}'
# Note for Pv switch out to: -ibam {input} -g {GFF} -max 500 | grep "genome\|chr\|NC_" \
	
rule summary_stats:
	input: 'aln/{sample}.recal.realn.bam'
	output: 'SumSTATsandQC/FlagStats/{sample}.samtools.flagstats'
	shell: 'samtools flagstat {input} > {output}'
# Provides alignment information
		
rule AlignSummaryMetrics:
	input: 'aln/{sample}.recal.realn.bam'
	output: 'SumSTATsandQC/AlignSummary/{sample}.AlignSummary.Metrics'
	shell: 'java -jar {PICARD} CollectAlignmentSummaryMetrics \
          R={REF} \
          I={input} \
          O={output}'
# Need this for read length averages, etc. 

###############################################
########   VALIDATE SAM FILE MODULE   #########
###############################################
rule ValidateSamFile:
	input:'aln/{sample}.realn.bam'
	output: 'SumSTATsandQC/ValidateBAM/{sample}.validatebams.txt'
	shell: 'java -jar {PICARD} ValidateSamFile \
	    MODE=SUMMARY \
		I={input} \
		OUTPUT={output}'
	# After this step you should run from the command line `cat 'SumSTATsandQC/ValidateBams.tab.txt | grep "ERROR"` -- if there are errors, STOP and figure out why	
		
######################################
########   FASTQC MODULE   #########
#####################################
rule fastqc:
	input: expand('symlinks/{fq}.fastq.qz', fq=FQS)
	output: outdir='symlinks/FASTQC/', log='symlinks/FASTQC/Log.txt',
	shell: 'fastqc {input} --outdir {output.outdir}; echo "Fastqc succesfully completed" > {output.log} '

	

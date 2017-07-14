###############################################################################
# Purpose:  SnakeMake File to take FASTQs and make BAMs
# Authors: Christian Parobek & Nick Brazeau                                        
#Given: FASTQ
#Return: BAM 
#Remarks: Note, you will have to change the paths and project specifics as well as the tools paths (depending on user) for each project/each new directory                               
#		  You will also need to make your own intervals/all_chrs.intervals file depending on your project and your reference sequence
#         You will need to add the BamMergeCommander to your path
#         Input for BamMergeCommander is a tab delimited file that first column contains a header and sample names or sample IDs (some unique n.alphanumeric for your sample)   
###############################################################################


####### Working Directory and Project Specifics ############
WRKDIR = '/proj/ideel/meshnick/users/NickB/Projects/DRC_Pf_Project/PhyloSNPsMIPs/HathawayDownload/malaria-gen/'
readWD = '/proj/ideel/meshnick/users/NickB/Projects/DRC_Pf_Project/PhyloSNPsMIPs/HathawayDownload/malaria-gen/'
DATEDSAMPS, = glob_wildcards(WRKDIR + 'fastq/{ds}_R1.fastq.gz')
MERGEDSAMPS, = glob_wildcards(WRKDIR + 'aln/{merge}.merged.bam')
MTDT - 

####### Turn on for Pv ##########
#REF = '/proj/ideel/meshnick/Genomes/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1_Genome.fasta'
#GFF = '/proj/ideel/meshnick/Genomes/PvSAL1_v13.0/PlasmoDB-13.0_PvivaxSal1.gff'

##### Turn on for Pf ##########
#REF = '/proj/ideel/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7_Genome.fasta'
REF2 = '/proj/ideel/meshnick/Genomes/Pf3D7_Pf3KDownload/Pfalciparum.genome.fasta'
GFF = '/proj/ideel/julianog/refs/Pf3D7_v13.0/PlasmoDB-13.0_Pfalciparum3D7.gff'



######## Tools to Call #########
PICARD = '/nas02/apps/picard-2.2.4/picard-tools-2.2.4/picard.jar'
GATK = '/nas02/home/n/f/nfb/.linuxbrew/Cellar/gatk/3.6/share/java/GenomeAnalysisTK.jar'
TMPDIR = '/pine/scr/n/f/nfb/PicardandGATKscratch'
#TMPDIR = '/lustre/nfb/tmp_for_picard/'
#pine for longleaf and lustre for killdevil

##########################################################################################

############################
####### Target #############
############################
#Rule all checks to see if file is the same and follows directions up to specified point
rule all:
#	input: expand('aln/{ds}.bam', ds = DATEDSAMPS)  
#    input: expand('aln/{ds}.sorted.bam', ds = DATEDSAMPS)
#	input: expand('aln/{ds}.dedup.bam', ds = DATEDSAMPS) 
#	input: 'merge.log.file'                                                 # Run to here and then check file
#    input: expand('aln/{merge}.merged.bam.bai', merge = MERGEDSAMPS) 
#   input: expand('aln/{merge}.realigner.intervals', merge = MERGEDSAMPS) 
#	input: expand('aln/{merge}.realn.bam', merge = MERGEDSAMPS) 

###############################################################################




############################
######## Alignment #########
############################
rule index_realigned: 
	input: 'aln/{merge}.realn.bam'
	output: 'aln/{merge}.realn.bam.bai'
	shell: 'java -jar {PICARD} BuildBamIndex INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'

rule realn_indels:
	input: bam = 'aln/{merge}.merged.bam', chrs = 'intervals/all_chrs.intervals', targets = 'aln/{merge}.realigner.intervals', 
	output: 'aln/{merge}.realn.bam'
	shell: 'java -jar {GATK2} -T IndelRealigner \
		-R {REF} -I {input.bam} \
		-L {input.chrs} -targetIntervals {input.targets} \
		-o {output}' 
		# all_chrs.intervals includes just chrs and mito -- it is similar to a bed file for GATK Caller
		# -fixMisencodedQuals must be added for SRA data
		# Interval file from CMP direcotry

rule find_indels:
	input: bam = 'aln/{merge}.merged.bam', index = 'aln/{merge}.merged.bam.bai', chrs = 'intervals/all_chrs.intervals'
	output: 'aln/{merge}.realigner.intervals'
	shell: 'java -jar {GATK2} -T RealignerTargetCreator \
		-R {REF} -I {input.bam} \
		-L {input.chrs} -o {output}'
		# all_chrs.intervals includes just  chrs and mito

rule index_merged: 
	input: 'aln/{merge}.merged.bam'
	output: 'aln/{merge}.merged.bam.bai'
	shell: 'java -jar {PICARD} BuildBamIndex INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'


rule merge_dedups:
	input: 'aln/CmdlineMerge_PvRelapseProject.sh'
	output: 'merge.log.file'
	shell: 'bash {input}; echo "The merge was completed using the custom Rscript BamMergeCommander -- good idea to check the CmdlineMerge table to confirm merge was appropriate" > {output}'

rule make_merge_commandLine:
	input:  metadata='{MTDT}',
	params: dedupbams='{WRKDIR}/aln', 
	    outdir = '{WRKDIR}/aln/',
	output: 'aln/CmdlineMerge_PvRelapseProject.sh'
	shell: 'BamMergeCommander --input {params.dedupbams} --metadata {input.metadata} --outdir {params.outdir} --output {output}'
# note outdir needs forward slash, dedupbams does not

rule mark_dups:
	input: 'aln/{ds}.sorted.bam'
	output:'aln/{ds}.dedup.bam','aln/{ds}.dedup.metrics'
	shell: 'java -jar {PICARD} MarkDuplicates \
		I={input} O={output[0]} \
		METRICS_FILE={output[1]} \
		TMP_DIR={TMPDIR} REMOVE_DUPLICATES=TRUE \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000'

rule sort_bam:
	input: 'aln/{ds}.matefixed.bam'
	output: 'aln/{ds}.sorted.bam'
	shell: 'java -jar {PICARD} SortSam \
		I={input} O={output} \
		SO=coordinate TMP_DIR={TMPDIR}'

rule matefix:
	input: 'aln/{samp}.bam'
	output: 'aln/{samp}.matefixed.bam'
	shell: 'java -jar {PICARD} FixMateInformation INPUT={input} OUTPUT={output} TMP_DIR={TMPDIR}'


rule fastq_to_bam:
	input: 'symlinks/{ds}_R1.fastq.gz', 'symlinks/{ds}_R2.fastq.gz'
	output: 'aln/{ds}.bam'
	shell: 'bwa mem {REF2} {readWD}{input[0]} {readWD}{input[1]} \
		-R "@RG\tID:bwa\tPL:illumina\tLB:{wildcards.ds}\tSM:{wildcards.ds[0]}{wildcards.ds[1]}{wildcards.ds[2]}{wildcards.ds[3]}{wildcards.ds[4]}" \
		 samtools view -Sb - > {output}'
		# calling the @RG ID: 'bwa' because this resolves a clash with @PG ID --> I updated this recently to make it more unique for MERGING
		# Can controls how long the read sample name is by wild cards for tSM, important if want to merge file later and need different library names but same sample names

# rule trim_illumina_Adaptors_fastqs:
# 	 input: 'symlinks/{ds}_R1.fastq.gz', 'symlinks/{ds}_R2.fastq.gz', 
# 	 output: 'symlinks/{ds}_R1.PAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/{ds}_R1.UNPAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/{ds}_R2.PAIREDtrimmomatictrimmed.fastq.gz', 'symlinks/{ds}_R2.UNPAIREDtrimmomatictrimmed.fastq.gz',  
# 	 shell: 'trimmomatic PE -threads 12 -trimlog symlinks/Relapse/trim_log.txt {input[0]} {input[1]} {output[0]} {output[1]} {output[2]} {output[3]} ILLUMINACLIP:/nas/longleaf/home/nfb/.linuxbrew/Cellar/trimmomatic/0.36/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
#     # Trimmomatic needed if illumina adpators are attached. The TRUE at the end keeps the paired end reads in R2
#     # Want to align the PAIRED trimmed
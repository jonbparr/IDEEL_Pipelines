# Amplicon Tutorial for IDEElers
## Session 3: Alignments

Here we will explore basic alignment strategies.

### Resources
1. Read in its entirety -- including frequently asked questions and review each ppt available -- GATK's Best Practices for [Germline Muations in WGS and Exome Sequences, Pre-Processing Page](https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS&p=1)
2. Review the articles and manuals for the most commonly used aligners: [bwa-mem](http://bio-bwa.sourceforge.net/) and [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
3. Review the format of a [Sam File](http://samtools.github.io/hts-specs/SAMv1.pdf)


### Exercises


Simulate Fasta and align it to a referent genome following the Rmarkdown file `ToyDataGenerator.Rmd`. Make sure to change the probability of the A,C,T,G calls and see how it affects the fasta file. Open the fasta file with a program like `Mega` or `AliView`. Paste several fasta files together and perform multiple pairwise alignment using `Muscle` or `CLUSTALW` (preloaded in the previously mentioned programs or in the R package `Ape`. 
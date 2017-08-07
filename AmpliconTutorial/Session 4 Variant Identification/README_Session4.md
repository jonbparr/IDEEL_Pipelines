# Amplicon Tutorial for IDEElers
## Session 4: Variant Discovery

Here we will explore basic alignment strategies.

### Resources
1. Read in its entirety -- including frequently asked questions and review each ppt available -- GATK's Best Practices for [Germline Muations in WGS and Exome Sequences, Variant Discovery Page](https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS&p=2) and [Germline Muations in WGS and Exome Sequences, Callset Refinement](https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS&p=3) 
2. **Read in its entirety** [What is a VCF and how should I interpret it? by GATK's Geraldine](https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it)

### Exercises

1. Finish the simulation from the Rmarkdown file `ToyDataGenerator.Rmd`. Explore the VCF that you created. Does it produce the mutation that you introduced in the `ToyDataGenerator.Rmd`? Introduce more mutations or change the frequency at which you expect to see the mutatin by concatenating several different fastas together in R or on the command line. 

2. Load the VCF file into R using the R package `VCFR`, you can dowload it following the directions [here](https://github.com/knausb/vcfR). Manipulate the VCF and look at GT and depth of coverage. Complete the `VCFR` tutorials (i.e Vignettes) located [here](https://cran.r-project.org/web/packages/vcfR/index.html). **This will help provide an issue and manageable way to refine your VCF discoveries to higher quality SNPs.**
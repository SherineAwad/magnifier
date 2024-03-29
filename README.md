
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


Magnifier 
============================================================

Magnifier is a snakemake  pipeline that help magnifies what's in your sample, check GC bias, coverage, etc. More ways of investigating your samples will be included as we go.
 
Update config files for different with you samples names and other parameters. 

## Plots
 
#### Mapping Percentages  

You will get a plot of percentages of primary and secondary alignments. In addition to percentage of unmapped reads and quality of mapped reads for each sample. 

   ![sample1_alignments.png](outputs/sample1.s_1.alignments.png)

#### Mapping Quality 

A density plot of the mapping quality of your  reads.

   ![sample1.s_1.mapq.png](outputs/sample1.s_1.mapq.png) 

#### GCBias  

   Investigates GCbias in your samples. It output several metrics (txt files) and wraps them up in a figure, as follows: 
    
   ![GCBias.png](outputs/GCbias.png)


#### Overall Coverage 

Outputs overall coverage plot: 

   ![sample1.s_1.coverage.png](outputs/sample1.s_1.coverage.png)

#### Coverage Histogram  

Outputs a histogram of coverage of your sample, a sample output example is:
  
   ![coveragehist.png](outputs/coveragehist.png) 

#### InsertSize

Outputs several metrics for insert size, and a wrapped up figure as follows: 

   ![insertsize.png](outputs/insertsize.png)


#### Check Contamination 

Outputs a nice plot as below to check contamination: 
  
   ![sample1_screen.png](outputs/sample1_screen.png)

### Wrapper all output to HTML Page 

   [sample1.s_1.html](outputs/sample1.s_1.html) 

## More Detailed Metrics  

The pipeline will also generate more text files with more detailed stats: 

#### More Detailed Alignments Metrics 

   [sample1.s_1.alignment_metrics.txt](/outputs/sample1.s_1.alignment_metrics.txt)
 

#### More Detailed GC Bias Metrics 

   [sample1.s_1.gc_bias_metrics.txt](/outputs/sample1.s_1.gc_bias_metrics.txt)

#### More Detailed Insert Size Metrics 

   [sample1.s_1.insert_size_metrics.txt](/outputs/sample1.s_1.insert_size_metrics.txt)


### Run Snakemake 

Use: 

   snakemake -jn 
 
to run the pipeline where n is the number of cores for example for 10 cores use. Snakemake has to be installed. 

### Use conda 

For less froodiness, use conda:


    snakemake -jn --use-conda 


For example, for 10 cores use: 

    snakemake -j10 --use-conda 

This will pull automatically the same versiosn of tools we used. Conda has to be installed in the system, in addition to snakemake. 
Conda as well has to be installed. 


### Dry Run


For a dry run use: 
  
  
    snakemake -j1 -n 


and to print command in dry run use: 

  
    snakemake -j1 -n -p 


#### References 

1. https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen
2. http://broadinstitute.github.io/picard/  

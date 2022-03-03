
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


Magnifier 
============================================================

Magnifier is a snakemake  pipeline that help magnifies what's in your sample, check GC bias, coverage, etc. More ways of investigating your samples will be included as we go.
 
Update config files for different with you samples names and other parameters. 


### Wrapper all output to HTML Page 

<h1> Contamination Check </h1> <br> </br>  <center><img src='sample1.s_1.r_1_screen.png' alt='Mapping to Multiple Genomes' class='cover' width='75%' height='75%' /></center><h1> Alignments </h1><h2> <a href='sample1.s_1.alignment_metrics.txt '> Alignments Metrics detailed here</a></h2><h2> Alignments summary </h2> <br> </br>  <img src='sample1.s_1.alignments.png' alt='Alignments' width='100% height='100%'/><h2> Mapping Quality </h2> <br> </br>  <img src='sample1.s_1.mapq.png' alt='Mapping Quality' width='100% height='100%'/><h1> Coverage</h1> <iframe src='sample1.s_1.coverage.histogram.txt#scrollbar=0&toolbar=0&scrolling=0' frameBorder='0' height='600' width='600'> </iframe> </center><h1> GC Bias</h1> <h2><a href='sample1.s_1.gc_bias_metrics.txt '> GC Bias Metrics detailed here</a></h2><h2> GC Bias plot </h2><center><iframe src='sample1.s_1.gc_bias_metrics.pdf#scrollbar=0&toolbar=0&scrolling=0' scrolling='no' seamless='seamless'frameBorder='0' height='600' width='600'> </iframe> </center><h1> Insert Size</h1> <h2> <a href='sample1.s_1.insert_size_metrics.txt'> Insert Size Metrics detailed here</a></h2><h2> Insert Size </h2><center><iframe src='sample1.s_1.insert_size_histogram.pdf#toolbar=0&scrollbar=0&scrolling=0' scrolling='no' seamless='seamless' frameBorder='0' height='600' width='600'> </iframe> </center></html><html>


### Run Snakemake 

    snakemake -jn 

where n is the number of cores for example for 10 cores use:


    snakemake -j10 

### Use conda 

For less froodiness, use conda:


    snakemake -jn --use-conda 


For example, for 10 cores use: 

    snakemake -j10 --use-conda 

This will pull automatically the same versiosn of tools we used. Conda has to be installed in the system, in addition to snakemake. 


### Dry Run


For a dry run use: 
  
  
    snakemake -j1 -n 


and to print command in dry run use: 

  
    snakemake -j1 -n -p 


#### References 

1. https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen
2. http://broadinstitute.github.io/picard/  

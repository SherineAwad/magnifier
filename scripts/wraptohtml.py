#! /usr/bin/env python
import argparse
import os
import sys
import math
import glob

def html_report(sample):
    outfile = sample+".html" 
    outtext = sample
    hfile = open(outfile, "w+") 
    message = """<html>
    """
    hfile.write("<h1> Contamination Check </h1> <br> </br>  ")
    contamination_image =sample+".r_1_screen.png"
    hfile.write("<center>")
    contamination_plots = "<img src='" + contamination_image+ "' alt='Mapping to Multiple Genomes' class='cover' width='75%' height='75%' />"
    hfile.write(contamination_plots)
    hfile.write("</center>")

    hfile.write("<h1> Alignments </h1>")
    alignments_txt= sample+".alignment_metrics.txt "
    alignments_link = "<a href='"+alignments_txt+"'> Alignments Metrics detailed here</a>"
    hfile.write("<h2> "+alignments_link+"</h2>")
    
    hfile.write("<h2> Alignments summary </h2> <br> </br>  ")
    alignments_image =sample+".alignments.png"
    alignments_plots = "<img src='" + alignments_image+ "' alt='Alignments' width='100% height='100%'/>"
    hfile.write(alignments_plots)

    hfile.write("<h2> Mapping Quality </h2> <br> </br>  ")
    mquality_image =sample+".mapq.png"
    mquality_plots = "<img src='" + mquality_image+ "' alt='Mapping Quality' width='100% height='100%'/>"
    hfile.write(mquality_plots)

    hfile.write("<h1> Coverage</h1> ")
    coverage_txt= sample+".coverage.histogram.txt"
    coverage_plots = "<iframe src='" +coverage_txt+"#scrollbar=0&toolbar=0&scrolling=0' frameBorder='0' height='600' width='600'> </iframe> " 
    hfile.write(coverage_plots)
    hfile.write("</center>")

         
    hfile.write("<h1> GC Bias</h1> ")
    gcbias_txt= sample+".gc_bias_metrics.txt "
    gcbias_link = "<a href='"+gcbias_txt+"'> GC Bias Metrics detailed here</a>"
    hfile.write("<h2>"+gcbias_link+"</h2>")   
    hfile.write("<h2> GC Bias plot </h2>")
    gcbias_image = sample+".gc_bias_metrics.pdf"
    hfile.write("<center>")
    gcbias_plot ="<iframe src='" +gcbias_image+"#scrollbar=0&toolbar=0&scrolling=0' scrolling='no' seamless='seamless'frameBorder='0' height='600' width='600'> </iframe> "
    hfile.write(gcbias_plot)
    hfile.write("</center>")

    hfile.write("<h1> Insert Size</h1> ")
    insertsize_txt= sample+".insert_size_metrics.txt"
    insertsize_link = "<a href='"+insertsize_txt+"'> Insert Size Metrics detailed here</a>"
    hfile.write("<h2> "+insertsize_link+"</h2>")
    hfile.write("<h2> Insert Size </h2>")
    insertsize_image = sample+".insert_size_histogram.pdf"
    hfile.write("<center>")
    insertsize_plot = "<iframe src='" +insertsize_image+"#toolbar=0&scrollbar=0&scrolling=0' scrolling='no' seamless='seamless' frameBorder='0' height='600' width='600'> </iframe> "
    hfile.write(insertsize_plot)
    hfile.write("</center>")

    hfile.write("</html>") 
    hfile.write(message)
    hfile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sample', default=False)
    args = parser.parse_args()
    html_report(args.sample) 
if __name__ == '__main__':
    main()




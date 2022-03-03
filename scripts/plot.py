#e! /usr/bin/env python
import argparse
import os
import sys
import math
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

def plot_stats(sample):
   stats = []  
   stats_file = sample+".alignments.txt"
   outfile = sample +".alignments.png"
   infile = open(stats_file, 'r')
   primary_reads = int(infile.readline()) 
   secondary_reads = int(infile.readline()) 
   unmapped_reads = int(infile.readline())
   total_reads = int(infile.readline())
   primary_reads = round(primary_reads / total_reads * 100, 2) 
   secondary_reads = round(secondary_reads / total_reads * 100, 2) 
   unmapped_reads = round(unmapped_reads / total_reads * 100, 2) 

   fig = go.Figure(data=[go.Table(
   header=dict(values=['Reads', 'Percentages'],
                line_color='grey',
                fill_color='cornflowerblue',
                align='left', font=dict(color='black', family="Open Sans", size=14)),
   cells=dict(values=[['Primary Reads', 'Secondary Reads', 'Unmapped Reads'], # 1st column
                       [str(primary_reads)+" %", str(secondary_reads)+" %", str(unmapped_reads)+" %"]], # 2nd column
               line_color='grey',
               fill_color='whitesmoke',
               align='left',font=dict(color='black', family="Open Sans", size=12)) ) ])
   fig.update_layout(height=200, margin=dict(t=5, b=5))
   fig.write_image(outfile) 
   infile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sample', default=False) 
    args = parser.parse_args()
    plot_stats(args.sample) 
if __name__ == '__main__':
    main()




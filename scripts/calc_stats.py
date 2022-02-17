#! /usr/bin/env python
import argparse
import os
import sys
import math
import matplotlib.pyplot as plt

def plot_stats(stats_file, outfile):
   stats = []  
   infile = open(stats_file, 'r')
   primary_reads = int(infile.readline())
   secondary_reads = int(infile.readline()) 
   unmapped_reads = int(infile.readline())
   total_reads = int(infile.readline())
   qbasel5 = int(infile.readline()) 
   qbaseg20 = int(infile.readline())
   infile.close()
   stats.append(primary_reads)
   stats.append(secondary_reads)
   stats.append(unmapped_reads)  
   stats.append(qbasel5)
   stats.append(qbaseg20)
   categories = ['Primary', 'Secondary', 'Unmapped reads', 'Qual <5','Qual >=20'] 
   for i in range(0, len(stats)):
       stats[i] = round(stats[i] / total_reads *100, 2)
   plt.bar(categories, stats, 0.75, color="blue")
   plt.title('Percentage of Mapped/Unmapped Reads', fontsize=10)
   plt.ylabel("Percentage to total Reads", fontsize=8)
   plt.xticks(fontsize=8)
   plt.ylim(0, 100)
   i = 1.0
   j = 2000 
   for i in range(len(categories)):
       plt.annotate(stats[i], (-0.1 + i, stats[i] + j)) 
   plt.savefig(outfile, format='png')
   plt.clf()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', default=False) 
    parser.add_argument('outfile', default=False)
    args = parser.parse_args()
    plot_stats(args.infile, args.outfile) 
if __name__ == '__main__':
    main()




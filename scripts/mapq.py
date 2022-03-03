#!/usr/bin/env python

import pysam
import seaborn as sns
import argparse
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description="Mapping Quality Density Plot")

parser.add_argument('--bam', required=True, help='path to input bam file', action='store')
parser.add_argument('--png', required=True, help='path to output png file', action='store')

args=parser.parse_args()

output = args.png
input = args.bam 

bam = pysam.AlignmentFile(input, "rb")

quals = [read.mapping_quality for read in bam.fetch()]

sns.set_style("darkgrid")
sns.set_context("paper")
plt.figure(figsize=(7,4))

ax = sns.kdeplot(quals, shade = True, color ="blue")
ax.set_ylabel('Density')
ax.set_xlabel('Mapping Quality Score')

plt.savefig(output)

#!/usr/bin/python3

import sys
import os
import time
import argparse
import pandas as pd
#----------------------------------------------------------------------------------------#
parser = argparse.ArgumentParser(description='Pipeline Usage')
parser.add_argument('1', metavar='<38 or 19>' ,help='Select Reference version')
parser.add_argument('2', metavar='<Core>' ,help='Set Core')
parser.add_argument('3', metavar='<Step>' ,help='All, Align, Annotation, SV, etc')
args = parser.parse_args()
#----------------------------------------------------------------------------------------#
Sample = pd.read_csv('SampleSheet.txt', sep='\t', header=None)
Name = Sample.iloc[0,0]
R1 = Sample.iloc[0,1]
R2 = Sample.iloc[0,2]
#----------------------------------------------------------------------------------------#
def SV(Name):
    command = f'delly call \
                -g /media/src/hg{sys.argv[1]}/hg{sys.argv[1]}.fa \
                -o {Name}.sv.bcf \
                -x /Bioinformatics/00.Tools/delly/excludeTemplates/human.hg{sys.argv[1]}.excl.tsv \
                {Name}.Real.bam'
    os.system(command)
#----------------------------------------------------------------------------------------#
if sys.argv[3] == 'All':
    SV(Name)
elif sys.argv[3] == 'FastQC':
    pass
elif sys.argv[3] == 'Align':
    pass
elif sys.argv[3] == 'SV':
    SV(Name)
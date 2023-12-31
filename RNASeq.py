#!/usr/bin/python3
#231101.ver

import os
import time
import argparse
import sys
import smtplib
import math
import pandas as pd
import numpy as np
import glob
from collections import defaultdict
from fpdf import FPDF
from datetime import datetime

Start_time = time.time()
dt = datetime.now()

parser = argparse.ArgumentParser(description='Pipeline Usage')
parser.add_argument('1', metavar='<38 or 19>' ,help='Select Reference version')
parser.add_argument('2', metavar='<Core>' ,help='Set Core')
parser.add_argument('3', metavar='<Step>' ,help='All, STAR, Mutation, QC')
args = parser.parse_args()
#----------------------------------------------------------------------------------------#
Sample = pd.read_csv('SampleSheet.txt', sep='\t', header=None)
Name = Sample.iloc[0,0]
R1 = Sample.iloc[0,1]
R2 = Sample.iloc[0,2]
print(R1, R2)
#----------------------------------------------------------------------------------------#
def PreQC(r1, r2):
    if os.path.isdir('00.PreQC'):
        pass
    else:
        command = 'mkdir 00.PreQC'
        os.system(command)

    command =f'fastqc -o 00.PreQC\
            -t {int(sys.argv[2])*2}\
            {r1}\
            {r2}'
    os.system(command)
#----------------------------------------------------------------------------------------#
def Trimming(name, r1, r2):
    if os.path.isdir('02.Trimmed'):
        pass
    else:
        command = 'mkdir 02.Trimmed'
        os.system(command)

    command = f'trim_galore --paired --gzip \
                -j {int(sys.argv[2])*2} \
                -o 02.Trimmed --basename {name} \
                {r1} {r2}'
    os.system(command)
#----------------------------------------------------------------------------------------#
def PostQC(name):
    if os.path.isdir('01.PostQC'):
        pass
    else:
        command = 'mkdir 01.PostQC'
        os.system(command)

    command = f'fastqc -o 01.PostQC \
                -t {int(sys.argv[2])*2} \
                02.Trimmed/{name}_val_1.fq.gz \
                02.Trimmed/{name}_val_2.fq.gz'
    os.system(command)
#----------------------------------------------------------------------------------------#
def Refindex():
    if os.path.isdir(f'/media/src/hg{sys.argv[1]}/00.RNA/Index/'):
        pass
    else:
        command = f'STAR --runThreadN {sys.argv[2]} --runMode genomeGenerate \
                    --genomeDir /media/src/hg{sys.argv[1]}/00.RNA/Index/ \
                    --sjdbOverhang 100 \
                    --sjdbGTFfile /media/src/hg{sys.argv[1]}/00.RNA/hg{sys.argv[1]}.GENCODE.v44.annotation.gtf \
                    --genomeFastaFiles /media/src/hg{sys.argv[1]}/hg{sys.argv[1]}.GENCODE.fa'
        os.system(command)
#----------------------------------------------------------------------------------------#
def STAR(name):
    command = f'STAR --runThreadN {int(sys.argv[2])*2} --genomeDir /media/src/hg{sys.argv[1]}/00.RNA/Index/ \
                --readFilesIn 02.Trimmed/{name}_val_1.fq.gz 02.Trimmed/{name}_val_2.fq.gz --readFilesCommand zcat \
                --outSAMtype BAM Unsorted \
                --twopassMode Basic \
                --quantMode GeneCounts \
                --outFilterMultimapNmax 1 --outFilterMismatchNmax 10 \
                --outFileNamePrefix 03.Output/{name}_'
    os.system(command)

    command = f'cp 03.Output/{name}_ReadsPerGene.out.tab ../Results/00.Reads/'
    os.system(command)
#----------------------------------------------------------------------------------------#
def QC(name, r1, r2):
    if os.path.isdir('04.QC/'):
        pass
    else:
        command = 'mkdir 04.QC'
        os.system(command)

    command = f'samtools sort 03.Output/{name}_Aligned.out.bam -o 03.Output/{name}_Sorted.out.bam'
    os.system(command)

    command = f'samtools index 03.Output/{name}_Sorted.out.bam'
    os.system(command)

    command = f'mv 03.Output/{name}_Log.final.out 04.QC/'
    os.system(command)

    command = f'samtools stats 03.Output/{name}_Aligned.out.bam > 04.QC/{name}.stats'
    os.system(command)

    command = f'samtools stats -t /Bioinformatics/01.Reference/Panel/rna.gene.bed 03.Output/{name}_Sorted.out.bam > 04.QC/{name}.Ontarget.stats'
    os.system(command)

    command = f'samtools flagstat 03.Output/{name}_Sorted.out.bam > 04.QC/{name}.Depth.txt'
    os.system(command)

    command = f'samtools depth 03.Output/{name}_Sorted.out.bam > 04.QC/{name}.Total.Depth.txt'
    os.system(command)

    command = f'samtools depth -b /Bioinformatics/01.Reference/Panel/BCR-ABL1.bed \
                03.Output/{name}_Sorted.out.bam > 04.QC/{name}.Target.Depth.txt'
    os.system(command)

    command = f'FxTools_Linux Fqtools fqcheck -i {r1} {r2} -o 04.QC/{name}.R1 04.QC/{name}.R2'
    os.system(command)
#----------------------------------------------------------------------------------------#
def QCPDF(name):
    if os.path.isdir(f'00.PreQC/{name}_R1_001_fastqc'):
        pass
    else:
        command = f'unzip 00.PreQC/{name}_R1_001_fastqc.zip \
                    -d 00.PreQC/ '
        os.system(command)

    if os.path.isdir(f'01.PostQC/{name}_val_1_fastqc'):
        pass
    else:
        command = f'unzip 01.PostQC/{name}_val_1_fastqc.zip \
                    -d 01.PostQC/'
        os.system(command)

    pdf = FPDF()
        
# Adding a page
    pdf.add_page()
        
# set style and size of font 
    pdf.set_font("helvetica", size = 15)
        
# create a cell
    pdf.set_fill_color(r= 100, g= 100, b= 100 )
    pdf.line(10,20,200,20)
    pdf.text(20, 11, txt = 'Imatinib Resistance')

#Sample Name
    pdf.set_font("helvetica", style = 'B', size = 10)
    pdf.text(20, 17, txt = f"Sample : {name}")

#Date
    pdf.set_font("helvetica", style = 'B', size = 10)
    pdf.text(170, 17, txt = dt.strftime('%Y-%m-%d'))

#Category 1
    pdf.set_font("helvetica", style = 'B', size = 12)
    pdf.text(20, 27, txt = '1. FASTQ Statistics - PreQC')

#Category1 Image
    pdf.set_font("helvetica", size = 9)

    pdf.text(26, 35, txt = "Per base quality : Untrimmed")
    pdf.image(f'00.PreQC/{name}_R1_001_fastqc/Images/per_base_quality.png', x = 20, y = 36, w = 50, h = 50, type = 'PNG')

    pdf.text(96, 35, txt = "Sequence length")
    pdf.image(f'00.PreQC/{name}_R1_001_fastqc/Images/per_sequence_quality.png', x = 80, y = 36, w = 50, h = 50, type = 'PNG')

    pdf.text(147, 35, txt = "Per base sequence content")
    pdf.image(f'00.PreQC/{name}_R1_001_fastqc/Images/per_base_sequence_content.png', x = 140, y = 36, w = 50, h = 50, type = 'PNG')

#Category 2
    pdf.set_font("helvetica",style = 'B', size = 12)
    pdf.text(20, 100, txt = '2. FASTQ Statistics - PostQC')

#Category2 Image
    pdf.set_font("helvetica", size = 9)

    pdf.text(26, 108, txt = "Per base quality : Trimmed")
    pdf.image(f'01.PostQC/{name}_val_1_fastqc/Images/per_base_quality.png', x = 20, y = 110, w = 50, h = 50, type = 'PNG')

    pdf.text(96, 108, txt = "Sequence length")
    pdf.image(f'01.PostQC/{name}_val_1_fastqc/Images/per_sequence_quality.png', x = 80, y = 110, w = 50, h = 50, type = 'PNG')

    pdf.text(147, 108, txt = "Per base sequence content")
    pdf.image(f'01.PostQC/{name}_val_1_fastqc/Images/per_base_sequence_content.png', x = 140, y = 110, w = 50, h = 50, type = 'PNG')

#Category 3 
    pdf.set_font("helvetica",style = 'B', size = 12)
    pdf.text(20, 170, txt = '3. Sequencing Statistics')

    Info = {}
    with open(f'04.QC/{name}.stats', 'r') as handle:
        for line in handle:
            if line.startswith('SN'):
                line = line.strip()
                splitted = line.split('\t')
                Info[splitted[1][:-1]] = splitted[2]
                
    Total_base = int(Info['bases mapped'])
    Total_read = int(Info['raw total sequences'])
    Mapped_read = int(Info['reads mapped'])
    Average_length = Info['average length']
    Average_qual = Info['average quality']
    Average_insert = Info['insert size average']

    Target_Info = {}
    with open(f'04.QC/{name}.Ontarget.stats', 'r') as handle:
        for line in handle:
            if line.startswith('SN'):
                line = line.strip()
                splitted = line.split('\t')
                Target_Info[splitted[1][:-1]] = splitted[2]
                
    Target_Mapped_read = int(Target_Info['reads mapped'])
    
    Percent = str(round(int(Mapped_read)/int(Total_read)*100))
    Ontarget = str(round(int(Target_Mapped_read)/int(Mapped_read)*100))

    with open(f'04.QC/{name}.Target.Depth.txt', 'r') as handle:
        Num = 0
        Total_Depth = 0
        for line in handle:
            line = line.strip()
            splitted = line.split('\t')
            Depth = int(splitted[2])
            Num += 1
            Total_Depth += Depth

        Total_AVG_Depth = str(round(Total_Depth/Num, 2))

    #R1_Phred_Score = os.popen(f"tail -n 1 04.QC/{name}.R1.fqcheck").read().strip().split('\t')[3]
    #R2_Phred_Score = os.popen(f"tail -n 1 04.QC/{name}.R2.fqcheck").read().strip().split('\t')[3]

    pdf.set_font("helvetica", size = 11)
    pdf.set_xy(20, 175)
    pdf.set_fill_color(r = 150, g = 150, b = 150)
    pdf.cell(57,10, txt = 'Total base', align = 'C', border=1, ln=0, fill = True)
    pdf.cell(57,10, txt = 'Total read', align = 'C', border=1, fill = True)
    pdf.cell(57,10, txt = 'Mapped read', align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(20, 185)
    pdf.cell(57,10, txt = str(format(Total_base, ',')), align = 'C', border=1)
    pdf.cell(57,10, txt = str(format(Total_read, ',')), align = 'C', border=1)
    pdf.cell(57,10, txt = "{:,} ({} %)".format(Mapped_read, Percent), align = 'C', border=1)

    pdf.set_xy(20, 195)
    pdf.set_fill_color(r = 150, g = 150, b = 150)
    pdf.cell(57,10, txt = 'Average Length', align = 'C', border=1, ln=0, fill = True)
    pdf.cell(57,10, txt = 'Average Insert size', align = 'C', border=1, fill = True)
    pdf.cell(57,10, txt = 'Average Quality', align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(20, 205)
    pdf.cell(57,10, txt = Average_length, align = 'C', border=1)
    pdf.cell(57,10, txt = Average_insert, align = 'C', border=1)
    pdf.cell(57,10, txt = Average_qual, align = 'C', border=1)

    pdf.set_xy(20, 215)
    pdf.set_fill_color(r = 150, g = 150, b = 150)
    pdf.cell(85.5,10, txt = 'Average Depth', align = 'C', border=1, ln=0, fill = True)
    pdf.cell(85.5,10, txt = 'Ontarget %', align = 'C', border=1, fill = True)
    #pdf.cell(57,10, txt = 'Phred score(Q>=30, R1 | R2)', align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(20, 225)
    pdf.cell(85.5,10, txt = Total_AVG_Depth, align = 'C', border=1)
    pdf.cell(85.5,10, txt = Ontarget + ' %', align = 'C', border=1)
    #pdf.cell(57,10, txt = R1_Phred_Score + ' | ' + R2_Phred_Score, align = 'C', border=1)

#Category 4 
    pdf.set_font("helvetica", style = 'B', size = 12)
    pdf.text(20, 250, txt = '4. STAR Alignment Statistics')

    Info = {}
    with open(f'04.QC/{name}_Log.final.out', 'r') as handle:
        for line in handle:
            line = line.strip()
            if '%' in line:
                line = line.replace('\t', '')
                splitted = line.split('|')
                Info[splitted[0]] = splitted[1]

    pdf.set_font("helvetica", size = 11)
    pdf.set_xy(20, 255)
    pdf.set_fill_color(r = 150, g = 150, b = 150)
    pdf.cell(57,10, txt = list(Info.keys())[0], align = 'C', border=1, ln=0, fill = True)
    pdf.cell(57,10, txt = list(Info.keys())[1].replace(',',''), align = 'C', border=1, ln=0, fill = True)
    pdf.cell(57,10, txt = list(Info.keys())[-1], align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(20, 265)
    pdf.cell(57,10, txt = Info[list(Info.keys())[0]], align = 'C', border=1)
    pdf.cell(57,10, txt = Info[list(Info.keys())[1]], align = 'C', border=1)
    pdf.cell(57,10, txt = Info[list(Info.keys())[-1]], align = 'C', border=1)

# save the pdf
    pdf.output(f"03.Output/{name}_QC.pdf")
    
    command = f'ln 03.Output/{name}_QC.pdf ../Results/{name}_QC.pdf'
    os.system(command)
#----------------------------------------------------------------------------------------#
#RUN Pipeline
if sys.argv[3] == 'All':
    PreQC(R1, R2)
    Trimming(Name, R1, R2)
    PostQC(Name)
    Refindex()
    STAR(Name)
    QC(Name, R1, R2)
    # QCPDF(Name)
    # Fusion(Name)
elif sys.argv[3] == 'FastQC':
    PreQC(R1, R2)
    Trimming(Name, R1, R2)
    PostQC(Name)
elif sys.argv[3] == 'Align':
    Trimming(Name, R1, R2)
    Refindex()
    STAR(Name)
elif sys.argv[3] == 'Fusion':
    Trimming(Name, R1, R2)
    STAR(Name)
    # Fusion(Name)
elif sys.argv[3] == 'Indexing':
    Refindex()

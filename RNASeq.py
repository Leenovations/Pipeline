#!/home/lab/anaconda3/envs/NGS/bin/python3
#240117.ver

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
from functools import reduce
#----------------------------------------------------------------------------------------#
Start_time = time.time()
dt = datetime.now()
#----------------------------------------------------------------------------------------#
parser = argparse.ArgumentParser(description='Pipeline Usage')
args = parser.parse_args()
#----------------------------------------------------------------------------------------#
Sample = pd.read_csv('SampleSheet.txt', sep='\t', header=None)
Name = Sample.iloc[0,0]
R1 = Sample.iloc[0,1]
R2 = Sample.iloc[0,2]
#----------------------------------------------------------------------------------------#
BATCH = {}
with open(f"{Name}.batch.config", "r") as batch:
    for line in batch:
        line = line.strip()
        splitted = line.split("=")
        Key = splitted[0]
        Value = splitted[1]
        BATCH[Key] = Value
#----------------------------------------------------------------------------------------#
def PreQC(r1, r2):
    if os.path.isdir('00.PreQC'):
        pass
    else:
        command = 'mkdir 00.PreQC'
        os.system(command)

    command =f'fastqc -o 00.PreQC \
            -t 2 \
            {r1} \
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
                -j {BATCH["CPU"]} \
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
                -t 2 \
                02.Trimmed/{name}_val_1.fq.gz \
                02.Trimmed/{name}_val_2.fq.gz'
    os.system(command)
#----------------------------------------------------------------------------------------#
def Refindex():
    if os.path.isdir(f'/media/src/hg{BATCH["Ref.ver"].split("g")[1]}/00.RNA/Index/'):
        pass
    else:
        command = f'STAR --runThreadN {BATCH["CPU"]} --runMode genomeGenerate \
                    --genomeDir /media/src/hg{BATCH["Ref.ver"].split("g")[1]}/00.RNA/Index/ \
                    --sjdbOverhang {BATCH["sjdbOverhang"]} \
                    --sjdbGTFfile /media/src/hg{BATCH["Ref.ver"].split("g")[1]}/00.RNA/hg{BATCH["Ref.ver"].split("g")[1]}.GENCODE.v44.annotation.gtf \
                    --genomeFastaFiles /media/src/hg{BATCH["Ref.ver"].split("g")[1]}/hg{BATCH["Ref.ver"].split("g")[1]}.GENCODE.fa'
        os.system(command)
#----------------------------------------------------------------------------------------#
def STAR(name):
    command = f'STAR --runThreadN {BATCH["CPU"]} --genomeDir /media/src/hg{BATCH["Ref.ver"].split("g")[1]}/00.RNA/Index/ \
                --readFilesIn 02.Trimmed/{name}_val_1.fq.gz 02.Trimmed/{name}_val_2.fq.gz --readFilesCommand zcat \
                --outSAMtype BAM Unsorted \
                --outFilterScoreMin {BATCH["FilterScoreMin"]} \
                --twopassMode {BATCH["twopassMode"]} \
                --quantMode {BATCH["quantMode"]} \
                --outFilterMultimapNmax {BATCH["FilterMultimapNmax"]} --outFilterMismatchNmax {BATCH["FilterMismatchNmax"]} \
                --outFileNamePrefix 03.Output/{name}_'
    os.system(command)

    Gene = pd.read_csv(f"/media/src/hg19/00.RNA/Index/geneInfo.tab", sep='\t', header=None)
    column_dict = {}
    for index, row in Gene.iterrows():
        column_dict[row[0]] = row[1]

    Genecount = pd.read_csv(f"03.Output/{name}_ReadsPerGene.out.tab", sep='\t', header=None, names=['ID', 0, 1, 2])
    Genecount = Genecount.drop(index=range(0, 4))
    Genecount['ID'] = [column_dict[key] for key in Genecount['ID']]
    Genecount = Genecount[['ID', int(BATCH["Stranded"])]]
    Genecount.columns = ['GeneSymbol', f"{name}"]
    Genecount.to_csv(f"03.Output/{name}.Genecount.txt", sep='\t', header='infer', index=False)

    names = BATCH['Sample.Name'].split(',')
    Dir = BATCH['Sample.Dir'].split(',')
    file_paths = [file + f"03.Output/{sample}.Genecount.txt" for file, sample in zip(Dir, names)]

    data_frames = [pd.read_csv(file, 
                            sep='\t', 
                            header='infer')
                            for file, sample in zip(file_paths, names)]
    DATA = reduce(lambda left, right: pd.merge(left, right, on=['GeneSymbol']), data_frames)
    DATA.to_csv(f"../Genecount/Total.Genecount.txt", sep='\t', header='infer', index=False)

    if BATCH["TPM"] == "Y":
        Length = pd.read_csv('/media/src/hg19/00.RNA/hg19.GENCODE.v44.GeneLength.txt',
                        sep='\t',
                        header='infer',
                        names =['ID', 'Type', 'GeneSymbol', 'Length'])
        LENGTH = dict(zip(Length['ID'].to_list(), Length['Length'].to_list()))
        GENE = dict(zip(Length['ID'].to_list(), Length['GeneSymbol'].to_list()))
        #--------------------------------------------------------------------------------#
        GeneCount = pd.read_csv('../Genecount/Total.Genecount.txt',
                                sep='\t',
                                header='infer')
        GeneCount.index = GeneCount['GeneSymbol']
        GeneCount = GeneCount.rename_axis(None)
        #--------------------------------------------------------------------------------#
        GeneCount['Length'] = GeneCount['GeneSymbol'].map(LENGTH)
        GeneCount = GeneCount.drop(columns=[GeneCount.columns[0]])
        Count = GeneCount.iloc[:, :GeneCount.columns.get_loc('Length')]
        Sum_Depth_Length = Count.sum(axis=0)
        #--------------------------------------------------------------------------------#
        Exon_length = np.array(GeneCount['Length'])
        Sum_Depth_Length = np.array(Sum_Depth_Length)
        result_matrix = np.outer(Exon_length, Sum_Depth_Length)
        NormFactor = pd.DataFrame(result_matrix)
        NormFactor.columns = list(Count.columns)
        NormFactor.index = GeneCount.index.to_list()
        TPM = Count.div(NormFactor)
        #--------------------------------------------------------------------------------#
        TPM['GeneSymbol'] = TPM.index
        TPM['GeneSymbol'] = TPM.index.map(GENE)
        last_column = TPM.pop('GeneSymbol')
        TPM.insert(0, 'GeneSymbol', last_column)
    if BATCH["FPKM"] == "Y":
        pass
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

    command = f'samtools stats 03.Output/{name}_Sorted.out.bam > 04.QC/{name}.stats'
    os.system(command)

    command = f'samtools stats -t /media/src/hg{BATCH["Ref.ver"].split("g")[1]}/04.cnv/NCBI.RefSeq.Selected.Gene.bed \
                03.Output/{name}_Sorted.out.bam > 04.QC/{name}.Ontarget.stats'
    os.system(command)

    command = f'samtools flagstat 03.Output/{name}_Sorted.out.bam > 04.QC/{name}.Depth.txt'
    os.system(command)

    command = f'samtools depth 03.Output/{name}_Sorted.out.bam > 04.QC/{name}.Total.Depth.txt'
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
    pdf.cell(171,10, txt = 'Ontarget %', align = 'C', border=1, fill = True)
    pdf.set_xy(20, 225)
    pdf.cell(171,10, txt = Ontarget + ' %', align = 'C', border=1)

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
    pdf.output(f"04.QC/{name}_QC.pdf")
#----------------------------------------------------------------------------------------#
if BATCH["Step"] == 'All':
    PreQC(R1, R2)
    Trimming(Name, R1, R2)
    PostQC(Name)
    Refindex()
    STAR(Name)
    QC(Name, R1, R2)
    QCPDF(Name)
    # Fusion(Name)
elif BATCH["Step"] == 'FastQC':
    PreQC(R1, R2)
    Trimming(Name, R1, R2)
    PostQC(Name)
elif BATCH["Step"] == 'Align':
    Trimming(Name, R1, R2)
    Refindex()
    STAR(Name)
elif BATCH["Step"] == 'Fusion':
    Trimming(Name, R1, R2)
    STAR(Name)
    # Fusion(Name)
elif BATCH["Step"] == 'Indexing':
    Refindex()

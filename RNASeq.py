#!/home/lab/anaconda3/envs/NGS/bin/python3
#240203.ver

import os
import time
import argparse
import pandas as pd
import numpy as np
from RNASeqNorm import *
from RNASeqQC import *
from datetime import datetime
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

    command = f'/media/src/Tools/TrimGalore-0.6.10/trim_galore \
                --paired --gzip \
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
    command = f'STAR --runThreadN {BATCH["CPU"]} \
                --genomeDir /media/src/hg{BATCH["Ref.ver"].split("g")[1]}/00.RNA/Index/ \
                --sjdbGTFfile /media/src/hg{BATCH["Ref.ver"].split("g")[1]}/00.RNA/hg{BATCH["Ref.ver"].split("g")[1]}.GENCODE.v44.annotation.gtf \
                --readFilesCommand zcat \
                --readFilesIn 02.Trimmed/{name}_val_1.fq.gz 02.Trimmed/{name}_val_2.fq.gz \
                --outSAMtype {BATCH["outSAMtype"]} \
                --twopassMode {BATCH["twopassMode"]} \
                --quantMode {BATCH["quantMode"]} \
                --outFilterMultimapNmax {BATCH["FilterMultimapNmax"]} \
                --outFilterMismatchNmax {BATCH["FilterMismatchNmax"]} \
                --winAnchorMultimapNmax {BATCH["AnchorMultimapNmax"]} \
                --outFilterScoreMin {BATCH["FilterScoreMin"]} \
                --outFileNamePrefix 03.Output/{name}_'
    os.system(command)
#----------------------------------------------------------------------------------------#
def Expression(name):
    ReadCount(name, f'{BATCH["Ref.ver"].split("g")[1]}', f'{BATCH["Stranded"]}')

    if BATCH["TPM"] == "Y":
        TPM(name)
    if BATCH["FPKM"] == "Y":
        FPKM(name)

    Merge(name,  BATCH['Sample.Name'].split(','), BATCH['Sample.Dir'].split(','), BATCH["TPM"], BATCH["FPKM"])
#----------------------------------------------------------------------------------------#
def QC(name):
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
def Results(name):
    QCPDF(name)
    pdfconverter(name)
#----------------------------------------------------------------------------------------#
if BATCH["Step"] == 'All':
    PreQC(R1, R2)
    Trimming(Name, R1, R2)
    PostQC(Name)
    Refindex()
    STAR(Name)
    Expression(Name)
    QC(Name)
    Results(Name)
elif BATCH["Step"] == 'FastQC':
    PreQC(R1, R2)
    Trimming(Name, R1, R2)
    PostQC(Name)
elif BATCH["Step"] == 'Trimming':
    Trimming(Name, R1, R2)
elif BATCH["Step"] == 'Align':
    Refindex()
    STAR(Name)
elif BATCH["Step"] == 'Expression':
    Expression(Name)
elif BATCH["Step"] == 'QC':
    QC(Name)
elif BATCH["Step"] == 'Fusion':
    Trimming(Name, R1, R2)
    STAR(Name)
elif BATCH["Step"] == 'Results':
    Results(Name)
elif BATCH["Step"] == 'Indexing':
    Refindex()
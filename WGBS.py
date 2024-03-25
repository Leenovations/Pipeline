#!/home/lab/anaconda3/envs/NGS/bin/python3

import os
import pandas as pd
import numpy as np
import argparse
from DepthofCoverage import *
from WGBSQC import *
#-----------------------------------------------------------------------------#
parser = argparse.ArgumentParser(description="Pipeline Usage")
args = parser.parse_args()
#-----------------------------------------------------------------------------#
Sample = pd.read_csv("SampleSheet.txt", sep="\t", header=None)
Name = Sample.iloc[0, 0]
R1 = Sample.iloc[0, 1]
R2 = Sample.iloc[0, 2]
#-----------------------------------------------------------------------------#
BATCH = {}
with open(f"{Name}.batch.config", "r") as batch:
    for line in batch:
        line = line.strip()
        splitted = line.split("=")
        Key = splitted[0]
        Value = splitted[1]
        BATCH[Key] = Value
#-----------------------------------------------------------------------------#
if BATCH["Bismark"] == "Y":
    def PreQC(r1, r2):
        if os.path.isdir("00.PreQC"):
            pass
        else:
            command = "mkdir 00.PreQC"
            os.system(command)

        command =f"fastqc -o 00.PreQC \
                -t 2 \
                {r1}\
                {r2}"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def Trimming(name, r1, r2):
        if os.path.isdir("02.Trimmed"):
            pass
        else:
            command = "mkdir 02.Trimmed"
            os.system(command)

        command = f"trim_galore \
                    --paired --gzip \
                    -j {BATCH['CPU']} \
                    -o 02.Trimmed --basename {name} \
                    {r1} {r2}"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def PostQC(name):
        if os.path.isdir("01.PostQC"):
            pass
        else:
            command = "mkdir 01.PostQC"
            os.system(command)

        command = f"fastqc -o 01.PostQC \
                    -t 2 \
                    02.Trimmed/{name}_val_1.fq.gz \
                    02.Trimmed/{name}_val_2.fq.gz"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def Index():
        if os.path.isdir(f"/Bioinformatics/01.Reference/{BATCH['Ref.ver']}/Methylation/Bisulfite_Genome/"):
            pass
        else:
            command = f"bismark_genome_preparation --path_to_aligner \
                        /Bioinformatics/00.Tools/bowtie2-2.4.5 \
                        --parallel {BATCH['CPU']} \
                        --verbose \
                        /Bioinformatics/01.Reference/{BATCH['Ref.ver']}/Methylation/"
            os.system(command)
#----------------------------------------------------------------------------------------#
    def Align(name):
        if os.path.isdir(f"03.Align/"):
            pass
        else:
            command = f"mkdir -p 03.Align/"
            os.system(command)

        command = f"bismark \
                    --multicore {BATCH['CPU']} --un --ambiguous --gzip --nucleotide_coverage \
                    -N {BATCH['AllowMismatch']} -L {BATCH['SeedLength']} \
                    --temp_dir TEMP \
                    -o 03.Align/ \
                    --genome /Bioinformatics/01.Reference/{BATCH['Ref.ver']}/Methylation/ \
                    -1 02.Trimmed/{name}_val_1.fq.gz -2 02.Trimmed/{name}_val_2.fq.gz"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def Dedup(name):
        command = f"deduplicate_bismark -p \
                    --output_dir 03.Align/ \
                    -o 03.Align/{name} \
                    03.Align/{name}_val_1_bismark_bt2_pe.bam"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def Lambda(name):
        command = f"bismark --gzip -o 03.Align --genome /Bioinformatics/01.Reference/lambda \
                    -1 03.Align/{name}_val_1.fq.gz_unmapped_reads_1.fq.gz \
                    -2 03.Align/{name}_val_2.fq.gz_unmapped_reads_2.fq.gz"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def bamflt(name):
        if os.path.isdir(f"03.Align/"):
            pass
        else:
            command = f"mkdir -p 03.Align/"
            os.system(command)

        command = f"samtools view -f 2 -bq {BATCH['bq']} 03.Align/{name}.deduplicated.bam -o 03.Align/{name}.flt.bam"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def Extract(name):
        if os.path.isdir(f"03.Align/"):
            pass
        else:
            command = f"mkdir -p 03.Align/"
            os.system(command)
    
        command = f"bismark_methylation_extractor \
                    -p --no_overlap --bedGraph --gzip --multicore {BATCH['CPU']} --cytosine_report -zero_based \
                    --genome_folder /Bioinformatics/01.Reference/{BATCH['Ref.ver']}/Methylation \
                    --comprehensive --merge_non_CpG \
                    -o 03.Align \
                    03.Align/{name}.flt.bam"
        os.system(command)

        command = f"samtools sort -@ {int(BATCH['CPU'])} 03.Align/{name}.flt.bam -o 03.Align/{name}.sorted.bam"
        os.system(command)

        command = f"samtools index -@ {int(BATCH['CPU'])} 03.Align/{name}.sorted.bam"
        os.system(command)

        command = f'samtools stats 03.Align/{name}.sorted.bam > 03.Align/{name}.stats'
        os.system(command)

        command = f"samtools view -H 03.Align/{name}.sorted.bam > 03.Align/{name}.header.txt"
        os.system(command)

        command = f"samtools depth -a 03.Align/{name}.sorted.bam > 03.Align/{name}.all.depth.txt"
        os.system(command)

        DepthOfCoverage(name)

        command = f"samtools coverage 03.Align/{name}.sorted.bam -o 03.Align/{name}.coverage.txt"
        os.system(command)

        command = f"Rscript /labmed/00.Code/Pipeline/NGSQC.R"
        os.system(command)

        command = f"methylation_consistency 03.Align/{name}.flt.bam"
        os.system(command)

        QCPDF(name)
        pdfconverter(name)
#----------------------------------------------------------------------------------------#
    def ChromosomalCNV(name):
        if os.path.isdir("04.ChromosomeCNV"):
            pass
        else:
            command = "mkdir 04.ChromosomeCNV"
            os.system(command)

        command = f"samtools bedcov \
                    /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/04.cnv/1MB.exclude.centromere.bed \
                    03.Align/{name}.sorted.bam > 04.ChromosomeCNV/{name}.bedcov"
        os.system(command)

        Chromosome = [str(i) for i in range(1,23)] + ['X', 'Y']
        Data = pd.read_csv(f"04.ChromosomeCNV/{name}.bedcov",
                        sep='\t',
                        header=None,
                        names = ['Chr', 'Start', 'End', 'Count'],
                        low_memory=False)
        
        Data['Length'] = Data['End'] - Data['Start']
        Data['count_per_length'] = Data['Count'] / Data['Length']
        Data['Norm'] = Data['Count'] / (Data['Length'] * Data['count_per_length'].sum())
        Data['Norm'] = np.log10(Data['Norm'] + 1)
        Median_Norm = Data['Norm'].median()
        Data['Norm'] = Data['Norm'] - Median_Norm
        
        Sorted = []
        for chromosome in Chromosome:
            Data_sub = Data[Data['Chr'] == chromosome]
            Sorted.append(Data_sub)
        Sorted_Data = pd.concat(Sorted)
        Sorted_Data['Order'] = range(1, len(Sorted_Data) + 1)
        
        Sorted_Data.to_csv(f"04.ChromosomeCNV/{name}.Chromosome.CNV.txt",
                            sep='\t',
                            header='infer',
                            index=False)
        
        command = f"Rscript /labmed/00.Code/Pipeline/ChromosomalCNV.R {name}"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def HTML(name):
        command = f"cd 03.Align; bismark2report; bismark2summary"
        os.system(command)
        
        # command = f"bismark2report --output 03.Align/{name}.html \
        #             --alignment_report 03.Align/{name}_val_1_bismark_bt2_PE_report.txt \
        #             --dedup_report 03.Align/{name}.deduplication_report.txt \
        #             --splitting_report 03.Align/{name}.flt_splitting_report.txt \
        #             --mbias_report 03.Align/{name}.flt.M-bias.txt \
        #             --nucleotide_report 03.Align/{name}_val_1_bismark_bt2_pe.nucleotide_stats.txt"
        # os.system(command)

        # command = f"bismark2summary --basename 03.Align/{name}.summary.html --title {name}"
        # os.system(command)
#----------------------------------------------------------------------------------------#
    if BATCH["Step"] == "All":
        PreQC(R1, R2)
        Trimming(Name, R1, R2)
        PostQC(Name)
        Index()
        Align(Name)
        Dedup(Name)
        Lambda(Name)
        bamflt(Name)
        Extract(Name)
        ChromosomalCNV(Name)
        HTML(Name)    
    elif BATCH["Step"] == "FastQC":
        PreQC(R1, R2)
        Trimming(Name, R1, R2)
        PostQC(Name)
    elif BATCH["Step"] == "Trimming":
        Trimming(Name, R1, R2)
    elif BATCH["Step"] == "Align":
        Trimming(Name, R1, R2)
        Index()
        Align(Name)
    elif BATCH["Step"] == "Dedup":
        Dedup(Name)
    elif BATCH["Step"] == "Bamflt":
        bamflt(Name)
    elif BATCH["Step"] == "Extract":
        Extract(Name)
        HTML(Name)
    elif BATCH["Step"] == "ChromosomeCNV":
        ChromosomalCNV(Name)
#-----------------------------------------------------------------------------#
elif BATCH["LAST"] == "Y":
    def LAST():
        command = f"/media/src/Tools/Bisulfighter/bsf-call/bsf-call \
                    -p {BATCH['CPU']} \
                    -o {Name}.CpG.txt \
                    /media/src/{BATCH['Ref.ver']}/02.Fasta/{BATCH['Ref.ver']}.fa \
                    {R1},{R2}"
        os.system(command)
#-----------------------------------------------------------------------------#
    if BATCH["Step"] == "All":
        LAST()

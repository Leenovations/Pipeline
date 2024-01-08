#!/home/lab/anaconda3/envs/NGS/bin/python3

import os
import pandas as pd
import argparse
import sys
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
                -t {BATCH['CPU']} \
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

        command = f"trim_galore --paired --gzip \
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
                    -t {BATCH['CPU']} \
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
                        --parallel 20 \
                        --verbose \
                        /Bioinformatics/01.Reference/{BATCH['Ref.ver']}/Methylation/"
            os.system(command)
#----------------------------------------------------------------------------------------#
    def Align(name):
        if os.path.isdir(f"03.Output/"):
            pass
        else:
            command = f"mkdir -p 03.Output/"
            os.system(command)

        command = f"bismark \
                    --multicore {BATCH['CPU']} --un --ambiguous --gzip --nucleotide_coverage \
                    --temp_dir TEMP \
                    -o 03.Output/ \
                    --genome /Bioinformatics/01.Reference/{BATCH['Ref.ver']}/Methylation/ \
                    -1 02.Trimmed/{name}_val_1.fq.gz -2 02.Trimmed/{name}_val_2.fq.gz"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def Dedup(name):
        command = f"deduplicate_bismark -p \
                    --output_dir 03.Output/ \
                    -o 03.Output/{name} \
                    03.Output/{name}_val_1_bismark_bt2_pe.bam"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def Lambda(name):
        command = f"bismark --gzip -o 03.Output --genome /Bioinformatics/01.Reference/lambda \
                    -1 03.Output/{name}_val_1.fq.gz_unmapped_reads_1.fq.gz \
                    -2 03.Output/{name}_val_2.fq.gz_unmapped_reads_2.fq.gz"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def Extract(name):
        command = f"bismark_methylation_extractor \
                    -p --no_overlap --bedGraph --gzip --multicore 20 --cytosine_report \
                    --genome_folder /Bioinformatics/01.Reference/{BATCH['Ref.ver']}/Methylation \
                    --comprehensive --merge_non_CpG \
                    -o 03.Output \
                    03.Output/{name}.deduplicated.bam"
        os.system(command)

        command = f"gunzip 03.Output/{name}.deduplicated.bismark.cov.gz"
        os.system(command)
        command = f"ln 03.Output/{name}.deduplicated.bismark.cov ../CpG/{name}.deduplicated.bismark.cov"
        os.system(command)
#----------------------------------------------------------------------------------------#
    def HTML(name):
        command = f"bismark2report --output 03.Output/{name}.html \
                    --alignment_report 03.Output/{name}_val_1_bismark_bt2_PE_report.txt \
                    --dedup_report 03.Output/{name}_val_1_bismark_bt2_pe.deduplication_report.txt \
                    --splitting_report 03.Output/{name}.deduplicated_splitting_report.txt \
                    --mbias_report 03.Output/{name}.deduplicated.M-bias.txt \
                    --nucleotide_report 03.Output/{name}_val_1_bismark_bt2_pe.nucleotide_stats.txt"
        os.system(command)
#----------------------------------------------------------------------------------------#
    if BATCH["Step"] == "All":
        PreQC(R1, R2)
        Trimming(Name, R1, R2)
        PostQC(Name)
        Index()
        Align(Name)
        Dedup(Name)
        Lambda(Name)
        Extract(Name)
        HTML(Name)
    elif BATCH["Step"] == "FastQC":
        PreQC(R1, R2)
        Trimming(Name, R1, R2)
        PostQC(Name)
    elif BATCH["Step"] == "Align":
        Trimming(Name, R1, R2)
        Index()
        Align(Name)
    elif BATCH["Step"] == "Dedup":
        Dedup(Name)
#-----------------------------------------------------------------------------#
elif BATCH["LAST"] == "Y":
    def LAST():
        command = f"/media/src/Tools/Bisulfighter/bsf-call/bsf-call \
                    -o {Name}.CpG.txt \
                    /media/src/{BATCH['Ref.ver']}/02.Fasta/{BATCH['Ref.ver']}.fa \
                    {R1} \
                    {R2}"
        os.system(command)
#-----------------------------------------------------------------------------#
    if BATCH["Step"] == "All":
        LAST()
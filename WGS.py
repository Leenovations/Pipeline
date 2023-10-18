#!/usr/bin/python3

import sys
import os
import time
import argparse
import pandas as pd

#----------------------------------------------------------------------------------------#
parser = argparse.ArgumentParser(description="Pipeline Usage")
parser.add_argument("1", metavar="<38 or {sys.argv[1]}>", help="Select Reference version")
parser.add_argument("2", metavar="<Core>", help="Set Core")
parser.add_argument("3", metavar="<Step>", help="All, Align, Annotation, SV, etc")
args = parser.parse_args()
#----------------------------------------------------------------------------------------#
Sample = pd.read_csv("SampleSheet.txt", sep='\t', header=None)
Name = Sample.iloc[0, 0]
R1 = Sample.iloc[0, 1]
R2 = Sample.iloc[0, 2]
#----------------------------------------------------------------------------------------#
def PreQC(r1, r2):
    if os.path.isdir('00.PreQC'):
        pass
    else:
        command = 'mkdir 00.PreQC'
        os.system(command)

    command =f'fastqc -o 00.PreQC \
            -t {int(sys.argv[2])*2} \
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
def bwaindex():
    if os.path.exists(f'/media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta.amb'):
        pass
    else:
        command = f'bwa index -a bwtsw /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta'
        os.system(command)

    if os.path.isfile(f'/media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta.fai'):
        pass
    else:
        command = f'samtools faidx /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta'
        os.system(command)
#----------------------------------------------------------------------------------------#
def bwa(name):
    if os.path.isdir('03.Align'):
        pass
    else:
        command = 'mkdir 03.Align'
        os.system(command)

    command = f'bwa mem \
                -M \
                -R "@RG\\tID:{name}\\tPL:Illumina\\tLB:NovaSeq\\tSM:{name}" \
                -v 1 -t {sys.argv[2]} /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                02.Trimmed/{name}_val_1.fq.gz 02.Trimmed/{name}_val_2.fq.gz > 03.Align/{name}.sam'
    os.system(command)
#----------------------------------------------------------------------------------------#
def markduplicate(name):
    command = f'/media/src/Tools/gatk-4.4.0.0/gatk \
                MarkDuplicatesSpark \
                -I 03.Align/{name}.sam \
                -O 03.Align/{name}.MarkDuplicat.bam \
                -M 03.Align/{name}.MarkDuplicatesSpark.metrics.txt \
                --remove-all-duplicates true'
    os.system(command)
#----------------------------------------------------------------------------------------#
def makedict():
    if os.path.isfile(f'/media/src/hg{sys.argv[1]}/02.Fasta/hg{sys.argv[1]}.dict'):
        pass
    else:
        command = f'java -jar /Bioinformatics/00.Tools/picard/build/libs/picard.jar \
                    CreateSequenceDictionary \
                    R=/media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                    O=/media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.dict'
        os.system(command) 
#----------------------------------------------------------------------------------------#
def baserecalibrator(name):
    command = f"/media/src/Tools/gatk-4.4.0.0/gatk \
                BaseRecalibratorSpark \
                -I 03.Align/{name}.MarkDuplicate.bam \
                -R /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                --known-sites /media/src/hg{sys.argv[1]}/03.db/Homo_sapiens_assembly{sys.argv[1]}.dbsnp138.vcf \
                --known-sites /media/src/hg{sys.argv[1]}/03.db/Homo_sapiens_assembly{sys.argv[1]}.known_indels.vcf \
                --known-sites /media/src/hg{sys.argv[1]}/03.db/Mills_and_1000G_gold_standard.indels.hg{sys.argv[1]}.sites.vcf \
                -O 03.Align/{name}.Recalibrator.table"
    os.system(command)
#----------------------------------------------------------------------------------------#
def applyBQSR(name):
    command = f"/media/src/Tools/gatk-4.4.0.0/gatk \
                ApplyBQSRSpark \
                -R /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                -I 03.Align/{name}.MarkDuplicate.bam \
                -bqsr 03.Align/{name}.Recalibrator.table \
                -O 03.Align/{name}.bam \
                -static-quantized-quals 10 \
                -static-quantized-quals 20 \
                -static-quantized-quals 30"
    os.system(command)
#----------------------------------------------------------------------------------------#
def haplotypecaller(name):
    command = f"/media/src/Tools/gatk-4.4.0.0/gatk \
                HaplotypeCallerSpark \
                -R /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                -I 03.Align/{name}.bam \
                -O 03.Align/{name}.vcf \
                -L /media/src/hg{sys.argv[1]}/03.db/Homo_sapiens_assembly{sys.argv[1]}.whole_genome.interval_list \
                --bam-output 03.Align/{name}.Haplotype.bam \
                -OVI true \
                --emit-ref-confidence GVCF"
    os.system(command)
#----------------------------------------------------------------------------------------#
def mutect2(name):
    command = f"/media/src/Tools/gatk-4.4.0.0/gatk \
                Mutect2 \
                -R /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                -I 03.Align/{name}.bam \
                -O 03.Align/{name}.vcf \
                -L /media/src/hg{sys.argv[1]}/03.db/Homo_sapiens_assembly{sys.argv[1]}.whole_genome.interval_list \
                --bam-output 03.Align/{name}.mutect2.bam \
                -OVI true \
                --germline-resource /media/src/hg{sys.argv[1]}/03.db/af-only-gnomad.raw.sites.vcf \
                --panel-of-normals /media/src/hg{sys.argv[1]}/03.db/Mutect2-WGS-panel-hg{sys.argv[1]}.vcf"
    os.system(command)
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
def SV(Name):
    command = f"delly call \
                -g /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                -o {Name}.sv.bcf \
                -x /Bioinformatics/00.Tools/delly/excludeTemplates/human.hg{sys.argv[1]}.excl.tsv \
                {Name}.Real.bam"
    os.system(command)
#----------------------------------------------------------------------------------------#
def ChromosomeCNV(Name):
    command = f"cnvkit.py batch \
                {Name}.Real.bam -n \
                -\ /media/src/hg{sys.argv[1]}/hg{sys.argv[1]}.bed \
                -f /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                --output-dir ./"
    os.system(command)

    command = f"cnvkit.py segment \
                {Name}.cnr \
                -o {Name}.cns"
    # os.system(command)

    command = f"cnvkit.py scatter \
                {Name}.markdup.cnr -s {Name}.Real.cns -o {Name}.Chromosome.CNV.png"
    # os.system(command)
#----------------------------------------------------------------------------------------#
if sys.argv[3] == "All":
    PreQC(R1, R2)
    Trimming(Name, R1, R2)
    PostQC(Name)
    bwaindex()
    bwa(Name)
    markduplicate(Name)
    makedict()
    baserecalibrator(Name)
    applyBQSR(Name)
    haplotypecaller(Name)
    mutect2(Name)
    # SV(Name)
    # ChromosomeCNV(Name)
elif sys.argv[3] == "FastQC":
    pass
elif sys.argv[3] == "Align":
    pass
elif sys.argv[3] == "SV":
    SV(Name)
elif sys.argv[3] == "ChromosomeCNV":
    ChromosomeCNV(Name)
#----------------------------------------------------------------------------------------#
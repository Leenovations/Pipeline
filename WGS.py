#!/usr/bin/python3

import sys
import os
import time
import glob
import argparse
import pandas as pd

#----------------------------------------------------------------------------------------#
parser = argparse.ArgumentParser(description="Pipeline Usage")
parser.add_argument("1", metavar="<38 or 19>", help="Select Reference version")
parser.add_argument("2", metavar="<Core>", help="Set Core")
parser.add_argument("3", metavar="<Step>", help="All, Align, Annotation, SV, etc")
args = parser.parse_args()
#----------------------------------------------------------------------------------------#
Sample = pd.read_csv("SampleSheet.txt", sep='\t', header=None)
Name = Sample.iloc[0, 0]
R1 = Sample.iloc[0, 1]
R2 = Sample.iloc[0, 2]

Current_Dir = pd.read_csv('PWD.conf', sep='\t', header=None)
PWD = Current_Dir.iloc[0,0]
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
                {PWD}/02.Trimmed/{name}_val_1.fq.gz {PWD}/02.Trimmed/{name}_val_2.fq.gz > 03.Align/{name}.sam'
    os.system(command)
#----------------------------------------------------------------------------------------#
def markduplicate(name):
    if os.path.isdir('03.Align'):
        pass
    else:
        command = 'mkdir 03.Align'
        os.system(command)

    command = f'/media/src/Tools/gatk-4.4.0.0/gatk \
                --java-options "-Xmx3G -Xms1G" \
                MarkDuplicatesSpark \
                --spark-master local[{sys.argv[2]}] \
                -I {PWD}/03.Align/{name}.sam \
                -O 03.Align/{name}.MarkDuplicate.bam \
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
    if os.path.isdir('03.Align'):
        pass
    else:
        command = 'mkdir 03.Align'
        os.system(command)

    command = f"/media/src/Tools/gatk-4.4.0.0/gatk \
                BaseRecalibratorSpark \
                --conf spark.executor.cores={sys.argv[2]} \
                -I 03.Align/{name}.MarkDuplicate.bam \
                -R /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                --known-sites /media/src/hg{sys.argv[1]}/03.db/Homo_sapiens_assembly{sys.argv[1]}.dbsnp138.vcf \
                --known-sites /media/src/hg{sys.argv[1]}/03.db/Homo_sapiens_assembly{sys.argv[1]}.known_indels.vcf \
                --known-sites /media/src/hg{sys.argv[1]}/03.db/Mills_and_1000G_gold_standard.indels.hg{sys.argv[1]}.sites.vcf \
                -O 03.Align/{name}.Recalibrator.table"
    os.system(command)
#----------------------------------------------------------------------------------------#
def applyBQSR(name):
    command = f'rm -rf 03.Align/{name}.sam'
    os.system(command)
    
    if os.path.isdir('03.Align'):
        pass
    else:
        command = 'mkdir 03.Align'
        os.system(command)
    
    command = f"/media/src/Tools/gatk-4.4.0.0/gatk \
                ApplyBQSRSpark \
                --conf spark.executor.cores={sys.argv[2]} \
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
    if os.path.isdir('03.Align'):
        pass
    else:
        command = 'mkdir 03.Align'
        os.system(command)

    command = f"/media/src/Tools/gatk-4.4.0.0/gatk \
                HaplotypeCallerSpark \
                --conf spark.executor.cores={sys.argv[2]} \
                -R /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                -I {PWD}/03.Align/{name}.bam \
                -O 03.Align/{name}.haplotype.vcf \
                -L /media/src/hg{sys.argv[1]}/03.db/Homo_sapiens_assembly{sys.argv[1]}.whole_genome.interval_list \
                -OVI true \
                --emit-ref-confidence GVCF"
    os.system(command)
#----------------------------------------------------------------------------------------#
def mutect2(name):
    if os.path.isdir('03.Align'):
        pass
    else:
        command = 'mkdir 03.Align'
        os.system(command)

    command = f"/media/src/Tools/gatk-4.4.0.0/gatk \
                Mutect2 \
                --conf spark.executor.cores={sys.argv[2]} \
                -R /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                -I {PWD}/03.Align/{name}.bam \
                -O 03.Align/{name}.mutect2.vcf \
                -L /media/src/hg{sys.argv[1]}/03.db/Homo_sapiens_assembly{sys.argv[1]}.whole_genome.interval_list \
                -tumor {name} \
                --bam-output 03.Align/{name}.mutect2.bam \
                -OVI true \
                --callable-depth 5 \
                --germline-resource /media/src/hg{sys.argv[1]}/03.db/af-only-gnomad.raw.sites.vcf \
                --panel-of-normals /media/src/hg{sys.argv[1]}/03.db/Mutect2-WGS-panel-hg{sys.argv[1]}.vcf"
    os.system(command)
#----------------------------------------------------------------------------------------#
def getpileupsummaries(name):
    if os.path.isdir('03.Align'):
        pass
    else:
        command = 'mkdir 03.Align'
        os.system(command)

    command = f"/media/src/Tools/gatk-4.4.0.0/gatk \
                GetPileupSummaries \
                --conf spark.executor.cores={sys.argv[2]} \
                -I 03.Align/{name}.bam \
                -V /media/src/hg{sys.argv[1]}/03.db/small_exac_common_3.vcf \
                -L /media/src/hg{sys.argv[1]}/03.db/Homo_sapiens_assembly{sys.argv[1]}.whole_genome.interval_list \
                -O 03.Align/{name}.getpileupsummaries.table"
    os.system(command)
#----------------------------------------------------------------------------------------#
def calculatecontamination(name):
    if os.path.isdir('03.Align'):
        pass
    else:
        command = 'mkdir 03.Align'
        os.system(command)

    command = f"/media/src/Tools/gatk-4.4.0.0/gatk \
                CalculateContamination \
                --conf spark.executor.cores={sys.argv[2]} \
                -I 03.Align/{name}.getpileupsummaries.table \
                -O 03.Align/{name}.contamination.table"
    os.system(command)
#----------------------------------------------------------------------------------------#
def filtermutectcall(name):
    if os.path.isdir('03.Align'):
        pass
    else:
        command = 'mkdir 03.Align'
        os.system(command)

    command = f"/media/src/Tools/gatk-4.4.0.0/gatk \
                FilterMutectCalls \
                --conf spark.executor.cores={sys.argv[2]} \
                -R /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                -V 03.Align/{name}.mutect2.vcf \
                --contamination-table 03.Align/{name}.contamination.table \
                --stats 03.Align/{name}.mutect2.vcf.stats \
                -O 03.Align/{name}.mutect2.filtered.vcf"
    os.system(command)
#----------------------------------------------------------------------------------------#
def varscan2(name):
    if os.path.isdir('03.Align'):
        pass
    else:
        command = 'mkdir 03.Align'
        os.system(command)

    command = f'samtools mpileup \
                -f /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                --max-depth 1000000 \
                03.Align/{name}.bam > 03.Align/{name}.mpileup'
    os.system(command)

    command = f'java -jar /Bioinformatics/00.Tools/varscan-2.4.5/VarScan.v2.4.1.jar \
                mpileup2cns 03.Align/{name}.mpileup \
                --min-avg-qual 20 --min-coverage 10 --min-reads2 3 \
                --min-var-freq 0.001 --variants \
                --output-vcf 1 > 03.Align/{name}.varscan2.vcf'
    os.system(command)
#----------------------------------------------------------------------------------------#
def SV(name):
    if os.path.isdir('04.SV'):
        pass
    else:
        command = 'mkdir 04.SV'
        os.system(command)

    command = f"delly call \
                -g /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                -o 04.SV/{Name}.sv.bcf \
                -x /Bioinformatics/00.Tools/delly/excludeTemplates/human.hg{sys.argv[1]}.excl.tsv \
                03.Align/{name}.bam"
    os.system(command)

    command = f"bcftools view 04.SV/{Name}.sv.bcf -Oz > 04.SV/{Name}.sv.vcf.gz"
    os.system(command)

    command = f"bcftools view \
               -i 'FILTER=\"PASS\" & INFO/SVTYPE!=\"DEL\" & INFO/SVTYPE!=\"INS\" & INFO/SVTYPE!=\"DUP\" & INFO/SVTYPE!=\"INV\"' \
               04.SV/{Name}.sv.vcf.gz > 04.SV/{Name}.sv.filtered.vcf"
    os.system(command)
#----------------------------------------------------------------------------------------#
def ChromosomeCNV(name):
    if os.path.isdir('04.SV'):
        pass
    else:
        command = 'mkdir 04.SV'
        os.system(command)

    if os.path.exists(f'/media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.bed'):
        pass
    else:
        command = f'cnvkit.py access \
                    /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                    -o /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.bed'
        os.system(command)

    #Run each sample & merge all cnn & time delay
    # command = f'cnvkit.py autobin \
    #             03.Align/{name}.bam \
    #             -t /media/src/hg{sys.argv[1]}/04.cnv/whole.exome.exon.bed \
    #             -g /media/src/hg{sys.argv[1]}/04.cnv/access.hg{sys.argv[1]}.bed \
    #             --target-output-bed 04.SV/CNV.target.bed \
    #             --antitarget-output-bed 04.SV/CNV.antitarget.bed'
    # os.system(command)

    # command = f'cnvkit.py coverage \
    #             03.Align/{name}.bam \
    #             04.SV/CNV.target.bed \
    #             -o 04.SV/{name}.targetcoverage.cnn'
    # os.system(command)

    # command = f'cnvkit.py coverage \
    #             03.Align/{name}.bam \
    #             04.SV/CNV.antitarget.bed \
    #             -o 04.SV/{name}.antitargetcoverage.cnn'
    # os.system(command)

    command = f'cnvkit.py reference \
                04.SV/{name}.targetcoverage.cnn \
                04.SV/{name}.antitargetcoverage.cnn \
                -f /media/src/hg{sys.argv[1]}/02.Fasta/Homo_sapiens_assembly{sys.argv[1]}.fasta \
                -o 04.SV/Reference.cnn'
    os.system(command)

    command = f'cnvkit.py fix \
                04.SV/{name}.targetcoverage.cnn \
                04.SV/{name}.antitargetcoverage.cnn \
                04.SV/Reference.cnn \
                -o 04.SV/{name}.cnr'
    os.system(command)

    command = f'cnvkit.py segment \
                04.SV/{name}.cnr \
                -o 04.SV/{name}.cns'
    os.system(command)

    command = f'cnvkit.py scatter \
                04.SV/{name}.cnr \
                -s 04.SV/{name}.cns \
                -o 04.SV/{name}.whole.pdf'
    os.system(command)            

    command = f'cnvkit.py scatter \
                -s 04.SV/{name}.cnr \
                -s 04.SV/{name}.cns \
                -o 04.SV/{name}.BCR.scatter.pdf \
                -g BCR'
    os.system(command)                
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
    getpileupsummaries(Name)
    calculatecontamination(Name)
    filtermutectcall(Name)
    varscan2(Name)
    SV(Name)
    ChromosomeCNV(Name)
elif sys.argv[3] == "FastQC":
    PreQC(R1, R2)
    Trimming(Name, R1, R2)
    PostQC(Name)
elif sys.argv[3] == "Align":
    bwaindex()
    bwa(Name)
elif sys.argv[3] == "Dedup":
    markduplicate(Name)
    baserecalibrator(Name)
    applyBQSR(Name)
elif sys.argv[3] == "Mutation":
    haplotypecaller(Name)
    mutect2(Name)
    getpileupsummaries(Name)
    calculatecontamination(Name)
    filtermutectcall(Name)
    varscan2(Name)
elif sys.argv[3] == "SV":
    SV(Name)
elif sys.argv[3] == "ChromosomeCNV":
    ChromosomeCNV(Name)
#----------------------------------------------------------------------------------------#
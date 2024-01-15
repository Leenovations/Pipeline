#!/home/lab/anaconda3/envs/NGS/bin/python3

import sys
import os
import time
import glob
import argparse
import pandas as pd
import numpy as np
import re
from collections import defaultdict
from functools import reduce
#----------------------------------------------------------------------------------------#
parser = argparse.ArgumentParser(description="Pipeline Usage")
args = parser.parse_args()
#----------------------------------------------------------------------------------------#
Sample = pd.read_csv("SampleSheet.txt", sep="\t", header=None)
Name = Sample.iloc[0, 0]
R1 = Sample.iloc[0, 1]
R2 = Sample.iloc[0, 2]
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
def bwaindex():
    if os.path.exists(f"/media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta.amb"):
        pass
    else:
        command = f"bwa index -a bwtsw /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta"
        os.system(command)

    if os.path.isfile(f"/media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta.fai"):
        pass
    else:
        command = f"samtools faidx /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta"
        os.system(command)
#----------------------------------------------------------------------------------------#
def bwa(name):
    if os.path.isdir("03.Align"):
        pass
    else:
        command = "mkdir 03.Align"
        os.system(command)

    command = f"bwa mem \
                -M \
                -R '@RG\\tID:{name}\\tPL:Illumina\\tLB:NovaSeq\\tSM:{name}' \
                -v 1 -t {BATCH['CPU']} /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                02.Trimmed/{name}_val_1.fq.gz 02.Trimmed/{name}_val_2.fq.gz | \
                samtools view -bS - > 03.Align/{name}.bwa.bam"
    os.system(command)
#----------------------------------------------------------------------------------------#
def AddOrReplaceReadGroups(name):
    if os.path.isdir("03.Align"):
        pass
    else:
        command = "mkdir 03.Align"
        os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /Bioinformatics/00.Tools/picard/build/libs/picard.jar \
                AddOrReplaceReadGroups \
                I=03.Align/{name}.bwa.bam \
                O=03.Align/{name}.sorted.bam \
                TMP_DIR=TEMP \
                RGLB=NGS \
                RGPL=Illumina \
                RGPU={name} \
                RGSM={name} \
                CREATE_INDEX=true \
                VALIDATION_STRINGENCY=LENIENT \
                SO=coordinate"
    os.system(command)
#----------------------------------------------------------------------------------------#
def markduplicate(name):
    if os.path.isdir("03.Align"):
        pass
    else:
        command = "mkdir 03.Align"
        os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /Bioinformatics/00.Tools/picard/build/libs/picard.jar \
                MarkDuplicates \
                I=03.Align/{name}.sorted.bam \
                O=03.Align/{name}.MarkDuplicate.bam \
                M=03.Align/{name}.MarkDuplicatesSpark.metrics.txt \
                TMP_DIR=TEMP \
                REMOVE_DUPLICATES=true \
                VALIDATION_STRINGENCY=LENIENT \
                AS=true"
    os.system(command)
#----------------------------------------------------------------------------------------#
def makedict():
    if os.path.isfile(f"/media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/hg{BATCH['Ref.ver'].split('g')[1]}.dict"):
        pass
    else:
        command = f"java -jar /Bioinformatics/00.Tools/picard/build/libs/picard.jar \
                    CreateSequenceDictionary \
                    R=/media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                    O=/media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.dict"
        os.system(command) 
#----------------------------------------------------------------------------------------#
def baserecalibrator(name):
    if os.path.isdir("03.Align"):
        pass
    else:
        command = "mkdir 03.Align"
        os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                BaseRecalibrator \
                -I 03.Align/{name}.MarkDuplicate.bam \
                -R /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                --known-sites /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/03.db/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.dbsnp138.vcf \
                --known-sites /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/03.db/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.known_indels.vcf \
                --known-sites /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/03.db/Mills_and_1000G_gold_standard.indels.hg{BATCH['Ref.ver'].split('g')[1]}.sites.vcf \
                -O 03.Align/{name}.Recalibrator.table"
    os.system(command)
#----------------------------------------------------------------------------------------#
def applyBQSR(name):
    if os.path.isdir("03.Align"):
        pass
    else:
        command = "mkdir 03.Align"
        os.system(command)
    
    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                ApplyBQSR \
                -R /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                -I 03.Align/{name}.MarkDuplicate.bam \
                -bqsr 03.Align/{name}.Recalibrator.table \
                -O 03.Align/{name}.bam"
    os.system(command)
#----------------------------------------------------------------------------------------#
def haplotypecaller(name):
    if os.path.isdir("03.Align"):
        pass
    else:
        command = "mkdir 03.Align"
        os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                HaplotypeCaller \
                -R /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                -I 03.Align/{name}.bam \
                -O 03.Align/{name}.haplotype.vcf \
                -L /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/03.db/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.whole_genome.interval_list \
                -ERC GVCF \
                --standard-min-confidence-threshold-for-calling 20" 
    os.system(command)
#----------------------------------------------------------------------------------------#
def Variantfilter(name):
    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                GenotypeGVCFs \
                -R /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                -V 03.Align/{name}.haplotype.vcf \
                -O 03.Align/{name}.haplotype.genotype.vcf"
    os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                SelectVariants \
                -R /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                -V 03.Align/{name}.haplotype.genotype.vcf \
                --select-type-to-include SNP \
                -O 03.Align/{name}.SNPs.vcf"
    os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                SelectVariants \
                -R /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                -V 03.Align/{name}.haplotype.genotype.vcf \
                --select-type-to-include INDEL \
                -O 03.Align/{name}.INDELs.vcf"
    os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                VariantFiltration \
                -V 03.Align/{name}.SNPs.vcf \
                -O 03.Align/{name}.SNPs.flt.vcf \
                --filter-expression 'QD < 2.0' --filter-name 'QD2' \
                --filter-expression 'QUAL < 30.0' --filter-name 'QUAL30' \
                --filter-expression 'SOR > 3.0' --filter-name 'SOR3' \
                --filter-expression 'FS > 60.0' --filter-name 'FS60' \
                --filter-expression 'MQ < 40.0' --filter-name 'MQ40' \
                --filter-expression 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' \
                --filter-expression 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8'" 
    os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                VariantFiltration \
                -R /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                -V 03.Align/{name}.INDELs.vcf \
                -O 03.Align/{name}.INDELs.flt.vcf \
                --filter-expression 'QD < 2.0' --filter-name 'QD2' \
                --filter-expression 'FS > 200.0' --filter-name 'FS200' \
                --filter-expression 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' \
                --filter-expression 'SOR > 10.0' --filter-name 'SOR10'"
    os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                SortVcf \
                -I 03.Align/{name}.SNPs.flt.vcf \
                -I 03.Align/{name}.INDELs.flt.vcf \
                -O 03.Align/{name}.vcf"
    os.system(command)

    command = f"egrep '^#|PASS' 03.Align/{name}.vcf > 03.Align/{name}.PASS.vcf"
    os.system(command)
#----------------------------------------------------------------------------------------#
def mutect2(name):
    if os.path.isdir("03.Align"):
        pass
    else:
        command = "mkdir 03.Align"
        os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                Mutect2 \
                -R /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                -I 03.Align/{name}.bam \
                -O 03.Align/{name}.mutect2.vcf \
                -L /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/03.db/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.whole_genome.interval_list \
                -tumor {name} \
                --bam-output 03.Align/{name}.mutect2.bam \
                -OVI true \
                --callable-depth 5 \
                --germline-resource /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/03.db/af-only-gnomad.raw.sites.vcf \
                --panel-of-normals /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/03.db/Mutect2-WGS-panel-hg{BATCH['Ref.ver'].split('g')[1]}.vcf"
    os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                GetPileupSummaries \
                -I 03.Align/{name}.bam \
                -V /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/03.db/small_exac_common_3.vcf \
                -L /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/03.db/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.whole_genome.interval_list \
                -O 03.Align/{name}.getpileupsummaries.table"
    os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                CalculateContamination \
                -I 03.Align/{name}.getpileupsummaries.table \
                -O 03.Align/{name}.contamination.table"
    os.system(command)

    command = f"java \
                -Xmx32G \
                -XX:ParallelGCThreads={str(2*int({BATCH['CPU']}))} \
                -jar /media/src/Tools/gatk-4.4.0.0/gatk-package-4.4.0.0-local.jar \
                FilterMutectCalls \
                -R /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                -V 03.Align/{name}.mutect2.vcf \
                --contamination-table 03.Align/{name}.contamination.table \
                --stats 03.Align/{name}.mutect2.vcf.stats \
                -O 03.Align/{name}.mutect2.final.vcf"
    os.system(command)

    command = f"egrep '^#|PASS' 03.Align/{name}.mutect2.final.vcf > 03.Align/{name}.PASS.vcf"
    os.system(command)
#----------------------------------------------------------------------------------------#
def varscan2(name):
    if os.path.isdir("03.Align"):
        pass
    else:
        command = "mkdir 03.Align"
        os.system(command)

    command = f"samtools mpileup \
                -f /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                --max-depth 1000000 \
                03.Align/{name}.bam > 03.Align/{name}.mpileup"
    os.system(command)

    command = f"java -jar /Bioinformatics/00.Tools/varscan-2.4.5/VarScan.v2.4.1.jar \
                mpileup2cns 03.Align/{name}.mpileup \
                --min-avg-qual 20 --min-coverage 10 --min-reads2 3 \
                --min-var-freq 0.001 --variants \
                --output-vcf 1 > 03.Align/{name}.varscan2.vcf"
    os.system(command)
#----------------------------------------------------------------------------------------#
def Annotation(name):
    if os.path.isdir("04.Annotation"):
        pass
    else:
        command = "mkdir 04.Annotation"
        os.system(command)

    command = f"java -jar /media/src/Tools/snpEff/snpEff.jar \
                -v hg{BATCH['Ref.ver'].split('g')[1]} \
                03.Align/{name}.PASS.vcf > 03.Align/{name}.snpeff.vcf"
    os.system(command)

    command = f"convert2annovar.pl -includeinfo -allsample -withfreq -format vcf4 03.Align/{name}.snpeff.vcf > 04.Annotation/{name}.avinput"
    os.system(command)

    command = f"annotate_variation.pl -geneanno -out 04.Annotation/{name}.hgvs -build hg{BATCH['Ref.ver'].split('g')[1]} \
                -dbtype refGene \
                -hgvs 04.Annotation/{name}.avinput /Bioinformatics/00.Tools/annovar/humandb"
    os.system(command)

    command = f"table_annovar.pl 04.Annotation/{name}.avinput /Bioinformatics/00.Tools/annovar/humandb \
                -buildver hg{BATCH['Ref.ver'].split('g')[1]} -out 04.Annotation/{name} \
                -remove -protocol \
                refGene,dbnsfp33a,cosmic70,snp138,snp138NonFlagged,popfreq_max_20150413,popfreq_all_20150413,dbscsnv11,exac03nontcga,avsnp147,clinvar_20160302,gnomad_exome,gnomad_genome \
                -operation g,f,f,f,f,f,f,f,f,f,f,f,f \
                -nastring . -otherinfo"
    os.system(command)
#----------------------------------------------------------------------------------------#
def SV(name):
    if os.path.isdir("04.SV"):
        pass
    else:
        command = "mkdir 04.SV"
        os.system(command)

    command = f"delly call \
                -g /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                -o 04.SV/{Name}.sv.bcf \
                -x /Bioinformatics/00.Tools/delly/excludeTemplates/human.hg{BATCH['Ref.ver'].split('g')[1]}.excl.tsv \
                03.Align/{name}.bam"
    os.system(command)

    command = f"bcftools view 04.SV/{Name}.sv.bcf -Oz > 04.SV/{Name}.sv.vcf.gz"
    os.system(command)

    command = f"bcftools view \
               -i 'FILTER=\'PASS\' & INFO/SVTYPE!=\'DEL\' & INFO/SVTYPE!=\'INS\' & INFO/SVTYPE!=\'DUP\' & INFO/SVTYPE!=\'INV\'' \
               04.SV/{Name}.sv.vcf.gz > 04.SV/{Name}.sv.filtered.vcf"
    os.system(command)
#----------------------------------------------------------------------------------------#
def ChromosomeCNV(name):
    if os.path.isdir("04.SV"):
        pass
    else:
        command = "mkdir 04.SV"
        os.system(command)

    if os.path.exists(f"/media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.bed"):
        pass
    else:
        command = f"cnvkit.py access \
                    /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                    -o /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.bed"
        os.system(command)

    #Run each sample & merge all cnn & time delay
    # command = f"cnvkit.py autobin \
    #             03.Align/{name}.bam \
    #             -t /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/04.cnv/whole.exome.exon.bed \
    #             -g /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/04.cnv/access.hg{BATCH['Ref.ver'].split('g')[1]}.bed \
    #             --target-output-bed 04.SV/CNV.target.bed \
    #             --antitarget-output-bed 04.SV/CNV.antitarget.bed"
    # os.system(command)

    # command = f"cnvkit.py coverage \
    #             03.Align/{name}.bam \
    #             04.SV/CNV.target.bed \
    #             -o 04.SV/{name}.targetcoverage.cnn"
    # os.system(command)

    # command = f"cnvkit.py coverage \
    #             03.Align/{name}.bam \
    #             04.SV/CNV.antitarget.bed \
    #             -o 04.SV/{name}.antitargetcoverage.cnn"
    # os.system(command)

    command = f"cnvkit.py reference \
                04.SV/{name}.targetcoverage.cnn \
                04.SV/{name}.antitargetcoverage.cnn \
                -f /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/02.Fasta/Homo_sapiens_assembly{BATCH['Ref.ver'].split('g')[1]}.fasta \
                -o 04.SV/Reference.cnn"
    os.system(command)

    command = f"cnvkit.py fix \
                04.SV/{name}.targetcoverage.cnn \
                04.SV/{name}.antitargetcoverage.cnn \
                04.SV/Reference.cnn \
                -o 04.SV/{name}.cnr"
    os.system(command)

    command = f"cnvkit.py segment \
                04.SV/{name}.cnr \
                -o 04.SV/{name}.cns"
    os.system(command)

    command = f"cnvkit.py scatter \
                04.SV/{name}.cnr \
                -s 04.SV/{name}.cns \
                -o 04.SV/{name}.whole.pdf"
    os.system(command)            

    command = f"cnvkit.py scatter \
                -s 04.SV/{name}.cnr \
                -s 04.SV/{name}.cns \
                -o 04.SV/{name}.BCR.scatter.pdf \
                -g BCR"
    os.system(command)                
#----------------------------------------------------------------------------------------#
def ChromosomalCNV(name):
    if os.path.isdir("05.SV/00.ChromosomeCNV"):
        pass
    else:
        command = "mkdir -p 05.SV/00.ChromosomeCNV"
        os.system(command)
    
    command = f"samtools bedcov \
                /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/01.Methylation/00.Bed/1MB.exclude.centromere.bed \
                03.Align/{name}.bam > 05.SV/00.ChromosomeCNV/{name}.bedcov"
    os.system(command)

    Chromosome = [str(i) for i in range(1,23)] + ['X', 'Y']
    Data = pd.read_csv(f"05.SV/00.ChromosomeCNV/{name}.bedcov",
                    sep='\t',
                    header=None,
                    names = ['Chr', 'Start', 'End', 'Count'],
                    low_memory=False)
    
    Data['Length'] = Data['End'] - Data['Start'] + 1
    Data['count_per_length'] = Data['Count'] / Data['Length']
    Data['TPM'] = Data['Count'] / Data['Length'] * Data['count_per_length']
    Data['TPM'] = np.log10(Data['TPM'] + 1)
    Median_TPM = Data['TPM'].median()
    Data['TPM'] = Data['TPM'] - Median_TPM
    
    Sorted = []
    for chromosome in Chromosome:
        Data_sub = Data[Data['Chr'] == chromosome]
        Sorted.append(Data_sub)
    Sorted_Data = pd.concat(Sorted)
    Sorted_Data['Order'] = range(1, len(Sorted_Data) + 1)
    
    Sorted_Data.to_csv(f"05.SV/00.ChromosomeCNV/{name}.Chromosome.CNV.txt",
                        sep='\t',
                        header='infer',
                        index=False)
    
    command = f"Rscript /labmed/00.Code/Pipeline/WGS.ChromosomalCNV.R {name}"
    os.system(command)
#----------------------------------------------------------------------------------------#
def GeneCNV(name):
    if os.path.isdir("05.SV/01.GeneCNV"):
        pass
    else:
        command = "mkdir -p 05.SV/01.GeneCNV"
        os.system(command)

    Gene = pd.read_csv(f"/media/src/hg{BATCH['Ref.ver'].split('g')[1]}/04.cnv/NCBI.RefSeq.Selected.Gene.bed",
                       sep='\t',
                       low_memory=False,
                       header=None)
    Gene = list(set(Gene.iloc[:, 3].to_list()))
    Gene.sort()

    for gene in Gene[0:1]:
        command = f"samtools bedcov \
                    /media/src/hg{BATCH['Ref.ver'].split('g')[1]}/04.cnv/{gene}.cnv.bed \
                    03.Align/{name}.bam > 05.SV/01.GeneCNV/{name}.{gene}.bedcov"
        os.system(command)

    names = BATCH['Sample.Name'].split(',')
    Dir = BATCH['Sample.Dir'].split(',')
    for gene in Gene[0:1]:
        file_paths = [file + f"05.SV/01.GeneCNV/{name}.{gene}.bedcov" for file, name in zip(Dir, names)]

        data_frames = [pd.read_csv(file, 
                                   sep='\t', 
                                   low_memory=False, 
                                   header=None,
                                   names=['Chr', 'Start', 'End', 'Gene', 'Exon', 'Strand', name])
                                   for file, name in zip(file_paths, names)]

        DATA = reduce(lambda left, right: pd.merge(left, right, on=['Chr', 'Start', 'End', 'Gene', 'Exon', 'Strand']), data_frames)
        Info = DATA.iloc[:, :6]
        Coverage = DATA.iloc[:, 6:]

        Exon_length = DATA['End'] - DATA['Start']
        Depth_Length = Coverage.apply(lambda count: (count / Exon_length))
        Sum_Depth_Length = Depth_Length.sum(axis=0)

        Exon_length = np.array(Exon_length)
        Sum_Depth_Length = np.array(Sum_Depth_Length)
        result_matrix = np.outer(Exon_length, Sum_Depth_Length)
        NormFactor = pd.DataFrame(result_matrix)
        NormFactor.columns = list(Coverage.columns)
        NormFactor.index = DATA.index.to_list()
        CNV_ALL = Coverage.div(NormFactor)
        CNV_ALL = CNV_ALL + 1
        CNV_ALL = CNV_ALL.applymap(np.log2)
        CNV_ALL = CNV_ALL.apply(lambda value : value - CNV_ALL.mean(axis=1))
        DATA = pd.concat([Info, CNV_ALL], axis=1)
        DATA = DATA[['Chr', 'Start', 'End', 'Gene', 'Exon', 'Strand', f'{name}']]
        DATA.to_csv(f"05.SV/01.GeneCNV/{name}.{gene}.Norm.CNV.txt",
                    sep='\t',
                    index=False,
                    header='infer')
        
        command = f"Rscript /labmed/00.Code/Pipeline/WGS.GeneCNV.VIZ.R {name} {gene}"
        os.system(command)
#----------------------------------------------------------------------------------------#
def Results(name):
    clinvar = pd.read_csv(f"/media/src/hg{BATCH['Ref.ver'].split('g')[1]}/a.clinvar.guideline.txt", sep='\t')
    clinvar['id'] = 'chr' + clinvar['id']
    clinvar['significance.submission'] = clinvar['significance.submission'].str.split('~').str[0]
    classification = dict(zip(clinvar['id'].to_list(), clinvar['significance.submission'].to_list()))

    NM = pd.read_csv(f"/media/src/hg{BATCH['Ref.ver'].split('g')[1]}/06.Annotation/UCSC.hg{BATCH['Ref.ver'].split('g')[1]}.Canonical.NMnumber.txt", 
                     sep='\t')
    Main_Accession = dict(zip(NM['GeneSymbol'].to_list(), NM['NM_number'].to_list()))

    VARIANTS = defaultdict(list)
    Anno = pd.read_csv(f"04.Annotation/{Name}.hg{BATCH['Ref.ver'].split('g')[1]}_multianno.txt", 
                       sep='\t',
                       low_memory=False)
    Header_Info = list(Anno.columns[10:])
    Other = Anno.iloc[:,10:]
    Flag = Anno.iloc[:,141]
    SnpEff = Anno.iloc[:,142]
    Format_info = Anno.iloc[:,143]
    Value_info = Anno.iloc[:,144]
    Gene_detail = Anno.iloc[:, 7]
    Gene_detail.replace('.', '', inplace=True)
    AAchange = Anno.iloc[:, 9]
    AAchange.replace('.', '', inplace=True)
    NM_info = Anno.iloc[:, 7] + Anno.iloc[:, 9]

    for i in range(Anno.shape[0]):
        Chr = str(Anno.iloc[i,0])
        Start = str(Anno.iloc[i,1])
        End = str(Anno.iloc[i,2])
        Ref = Anno.iloc[i,138]
        Alt = Anno.iloc[i,139]
        Class = '_'.join([Chr, Start, End, Ref, Alt])
        Region = Anno.iloc[i,5]
        Gene = Anno.iloc[i,6]
        Other = Anno.iloc[i,10:]
        Flag = Anno.iloc[i,141]
        GATK = Value_info.iloc[i]
        GATK = GATK.split(':')
        Format = Format_info.iloc[i]
        Format = Format.split(':')

        GT_GATK = str(GATK[Format.index('GT')])
        GT_GATK = re.split('/|\|', GT_GATK)
        GT_GATK = list(map(int, GT_GATK))
        if GT_GATK[1] != 1:
            continue
        AD_GATK = int(list(map(int, str(GATK[Format.index('AD')]).split(',')))[GT_GATK[1]])
        DP_GATK = int(sum(list(map(int, str(GATK[Format.index('AD')]).split(',')))))
        VAF_GATK = float(round(int(AD_GATK) / int(DP_GATK) * 100, 2))

        SnpEff_info = SnpEff[i].replace(';LOF', '|LOF')
        SnpEff_info = SnpEff_info.replace(';NMD', '|NMD')
        
        SnpEff_info = SnpEff_info.split(';')[-1].split('|')[0:10]
        SnpEff_Prediction = SnpEff_info[2]
        SnpEff_Variation = SnpEff_info[9]

        if Anno.iloc[i,5] == 'exonic':
            vtype = Anno.iloc[i,8]
            vinfo = NM_info[i]
            if vtype != 'unknown':
                for variant in vinfo.split(','):
                    variant = variant.split(':')
                    NM = variant[1]
                    Exon = variant[2]
                    HGVSc = variant[3]
                    HGVSp = variant[4]

                    if Class in classification.keys():
                        acmg = classification[Class]
                    else:
                        acmg = 'Uncertain significance'

                    if NM in Main_Accession.values():
                        Main = 'O'
                    else:
                        Main = '.'

                    VARIANTS[Chr + Start + Ref + Alt + Gene + HGVSp + NM].extend([
                            '',\
                            acmg,\
                            Main,\
                            'chr'+ Chr + ':' + Start + '-' + End,\
                            Region,\
                            vtype,\
                            Gene,\
                            NM,\
                            SnpEff_Variation,\
                            Ref,\
                            Alt,\
                            HGVSp,\
                            VAF_GATK,\
                            AD_GATK,\
                            DP_GATK,\
                            SnpEff_Prediction,\
                            Flag,\
                            Exon])

                    for remain in Other:
                        VARIANTS[Chr + Start + Ref + Alt + Gene + HGVSp + NM].append(remain)

            else:
                NM = '.'
                Exon = '.'
                HGVSc = '.'
                HGVSp = '.'
                vinfo = '.'
                Main = '.'
                VARIANTS[Chr + Start + Ref + Alt + Gene + HGVSp + NM].extend([
                    '',\
                    acmg,\
                    Main,\
                    'chr'+ Chr + ':' + Start + '-' + End,\
                    Region,\
                    vinfo,\
                    Gene,\
                    NM,\
                    SnpEff_Variation,\
                    Ref,\
                    Alt,\
                    HGVSp,\
                    VAF_GATK,\
                    AD_GATK,\
                    DP_GATK,\
                    SnpEff_Prediction,\
                    Flag,\
                    Exon])
                for remain in Other:
                        VARIANTS[Chr + Start + Ref + Alt + Gene + HGVSp + NM].append(remain)

        else:
            NM = '.'
            Exon = '.'
            HGVSc = '.'
            HGVSp = '.'
            vinfo = '.'
            Main = '.'

            if Class in classification.keys():
                acmg = classification[Class]
            else:
                acmg = 'Uncertain significance'

            VARIANTS[Chr + Start + Ref + Alt + Gene + HGVSp + NM].extend([
                    '',\
                    acmg,\
                    Main,\
                    'chr'+ Chr + ':' + Start + '-' + End,\
                    Region,\
                    vinfo,\
                    Gene,\
                    NM,\
                    SnpEff_Variation,\
                    Ref,\
                    Alt,\
                    HGVSp,\
                    VAF_GATK,\
                    AD_GATK,\
                    DP_GATK,\
                    SnpEff_Prediction,\
                    Flag,\
                    Exon])

            for remain in Other:
                VARIANTS[Chr + Start + Ref + Alt + Gene + HGVSp + NM].append(remain)

#Annotation Exel
    Data = pd.DataFrame(VARIANTS)
    Data = Data.transpose()

    COLUMNS = ['Select', 'ACMG', 'Main', 'Chromosome Position', 'Region', 'Variant Type', 'Gene', 'NM number', 'HGVSc', 'REF', 'ALT',\
                'HGVSp', 'VAF (%)', 'AD', 'DP', 'SnpEff' ,'Flag', 'Exon locus'] + Header_Info
    Data.columns = COLUMNS

    writer = pd.ExcelWriter(f"04.Annotation/{Name}.results.xlsx") 

    Data.to_excel(writer, sheet_name='Result', index=False, na_rep='NaN')
    writer.save()
#----------------------------------------------------------------------------------------#
if BATCH["Step"] == "All":
    # PreQC(R1, R2)
    # Trimming(Name, R1, R2)
    # PostQC(Name)
    # bwaindex()
    # bwa(Name)
    # AddOrReplaceReadGroups(Name)
    # markduplicate(Name)
    # makedict()
    # baserecalibrator(Name)
    # applyBQSR(Name)
    # haplotypecaller(Name)
    # Variantfilter(Name)
    # mutect2(Name)
    # varscan2(Name)
    # Annotation(Name)
    # SV(Name)
    # ChromosomeCNV(Name)
    # ChromosomalCNV(Name)
    GeneCNV(Name)
    # Results(Name)
elif BATCH["Step"] == "FastQC":
    PreQC(R1, R2)
    Trimming(Name, R1, R2)
    PostQC(Name)
elif BATCH["Step"] == "Align":
    bwaindex()
    bwa(Name)
elif BATCH["Step"] == "Dedup":
    AddOrReplaceReadGroups(Name)
    markduplicate(Name)
    baserecalibrator(Name)
    applyBQSR(Name)
elif BATCH["Step"] == "Mutation":
    haplotypecaller(Name)
    mutect2(Name)
    varscan2(Name)
elif BATCH["Step"] == "SV":
    SV(Name)
elif BATCH["Step"] == "ChromosomeCNV":
    ChromosomeCNV(Name)
#----------------------------------------------------------------------------------------#
#!/home/lab/anaconda3/envs/NGS/bin/python3

import os
import pandas as pd
import numpy as np
import time
from functools import reduce
#--------------------------------------------------------------------------------#
def ReadCount(sample, Ref, Col):
    Gene = pd.read_csv(f'/media/src/hg{Ref}/00.RNA/Index/geneInfo.tab', sep='\t', header=None)
    column_dict = {}
    for index, row in Gene.iterrows():
        column_dict[row[0]] = row[1]

    Genecount = pd.read_csv(f"03.Output/{sample}_ReadsPerGene.out.tab", sep='\t', header=None, names=['ID', 0, 1, 2])
    Genecount = Genecount.drop(index=range(0, 4))
    Genecount['GeneSymbol'] = [column_dict[key] for key in Genecount['ID']]
    Genecount = Genecount[['ID', 'GeneSymbol', int(Col)]]
    Genecount.columns = ['ID', 'GeneSymbol', f"{sample}"]
    Genecount.to_csv(f"03.Output/{sample}.Genecount.txt", sep='\t', header='infer', index=False)

def TPM(sample):
    Length = pd.read_csv('/media/src/hg19/00.RNA/hg19.GENCODE.v44.GeneLength.txt',
                sep='\t',
                header='infer',
                names =['ID', 'Type', 'GeneSymbol', 'Length'])
    LENGTH = dict(zip(Length['ID'].to_list(), Length['Length'].to_list()))
    GENE = dict(zip(Length['ID'].to_list(), Length['GeneSymbol'].to_list()))

    GeneCount = pd.read_csv(f"03.Output/{sample}.Genecount.txt",
                            sep='\t',
                            header='infer')

    GeneCount['Length'] = GeneCount['ID'].map(LENGTH)
    GeneCount['count_per_length'] = GeneCount[f"{sample}"] / GeneCount['Length']
    GeneCount['TPM'] = GeneCount[f"{sample}"] / (GeneCount['Length'] * GeneCount['count_per_length'].sum())
    GeneCount = GeneCount[['ID', 'GeneSymbol', 'TPM']]
    GeneCount.columns = ['ID', 'GeneSymbol', f"{sample}"]

    GeneCount.to_csv(f"03.Output/{sample}.TPM.txt",
                     sep='\t',
                     index=False,
                     header='infer')
#--------------------------------------------------------------------------------#
def FPKM(sample):
    Length = pd.read_csv('/media/src/hg19/00.RNA/hg19.GENCODE.v44.GeneLength.txt',
                sep='\t',
                header='infer',
                names =['ID', 'Type', 'GeneSymbol', 'Length'])
    LENGTH = dict(zip(Length['ID'].to_list(), Length['Length'].to_list()))
    GENE = dict(zip(Length['ID'].to_list(), Length['GeneSymbol'].to_list()))

    GeneCount = pd.read_csv(f"03.Output/{sample}.Genecount.txt",
                            sep='\t',
                            header='infer')

    GeneCount['Length'] = GeneCount['ID'].map(LENGTH)
    GeneCount['count_per_length'] = GeneCount[f"{sample}"] / GeneCount['Length']
    GeneCount['FPKM'] = (GeneCount[f"{sample}"] * 10**9) / (GeneCount['Length'] * GeneCount[f"{sample}"].sum())
    GeneCount = GeneCount[['ID', 'GeneSymbol', 'FPKM']]
    GeneCount.columns = ['ID', 'GeneSymbol', f"{sample}"]

    GeneCount.to_csv(f"03.Output/{sample}.FPKM.txt",
                     sep='\t',
                     index=False,
                     header='infer')
#--------------------------------------------------------------------------------#
def Merge(sample, samples, Dirs, TPM , FPKM):
    Samples = samples
    Dir = Dirs
    Genecount_paths = [file + f"03.Output/{sample}.Genecount.txt" for file, sample in zip(Dir, Samples)]
    TPM_paths = [file + f"03.Output/{sample}.TPM.txt" for file, sample in zip(Dir, Samples)]
    FPKM_paths = [file + f"03.Output/{sample}.FPKM.txt" for file, sample in zip(Dir, Samples)]

    while True:
        Check = 0
        for genecountfile in Genecount_paths:
            if os.path.exists(genecountfile):
                Check += 1
        if Check == len(Genecount_paths):
            break
        else:
            time.sleep(600)
    
    if sample == Samples[0]:
        data_frames = [pd.read_csv(file, sep='\t', header='infer') for file, sample in zip(Genecount_paths, Samples)]
        DATA = reduce(lambda left, right: pd.merge(left, right, on=['ID', 'GeneSymbol']), data_frames)
        DATA.to_csv(f"../Genecount/Total.Genecount.txt", sep='\t', header='infer', index=False)
        
        if TPM == 'Y':
            data_frames = [pd.read_csv(file, sep='\t', header='infer') for file, sample in zip(TPM_paths, Samples)]
            DATA = reduce(lambda left, right: pd.merge(left, right, on=['ID', 'GeneSymbol']), data_frames)
            DATA.to_csv(f"../Genecount/Total.TPM.txt", sep='\t', header='infer', index=False)

        if FPKM == 'Y':
            data_frames = [pd.read_csv(file, sep='\t', header='infer') for file, sample in zip(FPKM_paths, Samples)]
            DATA = reduce(lambda left, right: pd.merge(left, right, on=['ID', 'GeneSymbol']), data_frames)
            DATA.to_csv(f"../Genecount/Total.FPKM.txt", sep='\t', header='infer', index=False)
    else:
        pass
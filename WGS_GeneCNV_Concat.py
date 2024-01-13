#!/home/lab/anaconda3/envs/NGS/bin/python3

import pandas as pd
import glob
import numpy as np
#--------------------------------------------------#
def Merge(Name):
    CNV = glob.glob("../Intermediate/*txt")
    CNV.sort()

    DATA = None
    for data in CNV:
        Sample = data.split('/')[-1]
        Sample = Sample.split('.')[0]
        Data = pd.read_csv(data,
                        sep='\t',
                        header=None,
                        names=['Chr', 'Start', 'End', 'Gene', 'Strand', f'{Sample}'])
        
        if DATA is None:
            DATA = Data
        else:
            DATA = pd.merge(DATA, Data, on=['Chr', 'Start', 'End', 'Gene', 'Strand'])

    Info = DATA.iloc[:, :5]
    Coverage = DATA.iloc[:, 5:]

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
    CNV_ALL = CNV_ALL.applymap(np.log10)
    CNV_ALL = CNV_ALL.apply(lambda value : value - CNV_ALL.mean(axis=1))
    DATA = pd.concat([Info, CNV_ALL], axis=1)
    DATA = DATA[['Chr', 'Start', 'End', 'Gene', f'{Name}']]
    DATA.to_csv(f"05.SV/01.GeneCNV/{Name}.Norm.CNV.txt",
                sep='\t',
                index=False,
                header='infer')
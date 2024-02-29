import pandas as pd
#----------------------------------------------------------------------------------------#
Sample = pd.read_csv("SampleSheet.txt", sep="\t", header=None)
Name = Sample.iloc[0, 0]
R1 = Sample.iloc[0, 1]
R2 = Sample.iloc[0, 2]
Path = '/'.join(R1.split('/')[0:-1])
#----------------------------------------------------------------------------------------#
def DepthOfCoverage(name):
    DICT = {}
    with open(f'03.Align/{name}.all.depth.txt', 'r') as depth:
        for line in depth:
            line = line.strip()
            splitted = line.split('\t')
            DP = splitted[2]

            if DP not in DICT.keys():
                DICT[DP] = 1
            else:
                DICT[DP] += 1

    with open(f'03.Align/{name}.DepthofCoverage.txt', 'w') as note:
        for num in DICT.keys():
            note.write(str(num) + '\t' + str(DICT[num]) + '\n')
#----------------------------------------------------------------------------------------#
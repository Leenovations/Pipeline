import sys
#-------------------------------------------------------------------#
with open('/labmed/08.Imatinib_MachineLearning/node03/Y253H_03-1_S5/03.Output/Y253H_03-1_S5.haplotype.parsing.vcf', 'w') as note1:
    with open('/labmed/08.Imatinib_MachineLearning/node03/Y253H_03-1_S5/03.Output/Y253H_03-1_S5.haplotype.vcf', 'r') as haplo:
        for line in haplo:
            line = line.strip()
            if line.startswith('##'):
                note1.write(line + '\n')
            elif line.startswith('#CHROM'):
                splitted = line.split('\t')
                note1.write('\t'.join(splitted[0:7]) + '\t' + '\t'.join(splitted[8:10]) + '\n')
            else:
                splitted = line.split('\t')
                Info1 = splitted[7]
                Info2 = splitted[8]
                Info3 = splitted[9]
                note1.write('\t'.join(splitted[0:7]) + '\t' + Info2 + '\t' + Info3 + '\n')

with open('/labmed/08.Imatinib_MachineLearning/node03/Y253H_03-1_S5/03.Output/Y253H_03-1_S5.varscan.parsing.vcf', 'w') as note1:
    with open('/labmed/08.Imatinib_MachineLearning/node03/Y253H_03-1_S5/03.Output/Y253H_03-1_S5.varscan.vcf', 'r') as varscan:
        for line in varscan:
            line = line.strip()
            if line.startswith('##'):
                note1.write(line + '\n')
            elif line.startswith('#CHROM'):
                line = line.replace('Sample1', 'Y253H_03-1_S5')
                splitted = line.split('\t')
                note1.write('\t'.join(splitted[0:7]) + '\t' + '\t'.join(splitted[8:10]) + '\n')
            else:
                splitted = line.split('\t')
                Info1 = splitted[7]
                Info2 = splitted[8]
                Info3 = splitted[9]
                note1.write('\t'.join(splitted[0:7]) + '\t' + Info2 + '\t' + Info3 + '\n')
#!/home/ltk/anaconda3/envs/NGS/bin/python3

import subprocess
import os
#--------------------------------------------------------------------------------#
if os.path.exists("SRR622461.DepthPerGene.txt"):
    command = "rm -rf SRR622461.DepthPerGene.txt"
    os.system(command)
else:
    pass
#--------------------------------------------------------------------------------#
command = f"samtools depth ../00.Practice/SRR622461.bam > SRR622461.depth.txt"
os.system(command)

command = "awk '{ sum += $3 } END { if (NR > 0) print sum / NR }' SRR622461.depth.txt"
os.system(command)
#--------------------------------------------------------------------------------#
command = f"samtools coverage ../00.Practice/SRR622461.bam -o SRR622461.coverage.txt"
os.system(command)
#--------------------------------------------------------------------------------#
with open('/media/src/hg19/08.bed/whole.exome.gene.bed', 'r') as bed:
    for line in bed:
        line = line.strip()
        splitted = line.split('\t')
        Chr = splitted[0]
        Start = splitted[1]
        End = splitted[2]
        Gene = splitted[3]

        Position = Chr + ':' + Start + '-' + End
        
        result = subprocess.run(["samtools", "coverage", "-r", Position, "../00.Practice/SRR622461.bam"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        coverage = result.stdout.splitlines()[1].strip()

        with open("SRR622461.DepthPerGene.txt", "a") as note:
            note.write((f"{Gene}\t{coverage}\n"))
#--------------------------------------------------------------------------------#
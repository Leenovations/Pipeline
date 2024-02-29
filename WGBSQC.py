#!/home/lab/anaconda3/envs/NGS/bin/python3

import os
import pandas as pd
from fpdf import FPDF
from datetime import datetime
from pdf2image import convert_from_path
#----------------------------------------------------------------------------------------#
def QCPDF(name):
    Date = datetime.now()

    if os.path.isdir(f'00.PreQC/{name}_R1_001_fastqc'):
        pass
    else:
        command = f'unzip 00.PreQC/{name}_R1_001_fastqc.zip \
                    -d 00.PreQC/'
        os.system(command)

    if os.path.isdir(f'01.PostQC/{name}_val_1_fastqc'):
        pass
    else:
        command = f'unzip 01.PostQC/{name}_val_1_fastqc.zip \
                    -d 01.PostQC/'
        os.system(command)

    pdf = FPDF()
        
# Adding a page
    pdf.add_page()
        
# set style and size of font 
    pdf.set_font("Arial", size = 15)

# set center location of x (210 : A4)
    Center_x = pdf.w / 2 # 105
        
# create a cell
    pdf.set_fill_color(r=100, g=100, b=100)
    pdf.line(Center_x-95, 20, Center_x+95, 20)
    pdf.text(20, 11, txt = 'WGBS QC')

#Sample Name
    pdf.set_font("Arial", style = 'B', size = 10)
    pdf.text(20, 17, txt = f"Sample : {name}")

#Date
    pdf.set_font("Arial", style = 'B', size = 10)
    pdf.text(170, 17, txt = Date.strftime('%Y-%m-%d'))

#Category 1
    pdf.set_font("Arial", style = 'B', size = 12)
    pdf.text(20, 27, txt = '1. FASTQ Statistics - PreQC')

#Category1 Image
    pdf.set_font("Arial", size = 9)

    pdf.text(26, 35, txt = "Per base quality : Untrimmed")
    pdf.image(f'00.PreQC/{name}_R1_001_fastqc/Images/per_base_quality.png', x = 20, y = 36, w = 50, h = 50, type = 'PNG')

    pdf.text(96, 35, txt = "Sequence length")
    pdf.image(f'00.PreQC/{name}_R1_001_fastqc/Images/per_sequence_quality.png', x = 80, y = 36, w = 50, h = 50, type = 'PNG')

    pdf.text(147, 35, txt = "Per base sequence content")
    pdf.image(f'00.PreQC/{name}_R1_001_fastqc/Images/per_base_sequence_content.png', x = 140, y = 36, w = 50, h = 50, type = 'PNG')

#Category 2
    pdf.set_font("Arial",style = 'B', size = 12)
    pdf.text(20, 100, txt = '2. FASTQ Statistics - PostQC')

#Category2 Image
    pdf.set_font("Arial", size = 9)

    pdf.text(26, 108, txt = "Per base quality : Trimmed")
    pdf.image(f'01.PostQC/{name}_val_1_fastqc/Images/per_base_quality.png', x = 20, y = 110, w = 50, h = 50, type = 'PNG')

    pdf.text(96, 108, txt = "Sequence length")
    pdf.image(f'01.PostQC/{name}_val_1_fastqc/Images/per_sequence_quality.png', x = 80, y = 110, w = 50, h = 50, type = 'PNG')

    pdf.text(147, 108, txt = "Per base sequence content")
    pdf.image(f'01.PostQC/{name}_val_1_fastqc/Images/per_base_sequence_content.png', x = 140, y = 110, w = 50, h = 50, type = 'PNG')

#Category 3 
    pdf.set_font("Arial",style = 'B', size = 12)
    pdf.text(20, 170, txt = '3. Sequencing Statistics')

    Info = {}
    with open(f'04.QC/{name}.stats', 'r') as handle:
        for line in handle:
            if line.startswith('SN'):
                line = line.strip()
                splitted = line.split('\t')
                Info[splitted[1][:-1]] = splitted[2]
                
    Total_base = int(Info['bases mapped'])
    Total_read = int(Info['raw total sequences'])
    Mapped_read = int(Info['reads mapped'])
    Average_length = Info['average length']
    Average_qual = Info['average quality']
    Average_insert = Info['insert size average']

    Target_Info = {}
    with open(f'04.QC/{name}.Ontarget.stats', 'r') as handle:
        for line in handle:
            if line.startswith('SN'):
                line = line.strip()
                splitted = line.split('\t')
                Target_Info[splitted[1][:-1]] = splitted[2]
                
    Target_Mapped_read = int(Target_Info['reads mapped'])
    
    Percent = str(round(int(Mapped_read)/int(Total_read)*100))
    Ontarget = str(round(int(Target_Mapped_read)/int(Mapped_read)*100))

    pdf.set_font("Arial", size = 11)
    pdf.set_xy(20, 175)
    pdf.set_fill_color(r = 150, g = 150, b = 150)
    pdf.cell(57,10, txt = 'Total base', align = 'C', border=1, ln=0, fill = True)
    pdf.cell(57,10, txt = 'Total read', align = 'C', border=1, fill = True)
    pdf.cell(57,10, txt = 'Mapped read', align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(20, 185)
    pdf.cell(57,10, txt = str(format(Total_base, ',')), align = 'C', border=1)
    pdf.cell(57,10, txt = str(format(Total_read, ',')), align = 'C', border=1)
    pdf.cell(57,10, txt = "{:,} ({} %)".format(Mapped_read, Percent), align = 'C', border=1)

    pdf.set_xy(20, 195)
    pdf.set_fill_color(r = 150, g = 150, b = 150)
    pdf.cell(57,10, txt = 'Average Length', align = 'C', border=1, ln=0, fill = True)
    pdf.cell(57,10, txt = 'Average Insert size', align = 'C', border=1, fill = True)
    pdf.cell(57,10, txt = 'Average Quality', align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(20, 205)
    pdf.cell(57,10, txt = Average_length, align = 'C', border=1)
    pdf.cell(57,10, txt = Average_insert, align = 'C', border=1)
    pdf.cell(57,10, txt = Average_qual, align = 'C', border=1)

    # pdf.set_xy(20, 215)
    # pdf.set_fill_color(r = 150, g = 150, b = 150)
    # pdf.cell(171,10, txt = 'Ontarget %', align = 'C', border=1, fill = True)
    # pdf.set_xy(20, 225)
    # pdf.cell(171,10, txt = Ontarget + ' %', align = 'C', border=1)

#Category 4 
    pdf.set_font("Arial", style = 'B', size = 12)
    pdf.text(20, 225, txt = '4. STAR Alignment Statistics')

    Info = {}
    with open(f'04.QC/{name}_Log.final.out', 'r') as handle:
        for line in handle:
            line = line.strip()
            if '%' in line:
                line = line.replace('\t', '')
                splitted = line.split('|')
                Info[splitted[0]] = splitted[1]

    pdf.set_font("Arial", size = 11)
    pdf.set_xy(20, 230)
    pdf.set_fill_color(r = 150, g = 150, b = 150)
    pdf.cell(57,10, txt = list(Info.keys())[0], align = 'C', border=1, ln=0, fill = True)
    pdf.cell(57,10, txt = list(Info.keys())[1].replace(',',''), align = 'C', border=1, ln=0, fill = True)
    pdf.cell(57,10, txt = "Chimeric reads %", align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(20, 240)
    pdf.cell(57,10, txt = Info[list(Info.keys())[0]], align = 'C', border=1)
    pdf.cell(57,10, txt = Info[list(Info.keys())[1]], align = 'C', border=1)
    pdf.cell(57,10, txt = Info[list(Info.keys())[-1]], align = 'C', border=1)

    pdf.set_xy(20, 250)
    pdf.set_fill_color(r = 150, g = 150, b = 150)
    pdf.cell(57,10, txt = "Multiple loci mapped reads %", align = 'C', border=1, ln=0, fill = True)
    pdf.cell(57,10, txt = "Too many loci mapped reads %", align = 'C', border=1, ln=0, fill = True)
    pdf.cell(57,10, txt = "Too short unmapped reads %", align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(20, 260)
    pdf.cell(57,10, txt = Info[list(Info.keys())[4]], align = 'C', border=1)
    pdf.cell(57,10, txt = Info[list(Info.keys())[5]], align = 'C', border=1)
    pdf.cell(57,10, txt = Info[list(Info.keys())[6]], align = 'C', border=1)

# save the pdf
    pdf.output(f"04.QC/{name}.QC.pdf")
#----------------------------------------------------------------------------------------#
def pdfconverter(name):
    images = convert_from_path(f"04.QC/{name}.QC.pdf")

    for page, image in enumerate(images):
        image.save(f"04.QC/{name}.QC.page{page}.jpg", "JPEG")
#----------------------------------------------------------------------------------------#
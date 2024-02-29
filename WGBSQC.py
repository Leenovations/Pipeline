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

    pdf.image(f'03.Align/QCplot.png', x = 20, y = 180, w = 70, h = 70, type = 'PNG')

    Info = {}
    with open(f'03.Align/{name}.stats', 'r') as handle:
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
    Percent = str(round(int(Mapped_read)/int(Total_read)*100))

    pdf.set_font("Arial", size = 9)
    pdf.set_xy(100, 180)
    pdf.set_fill_color(r = 150, g = 150, b = 150)
    pdf.cell(30,7, txt = 'Total base', align = 'C', border=1, ln=0, fill = True)
    pdf.cell(30,7, txt = 'Total read', align = 'C', border=1, fill = True)
    pdf.cell(30,7, txt = 'Mapped read', align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(100, 187)
    pdf.cell(30,7, txt = str(format(Total_base, ',')), align = 'C', border=1)
    pdf.cell(30,7, txt = str(format(Total_read, ',')), align = 'C', border=1)
    pdf.cell(30,7, txt = "{:,} ({} %)".format(Mapped_read, Percent), align = 'C', border=1)

    pdf.set_xy(100, 194)
    pdf.set_fill_color(r = 150, g = 150, b = 150)
    pdf.cell(30,7, txt = 'Average Length', align = 'C', border=1, ln=0, fill = True)
    pdf.cell(30,7, txt = 'Average Insert size', align = 'C', border=1, fill = True)
    pdf.cell(30,7, txt = 'Average Quality', align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(100, 201)
    pdf.cell(30,7, txt = Average_length, align = 'C', border=1)
    pdf.cell(30,7, txt = Average_insert, align = 'C', border=1)
    pdf.cell(30,7, txt = Average_qual, align = 'C', border=1)

    pdf.set_xy(100, 208)
    pdf.set_fill_color(r = 112, g = 128, b = 144)
    pdf.cell(30,7, txt = 'All methylated', align = 'C', border=1, ln=0, fill = True)
    pdf.cell(30,7, txt = 'All unmethylated ', align = 'C', border=1, fill = True)
    pdf.cell(30,7, txt = 'Mixed methylation', align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(100, 215)
    pdf.cell(30,7, txt = '84.41 %', align = 'C', border=1)
    pdf.cell(30,7, txt = '0.00 %', align = 'C', border=1)
    pdf.cell(30,7, txt = '0.45 %', align = 'C', border=1)

    pdf.set_xy(100, 222)
    pdf.set_fill_color(r = 112, g = 128, b = 144)
    pdf.cell(90,7, txt = 'Duplication rate', align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(100, 229)
    pdf.cell(90,7, txt = '74.57 %', align = 'C', border=1)

    pdf.set_xy(100, 236)
    pdf.set_fill_color(r = 112, g = 128, b = 144)
    pdf.cell(45,7, txt = 'Mapping efficiency', align = 'C', border=1, ln=0, fill = True)
    pdf.cell(45,7, txt = 'Conversion rate', align = 'C', border=1, ln=0, fill = True)
    pdf.set_xy(100, 243)
    pdf.cell(45,7, txt = '99.04 %', align = 'C', border=1)
    pdf.cell(45,7, txt = '99 %', align = 'C', border=1)

    pdf.set_fill_color(r=100, g=100, b=100)
    pdf.line(Center_x-95, pdf.h-20, Center_x+95, pdf.h-20)


# save the pdf
    pdf.output(f"03.Align/{name}.QC.pdf")
#----------------------------------------------------------------------------------------#
def pdfconverter(name):
    images = convert_from_path(f"03.Align/{name}.QC.pdf")

    for page, image in enumerate(images):
        image.save(f"03.Align/{name}.QC.page{page}.jpg", "JPEG")
#----------------------------------------------------------------------------------------#
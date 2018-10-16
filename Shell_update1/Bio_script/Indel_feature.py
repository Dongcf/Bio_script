# coding: utf-8
#Programming assessment of vcf indel
#author:Dongchunfa 

#Import python library
import os
import re
import numpy as np
from optparse import OptionParser


usage = "python %prog -i vcf_file_path  -o outfile"
version = "0.0.1"
parser = OptionParser(usage = usage,version = version)
parser.add_option("-i","--inp",dest = "inp", help = "vcf file storage path")
parser.add_option("-o","--oup",dest = "oup",help = "outdir",default = ".")
(options, args) = parser.parse_args()
if options.inp == None:
    parser.error(" -i Parameters cannot be empty,must be a vcf file path")


#Get the file path
path = options.inp
out_file = options.oup

def Cont_Indel(file,file_type):
    Num_indel,Num_za,Num_chun,Num_dbsnp,Num_novel = np.zeros(5,dtype=np.int64)
    with open (path+'/'+file) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split('\t')
            if line[0] == 'Chr':
                line1 = line
                for i in line1:
                    if re.findall("snp",i):
                        snp_index = line1.index(i)
            if line[0][0:3] == 'chr':
                #if line.split()[-3].split(';')[0] == 'INDEL':
                Num_indel += 1
                if line[-1].split(" ")[-1].split(":")[0] in ('1/1','0/0'):
                    Num_chun += 1
                else:
                    Num_za += 1
                if line[snp_index].startswith('rs'):
                    Num_dbsnp += 1
                else:
                    Num_novel += 1
        dbSNP_per = float(Num_dbsnp)/Num_indel
        if re.findall("somatic",file):
            sample = '_'.join(file.split('.')[0].split("_")[:-2])
        else:
            sample = '_'.join(file.split('.')[0].split("_")[:-1])
        OUT_file.write(tplt.format(sample,Num_indel, Num_chun,Num_za,[Num_dbsnp,dbSNP_per],Num_novel))

OUT_file = open(out_file+"/indel.feature.xls",'w')
tplt1 = "{0:^}\t{1:^}\t{2:^}\t{3:^}\t{4:^}\t{5:^}\n"
tplt = "{0:^}\t{1:^,}\t{2:^,} \t{3:^,} \t{4[0]:^,}  ({4[1]:^.2%})\t{5:^} \n"
OUT_file.write(tplt1.format("Sample","Total","Homozygote","Heterozygote","dbSNP_percentage","Novel"))
for file in os.listdir(path):
    if file.split('.')[0].split('_')[-1]  == 'indel' and file.split(".")[-1] == 'xls':
        Type1 = 'germline'
        Cont_Indel(file,Type1)
    
    if  re.findall("somatic_indel",file) and file.split(".")[-1] == 'xls':
        Type2 = 'somatic'
        Cont_Indel(file,Type2)
            

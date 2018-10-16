# coding: utf-8
#Programming assessment of vcf indel
#author:Dongchunfa

#Import python library
import os
import re
import numpy as np
from optparse import OptionParser

usage = "python %prog -i vcf_file_path  -o outfile_path"
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
def Cont_Snp(file,file_type):
    
#Initialization parameters
    num_snp,num_chun,num_za,num_ts,num_tv,num_dbSNP,num_novel,novel_ts,novel_tv = np.zeros(9,dtype=np.int64)
    with open (path+'/'+file) as f:
        lines = f.readlines()
        for line in lines:
             line = line.strip().split('\t')
             if line[0] == 'Chr':
                line1 = line
                ref_index = line1.index("Ref")
                Alt_index = line1.index("Alt")
                for i in line1:
                    if re.findall("snp",i):
                        snp_index = line1.index(i)
             if line[0][0:3] == "chr":
                num_snp += 1
                if line[-1].split(" ")[-1].split(":")[0] in ('1/1','0/0'):
                    num_chun += 1
                else:
                    num_za += 1
                if (line[ref_index],line[Alt_index]) in (('A','G'),('G','A'),('C','T'),('T','C')):
                    num_ts += 1
                else:
                    num_tv += 1
                if not line[snp_index].startswith('rs'):
                    num_novel += 1
                    if (line[ref_index],line[Alt_index]) in (('A','G'),('G','A'),('C','T'),('T','C')):
                        novel_ts += 1
                    else:
                        novel_tv += 1 
                else:
                        num_dbSNP += 1
                          
        tstv_percent = float(num_ts)/num_tv
        dbSNP_percent = float(num_dbSNP)/num_snp
        novel_tstv_percent = float(novel_ts)/novel_tv
        if re.findall("somatic",file):
            sample = '_'.join(file.split('.')[0].split("_")[:-2])
        else:
            sample = '_'.join(file.split('.')[0].split("_")[:-1])

        OUT_file.write(tplt.format(sample,num_snp,num_chun,num_za,[num_dbSNP,dbSNP_percent],num_ts,num_tv,tstv_percent,num_novel,novel_ts,novel_tv,novel_tstv_percent))

OUT_file = open(out_file+"/snp.feature.xls",'w')
tplt1 = "{0:^}\t{1:^}\t{2:^}\t{3:^}\t{4:^}\t{5:^}\t{6:^}\t{7:^}\t{8:^}\t{9:^}\t{10:^}\t{11:^}\n"
tplt = "{0:^}\t{1:^,}\t{2:^,} \t{3:^,} \t{4[0]:,}  ({4[1]:.2%})\t{5:^,}\t{6:^,}\t{7:^.2f}\t{8:^,}\t{9:^,}\t{10:^,}\t{11:^.2f} \n"
OUT_file.write(tplt1.format("Sample","Total","Homozygote","Heterozygote","dbSNP_percentage","Ts","Tv","Ts/Tv","Novel","Novel_Ts","Novel_Tv","Novel_Ts/Tv"))
for file in os.listdir(path):
    if file.split('.')[0].split('_')[-1] == 'snp' and file.split(".")[-1] == 'xls':
        Type1 = 'germline'
        Cont_Snp(file,Type1)
    
    if  re.findall("somatic_snp",file) and file.split(".")[-1] == 'xls':
        Type2 = 'somatic'
        Cont_Snp(file,Type2)

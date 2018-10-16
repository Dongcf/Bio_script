# coding: utf-8
#Programming assessment of vcf snp
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
def Cont_Snvmodel(file,file_type):
    tplt = "{0:^}\t{1:^10,}\t{2[0]:,}   ({2[1]:.2%})\t{3[0]:,}   ({3[1]:.2%})\t{4[0]:,}   ({4[1]:.2%})\t{5[0]:,}   ({5[1]:.2%})\t{6[0]:,}   ({6[1]:.2%})\t{7[0]:,}   ({7[1]:.2%})\n"
    #Initialization parameters
    num_CAGT,num_CGGC,num_CTGA,num_TAAT,num_TCAG,num_TGAC,num_snp = np.zeros(7,dtype=np.int64)
    with open (path+'/'+file) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip().split('\t')
            if line[0] == 'Chr':
                line1 = line
                ref_index = line1.index("Ref")
                Alt_index = line1.index("Alt")
            if line[0][0:3] == "chr":
                num_snp += 1
                if (line[ref_index],line[Alt_index]) in (('C','A'),('G','T')):                                                                    
                    num_CAGT += 1
                elif (line[ref_index],line[Alt_index]) in (('C','G'),('G','C')):                              
                    num_CGGC += 1
                elif (line[ref_index],line[Alt_index]) in (('C','T'),('G','A')):
                    num_CTGA +=1
                elif (line[ref_index],line[Alt_index]) in (('T','A'),('A','T')):
                    num_TAAT += 1
                elif (line[ref_index],line[Alt_index]) in (('T','C'),('A','G')):
                    num_TCAG += 1
                else:
                    num_TGAC += 1
        data = np.array((num_CAGT,num_CGGC,num_CTGA,num_TAAT,num_TCAG,num_TGAC),dtype = np.float32)
        data = data/num_snp
        if re.findall("somatic",file):
            sample = '_'.join(file.split('.')[0].split("_")[:-2])
        else:
            sample = '_'.join(file.split('.')[0].split("_")[:-1])

        OUT_file.write(tplt.format(sample,num_snp,[num_CAGT,data[0]],[num_CGGC,data[1]],[num_CTGA,data[2]],[num_TAAT,data[3]],[num_TCAG,data[4]],[num_TGAC,data[5]]))
        

OUT_file = open(out_file+"/snp_mutation_spectrum.xls",'w')
tplt1 = "{0:^}\t{1:^}\t{2:^}\t{3:^}\t{4:^}\t{5:^}\t{6:^}\t{7:^}\n"
OUT_file.write(tplt1.format("Sample","Total","C>A/G>T","C>G/G>C", "C>T/G>A", "T>A/A>T"," T>C/A>G"," T>G/A>C"))
for file in os.listdir(path):
    if file.split('.')[0].split('_')[-1] == 'snp' and (not re.findall("somatic_snp",file)) and file.split(".")[-1] == 'xls' and re.findall("anno.vcf",file):
        Type1 = 'germline'
        Cont_Snvmodel(file,Type1)
    
    if  re.findall("somatic_snp",file) and file.split(".")[-1] == 'xls' and re.findall("anno.vcf",file):
        Type2 = 'somatic'
        Cont_Snvmodel(file,Type2)

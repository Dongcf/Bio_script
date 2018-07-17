#!/usr/bin/env python 
# -*- coding: utf-8 -*- 

import os 
import shutil
import sys
import re
from optparse import OptionParser 
import json
path = os.path.abspath(os.path.dirname(__file__))
#python script/qc_report.py -i json.list -t template/src/QC_report/QC_report_temple.HTML -o ./ -n 6 -s 老鼠 -n 6 -g 10
def parseCommand():
	usage = "python qc_report.py -i json.list -t qc_template -o outfile"
	version = "0.0.1"
	A = "wenger"
	parser = OptionParser(usage = usage,version = version)
	parser.add_option("-i","--inp",dest = "inp", help = "fastp outut dir")
	parser.add_option("-o","--oup",dest = "oup",help = "outdir",default = ".")
	parser.add_option("-t","--tem",dest = "tem",help = "template",default="%s/template/src/QC_report/QC_report_temple.HTML" % (path))
	parser.add_option("-s","--specise",dest = "species",help = "the report species",default = "人")
	parser.add_option("-l","--sampletype",dest = "sampletype",help = "sample type",default = "")
	parser.add_option("-n","--num",dest = "num",help = "sample num")
	parser.add_option("-g","--numg",dest = "numg",help = "Contract requirement volume of data",default=10,type="int")
	parser.add_option("-r","--reportdir",dest = "reportdir",default = "%s/template/src/QC_report/files/" % (path))
	parser.add_option("-x","--xm",dest="xm",default ="转录组",help = "Analysis type")
	parser.add_option("-R","--rawq",dest = "rawq",default = "N",choices = ['Y','N'],help = "output rawdata q20 q30 or not")
	return parser.parse_args()
(options, args) = parseCommand()

def logwrite(content,output):
	with open("%s" % (output+'/qc.log'),"a") as file:
		file.write(content)

def parse_json(infile,iid,samg,output):
	name = infile.split("/")[-1].split(".")[0]
	mean1_sum = 0
	mean2_sum = 0
	GC1_sum = 0
	GC2_sum = 0
	gc1 = 0
	gc2 = 0
	m1 = 0
	m2 = 0
	with open(infile,"r") as file:
		json_obj = json.load(file)
		logname = infile.split("/")[-1].split(".")[0]
		Raw_Reads =int(json_obj["summary"]["before_filtering"]["total_bases"])/(1000**3)
		if Raw_Reads < samg :
			logstr = "Sample %s less than %s\n" % (logname,samg)
			logwrite(logstr,output)
		Clean_Reads =int(json_obj["summary"]["after_filtering"]["total_bases"])/(1000**3)
		Q20 =float('%.2f' %( (int(json_obj["read1_after_filtering"]["q20_bases"])+int(json_obj["read2_after_filtering"]["q20_bases"]))/int(json_obj["summary"]["after_filtering"]["total_bases"])*100))
		Q30 =float('%.2f' %( (int(json_obj["read1_after_filtering"]["q30_bases"])+int(json_obj["read2_after_filtering"]["q30_bases"]))/int(json_obj["summary"]["after_filtering"]["total_bases"])*100))
		RQ20 = float('%.2f' %( (int(json_obj["read1_before_filtering"]["q20_bases"])+int(json_obj["read2_before_filtering"]["q20_bases"]))/int(json_obj["summary"]["before_filtering"]["total_bases"])*100))
		RQ30 =float('%.2f' %( (int(json_obj["read1_before_filtering"]["q30_bases"])+int(json_obj["read2_before_filtering"]["q30_bases"]))/int(json_obj["summary"]["before_filtering"]["total_bases"])*100))
		mean1 = json_obj["read1_after_filtering"]["quality_curves"]["mean"]
		for i in mean1:
			m1 += 1
			mean1_sum += float(i)
		#mean1_sum = mean1_sum/m1
		GC1 = json_obj["read1_after_filtering"]["content_curves"]["GC"]
		for i in GC1:
			gc1 +=1
			GC1_sum += float(i)
		GC1_sum  = GC1_sum/gc1
		mean2 = json_obj["read2_after_filtering"]["quality_curves"]["mean"]
		for i in mean2:
			m2 += 1
			mean2_sum += float(i)
		#mean2_sum = mean2_sum/m2
		GC2 = json_obj["read2_after_filtering"]["content_curves"]["GC"]
		for i in GC2:
			gc2 += 1
			GC2_sum += float(i)
		GC2_sum = GC2_sum/gc2
		mean = (mean1_sum+mean2_sum)/(len(mean1)+len(mean2))
		Error_rate = float('%.2f' % (10**(-mean/10)*100))
		GC_Content = float('%.2f' % ((GC1_sum+GC2_sum)/2*100) )
		effisive = float("%.2f" % ((Clean_Reads/Raw_Reads)*100))
		Raw_Reads =float("%.2f" % (Raw_Reads))
		Clean_Reads =float("%.2f" % (Clean_Reads))
		out = {'Sample':name,'Raw Data(G)':Raw_Reads,'Clean Data(G)':Clean_Reads,'Clean Q20(%)':Q20,'Clean Q30(%)':Q30,'Error rate(%)': Error_rate,"GC Content(%)":GC_Content}
		out2 = {'Sample':name,'Raw Data(G)':Raw_Reads,'Clean Data(G)':Clean_Reads,"Effective(%)":effisive,'Clean Q20(%)':Q20,'Clean Q30(%)':Q30,'Error rate(%)': Error_rate,"GC Content(%)":GC_Content}
		Rout2 ={'Sample':name,'Raw Data(G)':Raw_Reads,'Clean Data(G)':Clean_Reads,"Effective(%)":effisive,'Raw Q20(%)': RQ20,'Clean Q20(%)':Q20,'Raw Q30(%)': RQ30,'Clean Q30(%)':Q30,'Error rate(%)': Error_rate,"GC Content(%)":GC_Content}
		vardata =str("<tr class=\"nobgcolor\"><td><div id='qual_read")+str(iid)+str("' style='width:500px;height: 400px'></div><script type='text/javascript'>")+str("var data=[{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["quality_curves"]['A'])+str(",name: 'A',mode:'lines',line:{color:'rgba(128,128,0,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["quality_curves"]["T"])+str(",name: 'T',mode:'lines',line:{color:'rgba(128,0,128,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["quality_curves"]['C'])+str(",name: 'C',mode:'lines',line:{color:'rgba(0,255,0,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["quality_curves"]['G'])+str(",name: 'G',mode:'lines',line:{color:'rgba(0,0,255,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["quality_curves"]['mean'])+str(",name: 'mean',mode:'lines',line:{color:'rgba(20,20,20,1.0)', width:1}},")+str("];var layout={title:'")+str(name)+str(" Read1 过滤后的质量分布图', xaxis:{title:'cycles'}, yaxis:{title:'quality'}};Plotly.newPlot('qual_read")+str(iid)+str("', data, layout);</script></td><td>")+str("<div id='qual_readtwo")+str(iid)+str("' style='width:500px;height: 400px'></div><script type='text/javascript'>var data=[{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read2_after_filtering"]["quality_curves"]['A'])+str(",name: 'A',mode:'lines',line:{color:'rgba(128,128,0,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["quality_curves"]["T"])+str(",name: 'T',mode:'lines',line:{color:'rgba(128,0,128,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read2_after_filtering"]["quality_curves"]['C'])+str(",name: 'C',mode:'lines',line:{color:'rgba(0,255,0,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read2_after_filtering"]["quality_curves"]['G'])+str(",name: 'G',mode:'lines',line:{color:'rgba(0,0,255,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read2_after_filtering"]["quality_curves"]['mean'])+str(",name: 'mean',mode:'lines',line:{color:'rgba(20,20,20,1.0)', width:1}},")+str("];var layout={title:'")+str(name)+str(" Read2 过滤后的质量分布图', xaxis:{title:'cycles'}, yaxis:{title:'quality'}};Plotly.newPlot('qual_readtwo")+str(iid)+str("', data, layout);</script></td></tr>")
		vardata2 = str("<tr class=\"nobgcolor\"><td><div id='content_read")+str(iid)+str("' style='width:500px;height: 400px'></div><script type='text/javascript'>")+str("var data=[{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["content_curves"]["A"])+str(",name: 'A',mode:'lines',line:{color:'rgba(128,128,0,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["content_curves"]["T"])+str(",name: 'T',mode:'lines',line:{color:'rgba(128,0,128,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["content_curves"]["C"])+str(",name:'C',mode:'lines',line:{color:'rgba(0,255,0,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["content_curves"]["G"])+str(",name: 'G',mode:'lines',line:{color:'rgba(0,0,255,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["content_curves"]['N'])+str(",name: 'N',mode:'lines',line:{color:'rgba(255, 0, 0, 1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read1_after_filtering"]["content_curves"]['GC'])+str(",name: 'GC',mode:'lines',line:{color:'rgba(20,20,20,1.0)', width:1}},")+str("];var layout={title:'")+str(name)+str(" Read1 过滤后的含量分布图', xaxis:{title:'cycles'},yaxis:{title:'quality'}};Plotly.newPlot('content_read")+str(iid)+str("', data, layout);</script></td><td>")+str("<div id='content_readtwo")+str(iid)+str("' style='width:500px;height: 400px'></div><script type='text/javascript'>var data=[{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read2_after_filtering"]["content_curves"]['A'])+str(",name:'A',mode:'lines',line:{color:'rgba(128,128,0,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read2_after_filtering"]["content_curves"]["T"])+str(",name: 'T',mode:'lines',line:{color:'rgba(128,0,128,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read2_after_filtering"]["content_curves"]['C'])+str(",name: 'C',mode:'lines',line:{color:'rgba(0,255,0,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read2_after_filtering"]["content_curves"]['G'])+str(",name: 'G',mode:'lines',line:{color:'rgba(0,0,255,1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read2_after_filtering"]["content_curves"]["N"])+str(",name: 'N',mode:'lines',line:{color:'rgba(255, 0, 0, 1.0)', width:1}},")+str("{x:")+str(list(range(1,152)))+str(",y:")+str(json_obj["read2_after_filtering"]["content_curves"]['GC'])+str(",name: 'GC',mode:'lines',line:{color:'rgba(20,20,20,1.0)', width:1}},")+str("];var layout={title:'")+str(name)+str(" Read2 过滤后的含量分布图',xaxis:{title:'cycles'}, yaxis:{title:'quality'}};Plotly.newPlot('content_readtwo")+str(iid)+str("',data, layout);</script></td></tr>")
		return "%s,\n" % (out),"%s,\n" % (out2),"%s,\n" % (Rout2),vardata,vardata2

def inlist(indir,out,samg,output):
	oo1=[]
	oo2=[]
	oo3=[]
	vv1=[]
	vv2=[]
	html = []
	i = 0
	j = -1
	shtml = ""
	snum = 0
	a= os.listdir(indir)
	for f in a:
		if re.match(".*json$",f):
			filepath =indir+f
			i = i+1
			o1,o2,o3,v1,v2 = parse_json(filepath,i,samg,output)
			oo1.append(o1)
			oo2.append(o2)
			oo3.append(o3)
			vv1.append(v1)
			vv2.append(v2)
	for f in a:	
		if re.match(".*html$",f):
			filepath =indir+f
			os.system("cp %s %s" % (filepath,out))
			name = f.split(".")[0]
			samh = str("<td><a href = \"single/")+f+str("\">")+name+str("</a></td>")
			html.append(samh)
	for td in html:
		j +=1
		if j % 4 == 0:
			shtml = shtml+str("</tr><tr>")+td
		else:
			shtml +=str(td)	
	snum = i
	return oo1,oo2,oo3,vv1,vv2,snum,shtml

def template(temp,in1,in2,in3,out,snum,samhtml):
	text = open(temp,"r").read()
	text = text.replace("{{qcinfo}}",in1)
	text = text.replace("{{base1}}",in2)
	text = text.replace("{{base2}}",in3)
	text = text.replace("{{XM}}",options.xm)
	text = text.replace("{{XXX}}",options.species)
	text = text.replace("{{sampletype}}",options.sampletype)
	text = text.replace("{{samplenum}}",str(snum))
	text = text.replace("{{sampleg}}",str(options.numg))
	text = text.replace("{{samplehtml}}",samhtml)
	outp = open(out+'/QC_report/QC_report.html',"w")
	outp.write(text)
	outp.close()
	

def main():
	out = options.oup
	if os.path.exists(out+'/qc.log'):
		os.remove(out+'/qc.log')
	elif os.path.exists(out+'/QC_report'):
		shutil.rmtree(out+'/QC_report')
	os.makedirs(out+'/'+'QC_report/single')
	os.system("cp -r %s %s" % (options.reportdir,out+'/'+'/QC_report/'))
	oo1,oo2,oo3,vv1,vv2,snum,sht = inlist(options.inp+"/",options.oup+"/"+"QC_report/single",options.numg,options.oup)
	#print("".join(oo1))
	if options.num:
		snum = options.num
	str_oo2 = "".join(oo2)
	str_vv1 = "".join(vv1)
	str_vv2 = "".join(vv2)
	str_oo3 = "".join(oo3)
	if options.rawq in ['Y']:
		template(options.tem,str_oo3,str_vv1,str_vv2,options.oup,snum,sht)
	else:
		template(options.tem,str_oo2,str_vv1,str_vv2,options.oup,snum,sht)
	if os.path.exists(out+'/qc.log'):
		print("Warning :same Samples data size below required\nSee %s/qc.log" % (options.oup))


if __name__  == "__main__":
	main()

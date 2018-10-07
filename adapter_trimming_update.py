#!/usr/bin/python

#Code contributed by Yuhua Zhang, bash code part credit to Alex Tsoi

import sys
import pandas as pd
import re
import glob
import os
import itertools	
import multiprocessing

def parse_file(filename):
	#Generate patient ID and description
	#input: core_info_file
	#output: patient ID and description(cell type etc)
	core_info=pd.read_table(filename)
	#tmp=core_info[core_info['Lane']==1]
	tmp=core_info.drop_duplicates(subset='SampleID',keep="last").reset_index()
	#tmp=core_info.groupby(['Lane']).get_group(1)
	tmp=tmp.drop(tmp.index[len(tmp)-1])

	#patient ID
	patient_ID=tmp['Description']
	patient_ID=pd.Series(patient_ID).str.extract('(\A\d*\d)',expand=True)
	
	cell_type=tmp['Description']
	cell_type=pd.Series(cell_type).str.extract('(CD4[a-zA-Z]+|CD8[a-zA-Z]+|mDC)',expand=True)
	cell_type=cell_type.apply(lambda x:x.astype(str).str.upper())
	cell_type=cell_type.apply(lambda x:x.astype(str).str.replace('MDC','mDC',case=False))

	hour=tmp['Description']
	hour=pd.Series(hour).str.extract('(\d*h)',expand=True)
	
	operation=tmp['Description']
	operation=pd.Series(operation).str.extract('(s[a-zA-Z]*|c[a-zA-Z]*)',expand=True)

	#special_operation=tmp['Description']
	#special_operation=pd.Series(special_operation).str.extract('([a-zA-Z]*[a-zA-Z]\Z)',expand=True)

	#Generate the description
	descp=pd.concat([cell_type,hour,operation],axis=1)
	descp=descp.apply(lambda x:x.astype(str).str.cat(sep='_'),axis=1)
	descp=descp.apply(lambda x:re.sub('_nan','',x))
	return(patient_ID,descp)

def get_lanes_samples(filename):
	#get sample number and lanes
	#input: core_info_file
	#output: lanes, sample numbers and sampleID
	core_info=pd.read_table(filename)
	lanes=len(core_info.groupby(['Lane']).count())
	samples=len(core_info.groupby(['SampleID']).count())-1
	#tmp=core_info[core_info['Lane']==1]
	tmp=core_info.drop_duplicates(subset='SampleID',keep="last").reset_index()
	tmp=tmp.drop(tmp.index[len(tmp)-1])
	sampleID=tmp['SampleID']
	return(lanes,samples,sampleID)

def generate_sample_info(filename):
	#combine with the information from parse_file and generate sample_info
	#input: core_info_file by user
	#output: sample info
	(patient_ID,Description)=parse_file(filename)
	(lanes,samples,sampleID)=get_lanes_samples(filename)
	tmp=re.search('(Run_\d*\d)',filename).group(1)
	tmp=re.sub('Run_','',tmp)
	batch=re.search('(Batch\d*\d)',filename).group(1)
	sample_info_name=batch+'_Run'+tmp
	tmp=sample_info_name
	batch_run=pd.DataFrame(columns=['A'])
	for i in range(samples):
		batch_run=batch_run.append({'A':tmp},ignore_index=True)
	sample_info=pd.concat([sampleID,batch_run,Description,patient_ID],axis=1,ignore_index=True)
	sample_info_name='SampleInfo_'+sample_info_name
	return(sample_info,sample_info_name)

def generate_trimmed_file_para(param):
	pathway=param[0]
	i=param[1]
	ele=param[2]
	Flag2=param[3]
	if Flag2:
		os.system('python /net/fantasia/home/alextsoi/Software/atactk/scripts/trim_adapters '
			+pathway+'/elder/Sample_'+ele+'/*L00'+str(i+1)+'*R1*001.reads_trimmed.fastq.gz '
			+pathway+'/elder/Sample_'+ele+'/*L00'+str(i+1)+'*R2*001.reads_trimmed.fastq.gz')
	else:
		os.system('python /net/fantasia/home/alextsoi/Software/atactk/scripts/trim_adapters '
			+pathway+'/elder/Sample_'+ele+'/*L00'+str(i+1)+'*R1*001.fastq.gz '
			+pathway+'/elder/Sample_'+ele+'/*L00'+str(i+1)+'*R2*001.fastq.gz')

def check_sampleID(filename,pathway):
	#check whether samples in the folder is consistant with the core_info
	#input: core_info_file, pathway provoded by user
	#output: a bool variable; if passed, true, vice versa
	(lanes,samples,sampleID)=get_lanes_samples(filename)
	Flag=True
	Samples=len(glob.glob(pathway+'/elder/Sample*'))
	#First compare number of samples
	if samples!=Samples:
		print("Numder of Samples is not correct")
		Flag=False

	else:
		print("Number of Samples is correct")
		#Next, compare sampleID
		SampleID=[]
		for file in os.listdir(pathway+'/elder'):
			if (re.search('Sample*',file)):
				SampleID.append(re.sub('Sample_','',file))
		SampleID.sort()
		if set(sampleID)!=set(SampleID):
			print("The SampleID doesn't match")
			Flag=False

		else:
			print("The SampleID matches")
			for ele in SampleID:
				count_=0
				for file in os.listdir(pathway+'/elder/Sample_'+str(ele)+'/'):
					if (re.search('\d.fastq.gz',file)):
						count_+=1
				#Compare Lanes
				if count_!=2*lanes:
					print("Number of lanes in Sample"+str(ele)+" doesn't match")
					Flag=False
					break
				else: 
					print('Number of lanes in Sample'+str(ele)+" matches")
	return(Flag)


def trim_reads(filename,reads,direct,pathway):
	#input: core_info file, the length of the reads to maintain
	(lanes,samples,sampleID)=adapter_trimming_update.get_lanes_samples(filename)
	core_info=pd.read_table(filename)
	for ele in sampleID:
		#tmp=core_info[core_info['SampleID'].str.match(ele)]
		#tmp_lane=tmp['Lane']
		print('for i in '+ele+'; do for j in `ls '+pathway+'/elder/Sample_${i}/*R1*.fastq.gz`;'
			+'do k=`echo $j | cut -d"/" -f11 `;  gunzip -c $j | /home/alextsoi/Software/bin/fastx_trimmer -'
			+direct+' '+str(reads)+' -z -o '+pathway+'/elder/$k.reads_trimmed.fastq.gz  ;done; done')

def concat_func(command):
	#read in the file specified by user
	#first check if the files match with core_info
	#then generate sample_info and adapter-trimmed file
	
	print('Processing adapter trimming...')
	pathway=command[0]
	filename=command[1]
	process_input=command[2]
	entire_output=command[3]
	Flag2=command[4] #Whether to trim the reads
	reads=command[5]
	direct=command[6]

	#First check
	Flag=check_sampleID(filename,pathway)
	
	if Flag:
		if Flag2:
			trim_reads(filename,reads,direct,pathway)
		(sample_info,sample_info_name)=generate_sample_info(filename)
		sample_info.to_csv(entire_output+'/'+sample_info_name,header=None,index=None,sep='\t')
		#call trim_adapter in parallel
		(lanes,sample,SampleID)=get_lanes_samples(filename)
		lane=range(lanes)
		path=[pathway]
		Flag2=[Flag2]
		paramlist=list(itertools.product(path,lane,SampleID,Flag2))
		process_input=int(process_input)
		pool=multiprocessing.Pool(processes=process_input)
		pool.map(generate_trimmed_file_para,paramlist)
		print('Adapter trimming done')

if __name__=="__main__":
	concat_func(command)

#!/usr/bin/python

#Code created by Alex Tsoi in bash script, modified by Yuhua Zhang to python

import re
import os
import sys

def generate_bed_unfiltered(pathway,outdir):
	#input: pathway to .bam files
	Flag=False
	command=[]
	for bam_file in os.listdir(pathway):
		if (re.search('Sample_',bam_file)):
			Flag=True
			break
	if Flag:
		for file in os.listdir(pathway):
			for sub_dir in os.listdir(pathway+'/'+file+'/bams'):
				if(re.search('.recal.bam\Z',sub_dir)):
					filename=sub_dir
					command.append('bedtools bamtobed -i '+pathway+'/'+file+'/bams/'+sub_dir+' > '+outdir+'/'+filename+'.bed')
	else:
		for file in os.listdir(pathway+'/bams'):
			if (re.search('.recal.bam\Z',file)):
				filename=file
				command.append('bedtools bamtobed -i '+pathway+'/bams/'+file+' > '+outdir+'/'+filename+'.bed')
	return(command)

def generate_bed_filtered(pathway,outdir):
	#input: pathway to filtered .bam files
	command=[]
	for file in os.listdir(pathway):
		if (re.search('.noGL.bam\Z',file)):
			filename=file
			command.append('bedtools bamtobed -i '+pathway+'/'+file+' > '+outdir+'/'+filename+'.bed')
	return(command)

def remove_MT(outdir):
	#input: dir to store the bed files
	command=[]
	for file in os.listdir(outdir):
		if (re.search('recal.bam.bed',file)):
			tmp_file=re.sub('.bed','',file)
			command.append('grep -v ^[MGH] '+outdir+'/'+file+' > '+outdir+'/'+tmp_file+'.noMT.bed')
	return(command)

def multi_processing1(command):
	command_=command[0]
	os.system(command_)

def multi_processing2(command):
	command_=command[0]
	os.system(command_)

def multi_processing3(command):
	command_=command[0]
	os.system(command_)

def concat_func(command):
	pathway1=command[0]
	pathway2=command[1]
	outdir=command[2]
	job=command[3]
	process_input=int(job)
	pool=multiprocessing.Pool(processes=process_input)
	print('Creating .bed files...')
	print('Creating .bed files for unfiltered .bam files')
	command1=generate_bed_unfiltered(pathway1,outdir)
	paramlist=list(itertools.product(command1))
	pool.map(multi_processing1,paramlist)
	
	print('Creating .bed files for filtered .bam files')
	command2=generate_bed_filtered(pathway2,outdir)
	paramlist=list(itertools.product(command2))
	pool.map(multi_processing2,paramlist)
	
	print('Creating remove_MT .bed files for unfiltered .bam files')
	command3=remove_MT(outdir)
	paramlist=list(itertools.product(command3))
	pool.map(multi_processing3,paramlist)

	print('.bed files have been created.')

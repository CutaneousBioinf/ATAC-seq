#!/usr/bin/python

#This part of the code concatenate all the rest part together, contributed by Yuhua Zhang

import adapter_trimming_update
import alignment_gotcloud
import filtering
import genome_coverage
import re
import os
import sys
import getopt
import make_bed
import signal_to_noise
import summary

'''if len(sys.argv)!=2:
	print('Usage: python command.py config_file')
	print('Refer to https://github.com/YuhuaZhang1995/ATACseq_pipeline for more help')
	sys.exit(2)'''
conf_file=''
def generate_template():
	file=open('config_file_template','a')
	file.write('--core_info_file <pahtway/../core_info_file>\n')
	file.write('--seq_data <pathway/.../Run_XXXX>\n')
	file.write('--batch <e.g. 14>\n')
	file.write('--run <e.g. 1789>\n')
	file.write('--conf <pathway/../config_file_for_gotcloud>\n')
	file.write('--bedprofile <pathway/../bedprofile>\n')
	file.write('--entire_output *<pathway/to/store/all_files>\n')
	file.write('--specific_output *<pathway/../specific_output_file>\n')
	file.write('--job_AT *<e.g. 10>\n')
	file.write('--job_align *<e.g. 3>\n')
	file.write('--job_filter *<e.g. 5>\n')
	file.write('--intermediate_file *<Yes>\n')
	file.write('(parameters with \"*\" are not required, but optional, you can delete them if not needed)\n')
	file.close()

	file=open('specific_output_template','a')
	file.write('--out_sampleinfo <outdir/../sample_info>\n')
	file.write('--out_bam <outdir/../BAM>\n')
	file.write('--out_proc_bam <outdir/../BAMprocessed>\n')
	file.write('--out_plot <outdir/../QCs/Insertsize>\n')
	file.write('--out_coverage <outdir/../QCs/breadth>\n')
	file.write('--out_bed <outdir/../BED>\n')
	file.write('--out_clus <outdir/../QCs/SignaltoNoise/Tss>\n')
	file.write('--out_s2n <outdir/../QCs/SignaltoNoise/Tss/analysis>\n')
	file.write('--out_summary <outdir/../QCs/Summary>\n')
	file.close()

	file=open('bedprofile_template','a')
	file.write('win_size=10		#Size of window for counting, min=1\n')
	file.write('step_size=10	#Step size for sliding window, min=1\n')
	file.write('regions=/net/assembly/alextsoi/Researches/Psoriasis_ATAC_seq_preliminary2/QCs/SignaltoNoise/data/HOMER_hg19_TSS_nochr.bed\n')
	file.write('norm_flag=0		#Normalize to read depth: 0=false|1=true\n')
	file.write('read_type=atac	#Read type: single|paired|atac\n')
	file.write('avg_frag_len=0	#Average length of fragments generated during library prep, If set to 0 will use read length\n')
	file.write('count_type=ends	#Count type should be either: ends | middle | overlap\n')
	file.write('chr_sizes=/net/fantasia/home/alextsoi/db/UCSC/hg19.chrom.sizes_nochr\n')
	file.close()

	file=open('config_file_gotcloud_template','a')
	file.write('#BATCH_TYPE = mosix\n')
	file.write('#BATCH_OPTS = -j21,22,23,24,30\n')
	file.write('REF_DIR = /data/local/ref/gotcloud.ref\n')
	file.write('AS = NCBI37\n')
	file.write('REF = /net/1000g/mktrost/seqshop/gotcloud/gotcloud.ref/human.g1k.v37.fa\n')
	file.write('DBSNP_VCF = $(REF_DIR)/dbsnp_135.b37.vcf.gz\n')
	file.write('HM3_VCF = $(REF_DIR)/hapmap_3.3.b37.sites.vcf.gz\n')
	file.write('MAP_TYPE = BWA_MEM\n')
	file.write('BAMUTIL_THINNING = --phoneHomeThinning 0\n')
	file.close()

try:
	opts, args=getopt.getopt(sys.argv[1:],'c:t:h',['conf','template','help'])
except getopt.GetoptError:
	print('This pipeline is to automatically run the pipeline for the ATACseq process.') 
	print('See https://github.com/YuhuaZhang1995/ATACseq_pipeline for more help.')
	print('python ATACseq_pipeline.py --help (for help)')
	print('python ATACseq_pipeline.py --template (for template for config file to run the whole pipeline)')
	print('python ATACseq_pipeline.py -c <config_file> (to run the pipeline)')
	sys.exit(2)

for opt, arg in opts:
	if opt in ('-h','--help'):
		print('This pipeline is to automatically run the pipeline for the ATACseq process.') 
		print('============================')
		print('Usage: (to run the pipeline)')
		print('============================')
		print('python ATACseq_pipeline.py -c <config_file>')
		print('To run the whole pipeline, at least prepare the config files, prepare at least:') 
		print('1. config file for the pipeline;')
		print('2. config file for gotCloud;')
		print('3. bedprofile for bedprofilecounts')
		print('Use \"python ATACseq_pipeline.py --template\" to generate templates for required files for the whole process')
		print('==============================')
		print('Usage: (to run separate steps)')
		print('==============================')
		print('python ATACseq_pipeline.py -c <config_file>')
		print('Specify the parameters in the config files, \"##\" for steps, \"--\" for input parameters. For each step user need to the corresponding parameter as shon below. Write each parameter in a new line.')
		print('##adapter_trimming')
		print('--job_AT <#of_process_for_adaptor_trimming>')
		print('--seq_data <pathway_to_seq_data>')
		print('--core_info_file <core_info_file>')
		print('--out_sampleinfo <output_sample_info>')
		print('##alignment')
		print('--job_align <#_of_process_for_alignment>')
		print('--trimmed_file <pathway_to_adapter_trimmed_file>')
		print('--out_bam <output_directory_to_store_generated_bam_file>')
		print('--conf <conf/conf_file>')
		print('--batch <batch>')
		print('--run <run>')
		print('##filtering')
		print('--job_filter <#_of_process_for_filtering>')
		print('--in_bam <metagot_file>')
		print('--out_proc_bam <Output_dir_for_processed_bam>')
		print('--intermediate_file Yes')
		print('--batch <batch>')
		print('--run <run>')
		print('##plotting')
		print('--in_bam <metagot_file>')
		print('--in_bam_filtered <metagot_filtered_file>')
		print('--out_plot <pathway_to_store_the_output_plot.pdf>')
		print('--batch <batch>')
		print('--run <run>')
		print('##getting_coverage')
		print('--dir_filtered_bam <pathway_to_filtered_bam_files>')
		print('--out_coverage <directory_to_store_genome_coverage>')
		print('--batch <batch>')
		print('--run <run>')
		print('##making_bed')
		print('--dir_filtered_bam <pathway_to_filtered_bam_files>')
		print('--dir_bam <pathway_to_>')
		print('--out_bed <output_directory_to_store_bed_file>')
		print('--batch <batch>')
		print('--run <run>')
		print('##signal_to_noise')
		print('--bedprofile <bedprofile>')
		print('--dir_bed <directory_to_bed_files>')
		print('--out_clus <output_dir_to_store_clus_file>')
		print('--sample_info <sample_info_file>')
		print('--out_s2n <out_dir_to_store_s2n_files>')
		print('--batch <batch>')
		print('--run <run>')
		print('##summary')
		print('--batch <batch>')
		print('--run <run>')
		print('--dir_sampleinfo <directory_to_sampleinfo>')
		print('--dir_proc_bam <directory_to_processed_bams>')
		print('--dir_converge <directory_to_coverage_files>')
		print('--dir_s2n <dir_to_signal_to_noise>')
		print('--out_summary <directory_to_store_the_summary>')
		print('\n')
		print('See https://github.com/CutaneousBioinf/ATAC-seq for more help.')
		sys.exit(2)
	elif opt in ('-t','--template'):
		generate_template()
		print('Four templates: config_file, specific_output, bedprofile, config file for gotcloud are generated, pls revise the parameter correspondingly')
		print('config_file will be used as the input parameters for the pipeline.')
		print('specific_output will be used if you\' like to set the specific output for each intermediate file (not necessary for the pipeline to work)')
		print('bedprofile will be used in the signal to noise analysis, required by bedprofilecounts')
		print('gotcloud config file will be used by gotcloud in alignment, you will need to specify the REF genome etc,.')
		sys.exit(2)
	elif opt in ('-c','--conf'):
		conf_file=arg
	else:
		print('This pipeline is to automatically run the pipeline for the ATACseq process.') 
		print('See https://github.com/YuhuaZhang1995/ATACseq_pipeline for more help.')
		print('python ATACseq_pipeline.py --help (for help)')
		print('python ATACseq_pipeline.py --template (for template for config file for the whole process)')
		print('python ATACseq_pipeline.py -c <config_file> (to run the pipeline)')
		sys.exit(2)


config_file=open(conf_file)
commands={}
func_call=[]
for ele in config_file:
	if (re.search('--',ele)):
		ele=re.sub('--','',ele)
		tmp=ele.rstrip().split(' ')
		commands[tmp[0]]=tmp[1:]
	elif (re.search('##',ele)):
		ele=re.sub('##','',ele)
		tmp=ele.rstrip().split(' ')
		func_call.append(tmp[0])
config_file.close() 

#first deal with the case where all the function will be called
#if all(x in func_call for x in ['alignment','adapter_trimming','filtering','plotting']):
if len(func_call)==0:
	spec_out={}
	if 'specific_output' in commands:
		tmp_file=open(commands['specific_output'][0])
		for ele in tmp_file:
			if (re.search('--',ele)):
				ele=re.sub('--','',ele)
				tmp=ele.rstrip().split(' ')
				spec_out[tmp[0]]=tmp[1:]
		tmp_file.close()

	entire_output=os.path.dirname(os.path.realpath(__file__))
	if 'entire_output' in commands:
		entire_output=commands['entire_output'][0]
		if not os.path.exists(entire_output):
			os.makedirs(entire_output)

	#trim_adapter
	if all(x in commands for x in ['core_info_file','seq_data','batch','run','conf','bedprofile']):
		core_info_file=commands['core_info_file'][0]
		seq_data=commands['seq_data'][0]
		job_AT=9
		if 'job_AT' in commands:
			job_AT=int(commands['job_AT'][0])
		out_sampleinfo=entire_output+'/data'
		trim_read=False
		if 'trim_read' in commands:
			if commands['trim_read'][0] in ['Yes']:
				trim_read=True
		reads=38
		if 'reads' in commands:
			reads=int(commands['reads'][0])
		direct='l'
		if 'direct' in commands:
			direct=commands['direct'][0]
		if not os.path.exists(out_sampleinfo):
			os.makedirs(out_sampleinfo)
		command=[seq_data,core_info_file,job_AT,out_sampleinfo,trim_read,reads,direct]
		if ('specific_output' in commands) & ('out_sampleinfo' in spec_out):
			out_sampleinfo=spec_out['out_sampleinfo'][0]
			if not os.path.exists(out_sampleinfo):
				os.makedirs(out_sampleinfo)
			command=[seq_data,core_info_file,job_AT,out_sampleinfo,trim_read,reads,direct]
		adapter_trimming_update.concat_func(command)
		#print(command)

		#alignment through gotcloud
		trimmed_file=seq_data
		job_align=5
		if 'job_align' in commands:
			job_align=int(commands['job_align'][0])
		conf=commands['conf'][0]
		batch=commands['batch'][0]
		run=commands['run'][0]
		out_conf=conf.rsplit('/',1)[0]
		if 'specific_output' in commands:
			if 'out_bam' in spec_out:
				out_bam=spec_out['out_bam'][0]
				if not os.path.exists(out_bam):
					os.makedirs(out_bam)
				if not os.path.exists(out_bam+'/Batch'+batch+'_Run'+run):
					os.system('mkdir '+out_bam+'/Batch'+batch+'_Run'+run)
				out_bam=out_bam+'/Batch'+batch+'_Run'+run
		else:
			if not os.path.exists(entire_output+'/BAM'):
				os.makedirs(entire_output+'/BAM')
			if not os.path.exists(entire_output+'/BAM/Batch'+batch+'_Run'+run):
				os.makedirs(entire_output+'/BAM/Batch'+batch+'_Run'+run)
			out_bam=entire_output+'/BAM/Batch'+batch+'_Run'+run

		command=[out_bam,trimmed_file,out_conf,job_align,batch,run,conf]
		alignment_gotcloud.concat_func(command)
		#print(command)

		#filter the reads
		in_bam=out_conf+'/metagotCloudbamfiles_Batch'+batch+'_Run'+run
		job_filter=9
		if 'job_filter' in commands:
			job_filter=int(commands['job_filter'][0])

		if ('specific_output' in commands) and ('out_proc_bam' in spec_out):
			out_proc_bam=spec_out['out_proc_bam'][0]
			if not os.path.exists(out_proc_bam):
				os.makedirs(out_proc_bam)
			if not os.path.exists(out_proc_bam+'/Batch'+batch+'_Run'+run):
				os.system('mkdir '+out_proc_bam+'/Batch'+batch+'_Run'+run)
			out_proc_bam=out_proc_bam+'/Batch'+batch+'_Run'+run
		else:
			if not os.path.exists(entire_output+'/BAMprocessed'):
				os.makedirs(entire_output+'/BAMprocessed')
			if not os.path.exists(entire_output+'/BAMprocessed/Batch'+batch+'_Run'+run):
				os.makedirs(entire_output+'/BAMprocessed/Batch'+batch+'_Run'+run)
			out_proc_bam=entire_output+'/BAMprocessed/Batch'+batch+'_Run'+run

		intermediate_file=False
		if 'intermediate_file' in commands:
			if commands['intermediate_file'][0] in ['Yes']:
				intermediate_file=True
		command=[in_bam,out_proc_bam,job_filter,intermediate_file,batch,run]
		filtering.concat_func(command)
		#print(command)

		#plotting
		in_bam_filtered=out_proc_bam+'/metagotCloudbamfiles_filtered_Batch'+batch+'_Run'+run
		if ('specific_output' in commands) and ('out_plot' in spec_out):
			out_plot=spec_out['out_plot'][0]
			if not os.path.exists(out_plot):
				os.makedirs(out_plot)
			if not os.path.exists(out_plot+'/Batch'+batch+'_Run'+run):
				os.makedirs(out_plot+'/Batch'+batch+'_Run'+run)
			out_plot=out_plot+'/Batch'+batch+'_Run'+run
		else:
			if not os.path.exists(entire_output+'/QCs'):
				os.makedirs(entire_output+'/QCs')
			if not os.path.exists(entire_output+'/QCs/InsertSize'):
				os.makedirs(entire_output+'/QCs/InsertSize')
			if not os.path.exists(entire_output+'/QCs/InsertSize/Batch'+batch+'_Run'+run):
				os.makedirs(entire_output+'/QCs/InsertSize/Batch'+batch+'_Run'+run)
			out_plot=entire_output+'/QCs/InsertSize/Batch'+batch+'_Run'+run

		cur_ptw=os.path.dirname(os.path.realpath(__file__))
		os.system('Rscript '+cur_ptw+'/insertsizehist.R '+in_bam+' '+out_plot+'/metagotCloudbamfiles_Batch'+batch+'_Run'+run+'.pdf')
		os.system('Rscript '+cur_ptw+'/insertsizehist.R '+in_bam_filtered+' '+out_plot+'/metagotCloudbamfiles_filtered_Batch'+batch+'_Run'+run+'.pdf')

		#get genome_coverage
		if ('specific_output' in commands) and ('out_coverage' in spec_out):
			out_coverage=spec_out['out_coverage'][0]
			if not os.path.exists(out_coverage):
				os.makedirs(out_coverage)
		else:
			if not os.path.exists(entire_output+'/QCs'):
				os.makedirs(entire_output+'/QCs')
			if not os.path.exists(entire_output+'/QCs/Breadth'):
				os.makedirs(entire_output+'/QCs/Breadth')
			out_coverage=entire_output+'/QCs/Breadth'
		dir_filtered_bam=out_proc_bam
		command=[dir_filtered_bam,out_coverage,batch,run]
		genome_coverage.concat_func(command)
		#print(command)

		#generate bed files
		if ('specific_output' in commands) and ('out_bed' in spec_out):
			out_bed=spec_out['out_bed'][0]
			if not os.path.exists(out_bed):
				os.makedirs(out_bed)
			out_bed=out_bed+'/Batch'+batch+'_Run'+run
			if not os.path.exists(out_bed):
				os.makedirs(out_bed)
		else:
			if not os.path.exists(entire_output+'/BED'):
				os.makedirs(entire_output+'/BED')
			if not os.path.exists(entire_output+'/BED/Batch'+batch+'_Run'+run):
				os.makedirs(entire_output+'/BED/Batch'+batch+'_Run'+run)
			out_bed=entire_output+'/BED/Batch'+batch+'_Run'+run
		pathway1=out_bam
		pathway2=out_proc_bam
		command=[pathway1,pathway2,out_bed]
		make_bed.concat_func(command)
		#print(command)

		#signal to noise
		if ('specific_output' in commands) and ('out_clus' in spec_out):
			out_clus=spec_out['out_clus'][0]
			if not os.path.exists(out_clus):
				os.makedirs(out_clus)
			if not os.path.exists(out_clus+'/Batch'+batch+'_Run'+run):
				os.makedirs(out_clus+'/Batch'+batch+'_Run'+run)
		else:
			if not os.path.exists(entire_output+'/QCs'):
				os.makedirs(entire_output+'/QCs')
			if not os.path.exists(entire_output+'/QCs/SignaltoNoise'):
				os.makedirs(entire_output+'/QCs/SignaltoNoise')
			if not os.path.exists(entire_output+'/QCs/SignaltoNoise/TSS'):
				os.makedirs(entire_output+'/QCs/SignaltoNoise/TSS')
			if not os.path.exists(entire_output+'/QCs/SignaltoNoise/TSS/Batch'+batch+'_Run'+run):
				os.makedirs(entire_output+'/QCs/SignaltoNoise/TSS/Batch'+batch+'_Run'+run)
			out_clus=entire_output+'/QCs/SignaltoNoise/TSS/Batch'+batch+'_Run'+run
		bedprofile=commands['bedprofile'][0]
		sample_info=out_sampleinfo+'/SampleInfo_Batch'+batch+'_Run'+run

		if ('specific_output' in commands) and ('out_s2n' in spec_out):
			out_s2n=spec_out['out_s2n'][0]
			if not os.path.exists(out_s2n):
				os.makedirs(out_s2n)
		else:
			if not os.path.exists(entire_output+'/QCs'):
				os.makedirs(entire_output+'/QCs')
			if not os.path.exists(entire_output+'/QCs/SignaltoNoise'):
				os.makedirs(entire_output+'/QCs/SignaltoNoise')
			if not os.path.exists(entire_output+'/QCs/SignaltoNoise/analysis'):
				os.makedirs(entire_output+'/QCs/SignaltoNoise/analysis')
			out_s2n=entire_output+'/QCs/SignaltoNoise/analysis'
			dir_bed=out_bed
		command=[bedprofile,dir_bed,out_clus,sample_info,out_s2n,batch,run]
		signal_to_noise.concat_func(command)
		#print(command)

		#get the summary
		if ('specific_output' in commands) and ('out_summary' in spec_out):
			out_summary=spec_out['out_summary'][0]
			if not os.path.exists(out_summary):
				os.makedirs(out_s2n)
		else:
			if not os.path.exists(entire_output+'/QCs'):
				os.makedirs(entire_output+'/QCs')
			if not os.path.exists(entire_output+'/QCs/Summary'):
				os.makedirs(entire_output+'/QCs/Summary')
			out_summary=entire_output+'/QCs/Summary'
		command=[out_sampleinfo,out_proc_bam,out_coverage,out_s2n,batch,run,out_summary]
		summary.concat_func(command)
		#print(command)

	else:
		print('Please specify at least the --core_info_file --seq_data --batch --run --conf')

elif ('adapter_trimming' in func_call):
	if all(x in commands for x in ['core_info_file','seq_data']):
		core_info_file=commands['core_info_file'][0]
		seq_data=commands['seq_data'][0]
		job_AT=5
		if 'job_AT' in commands:
			job_AT=int(commands['job_AT'][0])
		out_sampleinfo=os.path.dirname(os.path.realpath(__file__))
		trim_read=False
		if 'trim_read' in commands:
			if commands['trim_read'][0] in ['Yes']:
				trim_read=True
		reads=38
		if 'reads' in commands:
			reads=int(commands['reads'][0])
		direct='l'
		if 'direct' in commands:
			direct=commands['direct'][0]
		if 'out_sampleinfo' in commands:
			out_sampleinfo=commands['out_sampleinfo'][0]
		command=[seq_data,core_info_file,job,out_sampleinfo,trim_read,reads,direct]
		adapter_trimming_update.concat_func(command)
	else:
		print('Please specify at least the --core_info_file --seq_data')

elif ('alignment' in func_call):
	if all(x in commands for x in ['trimmed_file','out_bam','batch','run','conf']):
		trimmed_file=commands['trimmed_file'][0]
		job_align=3
		if 'job_align' in commands:
			job_align=int(commands['job_align'][0])
		batch=commands['batch'][0]
		run=commands['run'][0]
		out_bam=commands['out_bam'][0]
		if not os.path.exists(out_bam):
			os.makedirs(out_bam)
		if not os.path.exists(out_bam+'/Batch'+batch+'_Run'+run):
			os.system('mkdir '+out_bam+'/Batch'+batch+'_Run'+run)
		out_bam=out_bam+'/Batch'+batch+'_Run'+run
		conf=commands['conf'][0]
		out_conf=conf.rsplit('/',1)[0]
		command=[out_bam,trimmed_file,out_conf,job_align,batch,run,conf]
		alignment_gotcloud.concat_func(command)
	else:
		print('Please at least specify the pathway to trimmed_file(--trimmed_file),')
		print('pathway to store config files(--out_conf), pathway to store output bam files(--out_bam)')
		print('batch number(--batch) and run number(--run)')

elif ('filtering' in func_call):
	if all(x in commands for x in ['out_proc_bam','in_bam','batch','run']):
		in_bam=commands['in_bam'][0]
		job_filter=5
		if 'job_filter' in commands:
			job_filter=int(commands['job_filter'][0])
		batch=commands['batch'][0]
		run=commands['run'][0]
		out_proc_bam=commands['out_proc_bam'][0]
		if not os.path.exists(out_proc_bam):
			os.makedirs(out_proc_bam)
		if not os.path.exists(out_proc_bam+'/Batch'+batch+'_Run'+run):
			os.system('mkdir '+out_proc_bam+'/Batch'+batch+'_Run'+run)
		out_proc_bam=out_proc_bam+'/Batch'+batch+'_Run'+run
		intermediate_file=False
		if 'intermediate_file' in commands:
			if commands['intermediate_file'][0] in ['Yes','yes','Y','y']:
				intermediate_file=True
		command=[in_bam,out_proc_bam,job_filter,intermediate_file,batch,run]
		filtering.concat_func(command)
	else:
		print('Please at least specify the --in_bam, --out_proc_bam, --batch, --run')

elif ('plotting' in func_call):
	if all(x in commands for x in ['in_bam','in_bam_filtered','out_plot','batch','run']):
		in_bam=commands['in_bam'][0]
		in_bam_filtered=commands['in_bam_filtered'][0]
		out_plot=commands['out_plot'][0]
		if not os.path.exists(out_plot):
			os.makedirs(out_plot)
		batch=commands['batch'][0]
		run=commands['run'][0]
		os.system('Rscript insertsizehist.R '+in_bam+' '+out_plot+'/metagotCloudbamfiles_Batch'+batch+'_Run'+run+'.pdf')
		os.system('Rscript insertsizehist.R '+in_bam_filtered+' '+out_plot+'/metagotCloudbamfiles_filtered_Batch'+batch+'_Run'+run+'.pdf')

elif ('getting_coverage' in func_call):
	if all(x in commands for x in ['dir_filtered_bam','out_coverage','batch','run']):
		dir_filtered_bam=commands['dir_filtered_bam'][0]
		out_coverage=commands['out_coverage'][0]
		batch=commands['batch'][0]
		run=commands['run'][0]
		out_coverage=out_coverage+'/Batch'+batch+'_Run'+run
		if not os.path.exists(out_coverage):
			os.makedirs(out_coverage)
		command=[dir_filtered_bam,out_coverage,batch,run]
		genome_coverage.concat_func(command)

elif ('making_bed' in func_call):
	if all(x in commands for x in ['dir_filtered_bam','dir_bam','out_bed','batch','run']):
		pathway2=commands['dir_filtered_bam'][0]
		pathway1=commands['dir_bam'][0]
		out_bed=commands['out_bed'][0]
		batch=commands['batch'][0]
		run=commands['run'][0]
		if not os.path.exists(out_bed):
			os.makedirs(out_bed)
		out_bed=out_bed+'/Batch'+batch+'_Run'+run
		if not os.path.exists(out_bed):
			os.makedirs(out_bed)
		command=[pathway1,pathway2,out_bed]
		make_bed.concat_func(command)

elif ('signal_to_noise' in func_call):
	if all(x in commands for x in ['bedprofile','dir_bed','out_clus','sample_info','out_s2n','batch','run']):
		bedprofile=commands['bedprofile'][0]
		dir_bed=commands['dir_bed'][0]
		batch=commands['batch'][0]
		run=commands['run'][0]
		out_clus=commands['out_clus'][0]
		if not os.path.exists(out_clus):
			os.makedirs(out_clus)
		out_clus=out_clus+'/Batch'+batch+'_Run'+run
		if not os.path.exists(out_clus):
			os.makedirs(out_clus)
		sample_info=commands['sample_info'][0]
		out_s2n=commands['out_s2n'][0]
		if not os.path.exists(out_s2n):
			os.makedirs(out_s2n)
		'''out_s2n=out_s2n+'/Batch'+batch+'_Run'+run
		if not os.path.exists(out_s2n):
			os.makedirs(out_s2n)'''
		command=[bedprofile,dir_bed,out_clus,sample_info,out_s2n]
		signal_to_noise.concat_func(command)

elif ('summary' in func_call):
	if all(x in commands for x in ['dir_sampleinfo','dir_proc_bam','dir_coverage','dir_s2n','batch','run','out_summary']):
		dir_sampleinfo=commands['dir_sampleinfo'][0]
		dir_proc_bam=commands['dir_proc_bam'][0]
		dir_coverage=commands['dir_coverage'][0]
		dir_s2n=commands['dir_s2n'][0]
		batch=commands['batch'][0]
		run=commands['run'][0]
		out_summary=commands['out_summary'][0]
		if not os.path.exists(out_summary):
			os.makedirs(out_summary)
		command=[dir_sampleinfo,dir_proc_bam,dir_coverage,dir_s2n,batch,run,out_summary]
		summary.concat_func(command)

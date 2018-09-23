# ATACseq_pipeline

This is an automatic pipeline of processing ATAC-seq data. The pipeline will need the user to provide the raw ATACseq sequence data and will do the following steps including adapter-trimming, alignment, filtering, plotting chromosome pattern, getting the coverage depth, generating .bed filea, comparing signal to noise and summarizing all the previous steps.

The pipeline was contributed by Yuhua Zhang and Alex Tsoi. Also thanks to Matthew Patrick for sharing his suggestions.

## Usage of the ATACseq_pipeline: 

**_Running the pipeline_**

python ATACseq_pipeline.py -c <config_file>

**_Generating template files to modify_**

python ATACseq_pipeline.py --template

**_For help_**

python ATACseq_pipeline.py --help

## Component of config_file:

1. Procedures to be involved, start with **##**, e.g. **_##alignment_**
2. Parameters for the corresponding procedure, start with **--**,e.g. **_--batch 14_**

## A quick start
In this case, the program will go through the whole procedure, and generate the intermediate files at each step and final results. The user need to specify **at least** the required parameters.

The user will need to provide several files include 

1. config files for the whole pipelines;
2. bedprofile for bedprofilecounts;
3. config files for gtcloud alignment;

\*4. specific output for intermediate files if required;

Note that the template for these files can be generated by **_python ATACseq_pipeline.py --template_**. The user can modify the parameters as they like.

- **Required parameter**

`--seq_data`: Pathway to the sequence data, end with Run_XXX, no '/' included at the end; the program will search for '/elder' folder and find the sequence data, **_e.g. --seq_data /pathway/to/data/Run_XXX_**

`--core_info_file`: The core information files, **_e.g. /pathway/to/file/Batch14_Run_1789_elder.txt_**;

`--batch`: The batch number, e.g. **_--batch 14_**;

`--run`: The run number, e.g. **_--run 1789_**;

`--conf`: The config file used by gotcloud. Please don't specify OUT_DIR and FASTQ_LIST in the config file. **_e.g. --conf /patheway/to/config_file_**

`--bedprofile`: The paramters specified for the bedprofilecounts. **_e.g. --bedprofile /pathway/../bedprofile_**

- **Optional parameter**

`--specific_output`: The output pathway file that contains the specific output for each step; **_e.g. --specific_output /pathway/to/file_**. The output pathway file should contain the specific pathway to store the intermediate files generated by the pipeline, specified by **_--out_sampleinfo --out_bam --out_proc_bam --out_plot --out_coverage --out_bed --out_clus --out_s2n --out_summary_**. If not specified every output pathway, the default will be the current pathway. The user can use **_python ATACseq_pipeline.py --template_** to generate a template for this file.

`--entire_output`: Where you would like to store the all the output(intermediate files and final results). The default is the current direrctory. Several directories will be created under the entire_output directoies, including /BAM, /BAMprocessed, /BED, /QCs.

`--job_AT`: Number of process for adapter trimming step, the default is 5;

`--job_align`: Number of process for alignment, the default is 3. Please pay attention to the **memory**, as the demand for the memory is considerably high;

`--job_filter`: Number of process for filtering, the default is 5;

`--intermediate_file`: Whether to keep the intermediate .bam files generated in the filtering step, the default is no. e.g. **_--intermediate_file Yes_**


- **A brief example of config_file**

--entire_output .

--specific_output ./output_file

--batch 14

--run 1789

--job_AT 5

--seq_data /pathway/to/seq/data

--core_info_file /core/info/file.txt

--job_align 5

--conf /pathway/to/config_file

--job_filter 5

--intermediate_file Yes

--bedprofile /pathway/../bedprofile

## Separate each step
If the user prefer to run each procedure step by step and specify the output pathway for the intermediate files, just make the config_file contains only one operation each time.

- **Adapter_trimming**

This step will generate the trimmed files under the same directory as the input sequence data, as well as the extracted information from the input core info file(sampleID, cell line description and patient ID). 

`--core_info_file`: The core information files, **_e.g. /pathway/to/file/Batch14_Run_1789_elder.txt_**;

`--seq_data`: Pathway to the sequence data, end with Run_XXX, no '/' included at the end; the program will search for '/elder' folder and find the sequence data, **_e.g. --seq_data /pathway/to/data/Run_XXX_**

`--job_AT`: Optional parameter. Number of process for adapter trimming step, the default is 5;

`--out_sampleinfo`: Optional parameter. Directory to store the generated sample info file, **_e.g. out_sampleinfo /pathway/.._**


**_Example_**

##adapter_trimming

--core_info_file /core/info.txt

--job_AT 5

--seq_data /pathway/to/seq/data

--output_AT_core_info /pathway/to/store


- **Alignment**

This step will generate the aligned .bam files as well as a metaCloudbamfiles which will be used for the filtering. The .bam files will be stored in the pathway either specified by the user or by default the current disrctory. The metaCouldbamfiles will be stored in the same pathway as the config file for gotcloud.

`--trimmed_file`: Pathway to the trimmed files, end with Run_XXX, no '/' included at the end; the program will search for '/elder' folder and find the sequence data;

`--out_conf`: Directory to store config and index files that would be used by gotcloud;

`--out_bam`: Directory to store the output bam files. A new directory named as _BatchXX_RunXXXX_ will be created to store the generated .bam files;

`--batch`: The batch number, e.g. **_--batch 14_**;

`--run`: The run number, e.g. **_--run 1789_**;

`--conf`: The config_file used by gotcloud, **_e.g. --conf /pathway/to/config_file_**;

`--job_align`: Optional parameter. Number of process for alignment, the default is 3. Please pay attention to the **memory**, as the demand for the memory is considerably high;


**_Example_**

##alignment

--job_align 5

--trimmed_file /pathway/to/adapter/trimmed/file

--out_conf /directory/to/store/config/index/files/for/gotcloud

--out_bam /output/directory/to/store/generated/bam/file

--batch 14

--run 1789

--conf /pathway/to/config_file


- **Filtering**

This step will generate the filtered .bam file as well as the metagotCloudbamfiles which will be used in the plotting step. Whether to keep the intermediate files is up to the users.

`--in_bam`: MetagotCloud files. e.g. **_--in_bam /pathway/metaXXXX_**;

`--out_proc_bam`: Output dir for processed .bam files; A new directory named as _BatchXX_RunXXXX_ will be created to store the processed .bam files

`--batch`: The batch number, e.g. **_--batch 14_**;

`--run`: The run number, e.g. **_--run 1789_**;

`--job_filter`: Optional parameter. Number of process for filtering, the default is 5;

`--intermediate_file`: Optional parameter. Whether to keep the intermediate .bam files generated in the filtering step, the default is no. e.g. **_--intermediate_file Yes_**



**_Example_**

##filtering

--job_filter 5

--out_proc_bam /pathway

--in_bam /pathway/metagotCloudbamfiles_Batch14_Run1789

--intermediate_file Yes

--batch 14

--run 1789

- **Plotting**

This step will generate the plots of the frequency of the insert sizes, compared the filtered reads to the unfiltered reads.

`--in_bam`: Directory to MetagotCloud files for unfiltered reads, the pipeline will search for the metagotcloud files based on the run and batch number provided by user. e.g. **_--in_bam /pathway/..**;

`--in_bam_filtered`:Directory to MetagotCloud files for filtered reads, the pipeline will search for the metagotcloud files based on the run and batch number provided by user. e.g. **_--in_bam /pathway/metaXXXX_**;

`--out_plot`: The directory to store the generated plots, A new directory named as _BatchXX_RunXXXX_ will be created to store the generated plots;

`--batch`: The batch number, e.g. **_--batch 14_**;

`--run`: The run number, e.g. **_--run 1789_**;

**_Example_**

##plotting

--in_bam /pathway/meta_file

--in_bam_filtered /pathway/meta_filtered_file

--out_plot /pathway/plot

--batch 14

--run 1789

- **getting_coverage**

This step will generate the summary table of the coverage depth in the whole genome. The coverage depth will check 1X, 2X and 5X correspondingly.

`--dir_filtered_bam`: The directory to the filtered .bam files (Please include BatchXX_RunXX) e.g. **_--dir_filtered_bam /pathway/../BatchXX_RunXXXX_**;

`--out_coverage`: The output directory where you would like yo store the output, e.g. **_--out_coverage /pathway/to/store/output_**;

`--batch`: The batch number, e.g. **_--batch 14_**;

`--run`: The run number, e.g. **_--run 1789_**;

**_Example_**

##getting_coverage

--dir_filtered_bam /path/way

--out_coverage /pathway/for/output

--batch 14

--run 1789

- **making_bed**

This step will generate the bed file for both filtered and unfiltered .bam files.

`--dir_filtered_bam`: The directory to the filtered .bam files, e.g. **_--dir_filtered_bam /path/way_**;

`--dir_bam`: The directory to the unfiltered .bam files, e.g. **_--dir_bam /path/way_**;

`--out_bed`: The output directory where you would like yo store the output, the pipeline will generate a directory named BatchXX_RunXXXX under the directory provided by the user. e.g. **_--out_bed /pathway/to/store/output_**;

`--batch`: The batch number, e.g. **_--batch 14_**;

`--run`: The run number, e.g. **_--run 1789_**;

**_Example_**

##making_bed

--dir_filtered_bam /path/way

--dir_bam /path/way

--out_bed /pathway/to/store/output

--batch 14

--run 1789


- **signal_to_noise**

This step will generate the plots of signal to noise;

`--bedprofile`: The bedprofile that will be used by bedprofilecount, the template can be generated by the pipeline. e.g. **_--bedprofile /pathway/../bedprofile_**;

`--dir_bed`: The directory to the .bed files, e.g. **_--dir_bed /path/way_**;

`--out_clus`: The output directory where you would like yo store the .clus file, e.g. **_--out_bed /pathway/to/store/output_**;

`--sample_info`: The previous generated sampleInfo file, e.g. **_--sample_info /pathway/../SampleInfo_BatchXX_RunXXXX_**;

`--out_s2n`: The output directory to store the plot of signal to noise, e.g. **_--out_s2n /pathway/to/store/plot_**

`--batch`: The batch number, e.g. **_--batch 14_**;

`--run`: The run number, e.g. **_--run 1789_**;

**_Example_**

##signal_to_noise

--bedprofile /pathway/../bedprofile

--dir_bed /path/way

--out_clus /pathway/to/store/output

--sample_info /pathway/../sample_info

--out_s2n /pathway/to/store/plot

--batch 14

--run 1789

- **summary**

This step will generate the summary of the quality control steps;

`--dir_sampleinfo`: The directory to the sampleinfo file generated before, the pipeline will search for _SampleInfo_BatchXX_RunXXXX_ under the directory provided. e.g. **_--dir_sampleinfo /pathway/..**;

`--dir_proc_bam`: The directory to the processed .bam files, the pipeline will search for _readcount_BatchXX_RunXXXX_ under the directory provided. e.g. **_--dir_proc_bam /pathway/..**;

`--dir_coverage`: The directory to the previous summary table of coverage depth, the pipeline will search for _*BatchXX_RunXXXX.table_ under the directory provided. e.g. **_--dir_coverage /pathway/..**;

`--dir_s2n`: The directory to the previous summary table of signal to noise, the pipeline will search for _*BatchXX_RunXXXX.signal2noise_ under the directory provided by user. e.g. **_--dir_s2n /pathway/..**;

`--out_summary`: The output directory to store the summary table, e.g. **_--out_summary /pathway/..**

`--batch`: The batch number, e.g. **_--batch 14_**;

`--run`: The run number, e.g. **_--run 1789_**;

**_Example_**

##summary

--dir_sampleinfo /pathway/..

--dir_proc_bam /pathway/..

--dir_coverage /pathway/..

--dir_s2n /pathway/..

--out_summary /pathway/..

--batch 14

--run 1789



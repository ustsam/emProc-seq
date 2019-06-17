# emProc-seq
Python Script for treating emProc-seq data

## Ownership
[Xuhui HUANG Lab at HKUST]  
[Jihuang WANG Lab at HKUST]

## Status
Active Development

## Introduction
Circular sequencing was initially developed by [Andino Group at UCSF](https://andino.ucsf.edu/) in 2014 for detecting low-frequency variants in the study of RNA virus by reducing sequencing errors [1]. The authors also released their data analysis software, namely CirSeq, which is open-access at https://andino.ucsf.edu/CirSeq.

However, the original was written in hard-code form for repeat and read length. The code is rewritten to perform >=2 repeat identification and variable readlength. The relocalization scripts are rewritten to fit the use of bwa-mem as mapper. Selection of the relocalized read is also changed to be more stringent parameter.

## Simulation data

Script for generating simulation data is located at directory named as simulation_script/. It is a matlab script created by Dr. Biaobin Jiang from Prof. Jiguang Wang's group at HKUST. AY1273.fastq.gz are provided as one set of demo data, typical run time would be around 20-30 minutes.

Demo data can be run as following:

1. Download the zipped script file.
2. Unzip the file.
3. Enter script directory
3. Compile the codes by typing `python setup_newreloc.py build_ext --inplace`.
4. `mkdir demo_data_real_time`
5. Call the function using `./run_noQsfilter_bwa.sh ./demo_data_real_time ./reference/rDNA1.fa ./ DUMMY 2 600 ./demo_data/AY1273.fastq.gz `


## Usage
1. Download the zipped script file.
2. Unzip the file.
3. Enter script directory
3. Compile the codes by typing `python setup_newreloc.py build_ext --inplace`.
4. Call the function using `./run_noQsfilter_bwa.sh {PATH of the output directory} {PATH of the reference file} {PATH of the script directory} DUMMY 2 ${twice of the max readlength} ${PATH of the data file in gzipped form}`. The data file are suggested to process with triming software for adaptor and low qulaity base before running this script.

## Figure plotting

To obtain figures in the manuscript, please run the following three main functions:

calling_transcription_errors.m  
plotFigure1.m  
plotFigure2.m  

It is a matlab script created by Dr. Biaobin Jiang(https://github.com/bbjiang) from Prof. Jiguang Wang's group at HKUST.

## Methodology:

Please check methodology.pdf for detail of the script and the meaning of all generated files. Paper is not yet publish, link will be updated as soon as possible.

## System requirements

The following OS had been used for running the script.

Ubuntu 16.04 LTS (Xenial Xerus)

CentOS release 6.6 (Final)


The following packages are prerequisites for using emProc-seq

1. Python (version 2.7.12)    
2. Cython (version 0.23.4)   
3. NumPy (version 1.11.0)     
4. SciPy (version 0.17.0)    
5. bwa (version 0.7.17-r1188)   
6. samtools (version 1.4.1)

NOTE 1: Cython requires a compiler. For OSX this may require installation of Xcode.

NOTE 2: bwa and samtools binaries must be in the PATH.

NOTE 3: Install time are estimated to be no more than 1.5 hours.


## Reference
[1] Acevedo, A., Brodsky, L., & Andino, R. (2014). Mutational and fitness landscapes of an RNA virus revealed through population sequencing. Nature, 505, pp.686-690.

## Contact
For technical questions, please contact TinHang via email: thchong@connect.ust.hk



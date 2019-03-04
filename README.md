# emProc-seq
Python Script for treating emProc-seq data

## Ownership
[Huang Lab at HKUST]

## Status
Active Development

## Introduction
Circular sequencing was initially developed by [Andino Group at UCSF](https://andino.ucsf.edu/) in 2014 for detecting low-frequency variants in the study of RNA virus by reducing sequencing errors [1]. The authors also released their data analysis software, namely CirSeq, which is open-access at https://andino.ucsf.edu/CirSeq.

However, the original was written in hard-code form for repeat and read length. The code is rewritten to perform >=2 repeat identification and variable readlength. The relocalization scripts are rewritten to fit the use of bwa-mem as mapper. Selection of the relocalized read is also changed to be more stringent parameter.

## Usage
1. Download the zipped script file.
2. Unzip the file.
3. Enter script directory
3. Compile the codes by typing `python setup_newreloc.py build_ext --inplace`.
4. Call the function using `./run_noQsfilter_bwa.sh {PATH of the output directory} {PATH of the reference file} {PATH of the script directory} DUMMY 2 ${twice of the max readlength} ${PATH of the data file in gzipped form}`. The data file are suggested to procecss with triming software for adaptor and low qulaity base

##Methodology:

Please check methodology.pdf for detail of the script. Paper is not yet publish, link will be updated as soon as possible.

## Reference
[1] Acevedo, A., Brodsky, L., & Andino, R. (2014). Mutational and fitness landscapes of an RNA virus revealed through population sequencing. Nature, 505, pp.686-690.

## Contact
For technical questions, please contact TinHang via email: thchong@connect.ust.hk


# ONTbyAHR
Generate putative ONT reads with the help of query\`s assemblies, query\`s PacBio HiFi reads and reference genome.

## Overview
purge_dups is designed to remove haplotigs and contig overlaps in a de novo assembly based on read depth.
You can follow the Usage part and use our pipeline to purge your assembly or go to the Pipeline Guide to build your own pipeline.

## Dependencies
1. minimap2
2. racon
3. ragtag.py
4. nucmer
7. awk
8. hifiasm
9. samtools 

## Installation
Run the following commands to intall purge_dups (required):
```
1. wget http://cbi.hzau.edu.cn/ppsPCP/files/ppsPCP.zip or git clone https://github.com/wzxie/ONTbyAHR.git
2. export PATH=/path/to/ppsPCP/bin/:$PATH
```

# Usage
echo "Usage: l3-4.sh [-r reference genome] [-q query genome] [-c ccs reads] [-o output directory] [-t number of threads] [-n number of polishing rounds]
Options:
    Input/Output:
    	-r		reference genome
    	-q		query genome
    	-c		ccs reads
    	-o		output directory
    Polish:
    	-n INT		Number of Polishing Rounds [1-5]
    Supp:
    	-t INT		number of threads [1]
    	-version	show version number
See more information at https://github.com/***.\n

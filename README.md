# ONTbyAHR
Generate putative ONT reads with the help of query\`s assemblies, query\`s PacBio HiFi reads and reference genome.

## Overview
ONTbyAHR is capable of producing fake ONT reads by utilizing the query's assembly and query`s PacBio HiFi reads, in conjunction with a reference genome. You can follow the Usage part and use our pipeline to do it.
![test](https://github.com/wzxie/ONTbyAHR/assets/42645873/412f6cf6-b312-4b8c-9877-758e204f4996)

## Dependencies
1. minimap2 (2.22-r1101)
2. racon (v1.4.20)
3. ragtag.py (v2.1.0)
4. nucmer (4.0.0rc1)
7. awk (GNU Awk 4.0.2)
8. hifiasm (0.19.4-r575)
9. samtools (Version: 1.11)

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

# ONTbyAHR
Generate putative ONT reads with the help of query\`s assemblies, query\`s PacBio HiFi reads and reference genome.

## Overview
ONTbyAHR is capable of producing fake ONT reads by utilizing the query's assembly and query`s PacBio HiFi reads, in conjunction with a reference genome. You can follow the Usage part and use our pipeline to do it.
![image](https://github.com/wzxie/ONTbyAHR/assets/42645873/524c6b5e-7373-493a-b6ab-325785a5ff99)
![66f7079853b07c21475ab242bb889b4](https://github.com/wzxie/ONTbyAHR/assets/42645873/9a56b41c-6d5b-462a-b521-c45594b0a764)

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

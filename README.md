# ONTbyAHR
Generate putative ONT reads with the help of query\`s assemblies, query\`s PacBio HiFi reads and reference genome.

## Overview
ONTbyAHR is capable of producing fake ONT reads by utilizing the query's assembly and query`s PacBio HiFi reads, in conjunction with a reference genome. You can follow the Usage part and use our pipeline to do it.
[基因组质量提升示意图.pdf](https://github.com/wzxie/ONTbyAHR/files/12923535/default.pdf)

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
1. wget http://cbi.gxu.edu.cn/ONTbyAHR/files/ONTbyAHR.zip or git clone https://github.com/wzxie/ONTbyAHR.git
2. export PATH=/path/to/ONTbyAHR/bin/:$PATH
```

## Usage
```
echo "Usage: l3-4.sh [-r reference genome] [-q query genome] [-c ccs reads] [-o output directory] [-t number of threads] [-n number of polishing rounds]
Options:
    Input/Output:
    	-r		    reference genome
    	-q		    query genome
    	-c		    ccs reads
    	-o		    output directory
    Polish:
    	-n INT		Number of Polishing Rounds [1-5]
    Supp:
    	-t INT		number of threads [1]
    	-version	show version number
```

## Inputs and Outputs
### Input files
The reference genome sequence file
```
>Chr01
ATCGATCGATCGATCGATCGATCG...
```
The query genome sequence file
```
>Contig1
ATCGATCGATCGATCGATCGATCG...
```
The query PacBio HiFi sequence file
```
@*/ccs
ATCGATCGATCGATCGATCGATCG...
+
~{m~pI~~e~P~c~~o~[o~q_~~...
@*/ccs
ATCGATCGATCGATCGATCGATCG...
+
=<G~~n~~~~V~~~~~~~}~~~~~...
```
### Output files
* The 'gapA_area.svg' allows for visual observation of the extracted regions from the reference.
* The 'ontAlt_qry.fa' contains the final synthetic ONT reads used to aid in subsequent assembly.
* The 'ontAlt_qry.svg' provides a simple statistics of the 'ontAlt_qry.fa'.

## Example
```
```

## Contact
We hope this pipeline could be helpful for the groups which focused on plants genome assembly, you can use the github page to report issues or email us with any suggestions.
* Zhao rupeng:   2247290650@qq.com

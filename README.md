# RAGA
Generate putative long reads with the help of query\`s assemblies, query\`s PacBio HiFi reads and reference genome.

## Overview
RAGA is capable of producing fake long reads by utilizing the query's assembly and query`s PacBio HiFi reads, in conjunction with a reference genome. You can follow the Usage part and use our pipeline to do it.

![示意图](https://github.com/wzxie/RAGA/blob/main/RAGA%E7%A4%BA%E6%84%8F%E5%9B%BE.jpg)

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
1. git clone https://github.com/wzxie/RAGA.git
2. export PATH=/path/to/RAGA/bin/:$PATH
```

## Usage
```
Usage: l3-8.sh [-r reference genome] [-q query genome] [-c ccs reads] [options]
Options:
    Input/Output:
    -r          reference genome
    -q          query genome
    -c          ccs reads
    -o          output directory

    Polish:
    -n INT      Number of Polishing Rounds [>=3], default 3

    Filter:
    -i FLOAT    Set the minimum alignment identity [0, 100], default 90
    -l INT      Set the minimum alignment length, default 20,000
    -p FLOAT    Extract the PacBio HiFi read which align length is >= *% of its own length [0-1], default 0.9
    -P FLOAT    Extract the query longAlt read which align length is >= *% of its own length [0-1), default 0.5

    Supp:
    -t INT      number of threads, default 1
    -v|-version show version number
    -h|-help    show help information

See more information at https://github.com/wzxie/RAGA.
```

```
Usage: l4-6.sh [-r reference genome] [-c ccs reads] [options]
Options:
    Input/Output:
    -r          reference genome
    -c          ccs reads
    -o          output directory

    Polish:
    -n INT      Number of Polishing Rounds [>=10], default 10

    Filter:
    -l INT      Set the minimum alignment length, default 10,000

    Supp:
    -t INT      number of threads, default 1
    -v|-version show version number
    -h|-help    show help information

See more information at https://github.com/wzxie/RAGA.
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
* The 'longAlt_qry.fa' contains the final synthetic long reads used to aid in subsequent assembly.
* The 'longAlt_qry_lenDis.svg' provides a simple statistics of the 'longAlt_qry.fa'.

## Example
```
```

## Contact
We hope this pipeline could be helpful for the groups which focused on plants genome assembly, you can use the github page to report issues or email us with any suggestions.
* Zhao rupeng:    2247290650@qq.com
* Song jiaming:   jmsong@gxu.edu.cn
* Chen lingling:  llchen@gxu.edu.cn

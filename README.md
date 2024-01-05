# RAGA
Generate putative long reads with the help of reference genome, source\`s assemblies, source\`s PacBio HiFi reads.

## Overview
RAGA is capable of producing fake long reads by utilizing the source's assembly and source`s PacBio HiFi reads, in conjunction with a reference genome. You can follow the Usage part and use our pipeline to do it.

![workflow](https://github.com/wzxie/RAGA/blob/main/workflow.jpg)

## Dependencies
1. minimap2 (2.22-r1101)
2. racon (v1.4.20)
3. ragtag.py (v2.1.0)
4. nucmer (4.0.0rc1)
5. awk (GNU Awk 4.0.2)
6. hifiasm (0.19.6-r595)
7. samtools (Version: 1.11)
8. bedtools (v2.26.0)

## Installation
Run the following commands to intall RAGA (required):
```
1. git clone https://github.com/wzxie/RAGA.git
2. chmod 755 /path/to/RAGA/bin/*
3. export PATH=/path/to/RAGA/bin/:$PATH
```

## Usage
```
Usage: RAGA-same.sh [-r reference genome] [-q source genome] [-c ccs reads] [options]
Options:
    Input/Output:
    -r          reference genome
    -q          source genome
    -c          ccs reads
    -o          output directory

    Polish:
    -n INT      Number of Polishing Rounds [>=3], default 3

    Filter:
    -i FLOAT    Set the minimum alignment identity [0, 100], default 90
    -l INT      Set the minimum alignment length, default 20,000
    -p FLOAT    Extract the PacBio HiFi read which align length is >= *% of its own length [0-1], default 0.9
    -P FLOAT    Extract the source longAlt read which align length is >= *% of its own length [0-1), default 0.5

    Supp:
    -t INT      number of threads, default 1
    -v|-version show version number
    -h|-help    show help information

See more information at https://github.com/wzxie/RAGA.
```

```
Usage: RAGA-diff.sh [-r reference genome] [-c ccs reads] [options]
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
The source genome sequence file
```
>Contig1
ATCGATCGATCGATCGATCGATCG...
```
The source PacBio HiFi sequence file
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
* The 'longAlt_sur.fa' contains the final synthetic long reads used to aid in subsequent assembly.
* The 'longAlt_sur_lenDis.svg' provides a simple statistics of the 'longAlt_qry.fa'.

## Example
```
```

## Note
* When the input is the reference genome of the same species, RAGA needs to compare the reference genome with the source contigs to determine the gaps location and alignment block information. Therefore, the higher the quality of the input source assembly, the more reliable the output of RAGA is.
* Because the source of the non-related high-quality reference genome usually has a large and complex genome, it may take a lot of time to carry out an assembly. Therefore, when designing the process of RAGA input for the reference genome of the same species, the whole operation process of RAGA does not need source assembly. Therefore, users can directly take the source HiFi reads and the reference genome fasta as input, and the long sequence output by RAGA directly participates in the source assembly.

## Contact
We hope this pipeline could be helpful for the groups which focused on plants genome assembly, you can use the github page to report issues or email us with any suggestions.
* Zhao rupeng:    2247290650@qq.com
* Song jiaming:   jmsong@gxu.edu.cn
* Chen lingling:  llchen@gxu.edu.cn

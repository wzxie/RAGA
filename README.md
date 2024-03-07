# RAGA
RAGA is a pipeline designed to generate alternative long reads using the reference genome, the source\'s assembly, and the source\'s PacBio HiFi reads.

## Overview
RAGA is capable of producing alternative long reads by utilizing the source\'s assembly and the source\`s PacBio HiFi reads, in conjunction with a reference genome. You can follow the "Usage" part and use our pipeline to do it.

![workflow](https://github.com/wzxie/RAGA/blob/main/workflow.jpg)

## Dependencies
1. minimap2 (2.22-r1101)
2. racon (v1.4.20)
3. ragtag.py (v2.1.0)
4. nucmer (4.0.0rc1)
5. awk (GNU Awk 4.0.2)
6. hifiasm (0.19.6-r595 or higher)
7. samtools (Version: 1.11)
8. bedtools (v2.26.0)
9. seqkit (Version: 2.5.1)

## Installation
Run the following commands to intall RAGA (required):
```
1. git clone https://github.com/wzxie/RAGA.git
2. chmod 755 /path/to/RAGA/bin/*
3. export PATH=/path/to/RAGA/bin/:$PATH
```

## Usage
### A. Same Species as Reference.
```
Usage: RAGA-same.sh [-r reference genome] [-q source assembly] [-c source PacBio HiFi reads] [options]
Options:
    Input/Output:
    -r          reference genome
    -q          source assembly
    -c          source PacBio HiFi reads
    -o          output directory

    Polish:
    -n INT      number of Polishing Rounds [>=3], default 3

    Filter:
    -i FLOAT    set the minimum alignment identity [0, 100], default 90
    -l INT      set the minimum alignment length, default 20,000
    -p FLOAT    extract the source PacBio HiFi read which align length is >= *% of its own length [0-1], default 0.9
    -P FLOAT    extract the source longAlt read which align length is >= *% of its own length [0-1), default 0.5

    Supp:
    -t INT      number of threads, default 1
    -v|-version show version number
    -h|-help    show help information

See more information at https://github.com/wzxie/RAGA.
```

### B. Different species as reference.
```
Usage: RAGA-diff.sh [-r reference genome] [-c source PacBio HiFi reads] [options]
Options:
    Input/Output:
    -r          reference genome
    -c          source PacBio HiFi reads
    -o          output directory

    Polish:
    -n INT      number of Polishing Rounds [>=10], default 10

    Filter:
    -l INT      set the minimum alignment length, default 10,000

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
The source assembly sequence file
```
>Contig1
ATCGATCGATCGATCGATCGATCG...
>Contig2
ATCGATCGATCGATCGATCGATCG...
```
The source PacBio HiFi sequence file
```
@*/ccs1
ATCGATCGATCGATCGATCGATCG...
+
~{m~pI~~e~P~c~~o~[o~q_~~...
@*/ccs2
ATCGATCGATCGATCGATCGATCG...
+
=<G~~n~~~~V~~~~~~~}~~~~~...
```
### Output files
* The 'gapA_area.svg' allows for visual observation of the extracted regions from the reference.
* The 'longAlt_sur.fa' contains the final alternative long reads used to aid in subsequent assembly.
* The 'longAlt_sur_lenDis.svg' provides a simple statistics of the 'longAlt_qry.fa'.

## Example
### 1. Download example data
```
wget https://cncb-gsa.obs.cn-north-4.myhuaweicloud.com/data/gsapub/CRA008584/CRR591673/CRR591673.fastq.gz
wget https://github.com/schatzlab/Col-CEN/blob/main/v1.2/Col-CEN_v1.2.fasta.gz
gzip -d CRR591673.fastq.gz
gzip -d Col-CEN_v1.2.fasta.gz
```
### 2. De novo assembly
```
hifiasm -o denovo -t 6 --primary CRR591673.fastq
awk '/^S/{print ">"$2;print $3}' denovo.p_ctg.gfa > denovo.fa
```
### 3. Run RAGA pipeline
```
RAGA-same.sh -r Col-CEN_v1.2.fasta -q denovo.fa -c CRR591673.fastq -o output_same -t 6 -n 3 -i 90 -l 20000 -p 0.9 -P 0.5 &> output_same.log
```
### 4. Re de novo assembly
```
hifiasm -o re-denovo -t 6 --primary --ul longAlt_sur.fa CRR591673.fastq
```

## Note
* RAGA aligns the reference with the source contigs to identify the alignment blocks near the gaps when the input is the reference of the same species. Hence, RAGA's output will be more reliable if the input source assembly is of higher quality.
* Because the source of the non-related high-quality reference genome usually has a large and complex genome, it may take a lot of time to carry out an assembly. Therefore, when designing the process of RAGA input for the reference genome of the same species, the whole operation process of RAGA does not need source assembly. Therefore, users can directly take the source PacBio HiFi reads and the reference genome fasta as input, and the alternative long sequences output by RAGA directly participates in the source assembly.

## Contact
We hope this pipeline could be helpful for the groups which focused on plants genome assembly, you can use the github page to report issues or email us with any suggestions.
* Zhao rupeng:    2247290650@qq.com
* Song jiaming:   jmsong@swu.edu.cn
* Chen lingling:  llchen@gxu.edu.cn

#!/bin/bash

######################################################################
## File Name: RAGA.sh (Tuesday August 13 10:30:00 2024);
######################################################################
set -o nounset
#bak=$IFS
#IFS=$'\n'


######################################################################
## How to Use RAGA: An Introduction;
######################################################################

# a. A function checks whether a file is followed by a parameter
check_required_file_params() {
	opt="$1"
	arg="$2"
	if [ -z "$arg" ]; then
		echo "ERROR: there is a parameter required for the option $opt"
		exit 1
	fi
}

# b. Parameter passing
ref=""
ccs=""
hic1=""
hic2=""
thr="1"
npr="3"
dfi="99"
dfl="20000"
per="0.9"
PER="0.5"
hep=""
ver=""
homo=""
hetero=""
haplotype=""
trio=""

while [[ $# -gt 0 ]]
do
	key="$1"
	case "$key" in
		-r)
		shift
		if [[ $# -gt 0 ]]; then
			arg="$1"
			check_required_file_params "$key" "$arg"
			ref="$arg"
		shift
		else
			echo "ERROR: there is a parameter required for the option $key"
			exit 1
		fi
		;;
		-c)
		shift
		if [[ $# -gt 0 ]]; then
			arg="$1"
			check_required_file_params "$key" "$arg"
			ccs="$arg"
		shift
		else
			echo "ERROR: there is a parameter required for the option $key"
			exit 1
		fi
		;;
		-hic1)
		shift
		if [[ $# -gt 0 ]]; then
			arg="$1"
			check_required_file_params "$key" "$arg"
			hic1="$arg"
		shift
		else
			echo "ERROR: there is a parameter required for the option $key"
			exit 1
		fi
		;;
		-hic2)
		shift
		if [[ $# -gt 0 ]]; then
			arg="$1"
			check_required_file_params "$key" "$arg"
			hic2="$arg"
		shift
		else
			echo "ERROR: there is a parameter required for the option $key"
			exit 1
		fi
		;;
		-t)
		thr="$2"
		shift
		shift
		;;
		-n)
		npr="$2"
		shift
		shift
		;;
		-i)
		dfi="$2"
		shift
		shift
		;;
		-l)
		dfl="$2"
		shift
		shift
		;;
		-p)
		per="$2"
		shift
		shift
		;;
		-P)
		PER="$2"
		shift
		shift
		;;
		-h|-help)
		hep="help"
		shift
		;;
		-v|-version)
		ver="version"
		shift
		;;
		-homo)
		homo="homo"
		shift
		;;
		-hetero)
		hetero="hetero"
		shift
		;;
		-haplotype)
		haplotype="haplotype"
		shift
		;;
		-trio)
		trio="trio"
		shift
		;;
		*)
		echo "Sorry, the parameter you provided does not exist."
		shift
		exit 1
		;;
	esac
done

# c. Result return
if [ "$hep" == "help" ]; then
	echo "Usage: RAGA.sh [-r reference genome] [-c target PacBio HiFi reads] [options]

Options:
	Input/Output:
	-r          reference genome
	-c          target PacBio HiFi reads
	-hic1       target Hi-C_r1 reads
	-hic2       target Hi-C_r2 reads

	Assembly type:
	-homo       assemble inbred/homozygous genomes
	-hetero     assemble heterozygous genomes
	-haplotype  generates a pair of haplotype-resolved assemblies with paired-end Hi-C reads

	Polish:
	-n INT      number of Polishing Rounds [>=3], default 3

	Filter:
	-i FLOAT    set the minimum alignment identity [0, 100], default 99
	-l INT      set the minimum alignment length, default 20,000
	-p FLOAT    extract the target PacBio HiFi read which align length is >= *% of its own length [0-1], default 0.9
	-P FLOAT    extract the target longAlt read which align length is >= *% of its own length [0-1), default 0.5

	Supp:
	-t INT      number of threads, default 1
	-v|-version show version number
	-h|-help    show help information

See more information at https://github.com/wzxie/RAGA.
	"
	exit 1

elif [ "$ver" == "version" ]; then
	echo "Version: 1.0.0"
	exit 1

else
	[[ $ref == "" ]] && echo -e "ERROR: path to reference genome not found, assign using -r." && exit 1
	[[ $ccs == "" ]] && echo -e "ERROR: path to target PacBio HiFi reads not found, assign using -c." && exit 1
	[[ $npr -lt 3 ]] && echo -e "ERROR: -n INT number of Polishing Rounds [>=3], default 3." && exit 1
	([[ $dfi -lt 0 ]] || [[ $dfi -gt 100 ]]) && echo -e "ERROR: -i FLOAT	set the minimum alignment identity [0, 100], default 99." && exit 1
fi


#=====================================================================
## step0.1: Check other scripts if they are found in the path;
#=====================================================================
echo -e "Step0: verifying the availability of related dependencies."

for scr in minimap2 racon ragtag.py nucmer delta-filter show-coords awk hifiasm samtools seqkit
do
	check=$(command -v $scr)
	if [ "$check" == "" ]; then
		echo -e "\tERROR: command $scr is NOT in you PATH. Please check."
		exit 1
	else
		echo -e "\t$scr is ok"
	fi
done

echo -e "\tAll dependencies have been checked.\n"


#=====================================================================
## step0.2: File name parser
#=====================================================================
#echo -e "Step1: file name parser."

ccsbase1=$(basename $ccs)
ccsbase2=${ccsbase1%.*}
if [ "$haplotype" == "haplotype" ]; then
	hic1base1=$(basename $hic1)
	hic2base1=$(basename $hic2)
	hic1base2=${hic1base1%.*}
	hic2base2=${hic2base1%.*}
fi


#=====================================================================
## step1: Initial Assembly
#=====================================================================
# echo -e "Step1: Initial Assembly."

# 1.1 ln -s
if [ -e $ccsbase1 ]; then
	cm=no
else
	ln -s $ccs $ccsbase1
	cm=yes
fi

if [ "$haplotype" == "haplotype" ]; then
	if [ -e $hic1base1 ]; then
		h1=no
	else
		ln -s $hic1 $hic1base1
		h1=yes
	fi
fi

if [ "$haplotype" == "haplotype" ]; then
	if [ -e $hic2base1 ]; then
		h2=no
	else
		ln -s $hic2 $hic2base1
		h2=yes
	fi
fi

# 1.2 initial assembly
mkdir Initial_assembly && cd Initial_assembly
if [ "$homo" == "homo" ]; then
	# a. Assemble inbred/homozygous genomes
	echo -e "Step1: Initial Assembly (assemble inbred/homozygous genome)."
	hifiasm -o contigs -l0 ../$ccsbase1 -t $thr
	[[ $? -eq 0 ]] || exit 1
	awk '/^S/{print ">"$2;print $3}' contigs.bp.p_ctg.gfa > contigs.bp.p_ctg.fa
	ln -s contigs.bp.p_ctg.fa initial.fa

elif [ "$hetero" == "hetero" ]; then
	# b. Assemble heterozygous genomes
	echo -e "Step1: Initial Assembly (assemble heterozygous genome)."
	hifiasm -o contigs ../$ccsbase1 -t $thr
	[[ $? -eq 0 ]] || exit 1
	awk '/^S/{print ">"$2;print $3}' contigs.bp.p_ctg.gfa > contigs.bp.p_ctg.fa
	ln -s contigs.bp.p_ctg.fa initial.fa

elif [ "$haplotype" == "haplotype" ]; then
	# c. generate a pair of haplotype-resolved assemblies with paired-end Hi-C reads
	echo -e "Step1: Initial Assembly (assemble haplotype-resolved genome with Hi-C reads)."
	hifiasm -o contigs --h1 ../$hic1base1 --h2 ../$hic2base1 ../$ccsbase1 -t $thr
	[[ $? -eq 0 ]] || exit 1
	awk '/^S/{print ">"$2;print $3}' contigs.hic.hap1.p_ctg.gfa > contigs.hic.hap1.p_ctg.fa
	awk '/^S/{print ">"$2;print $3}' contigs.hic.hap2.p_ctg.gfa > contigs.hic.hap2.p_ctg.fa
	cat contigs.hic.hap1.p_ctg.fa contigs.hic.hap2.p_ctg.fa > contigs.hic.hap.p_ctg.fa
	ln -s contigs.hic.hap.p_ctg.fa initial.fa

#elif [ "$trio" == "trio" ]; then
#	# d. Trio binning assembly
#	echo -e "Step1: Initial Assembly (assemble phased genome with Trio binning)."
#	yak count -b37 -t16 -o pat.yak <(cat pat_1.fq.gz pat_2.fq.gz) <(cat pat_1.fq.gz pat_2.fq.gz)
#	yak count -b37 -t16 -o mat.yak <(cat mat_1.fq.gz mat_2.fq.gz) <(cat mat_1.fq.gz mat_2.fq.gz)
#	hifiasm -o HG002.asm -t32 -1 pat.yak -2 mat.yak HG002-HiFi.fa.gz

else
	# e. default assembly
	echo -e "Step1: Initial Assembly (the default assembly mode 'homo' is used)."
	hifiasm -o contigs -l0 ../$ccsbase1 -t $thr
	[[ $? -eq 0 ]] || exit 1
	awk '/^S/{print ">"$2;print $3}' contigs.bp.p_ctg.gfa > contigs.bp.p_ctg.fa
	ln -s contigs.bp.p_ctg.fa initial.fa
fi

# 1.3 rm ln -s
if [ "$cm" = "yes" ]; then
	rm ../$ccsbase1
fi

if [ "$haplotype" == "haplotype" ]; then
	if [ "$h1" = "yes" ]; then
		rm ../$hic1base1
	fi
fi

if [ "$haplotype" == "haplotype" ]; then
	if [ "$h2" = "yes" ]; then
		rm ../$hic2base1
	fi
fi
unset cm h1 h2

# 1.4 $?
if [[ $? -eq 0 ]]; then
	echo -e "Step1 is done.\n"
else
	exit 1
fi
cd ..


#=====================================================================
## step2: Generate alternative long reads
#=====================================================================
echo -e "Step2: generate alternative long reads with a reference genome."

# 2.1
RAGA-same.sh -r $ref -q Initial_assembly/initial.fa -c $ccs -o Alternative_reads -n $npr -i $dfi -l $dfl -p $per -P $PER -t $thr -solo

# 2.2 $?
if [[ $? -eq 0 ]]; then
	echo -e "Step2 is done.\n"
else
	exit 1
fi


#=====================================================================
## step3: RAGA optimized assembly
#=====================================================================
#echo -e "Step3: RAGA optimized Assembly."

# 3.1 ln -s
if [ -e $ccsbase1 ]; then
	cm=no
else
	ln -s $ccs $ccsbase1
	cm=yes
fi

if [ "$haplotype" == "haplotype" ]; then
	if [ -e $hic1base1 ]; then
		h1=no
	else
		ln -s $hic1 $hic1base1
		h1=yes
	fi
fi

if [ "$haplotype" == "haplotype" ]; then
	if [ -e $hic2base1 ]; then
		h2=no
	else
		ln -s $hic2 $hic2base1
		h2=yes
	fi
fi

# 3.2 assembly with HiFi / Hi-C and longAlt_tgt.fa
mkdir Optimized_assembly && cd Optimized_assembly
if [ "$homo" == "homo" ]; then
	# a. Assemble inbred/homozygous genomes
	echo -e "Step3: Optimized Assembly (assemble inbred/homozygous genome)."
	hifiasm -o contigs --ul ../Alternative_reads/longAlt_tgt.fa ../$ccsbase1 -t $thr
	[[ $? -eq 0 ]] || exit 1
	awk '/^S/{print ">"$2;print $3}' contigs.bp.p_ctg.gfa > contigs.bp.p_ctg.fa
	ln -s contigs.bp.p_ctg.fa optimized.fa
#	ln -s contigs.bp.p_ctg.fa RAGA-optimized_contigs.fa

elif [ "$hetero" == "hetero" ]; then
	# b. Assemble heterozygous genomes
	echo -e "Step3: Optimized Assembly (assemble heterozygous genome)."
	hifiasm -o contigs --ul ../Alternative_reads/longAlt_tgt.fa ../$ccsbase1 -t $thr
	[[ $? -eq 0 ]] || exit 1
	awk '/^S/{print ">"$2;print $3}' contigs.bp.p_ctg.gfa > contigs.bp.p_ctg.fa
	ln -s contigs.bp.p_ctg.fa optimized.fa
#	ln -s contigs.bp.p_ctg.fa RAGA-optimized_contigs.fa

elif [ "$haplotype" == "haplotype" ]; then
	# c. Generate a pair of haplotype-resolved assemblies with paired-end Hi-C reads
	echo -e "Step3: Optimized Assembly (assemble haplotype-resolved genome with Hi-C reads)."
	hifiasm -o contigs --h1 ../$hic1base1 --h2 ../$hic2base1 --ul ../Alternative_reads/longAlt_tgt.fa ../$ccsbase1 -t $thr
	[[ $? -eq 0 ]] || exit 1
	awk '/^S/{print ">"$2;print $3}' contigs.hic.hap1.p_ctg.gfa > contigs.hic.hap1.p_ctg.fa
	awk '/^S/{print ">"$2;print $3}' contigs.hic.hap2.p_ctg.gfa > contigs.hic.hap2.p_ctg.fa
	ln -s contigs.hic.hap1.p_ctg.fa optimized_hap1.fa
	ln -s contigs.hic.hap2.p_ctg.fa optimized_hap2.fa
#	ln -s contigs.hic.hap1.p_ctg.fa RAGA-optimized_contigs_hap1.fa
#	ln -s contigs.hic.hap2.p_ctg.fa RAGA-optimized_contigs_hap2.fa

#elif [ "$trio" == "trio" ]; then
#	d. Trio binning assembly
#	echo -e "Step3: Optimized Assembly (assemble phased genome with Trio binning)."
#	yak count -b37 -t16 -o pat.yak <(cat pat_1.fq.gz pat_2.fq.gz) <(cat pat_1.fq.gz pat_2.fq.gz)
#	yak count -b37 -t16 -o mat.yak <(cat mat_1.fq.gz mat_2.fq.gz) <(cat mat_1.fq.gz mat_2.fq.gz)
#	hifiasm -o HG002.asm -t32 -1 pat.yak -2 mat.yak HG002-HiFi.fa.gz

else
	# e. default assembly
	echo -e "Step3: Optimized Assembly (the default assembly mode 'homo' is used)."
	hifiasm -o contigs --ul ../Alternative_reads/longAlt_tgt.fa ../$ccsbase1 -t $thr
	[[ $? -eq 0 ]] || exit 1
	awk '/^S/{print ">"$2;print $3}' contigs.bp.p_ctg.gfa > contigs.bp.p_ctg.fa
	ln -s contigs.bp.p_ctg.fa optimized.fa
#	ln -s contigs.bp.p_ctg.fa RAGA-optimized_contigs.fa
fi

# 3.3 chromosome scaffolding
#if [ "$homo" == "homo" ] || [ "$hetero" == "hetero" ]; then
#	ragtag.py scaffold ../$refbase1 RAGA-optimized_contigs.fa -t $thr -o .
#	[[ $? -eq 0 ]] || exit 1
#	ln -s ragtag.scaffold.fasta RAGA-optimized.fa

#elif [ "$haplotype" == "haplotype" ]; then
#	ragtag.py scaffold ../$refbase1 RAGA-optimized_contigs_hap1.fa -t $thr -o .
#	[[ $? -eq 0 ]] || exit 1
#	ln -s ragtag.scaffold.fasta RAGA-optimized_hap1.fa
#	ragtag.py scaffold ../$refbase1 RAGA-optimized_contigs_hap2.fa -t $thr -o .
#	ln -s ragtag.scaffold.fasta RAGA-optimized_hap2.fa
#fi

# 3.4 rm ln -s
if [ "$cm" = "yes" ]; then
	rm ../$ccsbase1
fi

if [ "$haplotype" == "haplotype" ]; then
	if [ "$h1" = "yes" ]; then
		rm ../$hic1base1
	fi
fi

if [ "$haplotype" == "haplotype" ]; then
	if [ "$h2" = "yes" ]; then
		rm ../$hic2base1
	fi
fi
unset cm h1 h2

# 3.5 $?
if [[ $? -eq 0 ]]; then
	echo -e "Step3 is done.\n"
else
	exit 1
fi
cd ..

[[ $? -eq 0 ]] && echo -e "Congratulations! The RAGA.sh is done.\n"


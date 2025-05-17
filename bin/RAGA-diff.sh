#!/bin/bash

######################################################################
# File Name: RAGA-diff.sh (Tuesday August 13 10:30:00 2024);
######################################################################
set -o nounset
#bak=$IFS
#IFS=$'\n'


######################################################################
## How to Use diff.sh: An Introduction;
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
out="output_diff"
thr="1"
npr="10"
dfl="10000"
hep=""
ver=""
solo=""

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
		-o)
		out="$2"
		shift
		shift
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
		-l)
		dfl="$2"
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
		-solo)
		solo="yes"
		shift
		;;
		*)
		echo "Sorry, the parameter you provided does not exist."
		shift
		exit
		;;
	esac
done

# c. Result return
if [ "$hep" == "help" ]; then
	echo "Usage: RAGA-diff.sh [-r reference genome] [-c target PacBio HiFi reads] [options]

Options:
	Input/Output:
	-r          reference genome
	-c          target PacBio HiFi reads
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
	"
	exit 1

elif [ "$ver" == "version" ]; then
	echo "Version: 1.0.0"
	exit 1

else
	[[ $ref == "" ]] && echo -e "ERROR: path to reference genome not found, assign using -r." && exit 1
	[[ $ccs == "" ]] && echo -e "ERROR: path to target PacBio HiFi reads not found, assign using -c." && exit 1
	[[ $npr -lt 10 ]] && echo -e "ERROR: number of Polishing Rounds [>=10], default 10." && exit 1
	[[ -d $out ]] && echo -e "ERROR: output directory already exists, please remove or set alternate output." && exit 1
	[[ -d $out ]] || mkdir $out
fi


#=====================================================================
## step0: Check other scripts if they are found in the path;
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Verifying the availability of related dependencies."

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
fi


#=====================================================================
## step1: File name parser
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step1: Create a symbolic link for the files to enable access."
else
	echo -e "2.1 Create a symbolic link for the files to enable access."
fi

refbase1=$(basename $ref)
ccsbase1=$(basename $ccs)
refbase2=${refbase1%.*}
ccsbase2=${ccsbase1%.*}

[[ ! -e ${refbase2}_racon0.fa ]] && ln -s $ref ${refbase2}_racon0.fa

if [ -e $ccsbase1 ]; then
	cm=no
else
	ln -s $ccs $ccsbase1
	cm=yes
fi

if [ "$solo" != "yes" ]; then
	echo -e "Step1 is done.\n"
else
	echo -e "2.1 is done.\n"
fi


#=====================================================================
## step2: polish
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step2: polish the $refbase1 $npr rounds using ${ccsbase1}."
else
	echo -e "2.2 polish the $refbase1 $npr rounds using ${ccsbase1}."
fi

cd $out
for i in $(seq 1 $npr)
do
	# racon
	flag=$(($i-1))
	if [ "$flag" == 0 ]; then
		minimap2 -ax map-hifi ../${refbase2}_racon${flag}.fa ../$ccsbase1 -t $thr > ${refbase2}_racon${flag}.sam
		racon ../$ccsbase1 ${refbase2}_racon${flag}.sam ../${refbase2}_racon${flag}.fa -t $thr > ${refbase2}_racon${i}.fa
	else
		minimap2 -ax map-hifi ${refbase2}_racon${flag}.fa ../$ccsbase1 -t $thr > ${refbase2}_racon${flag}.sam
		racon ../$ccsbase1 ${refbase2}_racon${flag}.sam ${refbase2}_racon${flag}.fa -t $thr > ${refbase2}_racon${i}.fa
	fi

	# $?
	if [[ $? -eq 0 ]]; then
		if [[ $i -eq 10 ]]; then
			echo -e "The round $i of polish is done."
		else
			echo -e "The round $i of polish is done.\n"
		fi
	else
		exit 1
	fi
done

seqkit fx2tab -l -n ${refbase2}_racon${npr}.fa | awk '{print $1"\t0\t"$5}' > ${refbase2}_racon${npr}_len.txt

# $?
if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step2 is done.\n"
	else
		echo -e "2.2 is done.\n"
	fi
else
	exit 1
fi
unset i flag


#=====================================================================
## step3: minimap2 tgtHiFi reads to ref, and get longAlt_ref.fa
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step3: minimap2 $ccsbase1 to ${refbase2}_racon${npr}.fa, then get longAlt_ref.fa."
else
	echo -e "2.3 minimap2 $ccsbase1 to ${refbase2}_racon${npr}.fa, then get longAlt_ref.fa."
fi

minimap2 -x map-hifi ${refbase2}_racon${npr}.fa ../$ccsbase1 -t $thr > ${refbase2}_racon${npr}To${ccsbase2}.paf
awk '$13~/tp:A:P/ && $10/$2>=0.99{print $6"\t"$8"\t"$9}' ${refbase2}_racon${npr}To${ccsbase2}.paf | sort -k1,1 -k2,2n > ${refbase2}_racon${npr}To${ccsbase2}.bed
bedtools merge -i ${refbase2}_racon${npr}To${ccsbase2}.bed | awk '$3-$2>="'$dfl'"{print $1"\t"$2"\t"$3}' > ${refbase2}_racon${npr}To${ccsbase2}_m.bed
[[ $? -eq 0 ]] || exit 1

mkdir longAlt_ref
declare -i flag=0
cat ${refbase2}_racon${npr}To${ccsbase2}_m.bed | while read id start end
do
	flag=$(($flag + 1))
	echo -e "$id\t$start\t$end" > longAlt_ref/block_${flag}.bed
	seqkit subseq --bed longAlt_ref/block_${flag}.bed ${refbase2}_racon${npr}.fa -j $thr | seqkit replace -p ":." -j $thr > longAlt_ref/block_${flag}.fa
	[[ $? -eq 0 ]] || exit 1
	rm longAlt_ref/block_${flag}.bed
done

if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step3 is done.\n"
	else
		echo -e "2.3 is done.\n"
	fi
else
	exit 1
fi
unset flag


#=====================================================================
## step4: minimap2 tgtHiFi reads to ref, and get ccsAlt_tgt.fq
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step4: minimap2 $ccsbase1 to ${refbase2}_racon${npr}.fa, then get ccsAlt_tgt.fq."
else
	echo -e "2.4 minimap2 $ccsbase1 to ${refbase2}_racon${npr}.fa, then get ccsAlt_tgt.fq."
fi

mkdir ccsAlt_tgt
declare -i flag=0
cat ${refbase2}_racon${npr}To${ccsbase2}_m.bed | while read id start end
do
	flag=$(($flag + 1))
	echo -e "$id\t$start\t$end" > ccsAlt_tgt/block_${flag}.bed
	awk 'ARGIND==1{a[$1];b[$1]=$2;c[$1]=$3} ARGIND==2{if($6 in a && $9>=b[$6] && $8<=c[$6] && $10/$2>=0.9 && $13~/tp:A:P/){print $1}}' ccsAlt_tgt/block_${flag}.bed ${refbase2}_racon${npr}To${ccsbase2}.paf | sort | uniq > ccsAlt_tgt/block_${flag}.id
	seqkit grep -f ccsAlt_tgt/block_${flag}.id ../$ccsbase1 -j $thr > ccsAlt_tgt/block_${flag}.fq
	[[ $? -eq 0 ]] || exit 1
	rm ccsAlt_tgt/block_${flag}.bed ccsAlt_tgt/block_${flag}.id
done

if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step4 is done.\n"
	else
		echo -e "2.4 is done.\n"
	fi
else
	exit 1
fi
unset flag


#=====================================================================
## step5: partial assembly
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step5: partial assembly with ccsAlt_tgt.fq and longAlt_ref.fa."
else
	echo -e "2.5 partial assembly with ccsAlt_tgt.fq and longAlt_ref.fa."
fi

mkdir partial_assembly
flag=$(wc -l ${refbase2}_racon${npr}To${ccsbase2}_m.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	if [[ -s longAlt_ref/block_${i}.fa ]] && [[ -s ccsAlt_tgt/block_${i}.fq ]]; then
		hifiasm --primary -o partial_assembly/block_${i} --ul longAlt_ref/block_${i}.fa ccsAlt_tgt/block_${i}.fq -t $thr
		if [[ $? -eq 0 ]]; then
			echo -e "The partial assembly using longAlt_ref/block_${i}.fa and ccsAlt_tgt/block_${i}.fq is done.\n"
		else
			echo -e "The partial assembly using longAlt_ref/block_${i}.fa and ccsAlt_tgt/block_${i}.fq is failed.\n"
			exit 1
		fi
	else
		echo -e "The longAlt_ref/block_${i}.fa or ccsAlt_tgt/block_${i}.fq is an empty file, skip.\n"
	fi
done

if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step5 is done.\n"
	else
		echo -e "2.5 is done.\n"
	fi
else
	exit 1
fi
unset flag i 


#=====================================================================
## step6: get final longAlt reads
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step6: Get final longAlt reads."
else
	echo -e "2.6 Get final longAlt reads."
fi

# get each block_${i}_pa.fa
mkdir longAlt_tgt
flag=$(wc -l ${refbase2}_racon${npr}To${ccsbase2}_m.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	if [[ -s longAlt_ref/block_${i}.fa ]] && [[ -s ccsAlt_tgt/block_${i}.fq ]]; then
		awk '/^S/{print ">"$2;print $3}' partial_assembly/block_${i}.p_ctg.gfa > longAlt_tgt/block_${i}.p_ctg.fa
		awk '/^S/{print ">"$2;print $3}' partial_assembly/block_${i}.a_ctg.gfa > longAlt_tgt/block_${i}.a_ctg.fa
		cat longAlt_tgt/block_${i}.p_ctg.fa longAlt_tgt/block_${i}.a_ctg.fa > longAlt_tgt/block_${i}_pa.fa
		seqkit seq -n longAlt_tgt/block_${i}_pa.fa -j $thr > longAlt_tgt/block_${i}_pa.id
		[[ $? -eq 0 ]] || exit 1
		rm longAlt_tgt/block_${i}.p_ctg.fa longAlt_tgt/block_${i}.a_ctg.fa
	fi
done
unset flag i

# filter longAlt_tgt.fa by reads depth
flag=$(wc -l ${refbase2}_racon${npr}To${ccsbase2}_m.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	if [[ -s longAlt_ref/block_${i}.fa ]] && [[ -s ccsAlt_tgt/block_${i}.fq ]]; then
		minimap2 -ax map-hifi longAlt_tgt/block_${i}_pa.fa ccsAlt_tgt/block_${i}.fq -t $thr | samtools view -bS - | samtools sort -o longAlt_tgt/block_${i}_s.bam -
		samtools depth -a longAlt_tgt/block_${i}_s.bam > longAlt_tgt/block_${i}_dep.txt
		awk '$3==0' longAlt_tgt/block_${i}_dep.txt | awk '{print $1}' | sort | uniq > longAlt_tgt/block_${i}_dep0.id
		awk 'ARGIND==1{a[$1]} ARGIND==2{if(!($1 in a)){print $1}}' longAlt_tgt/block_${i}_dep0.id longAlt_tgt/block_${i}_pa.id > longAlt_tgt/block_${i}_pa_f1.id
		seqkit grep -f longAlt_tgt/block_${i}_pa_f1.id longAlt_tgt/block_${i}_pa.fa >> longAlt_tgt/block_pa_f1.fa
		[[ $? -eq 0 ]] || exit 1
	fi
done
unset flag i

# filter longAlt_tgt.fa by reads length
seqkit seq -m $dfl longAlt_tgt/block_pa_f1.fa | seqkit fx2tab -l -j $thr | awk '{print ">longAlt_"NR; print $2}' > longAlt_tgt/block_pa_f2.fa
[[ $? -eq 0 ]] || exit 1

ln -s longAlt_tgt/block_pa_f2.fa longAlt_tgt.fa

if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step6 is done.\n"
	else
		echo -e "2.6 is done.\n"
	fi
else
	exit 1
fi


#=====================================================================
## step7: Visualization of result statistics
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step7: Visualization of result statistics."
else
	echo -e "2.7 Visualization of result statistics."
fi

pl_lenDis.pl -f longAlt_tgt.fa -split-length $dfl -o longAlt_tgt_lenDis.svg

if [[ $? -eq 0 ]]; then
	echo -e "The longAlt_tgt_lenDis.svg is done!"
else
	echo -e "The longAlt_tgt_lenDis.svg is failed."
#	exit 1
fi

# $?
if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step7 is done.\n"
	else
		echo -e "2.7 is done.\n"
	fi
else
	exit 1
fi


#=====================================================================
## step8: rm and mv result to ouput
#=====================================================================
rm ../${refbase2}_racon0.fa

if [ "$cm" = "yes" ]; then
	rm ../$ccsbase1
fi
unset cm
cd ..

#mv *racon* ccsAlt_* part* *longAlt_* $out

if [ "$solo" != "yes" ]; then
	[[ $? -eq 0 ]] && echo -e "Congratulations! The RAGA-diff.sh is done.\n"
fi


#!/bin/bash

######################################################################
## File Name: RAGA_ngs.sh (Saturday May 17 09:30:00 2025);
######################################################################
set -o nounset
#bak=$IFS
#IFS=$'\n'


######################################################################
## How to Use RAGA-ngs.sh: An Introduction;
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
ref=""							# ok
read1=""						# ok
read2=""						# ok
out="output_ngs"				# ok
thr="1"							# ok
npr="3"							# ok
hep=""							# ok
ver=""							# ok
solo=""							# ok

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
		-read1)
		shift
		if [[ $# -gt 0 ]]; then
			arg="$1"
			check_required_file_params "$key" "$arg"
			read1="$arg"
		shift
		else
			echo "ERROR: there is a parameter required for the option $key"
			exit 1
		fi
		;;
		-read2)
		shift
		if [[ $# -gt 0 ]]; then
			arg="$1"
			check_required_file_params "$key" "$arg"
			read2="$arg"
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
	echo "Usage: RAGA-ngs.sh [-r reference genome] [-read1 NGS_1.fq] [-read2 NGS_2.fq] [options]

Options:
	Input/Output:
	-r          reference genome
	-read1      files with #1 mates, paired with files in <m2>
	-read2      files with #2 mates, paired with files in <m1>
	-o          output directory

	Polish:
	-n INT      number of Polishing Rounds [>=3], default 3

	Filter:

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
	[[ $read1 == "" ]] && echo -e "ERROR: path to NGS_read1.fq not found, assign using -read1." && exit 1
	[[ $read2 == "" ]] && echo -e "ERROR: path to NGS_read2.fq not found, assign using -read2." && exit 1
	[[ $npr -lt 3 ]] && echo -e "ERROR: -n INT number of Polishing Rounds [>=3], default 3." && exit 1
	[[ -d $out ]] && echo -e "ERROR: output directory already exists, please remove or set alternate output." && exit 1
	[[ -d $out ]] || mkdir $out
fi


#=====================================================================
## step0: Check other scripts if they are found in the path;
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Verifying the availability of related dependencies."

	for scr in bwa samtools pilon seqkit grep spades.py
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
read1base1=$(basename $read1)
read2base1=$(basename $read2)
refbase2=${refbase1%.*}
read1base2=${read1base1%.*}
read2base2=${read2base1%.*}

[[ ! -e ${refbase2}_pilon0.fa ]] && ln -s $ref ${refbase2}_pilon0.fa

if [ -e $read1base1 ]; then
	m1=no
else
	ln -s $read1 $read1base1
	m1=yes
fi

if [ -e $read2base1 ]; then
	m2=no
else
	ln -s $read2 $read2base1
	m2=yes
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
	echo -e "Step2: polish the $refbase1 $npr rounds using ${read1base1} and ${read2base1}."
else
	echo -e "2.2 polish the $refbase1 $npr rounds using ${ccsbase1}."
fi

pilon=$(find ~ -name pilon-1.24.jar | awk 'NR==1{print $1}')

cd $out
for i in $(seq 1 $npr)
do
	# pilon
	flag=$(($i-1))
	if [ "$flag" == 0 ]; then
		bwa index ../${refbase2}_pilon${flag}.fa
		bwa mem -t $thr ../${refbase2}_pilon${flag}.fa ../$read1base1 ../$read2base1 > ${refbase2}_pilon${flag}.sam
		samtools view -S -b -h -q 20 -@ $thr ${refbase2}_pilon${flag}.sam | samtools sort -@ $thr > ${refbase2}_pilon${flag}_s.bam
		samtools index ${refbase2}_pilon${flag}_s.bam
		java -Xmx100G -jar $pilon --genome ../${refbase2}_pilon${flag}.fa --fix all --bam ${refbase2}_pilon${flag}_s.bam --output ${refbase2}_pilon${i} --outdir .
		mv ${refbase2}_pilon${i}.fasta ${refbase2}_pilon${i}.fa
#		rm *.sam *.bam *.bai ../*.bwt ../*.pac ../*.ann ../*.amb ../*.sa

	else
		bwa index ${refbase2}_pilon${flag}.fa
		bwa mem -t $thr ${refbase2}_pilon${flag}.fa ../$read1base1 ../$read2base1 > ${refbase2}_pilon${flag}.sam
		samtools view -S -b -h -q 20 -@ $thr ${refbase2}_pilon${flag}.sam | samtools sort -@ $thr > ${refbase2}_pilon${flag}_s.bam
		samtools index ${refbase2}_pilon${flag}_s.bam
		java -Xmx100G -jar $pilon --genome ${refbase2}_pilon${flag}.fa --fix all --bam ${refbase2}_pilon${flag}_s.bam --output ${refbase2}_pilon${i} --outdir .
		mv ${refbase2}_pilon${i}.fasta ${refbase2}_pilon${i}.fa
#		rm *.sam *.bam *.bai *.bwt *.pac *.ann *.amb *.sa
	fi

	# $?
	if [[ $? -eq 0 ]]; then
		if [[ $i -eq 3 ]]; then
			echo -e "The round $i of polish is done."
		else
			echo -e "The round $i of polish is done.\n"
		fi
	else
		exit 1
	fi
done

if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step2 is done.\n"
	else
		echo -e "2.3 is done.\n"
	fi
else
	exit 1
fi
unset i flag


#=====================================================================
## step3: get filter reads within each window by reads depth
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step3: get filtered read1.fq and read2.fq."
else
	echo -e "2.3 get filtered read1.fq and read2.fq."
fi

## get reads depth
bwa index ${refbase2}_pilon${npr}.fa
bwa mem -t $thr ${refbase2}_pilon${npr}.fa ../$read1base1 ../$read2base1 > ${refbase2}_pilon${npr}.sam
samtools view -S -b -@ $thr -q 20 ${refbase2}_pilon${npr}.sam | samtools sort -@ $thr > ${refbase2}_pilon${npr}_s.bam
samtools index ${refbase2}_pilon${npr}_s.bam

awk '$5==60{print $1"\t"$3"\t"$4"\t"$5}' ${refbase2}_pilon${npr}.sam > ${refbase2}_pilon${npr}_q60.txt

seqkit fx2tab -l -n ${refbase2}_pilon${npr}.fa | awk 'BEGIN{FS="\t"}{print $1}' > ${refbase2}_pilon${npr}.id
seqkit fx2tab -l -n ${refbase2}_pilon${npr}.fa | awk 'BEGIN{FS="\t"}{print $NF}' > ${refbase2}_pilon${npr}.len
paste ${refbase2}_pilon${npr}.id ${refbase2}_pilon${npr}.len > ${refbase2}_pilon${npr}_len.txt

## get window reads
for i in 10 20 30 40 50
do
	mkdir -p block/${i}K
	bedtools makewindows -g ${refbase2}_pilon${npr}_len.txt -w ${i}000 -s ${i}000 > block/${i}K/${refbase2}_pilon${npr}_w${i}.bed
	samtools bedcov block/${i}K/${refbase2}_pilon${npr}_w${i}.bed ${refbase2}_pilon${npr}_s.bam | awk '{print $1"\t"$2"\t"$3"\t"$4/"'${i}'000"}' | awk '$4>=30{print $1"\t"$2"\t"$3}' > block/${i}K/${refbase2}_pilon${npr}_w${i}_f.dep

	declare -i flag=1
	cat block/${i}K/${refbase2}_pilon${npr}_w${i}_f.dep | while read chr start end
	do
		awk '$2=="'$chr'" && $3>="'$start'" && $3<="'$end'"{print $1}' ${refbase2}_pilon${npr}_q60.txt | sort | uniq > block/${i}K/block${flag}.id
		cat block/${i}K/block${flag}.id | while read id
		do
			seqkit grep -nrp $id ../$read1base1 >> block/${i}K/block${flag}_1.fq
			seqkit grep -nrp $id ../$read2base1 >> block/${i}K/block${flag}_2.fq
#			less ../$read1base1 | awk '$1~"'$id'"{print; count=3; next} count>0 {print; count--}' >> block/${i}K/block${flag}_1.fq
#			less ../$read2base1 | awk '$1~"'$id'"{print; count=3; next} count>0 {print; count--}' >> block/${i}K/block${flag}_2.fq
		done
		flag=$flag+1
	done
done

## $?
if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step3 is done.\n"
	else
		echo -e "2.3 is done.\n"
	fi
else
	exit 1
fi
#unset i flag


#=====================================================================
## step4: partial assembly
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step4: partial assembly."
else
	echo -e "2.7 partial assembly with ccsAlt_sur.fq and longAlt_ref.fa."
fi

for i in 10 20 30 40 50
do
#	mkdir -p partial_assembly/${i}K
	for j in $(seq 1 $((flag-1)))
	do
		spades.py -1 block/${i}K/block${j}_1.fq -2 block/${i}K/block${j}_2.fq -o partial_assembly/${i}K -t $thr
		mv partial_assembly/${i}K/contigs.fasta partial_assembly/${i}K/block${j}.fa
	done
done

# $?
if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step4 is done.\n"
	else
		echo -e "2.7 is done.\n"
	fi
else
	exit 1
fi
#unset flag i


#=====================================================================
## step5: get final longAlt_ngs reads
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step5: Get final longAlt reads."
else
	echo -e "2.8 Get final longAlt reads."
fi

## filter1 by depth 0
for i in 10 20 30 40 50
do
	mkdir -p longAlt_ngs/${i}K
	for j in $(seq 1 $((flag-1)))
	do
		bwa index partial_assembly/${i}K/block${j}.fa
		seqkit seq -n partial_assembly/${i}K/block${j}.fa > longAlt_ngs/${i}K/block${j}.id
		bwa mem -t $thr partial_assembly/${i}K/block${j}.fa block/${i}K/block${j}_1.fq block/${i}K/block${j}_2.fq > longAlt_ngs/${i}K/block${j}.sam
		samtools view -S -b -@ $thr longAlt_ngs/${i}K/block${j}.sam | samtools sort -@ $thr -o longAlt_ngs/${i}K/block${j}_s.bam
		samtools index longAlt_ngs/${i}K/block${j}_s.bam
		samtools depth -a longAlt_ngs/${i}K/block${j}_s.bam | awk '$3==0{print $1}' | sort | uniq > longAlt_ngs/${i}K/block${j}_dep0.id
		grep -v -f longAlt_ngs/${i}K/block${j}_dep0.id longAlt_ngs/${i}K/block${j}.id > longAlt_ngs/${i}K/block${j}_f.id
		seqkit grep -f longAlt_ngs/${i}K/block${j}_f.id partial_assembly/${i}K/block${j}.fa > longAlt_ngs/${i}K/block${j}_f.fa
	done

	# merge fa by 10k 20k 30k 40k 50k
	for h in $(ls longAlt_ngs/${i}K/*_f.fa)
	do
		cat $h >> longAlt_ngs/${i}K/block_t.fa
	done

	seqkit rename longAlt_ngs/${i}K/block_t.fa | seqkit seq -m 1000 > longAlt_ngs/${i}K/block.fa
	cat longAlt_ngs/${i}K/block.fa >> longAlt_ngs/block.fa
	rm longAlt_ngs/${i}K/block_t.fa
done

## filter3
less longAlt_ngs/block.fa | awk 'BEGIN{Cnt=0}{if($0~/>/){Cnt=Cnt+1;tmp=">seq"Cnt; $0=tmp; print$0}else{print$0}}' > longAlt_ngs/block_t1.fa

minimap2 -x asm5 -t $thr ${refbase2}_pilon${npr}.fa longAlt_ngs/block_t1.fa | grep "tp:A:P" > longAlt_ngs/all.paf
awk '{print $1"\t"$3"\t"$4}' longAlt_ngs/all.paf > longAlt_ngs/all-map.bed
seqkit subseq --bed longAlt_ngs/all-map.bed longAlt_ngs/block_t1.fa | seqkit seq -m 10000 | awk 'BEGIN{Cnt=0}{if($0~/>/){Cnt=Cnt+1;tmp=">seq"Cnt; $0=tmp; print$0}else{print$0}}' > longAlt_ngs/block_t2.fa

minimap2 -x asm20 -t $thr ../${refbase2}_pilon0.fa longAlt_ngs/block_t2.fa | grep "tp:A:P" > longAlt_ngs/all-ref.paf
awk '{print $1"\t"$3"\t"$4}' longAlt_ngs/all-ref.paf > longAlt_ngs/all-ref-map.bed
seqkit subseq --bed longAlt_ngs/all-ref-map.bed longAlt_ngs/block_t2.fa | seqkit seq -m 10000 |awk 'BEGIN{Cnt=0}{if($0~/>/){Cnt=Cnt+1;tmp=">seq"Cnt; $0=tmp; print$0}else{print$0}}' > longAlt_ngs/block.fa

ln -s longAlt_ngs/block.fa longAlt_ngs.fa

# $?
if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step5 is done.\n"
	else
		echo -e "2.8 is done.\n"
	fi
else
	exit 1
fi


#=====================================================================
## step6: rm and mv result to ouput
#=====================================================================
rm ../${refbase2}_racon0.fa

rm ../${refbase2}_pilon0.fa

if [ "$m1" = "yes" ]; then
	rm ../$read1base1
fi

if [ "$m2" = "yes" ]; then
	rm ../$read2base1
fi
unset m1 m2

rm ../*.bwt ../*.pac ../*.ann ../*.amb ../*.sa 
rm *.sam *.bam *.bai *.bwt *.pac *.ann *.amb *.sa *.id *.len
for i in 10 20 30 40 50
do
	rm block/${i}K/
done


cd ..

if [ "$solo" != "yes" ]; then
	[[ $? -eq 0 ]] && echo -e "Congratulations! The RAGA-ngs.sh is done.\n"
fi


rm ${refbase2}_pilon${npr}.id ${refbase2}_pilon${npr}.len ${refbase2}_pilon${npr}_len.txt ${refbase2}_pilon${npr}_w10.bed
rm *.sam *.bam *.bai ../*.bwt ../*.pac ../*.ann ../*.amb ../*.sa





















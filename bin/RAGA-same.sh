#!/bin/bash

######################################################################
## File Name: RAGA-same.sh (Tuesday August 13 10:30:00 2024);
######################################################################
set -o nounset
#bak=$IFS
#IFS=$'\n'


######################################################################
## How to Use same.sh: An Introduction;
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
qry=""
ccs=""
out="output_same"
thr="1"
npr="3"
dfi="99"
dfl="20000"
per="0.9"
PER="0.5"
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
		-q)
		shift
		if [[ $# -gt 0 ]]; then
			arg="$1"
			check_required_file_params "$key" "$arg"
			qry="$arg"
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
	echo "Usage: RAGA-same.sh [-r reference genome] [-q target assembly] [-c target PacBio HiFi reads] [options]

Options:
	Input/Output:
	-r          reference genome
	-q          target assembly
	-c          target PacBio HiFi reads
	-o          output directory

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
	[[ $qry == "" ]] && echo -e "ERROR: path to target assembly not found, assign using -q." && exit 1
	[[ $ccs == "" ]] && echo -e "ERROR: path to target PacBio HiFi reads not found, assign using -c." && exit 1
	[[ $npr -lt 3 ]] && echo -e "ERROR: -n INT number of Polishing Rounds [>=3], default 3." && exit 1
	([[ $dfi -lt 0 ]] || [[ $dfi -gt 100 ]]) && echo -e "ERROR: -i FLOAT	set the minimum alignment identity [0, 100], default 99." && exit 1
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
qrybase1=$(basename $qry)
ccsbase1=$(basename $ccs)
refbase2=${refbase1%.*}
qrybase2=${qrybase1%.*}
ccsbase2=${ccsbase1%.*}

[[ ! -e ${refbase2}_racon0.fa ]] && ln -s $ref ${refbase2}_racon0.fa

if [ -e $qrybase1 ]; then
	qm=no
else
	ln -s $qry $qrybase1
	qm=yes
fi

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
## step2: filter target`s hifi reads, which just used to step3
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step2: filter the $ccsbase1, which just used to Step3."
else
	echo -e "2.2 filter the $ccsbase1, which just used to 2.3."
fi

cd $out
minimap2 -x map-hifi ../$qrybase1 ../$ccsbase1 -t $thr > ${qrybase2}_${ccsbase2}.paf
less ${qrybase2}_${ccsbase2}.paf | awk '$10/$2>=0.90{print $1}' | sort | uniq > ${ccsbase2}_f.id
seqkit grep -f ${ccsbase2}_f.id ../$ccsbase1 > ${ccsbase2}_f.fq

if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step2 is done.\n"
	else
		echo -e "2.2 is done.\n"
	fi
else
	exit 1
fi


#=====================================================================
## step3: polish
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step3: polish the $refbase1 $npr rounds using ${ccsbase1}."
else
	echo -e "2.3 polish the $refbase1 $npr rounds using ${ccsbase1}."
fi

for i in $(seq 1 $npr)
do
	# racon
	flag=$(($i-1))
	if [ "$flag" == 0 ]; then
		minimap2 -ax map-hifi ../${refbase2}_racon${flag}.fa ${ccsbase2}_f.fq -t $thr > ${refbase2}_racon${flag}.sam
		racon ${ccsbase2}_f.fq ${refbase2}_racon${flag}.sam ../${refbase2}_racon${flag}.fa -t $thr > ${refbase2}_racon${i}.fa
	else
		minimap2 -ax map-hifi ${refbase2}_racon${flag}.fa ${ccsbase2}_f.fq -t $thr > ${refbase2}_racon${flag}.sam
		racon ${ccsbase2}_f.fq ${refbase2}_racon${flag}.sam ${refbase2}_racon${flag}.fa -t $thr > ${refbase2}_racon${i}.fa
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

seqkit fx2tab -l -n ${refbase2}_racon${npr}.fa | awk '{print $1"\t0\t"$5}' > ${refbase2}_racon${npr}_len.txt

if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step3 is done.\n"
	else
		echo -e "2.3 is done.\n"
	fi
else
	exit 1
fi
	
unset i flag


#=====================================================================
## step4: scaffolding, locate expanded gap area
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step4: homology-based scaffolding ${qrybase1} by ${refbase2}_racon${npr}.fa, then locate expanded gap area."
else
	echo -e "2.4 homology-based scaffolding ${qrybase1} by ${refbase2}_racon${npr}.fa, then locate expanded gap area."
fi

# scaffolding
ragtag.py scaffold ${refbase2}_racon${npr}.fa ../$qrybase1 -t $thr -o .
[[ $? -eq 0 ]] || exit 1
ln -s ragtag.scaffold.fasta ${qrybase2}_ragtag.fa

# get gap between ref_racon${npr}.fa and tgt_ragtag.fa
awk '$5=="U"{if($2-1000000>=0){print $1"\t"$2-1000000"\t"$3+1000000} else {print $1"\t0\t"$3+1000000}}' ragtag.scaffold.agp > gap_around.bed
if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step4 is done.\n"
	else
		echo -e "2.4 is done.\n"
	fi
else
	exit 1
fi

# get not ragtag contig as alternative long reads
#seqkit fx2tab -l -n ${qrybase2}_ragtag.fa -j $thr | awk '$1!~/_RagTag/ && $2>=20000{print $1}' > ${qrybase2}_ragtagN.id
#seqkit grep -f ${qrybase2}_ragtagN.id ${qrybase2}_ragtag.fa -j $thr > longAlt_surOwn.fa
#[[ $? -eq 0 ]] && echo -e "Get longAlt_surOwn.fa is done!\n"


#=====================================================================
## step5: locating the gap at ref, and get longAlt_ref.fa
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step5: locate expanded gap of ${refbase2}_racon${npr}.fa, then get longAlt_ref.fa."
else
	echo -e "2.5 locate expanded gap of ${refbase2}_racon${npr}.fa, then get longAlt_ref.fa."
fi

# get algin block
nucmer -p ${refbase2}_racon${npr}To${qrybase2}_ragtag ${refbase2}_racon${npr}.fa ${qrybase2}_ragtag.fa -t $thr
delta-filter -i $dfi -l $dfl ${refbase2}_racon${npr}To${qrybase2}_ragtag.delta > ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.delta
show-coords -b -B ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.delta | awk '{if($9<$10){print $1"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12} else{print $1"\t"$8"\t"$10"\t"$9"\t"$11"\t"$12}}' > ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.coords
[[ $? -eq 0 ]] || exit 1

# get gap.bed of ref; and get alternative long reads
awk 'ARGIND==1{a[$1];b[$1]=$2;c[$1]=$3} ARGIND==2{if($1 in a && $2"_RagTag" in a && $4>=b[$1] && $3<=c[$1]){print $2"\t"$5"\t"$6}}' gap_around.bed ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.coords > ${refbase2}_racon${npr}_gapA.bed
[[ $? -eq 0 ]] || exit 1

# longAlt_ref.fa
mkdir longAlt_ref
declare -i flag=0
cat ${refbase2}_racon${npr}_gapA.bed | while read id start end
do
	flag=$(($flag + 1))
	echo -e "$id\t$start\t$end" > longAlt_ref/gapA_${flag}.bed
	seqkit subseq --bed longAlt_ref/gapA_${flag}.bed ${refbase2}_racon${npr}.fa -j $thr | seqkit replace -p ":." -j $thr > longAlt_ref/gapA_${flag}.fa
	[[ $? -eq 0 ]] || exit 1
	rm longAlt_ref/gapA_${flag}.bed
done
unset flag

# visualization the location of gapA
pl_locate.pl -f ${refbase2}_racon${npr}_gapA.bed -l ${refbase2}_racon${npr}_len.txt -o ${refbase2}_racon${npr}_gapA.svg

if [[ $? -eq 0 ]]; then
	echo -e "The ${refbase2}_racon${npr}_gapA.svg is done!"
else
	echo -e "The ${refbase2}_racon${npr}_gapA.svg is failed."
#	exit 1
fi

# $?
if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step5 is done.\n"
	else
		echo -e "2.5 is done.\n"
	fi
else
	exit 1
fi


#=====================================================================
## step6: minimap2 tgtHiFi reads to ref, and get ccsAlt_tgt.fq
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step6: minimap2 $ccsbase1 to ${refbase2}_racon${npr}.fa, then get ccsAlt_tgt.fq."
else
	echo -e "2.6 minimap2 $ccsbase1 to ${refbase2}_racon${npr}.fa, then get ccsAlt_tgt.fq."
fi

minimap2 -x map-hifi ${refbase2}_racon${npr}.fa ../$ccsbase1 -t $thr > ${refbase2}_racon${npr}To${ccsbase2}.paf
[[ $? -eq 0 ]] || exit 1

# ccsAlt_tgt.fa
mkdir ccsAlt_tgt
declare -i flag=0
cat ${refbase2}_racon${npr}_gapA.bed | while read id start end
do
	flag=$(($flag + 1))
	echo -e "$id\t$start\t$end" > ccsAlt_tgt/gapA_${flag}.bed
	awk 'ARGIND==1{a[$1];b[$1]=$2;c[$1]=$3} ARGIND==2{if($6 in a && $9>=b[$6] && $8<=c[$6] && $10/$2>="'$per'" && $13~/tp:A:P/){print $1}}' ccsAlt_tgt/gapA_${flag}.bed ${refbase2}_racon${npr}To${ccsbase2}.paf | sort | uniq > ccsAlt_tgt/gapA_${flag}.id
	seqkit grep -f ccsAlt_tgt/gapA_${flag}.id ../$ccsbase1 -j $thr > ccsAlt_tgt/gapA_${flag}.fq
	[[ $? -eq 0 ]] || exit 1
	rm ccsAlt_tgt/gapA_${flag}.bed ccsAlt_tgt/gapA_${flag}.id
done

# $?
if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step6 is done.\n"
	else
		echo -e "2.6 is done.\n"
	fi
else
	exit 1
fi

unset flag


#=====================================================================
## step7: partial assembly
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step7: partial assembly with ccsAlt_tgt.fq and longAlt_ref.fa."
else
	echo -e "2.7 partial assembly with ccsAlt_tgt.fq and longAlt_ref.fa."
fi

# partial_assembly.fa
mkdir partial_assembly
flag=$(wc -l ${refbase2}_racon${npr}_gapA.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	if [[ -s longAlt_ref/gapA_${i}.fa ]] && [[ -s ccsAlt_tgt/gapA_${i}.fq ]]; then
		hifiasm --primary -o partial_assembly/gapA_${i} --ul longAlt_ref/gapA_${i}.fa ccsAlt_tgt/gapA_${i}.fq -t $thr

		if [[ $? -eq 0 ]]; then
			echo -e "The partial assembly using longAlt_ref/gapA_${i}.fa and ccsAlt_tgt/gapA_${i}.fq is done.\n"

		else
			echo -e "The partial assembly using longAlt_ref/gapA_${i}.fa and ccsAlt_tgt/gapA_${i}.fq is failed.\n"
			exit 1
		fi

	else
		echo -e "The longAlt_ref/gapA_${i}.fa or ccsAlt_tgt/gapA_${i}.fq is an empty file, skip.\n"
	fi
done

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
unset flag i


#=====================================================================
## step8: get final longAlt reads
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step8: Get final longAlt reads."
else
	echo -e "2.8 Get final longAlt reads."
fi

# longAlt_tgt.fa
mkdir longAlt_tgt
flag=$(wc -l ${refbase2}_racon${npr}_gapA.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	if [[ -s longAlt_ref/gapA_${i}.fa ]] && [[ -s ccsAlt_tgt/gapA_${i}.fq ]]; then
		awk '/^S/{print ">"$2;print $3}' partial_assembly/gapA_${i}.p_ctg.gfa > longAlt_tgt/gapA_${i}.p_ctg.fa
		awk '/^S/{print ">"$2;print $3}' partial_assembly/gapA_${i}.a_ctg.gfa > longAlt_tgt/gapA_${i}.a_ctg.fa
		cat longAlt_tgt/gapA_${i}.p_ctg.fa longAlt_tgt/gapA_${i}.a_ctg.fa > longAlt_tgt/gapA_${i}_pa.fa
		seqkit seq -n longAlt_tgt/gapA_${i}_pa.fa -j $thr > longAlt_tgt/gapA_${i}_pa.id
		[[ $? -eq 0 ]] || exit 1
		rm longAlt_tgt/gapA_${i}.p_ctg.fa longAlt_tgt/gapA_${i}.a_ctg.fa
	fi
done
unset flag i

# filter longAlt_tgt.fa by reads depth
flag=$(wc -l ${refbase2}_racon${npr}_gapA.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	if [[ -s longAlt_ref/gapA_${i}.fa ]] && [[ -s ccsAlt_tgt/gapA_${i}.fq ]]; then
		minimap2 -ax map-hifi longAlt_tgt/gapA_${i}_pa.fa ccsAlt_tgt/gapA_${i}.fq -t $thr | samtools view -bS - | samtools sort -o longAlt_tgt/gapA_${i}_s.bam -
		samtools depth -a longAlt_tgt/gapA_${i}_s.bam > longAlt_tgt/gapA_${i}_dep.txt
		awk '$3==0' longAlt_tgt/gapA_${i}_dep.txt | awk '{print $1}' | sort | uniq > longAlt_tgt/gapA_${i}_dep0.id
		awk 'ARGIND==1{a[$1]} ARGIND==2{if(!($1 in a)){print $1}}' longAlt_tgt/gapA_${i}_dep0.id longAlt_tgt/gapA_${i}_pa.id > longAlt_tgt/gapA_${i}_pa_f1.id
		seqkit grep -f longAlt_tgt/gapA_${i}_pa_f1.id longAlt_tgt/gapA_${i}_pa.fa >> longAlt_tgt/gapA_pa_f1_t.fa
		[[ $? -eq 0 ]] || exit 1
	fi
done
seqkit rename longAlt_tgt/gapA_pa_f1_t.fa > longAlt_tgt/gapA_pa_f1.fa
rm longAlt_tgt/gapA_pa_f1_t.fa
unset flag i

# filter longAlt_tgt.fa by reads coverage
minimap2 -x asm5 ../$qrybase1 longAlt_tgt/gapA_pa_f1.fa -t $thr > longAlt_tgt/${qrybase2}TogapA_pa_f1.paf
awk '$10/$2>="'$PER'" && $10/$2<1{print $1}' longAlt_tgt/${qrybase2}TogapA_pa_f1.paf | sort | uniq > longAlt_tgt/gapA_pa_f.id
seqkit grep -f longAlt_tgt/gapA_pa_f.id longAlt_tgt/gapA_pa_f1.fa > longAlt_tgt/gapA_pa_f2.fa
[[ $? -eq 0 ]] || exit 1

# filter longAlt_tgt.fa by reads length
seqkit seq -m $dfl -j $thr longAlt_tgt/gapA_pa_f2.fa | seqkit fx2tab -l -j $thr | awk '{print ">longAlt_"NR; print $2}' > longAlt_tgt/gapA_pa_f3.fa
[[ $? -eq 0 ]] || exit 1
ln -s longAlt_tgt/gapA_pa_f3.fa longAlt_tgt.fa

# $?
if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step8 is done.\n"
	else
		echo -e "2.8 is done.\n"
	fi
else
	exit 1
fi


#=====================================================================
## step9: Visualization of result statistics
#=====================================================================
if [ "$solo" != "yes" ]; then
	echo -e "Step9: Visualization of result statistics."
else
	echo -e "2.9 Visualization of result statistics."
fi

pl_lenDis.pl -f longAlt_tgt.fa -split-length $dfl -o longAlt_tgt_lenDis.svg

if [[ $? -eq 0 ]]; then
	echo -e "The longAlt_tgt_lenDis.svg is done."
else
	echo -e "The longAlt_tgt_lenDis.svg is failed."
#	exit 1
fi

# $?
if [[ $? -eq 0 ]]; then
	if [ "$solo" != "yes" ]; then
		echo -e "Step9 is done.\n"
	else
		echo -e "2.9 is done."
	fi
else
	exit 1
fi


#=====================================================================
## step10: rm and mv result to ouput
#=====================================================================
rm ../${refbase2}_racon0.fa

if [ "$qm" = "yes" ]; then
	rm ../$qrybase1
fi

if [ "$cm" = "yes" ]; then
	rm ../$ccsbase1
fi
unset qm cm

#mv ${ccsbase2}_f.id ${ccsbase2}_f.fq *racon* ragtag* ${qrybase2}.fa.fai ${qrybase2}_* part* longAlt_* ccsAlt_* gap* $out
mv ../${qrybase2}.fa.fai .
cd ..

if [ "$solo" != "yes" ]; then
	[[ $? -eq 0 ]] && echo -e "Congratulations! The RAGA-same.sh is done.\n"
fi


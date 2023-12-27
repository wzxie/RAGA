#!/bin/bash

######################################################################
## File Name: l3-8.sh (Thu Dec 25 10:30:00 2023);
######################################################################
set -o nounset
#bak=$IFS
#IFS=$'\n'


######################################################################
## Parameter;
######################################################################
ref=""
qry=""
ccs=""
out="output1"
thr="1"
npr="3"
dfi="90"
dfl="20000"
per="0.9"
PER="0.5"
hep=""
ver=""

while [[ $# -gt 0 ]]
do
	key="$1"
	case "$key" in
		-r)
		ref="$2"
		shift
		shift
		;;
		-q)
		qry="$2"
		shift
		shift
		;;
		-c)
		ccs="$2"
		shift
		shift
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
		*)
		echo "Sorry, the parameter you provided does not exist."
		shift
		exit
		;;
	esac
done


######################################################################
## An introduction to how to use this software;
######################################################################
if [ "$hep" == "help" ]; then
	echo "Usage: l3-8.sh [-r reference genome] [-q query genome] [-c ccs reads] [options]
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
	"
	exit

elif [ "$ver" == "version" ]; then
	echo "Version: 1.3.8"
	exit

else
	[[ $ref == "" ]] && echo -e "ERROR: Path to reference genome not found, assign using -r." && exit
	[[ $qry == "" ]] && echo -e "ERROR: Path to query genome not found, assign using -q." && exit
	[[ $ccs == "" ]] && echo -e "ERROR: PacBio HiFi reads of query not found, assign using -c." && exit
	[[ $npr -lt 3 ]] && echo -e "ERROR: -n INT	Number of Polishing Rounds [>=3], default 3." && exit
	([[ $dfi -lt 0 ]] || [[ $dfi -gt 100 ]]) && echo -e "ERROR: -i FLOAT	Set the minimum alignment identity [0, 100], default 90." && exit
	[[ -d $out ]] && echo -e "ERROR: Output directory already exists, please remove or set alternate output." && exit
	[[ -d $out ]] || mkdir $out
fi


######################################################################
## Check other scripts if they are found in path;
######################################################################
echo -e "Verifying the availability of related dependencies!"
for scr in minimap2 racon ragtag.py nucmer delta-filter show-coords awk hifiasm samtools
do
	check=$(command -v $scr)
	if [ "$check" == "" ]; then
		echo -e "\tERROR: Command $scr is NOT in you PATH. Please check."
		exit 1
	else
		echo -e "\t$scr is ok"
	fi
done
echo -e "All dependencies have been checked.\n"


#=====================================================================
## step0: file name parser
#=====================================================================
refbase1=$(basename $ref)
qrybase1=$(basename $qry)
ccsbase1=$(basename $ccs)
refbase2=${refbase1%.*}
qrybase2=${qrybase1%.*}
ccsbase2=${ccsbase1%.*}

[[ ! -e ${refbase2}_racon0.fa ]] && ln -s $ref ${refbase2}_racon0.fa
#[[ ! -e $qrybase1 ]] && ln -s $qry $qrybase1
#[[ ! -e $ccsbase1 ]] && ln -s $ccs $ccsbase1

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


#=====================================================================
## step1: filter query`s hifi reads, which just used to step2
#=====================================================================
echo -e "Step1: Filter the $ccsbase1, which just used to Step2!"
minimap2 -x map-hifi ${qrybase1} $ccsbase1 -t $thr > ${qrybase2}_${ccsbase2}.paf
less ${qrybase2}_${ccsbase2}.paf | awk '$10/$2>=0.90{print $1}' | sort | uniq > ${ccsbase2}_f.id
seqkit grep -f ${ccsbase2}_f.id $ccsbase1 > ${ccsbase2}_f.fq
[[ $? -eq 0 ]] && echo -e "Step1 is done!\n"


#=====================================================================
## step2: polish
#=====================================================================
echo -e "Step2: polish the $refbase1 $npr rounds using ${ccsbase1}!"
for i in $(seq 1 $npr)
do
	# racon
	flag=$(($i-1))
	minimap2 -ax map-hifi ${refbase2}_racon${flag}.fa ${ccsbase2}_f.fq -t $thr > ${refbase2}_racon${flag}.sam
	racon ${ccsbase2}_f.fq ${refbase2}_racon${flag}.sam ${refbase2}_racon${flag}.fa -t $thr > ${refbase2}_racon${i}.fa

	# $?
#	if [[ $? -eq 0 ]]; then
#		[[ $i == 1 ]] && echo -e "The round 1 of polish is done.\n"
#		[[ $i == 2 ]] && echo -e "The round 2 of polish is done.\n"
#		[[ $i == 3 ]] && echo -e "The round 3 of polish is done.\n"
#		[[ $i == 4 ]] && echo -e "The round 4 of polish is done.\n"
#		[[ $i == 5 ]] && echo -e "The round 5 of polish is done.\n"
#	fi

    # $?
    if [[ $? -eq 0 ]]; then
		echo -e "The round $i of polish is done.\n"
	fi
done
seqkit fx2tab -l -n ${refbase2}_racon${npr}.fa | awk '{print $1"\t0\t"$5}' > ${refbase2}_racon${npr}_len.txt
[[ $? -eq 0 ]] && echo -e "Step2 is done!\n"
unset i flag


#=====================================================================
## step3: scaffolding, locate expanded gap area
#=====================================================================
echo -e "Step3: homology-based scaffolding ${qrybase1} by ${refbase2}_racon${npr}.fa, then locate expanded gap area!"
# 3.1 scaffolding
ragtag.py scaffold ${refbase2}_racon${npr}.fa $qrybase1 -t $thr -o .
ln -s ragtag.scaffold.fasta ${qrybase2}_ragtag.fa

# 3.2 get gap between ref_racon${npr}.fa and qry_ragtag.fa
awk '$5=="U"{if($2-1000000>=0){print $1"\t"$2-1000000"\t"$3+1000000} else {print $1"\t0\t"$3+1000000}}' ragtag.scaffold.agp > gap_around.bed
#[[ $? -eq 0 ]] && echo -e "Scaffold and locate expanded gap of $qrybase1 is done!\n"
[[ $? -eq 0 ]] && echo -e "Step3 is done!\n"

# 3.3 get not ragtag contig as alternative long reads
#seqkit fx2tab -l -n ${qrybase2}_ragtag.fa -j $thr | awk '$1!~/_RagTag/ && $2>=20000{print $1}' > ${qrybase2}_ragtagN.id
#seqkit grep -f ${qrybase2}_ragtagN.id ${qrybase2}_ragtag.fa -j $thr > longAlt_qryOwn.fa
#[[ $? -eq 0 ]] && echo -e "Get longAlt_qryOwn.fa is done!\n"


#=====================================================================
## step4: locating the gap at ref, and get longAlt_ref.fa
#=====================================================================
echo -e "Step4: locate expanded gap of ${refbase2}_racon${npr}.fa, then get longAlt_ref.fa!"
# 4.1 get algin block
nucmer -p ${refbase2}_racon${npr}To${qrybase2}_ragtag ${refbase2}_racon${npr}.fa ${qrybase2}_ragtag.fa -t $thr
delta-filter -i $dfi -l $dfl ${refbase2}_racon${npr}To${qrybase2}_ragtag.delta > ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.delta
show-coords -b -B ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.delta | awk '{if($9<$10){print $1"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12} else{print $1"\t"$8"\t"$10"\t"$9"\t"$11"\t"$12}}' > ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.coords

# 4.2 get gap.bed of ref; and get alternative long reads
awk 'ARGIND==1{a[$1];b[$1]=$2;c[$1]=$3} ARGIND==2{if($1 in a && $2"_RagTag" in a && $4>=b[$1] && $3<=c[$1]){print $2"\t"$5"\t"$6}}' gap_around.bed ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.coords > ${refbase2}_racon${npr}_gapA.bed

mkdir longAlt_ref
declare -i flag=0
cat ${refbase2}_racon${npr}_gapA.bed | while read id start end
do
	flag=$(($flag + 1))
	echo -e "$id\t$start\t$end" > longAlt_ref/gapA_${flag}.bed
	seqkit subseq --bed longAlt_ref/gapA_${flag}.bed ${refbase2}_racon${npr}.fa -j $thr | seqkit replace -p ":." -j $thr > longAlt_ref/gapA_${flag}.fa
	rm longAlt_ref/gapA_${flag}.bed
done
unset flag

# 4.3 visualization the location of gapA
pl_locate.pl -f ${refbase2}_racon${npr}_gapA.bed -l ${refbase2}_racon${npr}_len.txt -o ${refbase2}_racon${npr}_gapA.svg
#[[ $? -eq 0 ]] && echo -e "Get longAlt_ref.fa is done!\n"
[[ $? -eq 0 ]] && echo -e "Step4 is done!\n"


#=====================================================================
## step5: minimap2 qryHifi reads to ref, and get ccsAlt_qry.fq
#=====================================================================
echo -e "Step5: minimap2 $ccsbase1 to ${refbase2}_racon${npr}.fa, then get ccsAlt_qry.fq!"
minimap2 -x map-hifi ${refbase2}_racon${npr}.fa $ccsbase1 -t $thr > ${refbase2}_racon${npr}To${ccsbase2}.paf

mkdir ccsAlt_qry
declare -i flag=0
cat ${refbase2}_racon${npr}_gapA.bed | while read id start end
do
	flag=$(($flag + 1))
	echo -e "$id\t$start\t$end" > ccsAlt_qry/gapA_${flag}.bed
	awk 'ARGIND==1{a[$1];b[$1]=$2;c[$1]=$3} ARGIND==2{if($6 in a && $9>=b[$6] && $8<=c[$6] && $10/$2>="'$per'" && $13~/tp:A:P/){print $1}}' ccsAlt_qry/gapA_${flag}.bed ${refbase2}_racon${npr}To${ccsbase2}.paf | sort | uniq > ccsAlt_qry/gapA_${flag}.id
	seqkit grep -f ccsAlt_qry/gapA_${flag}.id $ccsbase1 -j $thr > ccsAlt_qry/gapA_${flag}.fq
	rm ccsAlt_qry/gapA_${flag}.bed ccsAlt_qry/gapA_${flag}.id
done
#[[ $? -eq 0 ]] && echo -e "Get ccsAlt_qry.fq is done!\n"
[[ $? -eq 0 ]] && echo -e "Step5 is done!\n"
unset flag


#=====================================================================
## step6: partial assembly
#=====================================================================
echo -e "Step6: partial assembly with ccsAlt_qry.fq and longAlt_ref.fa!"
mkdir partial_assembly
flag=$(wc -l ${refbase2}_racon${npr}_gapA.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	hifiasm --primary -o partial_assembly/gapA_${i} --ul longAlt_ref/gapA_${i}.fa ccsAlt_qry/gapA_${i}.fq -t $thr
done
#[[ $? -eq 0 ]] && echo -e "Partial assembly is done!\n"
[[ $? -eq 0 ]] && echo -e "Step6 is done!\n"
unset flag i


#=====================================================================
## step7: get final longAlt reads
#=====================================================================
echo -e "Step7: Get final longAlt reads!"
# 7.1 get each gapA_${i}_pa.fa
mkdir longAlt_qry
flag=$(wc -l ${refbase2}_racon${npr}_gapA.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	awk '/^S/{print ">"$2;print $3}' partial_assembly/gapA_${i}.p_ctg.gfa > longAlt_qry/gapA_${i}.p_ctg.fa
	awk '/^S/{print ">"$2;print $3}' partial_assembly/gapA_${i}.a_ctg.gfa > longAlt_qry/gapA_${i}.a_ctg.fa
	cat longAlt_qry/gapA_${i}.p_ctg.fa longAlt_qry/gapA_${i}.a_ctg.fa > longAlt_qry/gapA_${i}_pa.fa
	seqkit seq -n longAlt_qry/gapA_${i}_pa.fa -j $thr > longAlt_qry/gapA_${i}_pa.id
	rm longAlt_qry/gapA_${i}.p_ctg.fa longAlt_qry/gapA_${i}.a_ctg.fa
done
unset flag i

# 7.2 filter ontAlt_qry.fa by reads depth
# module load samtools/1.11
flag=$(wc -l ${refbase2}_racon${npr}_gapA.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	minimap2 -ax map-hifi longAlt_qry/gapA_${i}_pa.fa ccsAlt_qry/gapA_${i}.fq -t $thr | samtools view -bS -@ $thr - | samtools sort -o longAlt_qry/gapA_${i}_s.bam -@ $thr -
	samtools depth -a longAlt_qry/gapA_${i}_s.bam > longAlt_qry/gapA_${i}_dep.txt
	awk '$3==0' longAlt_qry/gapA_${i}_dep.txt | awk '{print $1}' | sort | uniq > longAlt_qry/gapA_${i}_dep0.id
	awk 'ARGIND==1{a[$1]} ARGIND==2{if(!($1 in a)){print $1}}' longAlt_qry/gapA_${i}_dep0.id longAlt_qry/gapA_${i}_pa.id > longAlt_qry/gapA_${i}_pa_f1.id
	seqkit grep -f longAlt_qry/gapA_${i}_pa_f1.id longAlt_qry/gapA_${i}_pa.fa -j $thr >> longAlt_qry/gapA_pa_f1.fa
done
unset flag i

# 7.3 filter ontAlt_qry.fa by reads length
seqkit rename longAlt_qry/gapA_pa_f1.fa -j $thr | seqkit seq -m $dfl -j $thr > longAlt_qry/gapA_pa_f2.fa

# 7.4 filter ontAlt_qry.fa by reads coverage
minimap2 -x asm5 $qrybase1 longAlt_qry/gapA_pa_f2.fa -t $thr > longAlt_qry/${qrybase2}TogapA_pa_f2.paf
awk '$10/$2>="'$PER'" && $10/$2<1{print $1}' longAlt_qry/${qrybase2}TogapA_pa_f2.paf | sort | uniq > longAlt_qry/gapA_pa_f.id
seqkit grep -f longAlt_qry/gapA_pa_f.id longAlt_qry/gapA_pa_f2.fa -j $thr > longAlt_qry/gapA_pa_f3.fa

ln -s longAlt_qry/gapA_pa_f3.fa longAlt_qry.fa
#[[ $? -eq 0 ]] && echo -e "Get final ontAlt reads is done!\n"
[[ $? -eq 0 ]] && echo -e "Step7 is done!\n"


#=====================================================================
## step8: Visualization of result statistics
#=====================================================================
echo -e "Step8: Visualization of result statistics!"
pl_lenDis.pl -f longAlt_qry.fa -split-length $dfl -o longAlt_qry_lenDis.svg
#[[ $? -eq 0 ]] && echo -e "The Visualization is done!\n"
[[ $? -eq 0 ]] && echo -e "Step8 is done!\n"


#=====================================================================
## step9: rm and mv result to ouput
#=====================================================================
if [ "$qm" = "yes" ] && [ "$cm" = "yes" ]; then
	rm ${refbase2}_racon0.fa $qrybase1 $ccsbase1
fi
unset qm cm

mv *racon* ragtag* ${qrybase2}* part* longAlt_* ccsAlt_* gap* $out
[[ $? -eq 0 ]] && echo -e "Congratulations! The l3-8.sh is done!\n"


#!/bin/bash

######################################################################
# File Name: l3-6.sh (Thu Nov 14 20:00:00 2023);
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
out="output"
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
		-h)
		hep="help"
		shift
		;;
		-version)
		ver="version"
		shift
		;;
	esac
done


######################################################################
## creates error message and exits if there is an error;
######################################################################
if [ "$hep" == "help" ]; then
	echo "Usage: l3-6.sh [-r reference genome] [-q query genome] [-c ccs reads] [options]
Options:
	Input/Output:
	-r			reference genome
	-q			query genome
	-c			ccs reads
	-o			output directory

	Polish:
	-n INT		Number of Polishing Rounds [1-5], default 3

	Filter:
	-i FLOAT	Set the minimum alignment identity [0, 100], default 90
	-l INT		Set the minimum alignment length, default 20000
	-p FLOAT	Extract the PacBio HiFi read which align length is >= *% of its own length [0-1], default 0.9
	-P FLOAT	Extract the query ontAlt read which align length is >= *% of its own length [0-1), default 0.5

	Supp:
	-t INT		number of threads, default 1
	-version	show version number

See more information at https://github.com/wzxie/ONTbyAHR.
	"
	exit

elif [ "$ver" == "version" ]; then
	echo "Version: 1.3.6"
	exit

else
	[[ $ref == "" ]] && echo -e "ERROR: Path to reference genome not found, assign using -r." && exit
	[[ $qry == "" ]] && echo -e "ERROR: Path to query genome not found, assign using -q." && exit
	[[ $ccs == "" ]] && echo -e "ERROR: PacBio HiFi reads of query not found, assign using -c." && exit
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
		echo -e "\tError: Command $scr is NOT in you PATH. Please check."
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
[[ ! -e $qrybase1 ]] && ln -s $qry $qrybase1
[[ ! -e $ccsbase1 ]] && ln -s $ccs $ccsbase1


#=====================================================================
## step1: filter query`s hifi reads, which just used to step2
#=====================================================================
echo -e "Step1: Filter the $ccsbase1, which just used to Step2!"
minimap2 -x map-hifi ${qrybase1} $ccsbase1 -t $thr > ${qrybase2}_${ccsbase2}.paf
less ${qrybase2}_${ccsbase2}.paf | awk '$10/$2>=0.90{print $1}' | sort | uniq > ${ccsbase2}_f.id
seqkit grep -f ${ccsbase2}_f.id $ccsbase1 > ${ccsbase2}_f.fq
[[ $? -eq 0 ]] && echo -e "Filter is done!\n"


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
[[ $? -eq 0 ]] && echo -e "Scaffold and locate expanded gap of $qrybase1 is done!\n"

# 3.3 get not ragtag contig as alternative ONT reads
#seqkit fx2tab -l -n ${qrybase2}_ragtag.fa -j $thr | awk '$1!~/_RagTag/ && $2>=20000{print $1}' > ${qrybase2}_ragtagN.id
#seqkit grep -f ${qrybase2}_ragtagN.id ${qrybase2}_ragtag.fa -j $thr > ontAlt_qryOwn.fa
#[[ $? -eq 0 ]] && echo -e "Get ontAlt_qryOwn.fa is done!\n"


#=====================================================================
## step4: locating the gap at ref, and get ontAlt_ref.fa
#=====================================================================
echo -e "Step4: locate expanded gap of ${refbase2}_racon${npr}.fa, then get ontAlt_ref.fa!"
# 4.1 get algin block
nucmer -p ${refbase2}_racon${npr}To${qrybase2}_ragtag ${refbase2}_racon${npr}.fa ${qrybase2}_ragtag.fa -t $thr
delta-filter -i $dfi -l $dfl ${refbase2}_racon${npr}To${qrybase2}_ragtag.delta > ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.delta
show-coords -b -B ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.delta | awk '{if($9<$10){print $1"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12} else{print $1"\t"$8"\t"$10"\t"$9"\t"$11"\t"$12}}' > ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.coords

# 4.2 get gap.bed of ref; and get alternative ONT reads
awk 'ARGIND==1{a[$1];b[$1]=$2;c[$1]=$3} ARGIND==2{if($1 in a && $2"_RagTag" in a && $4>=b[$1] && $3<=c[$1]){print $2"\t"$5"\t"$6}}' gap_around.bed ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.coords > ${refbase2}_racon${npr}_gapA.bed

mkdir ontAlt_ref
declare -i flag=0
cat ${refbase2}_racon${npr}_gapA.bed | while read id start end
do
	flag=$(($flag + 1))
	echo -e "$id\t$start\t$end" > ontAlt_ref/gapA_${flag}.bed
	seqkit subseq --bed ontAlt_ref/gapA_${flag}.bed ${refbase2}_racon${npr}.fa -j $thr | seqkit replace -p ":." -j $thr > ontAlt_ref/gapA_${flag}.fa
	rm ontAlt_ref/gapA_${flag}.bed
done
unset flag

# 4.3 visualization the location of gapA
pl_locate.pl -f ${refbase2}_racon${npr}_gapA.bed -l ${refbase2}_racon${npr}_len.txt -o ${refbase2}_racon${npr}_gapA.svg
[[ $? -eq 0 ]] && echo -e "Get ontAlt_ref.fa is done!\n"


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
[[ $? -eq 0 ]] && echo -e "Get ccsAlt_qry.fq is done!\n"
unset flag


#=====================================================================
## step6: partial assembly
#=====================================================================
echo -e "Step6: partial assembly with ccsAlt_qry.fq and ontAlt_ref.fa!"
mkdir partial_assembly
flag=$(wc -l ${refbase2}_racon${npr}_gapA.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	hifiasm --primary -o partial_assembly/gapA_${i} --ul ontAlt_ref/gapA_${i}.fa ccsAlt_qry/gapA_${i}.fq -t $thr
done
[[ $? -eq 0 ]] && echo -e "Partial assembly is done!\n"
unset flag i


#=====================================================================
## step7: get final ontAlt reads
#=====================================================================
echo -e "Step7: Get final ontAlt reads!"
# 7.1 get each gapA_${i}_pa.fa
mkdir ontAlt_qry
flag=$(wc -l ${refbase2}_racon${npr}_gapA.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	awk '/^S/{print ">"$2;print $3}' partial_assembly/gapA_${i}.p_ctg.gfa > ontAlt_qry/gapA_${i}.p_ctg.fa
	awk '/^S/{print ">"$2;print $3}' partial_assembly/gapA_${i}.a_ctg.gfa > ontAlt_qry/gapA_${i}.a_ctg.fa
	cat ontAlt_qry/gapA_${i}.p_ctg.fa ontAlt_qry/gapA_${i}.a_ctg.fa > ontAlt_qry/gapA_${i}_pa.fa
	seqkit seq -n ontAlt_qry/gapA_${i}_pa.fa -j $thr > ontAlt_qry/gapA_${i}_pa.id
	rm ontAlt_qry/gapA_${i}.p_ctg.fa ontAlt_qry/gapA_${i}.a_ctg.fa
done
unset flag i

# 7.2 filter ontAlt_qry.fa by reads depth
# module load samtools/1.11
flag=$(wc -l ${refbase2}_racon${npr}_gapA.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	minimap2 -ax map-hifi ontAlt_qry/gapA_${i}_pa.fa ccsAlt_qry/gapA_${i}.fq -t $thr | samtools view -bS -@ $thr - | samtools sort -o ontAlt_qry/gapA_${i}_s.bam -@ $thr -
	samtools depth -a ontAlt_qry/gapA_${i}_s.bam > ontAlt_qry/gapA_${i}_dep.txt
	awk '$3==0' ontAlt_qry/gapA_${i}_dep.txt | awk '{print $1}' | sort | uniq > ontAlt_qry/gapA_${i}_dep0.id
	awk 'ARGIND==1{a[$1]} ARGIND==2{if(!($1 in a)){print $1}}' ontAlt_qry/gapA_${i}_dep0.id ontAlt_qry/gapA_${i}_pa.id > ontAlt_qry/gapA_${i}_pa_f1.id
	seqkit grep -f ontAlt_qry/gapA_${i}_pa_f1.id ontAlt_qry/gapA_${i}_pa.fa -j $thr >> ontAlt_qry/gapA_pa_f1.fa
done
unset flag i

# 7.3 filter ontAlt_qry.fa by reads length
seqkit rename ontAlt_qry/gapA_pa_f1.fa -j $thr | seqkit seq -m 20000 -j $thr > ontAlt_qry/gapA_pa_f2.fa

# 7.4 filter ontAlt_qry.fa by
minimap2 -x asm5 $qrybase1 ontAlt_qry/gapA_pa_f2.fa -t $thr > ontAlt_qry/${qrybase2}TogapA_pa_f2.paf
awk '$10/$2>="'$PER'" && $10/$2<1{print $1}' ontAlt_qry/${qrybase2}TogapA_pa_f2.paf | sort | uniq > ontAlt_qry/gapA_pa_f.id
seqkit grep -f ontAlt_qry/gapA_pa_f.id ontAlt_qry/gapA_pa_f2.fa -j $thr > ontAlt_qry/gapA_pa_f3.fa

ln -s ontAlt_qry/gapA_pa_f3.fa ontAlt_qry.fa
[[ $? -eq 0 ]] && echo -e "Get final ontAlt reads is done!\n"


#=====================================================================
## step8: Visualization of result statistics
#=====================================================================
echo -e "Step8: Visualization of result statistics!"
pl_lenDis.pl -f ontAlt_qry.fa -split-length 20000 -o ontAlt_qry_lenDis.svg
[[ $? -eq 0 ]] && echo -e "The Visualization is done!\n"


#=====================================================================
## step9: rm and mv result to ouput
#=====================================================================
rm ${refbase2}_racon0.fa $qrybase1 $ccsbase1
mv ref_racon* ragtag* query* part* ontAlt_* ccsAlt_* gap* $out


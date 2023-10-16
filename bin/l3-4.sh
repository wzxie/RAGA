#!/bin/bash

######################################################################
# File Name: hello.sh (Thu Sep 07 18:39:00 2023);
######################################################################
set -o nounset
#bak=$IFS
#IFS=$'\n'


######################################################################
## Check other scripts if they are found in path;
######################################################################
for scr in minimap2 racon ragtag.py nucmer delta-filter show-coords awk hifiasm samtools
do
	check=$(command -v $scr)
	if [ "$check" == "" ]; then
		echo -e "Error: Command $scr is NOT in you PATH. Please check."
		exit 1
	fi
done
echo -e "All dependencies have been checked.\n"


######################################################################
## parameter;
######################################################################
ref=""
qry=""
ccs=""
out="output"
thr="1"
npr="3"
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
	echo "Usage: l3-3.sh [-r reference genome] [-q query genome] [-c ccs reads] [-o output directory] [-t number of threads] [-n number of polishing rounds]
Options:
	Input/Output:
	-r	reference genome
	-q	query genome
	-c	ccs reads
	-o	output directory
	Polish:
	-n INT		Number of Polishing Rounds [1-5]
	Mapping:
#	-f FLOAT	filter out top FLOAT fraction of repetitive minimizers [0.0002]
#	-g NUM		stop chain enlongation if there are no minimizers in INT-bp [5000]
#	-G NUM		max intron length (effective with -xsplice; changing -r) [200k]
#	-F NUM		max fragment length (effective with -xsr or in the fragment mode) [800]
#	-r NUM		chaining/alignment bandwidth and long-join bandwidth [500,20000]
#	-n INT		minimal number of minimizers on a chain [3]
#	-m INT		minimal chaining score (matching bases minus log gap penalty) [40]
#	-X			skip self and dual mappings (for the all-vs-all mode)
#	-p FLOAT	min secondary-to-primary score ratio [0.8]
#	-N INT		retain at most INT secondary alignments [5]
	Alignment:
#	-A INT		matching score [2]
#	-B INT		mismatch penalty (larger value for lower divergence) [4]
#	-O INT		gap open penalty [4,24]
#	-E INT		gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [2,1]
#	-z INT		Z-drop score and inversion Z-drop score [400,200]
#	-s INT		minimal peak DP alignment score [80]
#	-u CHAR		how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]
	Supp:
	-t INT		number of threads [1]
	-version	show version number

See more information at https://github.com/***.\n
	"
	exit

elif [ "$ver" == "version" ]; then
	echo "Version: 1.3.3"
	exit

else
	[[ $ref == "" ]] && echo -e "ERROR: Path to reference genome not found, assign using -r." && exit
	[[ $qry == "" ]] && echo -e "ERROR: Path to query genome not found, assign using -q." && exit
	[[ $ccs == "" ]] && echo -e "ERROR: PacBio HiFi reads of query not found, assign using -c." && exit
	[[ -d $out ]] && echo -e "ERROR: Output directory already exists, please remove or set alternate output." && exit
	[[ -d $out ]] || mkdir $out
fi


#=====================================================================
## step0: file name parser
#=====================================================================
refbase1=$(basename $ref)
qrybase1=$(basename $qry)
ccsbase1=$(basename $ccs)
refbase2=${refbase1%.*}
qrybase2=${qrybase1%.*}
[[ ! -e ${refbase2}_racon0.fa ]] && ln -s $ref ${refbase2}_racon0.fa
[[ ! -e $qrybase1 ]] && ln -s $qry $qrybase1
[[ ! -e $ccsbase1 ]] && ln -s $ccs $ccsbase1


#=====================================================================
## step1: polish
#=====================================================================
echo -e "Step1: polish $refbase1 $npr rounds by racon!"
for i in $(seq 1 $npr)
do
	# racon
	flag=$(($i-1))
	minimap2 -ax map-hifi ${refbase2}_racon${flag}.fa $ccsbase1 -t $thr > ${refbase2}_racon${flag}.sam
	racon $ccsbase1 ${refbase2}_racon${flag}.sam ${refbase2}_racon${flag}.fa -t $thr > ${refbase2}_racon${i}.fa

	# $?
	if [[ $? -eq 0 ]]; then
		[[ $i == 1 ]] && echo -e "The round 1 of polish is done.\n"
		[[ $i == 2 ]] && echo -e "The round 2 of polish is done.\n"
		[[ $i == 3 ]] && echo -e "The round 3 of polish is done.\n"
		[[ $i == 4 ]] && echo -e "The round 4 of polish is done.\n"
		[[ $i == 5 ]] && echo -e "The round 5 of polish is done.\n"
	fi
done
seqkit fx2tab -l -n ${refbase2}_racon${npr}.fa | awk '{print $1"\t0\t"$5}' > ${refbase2}_racon${npr}_len.txt
unset i flag


#=====================================================================
## step2: scaffolding, locate expanded gap area
#=====================================================================
echo -e "Step2: Homology-based scaffolding ${refbase2}_racon${npr}.fa by RagTag!"
# 2.1 scaffolding
ragtag.py scaffold ${refbase2}_racon${npr}.fa $qrybase1 -t $thr -o .
[[ $? -eq 0 ]] && echo -e "Scaffolding is done!\n"
ln -s ragtag.scaffold.fasta ${qrybase2}_ragtag.fa

# 2.2 get gap between ref_racon${npr}.fa and qry_ragtag.fa
awk '$5=="U"{if($2-1000000>=0){print $1"\t"$2-1000000"\t"$3+1000000} else {print $1"\t0\t"$3+1000000}}' ragtag.scaffold.agp > gap_around.bed
[[ $? -eq 0 ]] && echo -e "Locate expanded gap of $qrybase1 is done!\n"

# 2.3 get not ragtag contig as alternative ONT reads
#seqkit fx2tab -l -n ${qrybase2}_ragtag.fa -j $thr | awk '$1!~/_RagTag/ && $2>=20000{print $1}' > ${qrybase2}_ragtagN.id
#seqkit grep -f ${qrybase2}_ragtagN.id ${qrybase2}_ragtag.fa -j $thr > ontAlt_qryOwn.fa
#[[ $? -eq 0 ]] && echo -e "Get ontAlt_qryOwn.fa is done!\n"


#=====================================================================
## step3: locating the gap at ref, and get ontAlt_ref.fa
#=====================================================================
echo -e "Step3: locate expanded gap of ${refbase2}_racon${npr}.fa, then get ontAlt_ref.fa!"
# 3.1 get algin block
nucmer -p ${refbase2}_racon${npr}To${qrybase2}_ragtag ${refbase2}_racon${npr}.fa ${qrybase2}_ragtag.fa -t $thr
delta-filter -i 90 -l 20000 ${refbase2}_racon${npr}To${qrybase2}_ragtag.delta > ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.delta
show-coords -b -B ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.delta | awk '{if($9<$10){print $1"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12} else{print $1"\t"$8"\t"$10"\t"$9"\t"$11"\t"$12}}' > ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.coords

# 3.2 get gap.bed of ref; and get alternative ONT reads
awk 'ARGIND==1{a[$1];b[$1]=$2;c[$1]=$3} ARGIND==2{if($1 in a && $2"_RagTag" in a && $4>=b[$1] && $3<=c[$1]){print $2"\t"$5"\t"$6}}' gap_around.bed ${refbase2}_racon${npr}To${qrybase2}_ragtag_f.coords > ${refbase2}_racon${npr}_gapA.bed
seqkit subseq --bed ${refbase2}_racon${npr}_gapA.bed ${refbase2}_racon${npr}.fa -j $thr | seqkit replace -p ":." -j $thr > ontAlt_ref.fa

# 3.3 visualization the location of gapA
pl_locate.pl -f ${refbase2}_racon${npr}_gapA.bed -l ${refbase2}_racon${npr}_len.txt -o gapA_area.svg
[[ $? -eq 0 ]] && echo -e "Get ontAlt_ref.fa is done!\n"


#=====================================================================
## step4: minimap2 qryHifi reads to ref, and get ccsAlt_qry.fq
#=====================================================================
echo -e "Step4: minimap2 $ccsbase1 to ${refbase2}_racon${npr}.fa, then get ccsAlt_qry.fq!"
minimap2 -x map-hifi ${refbase2}_racon${npr}.fa $ccsbase1 -t $thr > ${refbase2}_racon${npr}Toqry_hifi.paf
awk 'ARGIND==1{a[$1];b[$1]=$2;c[$1]=$3} ARGIND==2{if($6 in a && $9>=b[$6] && $8<=c[$6] && $10/$2>=0.9){print $1}}' ${refbase2}_racon${npr}_gapA.bed ${refbase2}_racon${npr}Toqry_hifi.paf | sort | uniq > ${refbase2}_racon${npr}Toqry_hifi_f.id
seqkit grep -f ${refbase2}_racon${npr}Toqry_hifi_f.id $ccsbase1 -j $thr > ccsAlt_qry.fq
[[ $? -eq 0 ]] && echo -e "Get ccsAlt_qry.fq is done!\n"


#=====================================================================
## step5: partial assembly
#=====================================================================
echo -e "Step5: partial assembly with ccsAlt_qry.fq and ontAlt_ref.fa!"
/home/rpzhao/soft/hifiasm-0.19.4/hifiasm --primary -o part --ul ontAlt_ref.fa ccsAlt_qry.fq -t $thr
[[ $? -eq 0 ]] && echo -e "Partial assembly is done!\n"


#=====================================================================
## step6: get final ontAlt reads
#=====================================================================
echo -e "Step6: Get final ontAlt reads!"
# 6.1 get raw ontAlt_qry.fa
awk '/^S/{print ">"$2;print $3}' part.p_ctg.gfa > part.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' part.a_ctg.gfa > part.a_ctg.fa
cat part.p_ctg.fa part.a_ctg.fa > part_pa.fa
seqkit seq -n part_pa.fa -j $thr > part_pa.id

# 6.2 filter ontAlt_qry.fa by reads depth
# module load samtools/1.11
minimap2 -ax map-hifi part_pa.fa ccsAlt_qry.fq -t $thr | samtools view -bS -@ $thr - | samtools sort -o part_ccsAlt_s.bam -@ $thr -
samtools depth -a part_ccsAlt_s.bam > part_ccsAlt_dep.txt
awk '$3==0' part_ccsAlt_dep.txt | awk '{print $1}' | sort | uniq > part_ccsAlt_dep0.id
awk 'ARGIND==1{a[$1]} ARGIND==2{if(!($1 in a)){print $1}}' part_ccsAlt_dep0.id part_pa.id > part_pa_f.id
seqkit grep -f part_pa_f.id part_pa.fa -j $thr | seqkit rename -j $thr | seqkit seq -m 20000 -j $thr > part_pa_f.fa
ln -s part_pa_f.fa ontAlt_qry.fa

# 6.3 get ontAlt_qry.fa
#cat ontAlt_qryOwn.fa ontAlt_qryReb.fa > ontAlt_qry.fa

# 6.4 filter ontAlt_qry.fa
minimap2 -x asm5 $qrybase1 ontAlt_qry.fa -t $thr > ${qrybase2}ToontAlt_qry.paf
awk '$10/$2>=0.5 && $10/$2<1{print $1}' ${qrybase2}ToontAlt_qry.paf | sort | uniq > ontAlt_qry_f.id
seqkit grep -f ontAlt_qry_f.id ontAlt_qry.fa -j $thr > ontAlt_qry_f.fa
[[ $? -eq 0 ]] && echo -e "Get final ontAlt reads is done!\n"


#=====================================================================
## step7: Visualization of result statistics
#=====================================================================
echo -e "Step7: Visualization of result statistics!"
pl_lenDis.pl -f ontAlt_qry_f.fa -split-length 20000 -o ontAlt_qry_lenDis.svg
[[ $? -eq 0 ]] && echo -e "The Visualization is done!\n"


#=====================================================================
## step8: rm and mv result to ouput
#=====================================================================
#rm ${refbase2}_racon0.fa $qrybase1 $ccsbase1
#mv ref* ragtag* lal* query* gap* ontAlt* ccs* $out


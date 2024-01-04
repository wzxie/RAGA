#!/bin/bash

######################################################################
# File Name: RAGA-diff.sh (Thu Jan 01 09:00:00 2024);
######################################################################
set -o nounset
#bak=$IFS
#IFS=$'\n'


######################################################################
## Parameter;
######################################################################
ref=""
ccs=""
out="output2"
thr="1"
npr="10"
dfl="10000"
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
	echo "Usage: RAGA-diff.sh [-r reference genome] [-c ccs reads] [options]
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
	"
	exit

elif [ "$ver" == "version" ]; then
	echo "Version: 1.0.0"
	exit

else
	[[ $ref == "" ]] && echo -e "ERROR: Path to reference genome not found, assign using -r." && exit
	[[ $ccs == "" ]] && echo -e "ERROR: PacBio HiFi reads of source assembly not found, assign using -c." && exit
	[[ $npr -lt 10 ]] && echo -e "ERROR: Number of Polishing Rounds [>=10], default 10." && exit
	[[ -d $out ]] && echo -e "ERROR: Output directory already exists, please remove or set alternate output." && exit
	[[ -d $out ]] || mkdir $out
fi


######################################################################
## Check other scripts if they are found in path;
######################################################################
echo -e "Verifying the availability of related dependencies!"
for scr in minimap2 racon awk bedtools seqkit hifiasm
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


#=====================================================================
## step1: polish
#=====================================================================
echo -e "Step1: polish the $refbase1 $npr rounds using ${ccsbase1}!"
for i in $(seq 1 $npr)
do
    # racon
    flag=$(($i-1))
    minimap2 -ax map-hifi ${refbase2}_racon${flag}.fa $ccsbase1 -t $thr > ${refbase2}_racon${flag}.sam
    racon $ccsbase1 ${refbase2}_racon${flag}.sam ${refbase2}_racon${flag}.fa -t $thr > ${refbase2}_racon${i}.fa

    # $?
    if [[ $? -eq 0 ]]; then
		echo -e "The round $i of polish is done.\n"
    fi
done
[[ $? -eq 0 ]] && echo -e "Step1 is done!\n"
unset flag i


#=====================================================================
## step2: minimap2 surHifi reads to ref, and get longAlt_ref.fa
#=====================================================================
echo -e "Step2: minimap2 $ccsbase1 to ${refbase2}_racon${npr}.fa, then get longAlt_ref.fa!"
minimap2 -x map-hifi ${refbase2}_racon${npr}.fa $ccsbase1 -t $thr > ${refbase2}_racon${npr}To${ccsbase2}.paf
awk '$13~/tp:A:P/ && $10/$2>=0.99{print $6"\t"$8"\t"$9}' ${refbase2}_racon${npr}To${ccsbase2}.paf | sort -k1,1 -k2,2n > ${refbase2}_racon${npr}To${ccsbase2}.bed
bedtools merge -i ${refbase2}_racon${npr}To${ccsbase2}.bed | awk '$3-$2>="'$dfl'"{print $1"\t"$2"\t"$3}' > ${refbase2}_racon${npr}To${ccsbase2}_m.bed

mkdir longAlt_ref
declare -i flag=0
cat ${refbase2}_racon${npr}To${ccsbase2}_m.bed | while read id start end
do
    flag=$(($flag + 1))
    echo -e "$id\t$start\t$end" > longAlt_ref/block_${flag}.bed
    seqkit subseq --bed longAlt_ref/block_${flag}.bed ${refbase2}_racon${npr}.fa -j $thr | seqkit replace -p ":." -j $thr > longAlt_ref/block_${flag}.fa
    rm longAlt_ref/block_${flag}.bed
done
#[[ $? -eq 0 ]] && echo -e "Get longAlt_ref.fa is done!\n"
[[ $? -eq 0 ]] && echo -e "Step2 is done!\n"
unset flag


#=====================================================================
## step3: minimap2 surHifi reads to ref, and get ccsAlt_sur.fq
#=====================================================================
echo -e "Step3: minimap2 $ccsbase1 to ${refbase2}_racon${npr}.fa, then get ccsAlt_sur.fq!"
mkdir ccsAlt_sur
declare -i flag=0
cat ${refbase2}_racon${npr}To${ccsbase2}_m.bed | while read id start end
do
	flag=$(($flag + 1))
	echo -e "$id\t$start\t$end" > ccsAlt_sur/block_${flag}.bed
	awk 'ARGIND==1{a[$1];b[$1]=$2;c[$1]=$3} ARGIND==2{if($6 in a && $9>=b[$6] && $8<=c[$6] && $10/$2>=0.9 && $13~/tp:A:P/){print $1}}' ccsAlt_sur/block_${flag}.bed ${refbase2}_racon${npr}To${ccsbase2}.paf | sort | uniq > ccsAlt_sur/block_${flag}.id
	seqkit grep -f ccsAlt_sur/block_${flag}.id $ccsbase1 -j $thr > ccsAlt_sur/block_${flag}.fq
	rm ccsAlt_sur/block_${flag}.bed ccsAlt_sur/block_${flag}.id
done
#[[ $? -eq 0 ]] && echo -e "Get ccsAlt_sur.fq is done!\n"
[[ $? -eq 0 ]] && echo -e "Step3 is done!\n"
unset flag


#=====================================================================
## step4: partial assembly
#=====================================================================
echo -e "Step4: partial assembly with ccsAlt_sur.fq and longAlt_ref.fa!"
mkdir partial_assembly
flag=$(wc -l ${refbase2}_racon${npr}To${ccsbase2}_m.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
    hifiasm --primary -o partial_assembly/block_${i} --ul longAlt_ref/block_${i}.fa ccsAlt_sur/block_${i}.fq -t $thr
done
#[[ $? -eq 0 ]] && echo -e "Partial assembly is done!\n"
[[ $? -eq 0 ]] && echo -e "Step4 is done!\n"
unset flag i 


#=====================================================================
## step5: get final longAlt reads
#=====================================================================
echo -e "Step5: Get final longAlt reads!"
# 5.1 get each block_${i}_pa.fa
mkdir longAlt_sur
flag=$(wc -l ${refbase2}_racon${npr}To${ccsbase2}_m.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	awk '/^S/{print ">"$2;print $3}' partial_assembly/block_${i}.p_ctg.gfa > longAlt_sur/block_${i}.p_ctg.fa
	awk '/^S/{print ">"$2;print $3}' partial_assembly/block_${i}.a_ctg.gfa > longAlt_sur/block_${i}.a_ctg.fa
	cat longAlt_sur/block_${i}.p_ctg.fa longAlt_sur/block_${i}.a_ctg.fa > longAlt_sur/block_${i}_pa.fa
	seqkit seq -n longAlt_sur/block_${i}_pa.fa -j $thr > longAlt_sur/block_${i}_pa.id
	rm longAlt_sur/block_${i}.p_ctg.fa longAlt_sur/block_${i}.a_ctg.fa
done
unset flag i

# 5.2 filter longAlt_sur.fa by reads depth
# module load samtools/1.11
flag=$(wc -l ${refbase2}_racon${npr}To${ccsbase2}_m.bed | awk '{print $1}')
for i in $(seq 1 $flag)
do
	minimap2 -ax map-hifi longAlt_sur/block_${i}_pa.fa ccsAlt_sur/block_${i}.fq -t $thr | samtools view -bS -@ $thr - | samtools sort -o longAlt_sur/block_${i}_s.bam -@ $thr -
	samtools depth -a longAlt_sur/block_${i}_s.bam > longAlt_sur/block_${i}_dep.txt
	awk '$3==0' longAlt_sur/block_${i}_dep.txt | awk '{print $1}' | sort | uniq > longAlt_sur/block_${i}_dep0.id
	awk 'ARGIND==1{a[$1]} ARGIND==2{if(!($1 in a)){print $1}}' longAlt_sur/block_${i}_dep0.id longAlt_sur/block_${i}_pa.id > longAlt_sur/block_${i}_pa_f1.id
	seqkit grep -f longAlt_sur/block_${i}_pa_f1.id longAlt_sur/block_${i}_pa.fa -j $thr >> longAlt_sur/block_pa_f1.fa
done
unset flag i

# 5.3 filter longAlt_sur.fa by reads length
seqkit seq -m $dfl -j $thr longAlt_sur/block_pa_f1.fa | seqkit fx2tab -l -j $thr | awk '{print ">longAlt_"NR; print $2}' > longAlt_sur/block_pa_f2.fa

ln -s longAlt_sur/block_pa_f2.fa longAlt_sur.fa
#[[ $? -eq 0 ]] && echo -e "Get final longAlt reads is done!\n"
[[ $? -eq 0 ]] && echo -e "Step5 is done!\n"


#=====================================================================
## step6: rm and mv result to ouput
#=====================================================================
if [ "$cm" = "yes" ]; then
	rm ${refbase2}_racon0.fa $ccsbase1
fi
unset cm

mv *racon* ccsAlt_* part* *longAlt_* $out
[[ $? -eq 0 ]] && echo -e "Congratulations! The RAGA-diff.sh is done!\n"


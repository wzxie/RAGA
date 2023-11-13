#!/bin/bash

######################################################################
# File Name: l4-1.sh (Thu Nov 10 09:30:00 2023);
######################################################################
set -o nounset
#bak=$IFS
#IFS=$'\n'


######################################################################
## Parameter;
######################################################################
ref=""
ccs=""
out="out2"
thr="1"
npr="10"
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
	echo "Usage: l4-1.sh [-r reference genome] [-c ccs reads] [options]
Options:
	Input/Output:
	-r		  reference genome
	-c		  ccs reads
	-o		  output directory

	Polish:
	-n INT	  Number of Polishing Rounds [>=10], default 10

	Supp:
	-t INT	  number of threads, default 1
	-version	show version number

See more information at https://github.com/wzxie/ONTbyAHR.
	"
	exit

elif [ "$ver" == "version" ]; then
	echo "Version: 1.4.1"
	exit

else
	[[ $ref == "" ]] && echo -e "ERROR: Path to reference genome not found, assign using -r." && exit
	[[ $ccs == "" ]] && echo -e "ERROR: PacBio HiFi reads of query not found, assign using -c." && exit
	[[ -d $out ]] && echo -e "ERROR: Output directory already exists, please remove or set alternate output." && exit
	[[ -d $out ]] || mkdir $out
fi


######################################################################
## Check other scripts if they are found in path;
######################################################################
echo -e "Verifying the availability of related dependencies!"
for scr in minimap2 racon awk bedtools seqkit
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
ccsbase1=$(basename $ccs)
refbase2=${refbase1%.*}
ccsbase2=${ccsbase1%.*}
[[ ! -e ${refbase2}_racon0.fa ]] && ln -s $ref ${refbase2}_racon0.fa
[[ ! -e $ccsbase1 ]] && ln -s $ccs $ccsbase1


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
unset flag i


#=====================================================================
## step2: get longAlt_ref.fa
#=====================================================================
###do minimap2
minimap2 -x map-hifi ${refbase2}_racon${npr}.fa $ccsbase1 -t $thr > ${refbase2}_racon${npr}To${ccsbase2}.paf
awk '$13~/tp:A:P/ && $10/$2>=0.99{print $6"\t"$8"\t"$9}' ${refbase2}_racon${npr}To${ccsbase2}.paf | sort -k1,1 -k2,2n > ${refbase2}_racon${npr}To${ccsbase2}.bed
bedtools merge -i ${refbase2}_racon${npr}To${ccsbase2}.bed > ${refbase2}_racon${npr}To${ccsbase2}_m.bed
seqkit subseq --bed ${refbase2}_racon${npr}To${ccsbase2}_m.bed ${refbase2}_racon${npr}.fa | awk 'BEGIN{Cnt=0}{if($0~/>/){Cnt=Cnt+1;tmp=">ref"Cnt"l"; $0=tmp; print$0}else{print$0}}' | seqkit seq -w 0 -m 10000 > ${refbase2}_racon${npr}To${ccsbase2}.fa
ln -s ${refbase2}_racon${npr}To${ccsbase2}.fa longAlt_${refbase2}.fa
[[ $? -eq 0 ]] && echo -e "Get final longAlt reads is done!\n"


#=====================================================================
## step3: rm and mv result to ouput
#=====================================================================
rm ${refbase2}_racon0.fa $ccsbase1
mv *racon* *longAlt* $out
[[ $? -eq 0 ]] && echo -e "The l4-1.sh is done!\n"


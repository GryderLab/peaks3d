#!/bin/bash -l

#AUTHOR: Hyunmin Kim (human.gim@gmail.com)

BASEDIR=$(dirname $(realpath "$0"));
for f in $BASEDIR/../src/*.sh;do
        . $f
done
export PATH=$BASEDIR:$PATH
intro="
LastUpdate: 3/06/22
By	: Hyunmin Kim (hxk728@case.edu)
Summary :
	This pipeline generates the clustered and not-clustered filtered bedpe files
	using ChIP-seq and HiChIP data in the output directory
	

<config.txt> example:

	tad=/home/ubuntu/shared/hg38/ref/TAD_goldsorted_span.hg38.bed
	peaks=(
		peaks/AR_master.bed
		peaks/AR_master.bed
		peaks/AR_master.bed
	)
	samples=(
		LNCaPXIP_AR_dupin
		LNCaPXIP_NT_low_AR
		LNCaPXIP_NT_med_AR
	)

	## function parameters
	binsize=5000  # binsize for annotate_loops.sh
	mindist=25k   # min distance
	maxdist=3m    # max distance to make fitered.bedpe files
	minscore=0.05  # minimum annotated.bedpe score (AQuA CPM)
	mincluster=2

	APAbinsize=1000
	APA_options=\"--cpml\"

	Annotate_options=\" -Q TRUE -T <tf|broad> (default tf) \"

	## output directory
	output=AR_GRACE_hmk_v2
"

usage="
	Usage: $BASH_SOURCE <config.txt>

$intro
"
if [ $# -lt 1 ];then echo "$usage"; exit 1; fi

## read input parameters
source $1
## default function parameters unless defined in the configure file
binsize=${binsize:-5000}
mindist=${mindist:-25k}
maxdist=${maxdist:-3m}  
minscore=${minscore:-0.1}  
mincluster=${mincluster:-2}
APAbinsize=${APAbinsize:-1000}
APA_options="$APA_options  --bin_size $APAbinsize"
Annotate_options=${Annotate_options:-""};
output=${output:-out_peak3d}

if [ ${#peaks[@]} -ne ${#samples[@]} ];then
	echo "Error: the number of peaks and samples are not equal" >&2
	exit 1;
fi



## remind input files and parameters
echo "
#Data:
	tad=$tad
#Options:
	binsize=$binsize;
	mindist=$mindist;
	maxdist=$maxdist;
	minscore=$minscore;
	mincluster=$mincluster;
	APAbinsize=$APAbinsize;
	APA_options=$APA_options;
	Annotate_options=$Annotate_options;
output=$output
"

## read functions
#source $BABY_HOME/inst/bin/src_bedpe.sh
#source $BABY_HOME/inst/bin/src_util.sh
#echo "tools available:"
#declare -F 



bedpe-unmatch(){
perl -e 'use strict; my ($f1, $f2)=("'$1'","'$2'");
        my %r=();
        sub order{ my ($xx,$yy)=@_;
                my @x=split/\t/,$xx;
                my @y=split/\t/,$yy;
                if( $x[1] + $x[2] > $y[2] + $y[2] ){
                        return($yy,$xx);
                }
                return ($xx,$yy);
        }
        open(my $fh,"<",$f2) or die "$!";
        while(<$fh>){chomp; my @d=split/\t/,$_;
                my ($x,$y)=(join("\t",@d[0..2]),join("\t",@d[3..5]));
                $r{$x}{$y} ++;
                $r{$y}{$x} ++;
        }
        close($fh);
        open(my $fh,"<",$f1) or die "$!";
        while(<$fh>){chomp; my @d=split/\t/,$_;
                my ($x,$y)=(join("\t",@d[0..2]),join("\t",@d[3..5]));
                if( !defined $r{$x}{$y} ){
                        print $_,"\n";
                }
        }
'
}
sz(){
	size=`wc -l $1 | cut -d" " -f 1`;
	echo "	=> $1 (n=$size)"
}

#myapa=$BABY_HOME/aqua_tools/plot_APA.sh

if [ -d $output ];then
	echo "Warn: overwriting into $output .. ">&2;
fi
mkdir -p $output


for (( i=0; i < ${#samples[@]}; i++ )); do
##STEP 1
	echo "STEP 1: build bedpe "

	sample=${samples[$i]};
	peak=${peaks[$i]};
	n=${peak##*/};n=$output/${n%.bed}_${sample};
	o=$n.bedpe
	if [ `nz $o` ];then 
		echo "	skip making $o .. " 
	else
		echo "	=> $o .. " 
		build_bedpe.sh -A $peak -B $peak -T $tad > $o
	fi

##STEP 2
	echo "STEP 2: annotate bedpe"
	o=$n.annotated.bedpe
	if [ `nz $o` ];then 
		echo "	skip making $o"
	else
		annotate_loops.sh -G hg38 -P $n.bedpe -A $sample -R $binsize $Annotate_options > $o
		#. $BABY_HOME/inst/bin/src_hic.sh;
		#annotate-loops -G hg38 -P $n.bedpe -A $sample -R $binsize $Annotate_options > $o
		echo " 	=> $o"
	fi
	if [ `wc -l $o | cut -d" " -f 1` -lt 10 ];then
	echo "
	Error: something is wrong!

	please rerun this :

	annotate_loops.sh -G hg38 -P $n.bedpe -A $sample -R $binsize $Annotate_options > $o
	head $o 

	"
	head $o
	echo " --- Stop here: report this issue @ Slack::tinker channel"
	exit 1
	fi

	echo "STEP 2.1 QC plot"
	o=$n.annotated.bedpe.qc
	if [ `nz $o.png` ];then
		echo "	skip making $o.png"
	else
		echo "	=> plot QC"
		qc-annotated-bedpe.sh $n.annotated.bedpe $o &> /dev/null
		qc-annotated-bedpe.sh <( bedpe-filter $n.annotated.bedpe 0 $mindist 0) ${o/.qc/.lt_$mindist.qc}  &>/dev/null
		qc-annotated-bedpe.sh <( bedpe-filter $n.annotated.bedpe $mindist $maxdist 0) ${o/.qc/.gt_$mindist.qc} &>/dev/null
	fi
## STEP 3
	echo "STEP 3: filtering"
	tmp=$minscore;tmp=${tmp/\./p};
	n2=${n}_${mindist}_${maxdist}_${tmp};
	o=$n2.filtered.bedpe
	bedpe-filter $n.annotated.bedpe $mindist $maxdist $minscore > $n2.filtered.bedpe
	sz $n2.filtered.bedpe
## STEP 4
	echo "STEP 4: clustering"
	bedpe-cluster $n2.filtered.bedpe  $mincluster > $n2.clustered.bedpe
	sz $n2.clustered.bedpe

	cat $n2.filtered.bedpe | bedpe-uniq - | perl -e 'use strict; my $mincluster='$mincluster'; my %r=();
	while(<STDIN>){chomp; my@d=split/\t/,$_;
		$r{ join("\t",@d[0..2]) } ++ ;
		$r{ join("\t",@d[3..5] )} ++ ;
	}
	map { print $_,"\t",$r{$_},"\n"; } grep { $r{$_} >= $mincluster } keys %r;
	' > $n2.clustered.bedGraph
	sz $n2.clustered.bedGraph

	groupBy -i $n2.clustered.bedpe -g 8 -c 1,2,6,7,8 -o distinct,min,max,distinct,distinct | cut -f 2- > $n2.clustered.bed 
	sz $n2.clustered.bed

	bedpe-unmatch $n2.filtered.bedpe $n2.clustered.bedpe > $n2.notclustered.bedpe 
	sz $n2.notclustered.bedpe


## STEP 5
	echo "STE 5: APA plotting"
	tmp=${n2#*\/};

	if [ ! -f $output/notclustered/${tmp}/$sample/APA_cpm.pdf ];then
		plot_APA.sh -O $output/notclustered -G hg38 -P $n2.notclustered.bedpe -A $sample $APA_options &> /dev/null;
		echo "	=> $output/notclustered/${tmp}/$sample/APA_cpm.pdf"
	else
		echo "	skip : exists $output/notclustered/${tmp}/$sample/APA_cpm.pdf"
	fi
	if [ ! -f $output/clustered/${tmp}/$sample/APA_cpm.pdf ];then
		plot_APA.sh -O $output/clustered -G hg38 -P $n2.clustered.bedpe -A $sample $APA_options  &> /dev/null;
		echo "	=> $output/clustered/${tmp}/$sample/APA_cpm.pdf"
	else
		echo "	skip : exists $output/clustered/${tmp}/$sample/APA_cpm.pdf"
	fi



## STEP 6
	echo "STEP 6: Generate unlooped features .. "
	## merge 02/14/22 
	bedpe-to-bed.sh $n2.notclustered.bedpe | sort -k1,1 -k2,3n | bedtools merge -i stdin  > $n2.notclustered.bed 

       	## fix peak columns 
	bedtools intersect -a <( cat $peak |cut -f 1-3 ) -b $n2.notclustered.bed -v > $n2.unlooped.bed
	bedtools intersect -a $n2.unlooped.bed -b $n2.clustered.bed -v > $n2.unlooped.nocluster.bed
	summarize_interval.sh -G hg38 -I $n2.unlooped.nocluster.bed -A $sample > $n2.unloop.density.bed
	summarize_interval.sh -G hg38 -I $n2.notclustered.bed -A $sample > $n2.notclustered.density.bed
	summarize_interval.sh -G hg38 -I $n2.clustered.bed -A $sample > $n2.clustered.density.bed
	for f in $n2.unlooped.bed $n2.unlooped.nocluster.bed $n2.unloop.density.bed $n2.notclustered.density.bed $n2.clustered.density.bed;do
		sz $f
	done

	mkdir -p $output/GRACE_plots
	loopsname=${peak##*/}; loopsname=${loopsname%.bed}
	peaks3Dname=${n2#*/};
	cat $BASEDIR/../src/grace_plots.r |\
		 sed  "s#%peaks3Dname%#$peaks3Dname#" |\
		 sed  "s#%loops.name%#$loopsname#" |\
		 sed  "s#%output%#$output#" > $output/GRACE_plots.R 
	R -f $output/GRACE_plots.R &>/dev/null
	for f in  $output/GRACE_plots/*.pdf;do
		echo "	=> $f";
	done


done


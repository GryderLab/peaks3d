#BASEDIR=$(dirname $(realpath "$0"));

km(){
        echo $1 | sed s/[mM]/000k/ | sed s/[kK]/000/
}


_hic2dev(){
local b=`km $2`
        cat<<-eof | Rscript  -
                hic="$1"; b=$b; c1="$3";c2="${4:-$3}";
                library(strawr); is.chr=FALSE;
                x=readHicChroms(hic)
                is.chr=ifelse( length(grep("chr",x[,1])) > 0,T,F);
                ## assume input has chr prefix
                if(!is.chr){
                        c1=gsub("chr","",c1); c2=gsub("chr","",c2);
                }
                c1h=strsplit(c1,":")[[1]][1]; c2h=strsplit(c2,":")[[1]][1];
                if( which(x[,1] == c1h) > which(x[,1] == c2h)){
                        tmp=c1h;c1h=c2h;c2h=tmp;
                        ## c1 c2 swapping is not needed
                }
        ff=function(x){
                x=strsplit(x,":")[[1]];
                if(length(x) == 3){
                        start=as.integer(x[2]);end=as.integer(x[3]);
                        end=ifelse(end-b < start, start, end-b)
                        return( paste(x[1],start,end,sep=":"))
                }
                return(x);

        }
                c1=ff(c1);
                c2=ff(c2);

                res=straw(norm="NONE",fname=hic,chr1loc=c1, chr2loc=c2,  unit="BP", binsize=b);
                if(nrow(res)>0){
                        if(!is.chr){ c1h=paste0("chr",c1h);c2h=paste0("chr",c2h);}
                        try({
                                dummy=apply(res,1,function(x){
                                        cat(paste(c1h,x[1],c2h,x[2],x[3],sep="\t"),"\n",sep="")
                                });
                        },silent=T);
                }

eof
}


hic2dev(){
## dev format : https://github.com/aidenlab/juicer/wiki/Pre
usage="$FUNCNAME <hic> <bsize> <cloc1> [<cloc2>]"
if [ $# -lt 3 ];then echo "$usage"; return;fi
        if [ -f $3 ];then
                cat $3 | while read -r line;do
                        cloc1=`echo "$line" | cut -f 1-3 | tr "\t" ":"`;
                        cloc2=`echo "$line" | cut -f 4-6 | tr "\t" ":"`;
                        _hic2dev ${@:1:2} $cloc1 $cloc2
                done
        else
                _hic2dev $@
        fi
}


hic2bedpe(){
	hic2dev $@ | awk -v OFS="\t" -v B=$2 '{ print $1,$2,$2+B,$3,$4,$4+B,$5;}'
}


norm-factor(){
        local f=~/lab-data/hg38/$1/mergeStats.txt
        hm_total1=`grep valid_interaction_rmdup $f | cut -f 2`
        mm_total1=`grep valid_interaction_rmdup $f | cut -f 3`
        total1=$(( $hm_total1 + $mm_total1 ))
        norm_factor1=`perl -e "print 1000000 / $total1"`;
        aqua_factor1=`perl -e "print $hm_total1 / $mm_total1"`;
        echo -e "$norm_factor1\t$aqua_factor1";
}

bedpe-norm(){
        local bedpe=$1;
        local sample=$2;
        local norm=$3;
        local FC=1;
        if [ $norm == "aqua" ];then
                FC=$( norm-factor $sample | perl -ne 'chomp;my ($x,$y)=split/\s+/,$_; print $x*$y,"\n";' );
        else
                FC=$( norm-factor $sample | perl -ne 'chomp;my ($x,$y)=split/\s+/,$_; print $x,"\n";' );
        fi
        cat $1 | awk -v OFS="\t" -v FC=$FC '{ print $1,$2,$3,$4,$5,$6,$7*FC;}'
}

bedpe-count(){
        intersectBed -a ${1/-/stdin}  -b $2 -wa -wb  |\
        perl -e 'use strict; my %r=();
        while(<STDIN>){ chomp; my @d=split/\t/,$_; my $n=7; # a columns
                #print "$d[3]==$d[10] && $d[4] < $d[12] && $d[5] > $d[11]\n";
                if( $d[3]==$d[$n+3] && $d[4] < $d[$n+5] && $d[5] > $d[$n+4] ){ # intersect
                        $r{ join("\t",@d[7..$#d]) }{s} += $d[6];  ## sum
                        $r{ join("\t",@d[7..$#d]) }{n} ++;
                }
        }
        #map { print $_,"\t",$r{$_}{s}/$r{$_}{n},"\n"; } keys %r;
        map { print $_,"\t",$r{$_}{s},"\n"; } keys %r;
        '
}



annotate-loops(){
usage="$FUNCNAME <hic> <bin_size> <bedpe> [options]
 [options]:
 	--aqua_norm
"
if [ $# -lt 3 ];then echo "$usage";return;fi
local hic=$( realpath $1);
local B=$2;
local P=$3;
	an=0;
	if [ ! -f $hic ];then
	 	hic=~/lab-data/hg38/$1/$1.allValidPairs.hic
 	fi

	if [ $4=="--aqua_norm" ];then
		an=1;
	fi
	cut -f 1 $P | sort -u | while read -r c;do
		hic2dev $hic $B $c > a 
		return
	done
#	{
#		if [ ${Q:0:1} == "T" ];then bedpe-norm - $A aqua;
#		else bedpe-norm - $A norm; fi
#	}
}



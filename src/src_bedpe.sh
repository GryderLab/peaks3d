#!/bin/bash -l
bedpe-count(){
        intersectBed -a $1 -b $2 -wa -wb  |\
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
#echo	"chr1	100	200	chr1	300	400	1
#chr1	100	200	chr1	400	500	2"	>	tmp.a
#echo	"chr1	100	200	chr1	300	500"	>	tmp.b
#bedpe-count  tmp.a tmp.b

bedpe-count-bed(){
        cat $1 | perl -ne 'chomp;my@d=split/\t/,$_;
                print join("\t",@d[0..2],"L",join("\t",@d)),"\n";
                print join("\t",@d[3..5],"R",join("\t",@d)),"\n";
                ' | intersectBed -a stdin -b $2 -wa -wb  | perl -e 'use strict; my %r=();
                while(<STDIN>){ chomp;my @d=split/\t/,$_;
                        my $k=join("\t",@d[4..9]),"\n";
                        my $p=join("\t",@d[11..13]); ## ignore other columns
                        $r{$k}{$d[3]}{$p} ++;
                        $r{$k}{V} = $d[10];
                }
                my %r2=();
                foreach my $k (keys %r){
                foreach my $l (keys %{$r{$k}{L}}){
                foreach my $r (keys %{$r{$k}{R}}){
                        $r2{ join("\t",$l,$r) } += $r{$k}{V};
                }}}

                foreach my $k (keys %r2){
                        print $k,"\t",$r2{$k},"\n";
                }

                '
}



bedpe-join(){
usage="$FUNCNAME <bedpe> [<bedpe>..] "
if [ $# -lt 2 ];then echo "$usage"; return;fi
perl -e 'use strict;
        my %r=();
        my $fid=-1;
        foreach my $f (@ARGV){
                $fid++;
                open(my $fh,"<",$f) or die "$!";
                while(<$fh>){chomp; my @d=split/\t/,$_;
			if($d[1]+$d[2] > $d[4] + $d[5] ){
				$r{ join("\t",@d[3..5],@d[0..2]) }{$fid}=join(",",@d[6..$#d]);
			}else{
				$r{ join("\t",@d[0..5]) }{$fid}=join(",",@d[6..$#d]);
			}
                }
                close($fh);
        }
        foreach my $k (keys %r){
                print "$k\t";
                print join("\t",map { defined $r{$k}{$_} ? $r{$k}{$_} : "NULL" } 0..$fid);
                print "\n";
        }
' $@
}


# correct end ordering and print uniq pair-ends 
bedpe-uniq(){
usage="$FUNCNAME <bedpe>"
if [ $# -lt 1 ];then echo "$usage";return; fi
        cat $1 | perl -ne 'chomp;my@d=split/\t/,$_;
		if($d[0] eq $d[3] && $d[1]+$d[2] > $d[4]+$d[5] ){
                        print join("\t",@d[3..5],@d[0..2],@d[6..$#d] ),"\n";
		}else{
                        print join("\t",@d),"\n";
                }
        '| sort -u
}


build-bedpe(){
usage="Usage: $FUNCNAME <bed> <bed> <tad.bed> "
if [ $# -lt 1 ];then echo "$usage"; return; fi
	local tmp=$( mktemp -d );
	cut -f1-3 $3 > $tmp/a
        {
                awk -v OFS="\t" '{ print $1,$2,$3,"A";}' $1;
                awk -v OFS="\t" '{ print $1,$2,$3,"B";}' $2;
        } |  intersectBed -a $tmp/a -b stdin  -wa -wb  |\
        groupBy -i stdin -g 1,2,3 -c 5,6,7 -o collapse,collapse,collapse |\
        perl -ne 'chomp;my @d=split/\t/,$_;
                my @st=split/,/,$d[3];
                my @ed=split/,/,$d[4];
                my @ab=split/,/,$d[5];
                my @z= sort { $a->[0]<=>$b->[0] } map {[$st[$_],$ed[$_],$ab[$_]]} 0..$#st;
                for(my$i=0;$i<$#z;$i++){
                        for(my$j=$i+1;$j<=$#z;$j++){
                                if($z[$i]->[0] != $z[$j]->[0] &&  $z[$i]->[2] ne $z[$j]->[2]){
                                        print join("\t",$d[0],$z[$i]->[0],$z[$i]->[1],$d[0],$z[$j]->[0],$z[$j]->[1]),"\n";
                                }
                        }
                }

        '
}


## example
#cd /home/ubuntu/shared/projects/aria
#tad=/home/ubuntu/shared/hg38/ref/TAD_goldstandard.hg38.bed
#build-bedpe \
#        <( head -n 1000 peaks/LNCaPXIP_AR_122121_SEK_NovoG_MACS_p-7.nobl.bed ) \
#        <( head -n 1000 peaks/LNCaPXIP_AR_122121_SEK_NovoG_MACS_p-7.nobl.bed ) \
#        $tad


bedpe-order(){
## not implemented
	cat $1 | perl -e '
                my @chroms= (1..22,"X","Y");
                my %ch=map { $chroms[$_] => $_ + 1 } 0..$#chroms;
                #map { print "$_:$ch{$_}\n"} keys %ch;
                sub _order{ my ($c1,$s1,$e1,$c2,$s2,$e2)=@_;
                        if($c1 eq $c2 && $s1+$e1 > $s2+$e2 || $ch{$c1} > $ch{$c2} ){
                                return($c2,$s2,$e2,$c1,$s1,$e1);
                        }
                        return($c1,$s1,$e1,$c2,$s2,$e2);
                }'
}
bedpe-in-bed(){
usage="
        Usage: $FUNCNAME <bedpe> <bed> [intersectBed options]
        [intersectBed options ] : default -f 1 -wa -u
				: use -f 1 -v to inverse
"
if [ $# -lt 2 ];then echo "$usage"; return; fi

local btop=" -f 1 -wa -u "
if [ $# -gt 2 ];then btop=" ${@:3}";fi

        cat $1 | perl -ne 'chomp;my @d=split/\t/,$_;
                my ($s,$e)=($d[1],$d[5]);
                if($s > $e){ my $tmp=$e;$e=$s;$s=$tmp;}
                print $d[0],"\t",$s,"\t",$e,"\t",join("@",@d),"\n";
        ' | intersectBed -a stdin -b $2 $btop  |\
        cut -f 4 | tr "@" "\t"
}


bedpe-filter(){
	usage="
	Usage: $FUNCNAME <bedpe> <mindist> <maxdist> <minvalue>
		mindist : default 5k
		maxdist : default 5m
		minvalue : default 0.000000001

	"
	if [ $# -lt 1 ];then echo "$usage"; return; fi

	km(){ echo $1 | sed s/[mM]/000k/ | sed s/[kK]/000/; }
	local mind=`km ${2:-5k}`;
	local maxd=`km ${3:-3m}`;
	local minv=${4:-0.000000001};
	cat $1 | perl -ne '
		my $mind='$mind';
		my $maxd='$maxd';
		my $minv='$minv';
		chomp; my ($c1,$s1,$e1,$c2,$s2,$e2,$v) = split /\t/,$_;
		my $d= abs($s1 + $e1 - $s2 - $e2) * 0.5;
		next if( $d > $maxd || $d < $mind || defined $v && $v < $minv);
		print $_,"\n";
	'
}



km(){
	echo $1 | sed s/[mM]/000k/ | sed s/[kK]/000/
}

bedpe-1d-cluster(){
usage="
	Usage: $FUNCNAME <bedpe> 
	generate bedGraph file of contact intervals
	
";

if [ $# -lt 1 ];then echo "$usage";return;fi
	cat $1 | perl -e 'use strict;
		my @chroms= (1..22,"X","Y");
		my %ch=map { $chroms[$_] => $_ + 1 } 0..$#chroms;
		#map { print "$_:$ch{$_}\n"} keys %ch;

		sub _order{ my ($c1,$s1,$e1,$c2,$s2,$e2)=@_;
			if($c1 eq $c2 && $s1+$e1 > $s2+$e2 || $ch{$c1} > $ch{$c2} ){
				return($c2,$s2,$e2,$c1,$s1,$e1);
			}
			return($c1,$s1,$e1,$c2,$s2,$e2);
		}
		my %r=();
		while(<STDIN>){ chomp; my @d=split/\t/,$_;
			## [    ]----[    ]
			## [              ]
			my @x=sort {$a<=>$b} ($d[1],$d[5]);
			$r{$d[0]}{$x[0]} ++;
			$r{$d[0]}{$x[$#x]} --;
		}
		foreach my $c (keys %r){
			my @x=sort {$a<=>$b} keys %{$r{$c}};
			my $n=0;
			map { 
				$n = $n + $r{$c}{$x[$_]};
				if($n > 0){
					print join("\t",$c,$x[$_],$x[$_+1],$n),"\n";	
				}
			} 0..($#x - 1);
		}

	'	
}




bedpe-cluster(){
usage="Usage: $0 <bedpe> [<min_loops>]
<min_loops> : minimum number of loops (defaut=2);

"
if [ $# -lt 1 ];then echo "$usage";return;fi

        cat $1 | perl -e 'use strict; my $thre='${2:-2}';
        sub order{ my ($xx,$yy)=@_;
                my @x=split/\t/,$xx;
                my @y=split/\t/,$yy;
                if( $x[1] + $x[2] > $y[2] + $y[2] ){
                        return($yy,$xx);
                }
                return ($xx,$yy);
        }
        #print order("a\t3\t4","a\t1\t2"),"\n";
        #exit(1);

        sub sear{ my ($t,$k,$v,$p)=@_;
                return if defined $v->{$k};
                $v->{$k} ++;
                foreach my $x (keys %{$t->{$k}}){
                        my ($le,$ri) = order($x,$k);
                        $p->{ "$le\t$ri" } ++;
                        sear($t,$x,$v,$p);
                }

        }

        my %r=();
        my %cpm=();
        while(<STDIN>){chomp; my @d=split/\s+/,$_;
		next if($#d < 6);
                my ($x,$y) = ( join("\t", @d[0..2]), join("\t",@d[3..5]));
                $r{$x}{$y}++;
                $r{$y}{$x}++;
                $cpm{$x}{$y}=$d[6];
                $cpm{$y}{$x}=$d[6];
        }

        my %v=(); ## visited
        my $gid=0;
        foreach my $x (keys %r){
                my @yy=keys %{$r{$x}};
                if ( $#yy + 1 >= $thre ){
                        my %p1=();
                        sear(\%r,$x, \%v,\%p1);
                        my @p=keys %p1;
                        if( $#p + 1 >= $thre){
                                $gid ++;
                                print join("\n", map {"$_\t".($#p+1)."\tcluster$gid"} @p),"\n";
                        }
                }
        }
'
}

bedpe-cluster-test(){
echo "chr1      100     200     chr2    300     400     2
chr2    300     400     chr3    500     600     2
chrX    300     400     chrX    500     600     2
chrX    300     400     chrX    600     700     2
" | bedpe-cluster - 1
}


bedpe-touch-bed_old(){
usage="
Usage: $FUNCNAME <bedpe> <bed> 
	-v : reverse meaing no touch
"

        if [ $# -lt 2 ];then echo "$usage";exit 1; fi

	cat $1 | perl -ne 'chomp; my@d=split/\t/,$_;
	next if($#d < 5);
	print join("\t",@d[0..2], join("@",@d)),"\n";
	print join("\t",@d[3..5], join("@",@d)),"\n";
	' |\
	intersectBed  -a stdin -b $2 -wa -wb 
}

bedpe-touch-bed(){
usage="
Usage: $FUNCNAME <bedpe> <bed>
"

        if [ $# -lt 2 ];then echo "$usage";exit 1; fi

        ## merge clusters
        cat $1 | perl -ne 'chomp; my@d=split/\t/,$_;
        next if($#d < 5);
        print join("\t",@d[0..2], join("@",@d[0..6])),"\n";
        print join("\t",@d[3..5], join("@",@d[0..6])),"\n";
        ' |\
        intersectBed  -a stdin -b $2 -wa -wb  |\
        perl -ne 'use strict;
                my %gg=();
                my %gp=();
                while(<STDIN>){chomp; my @d=split/\t/,$_;
                        my $k=join("\t",@d[0..2]);
                        my $p=join("\t",split/@/,$d[3]);
                        my $g=join(":",@d[4..7]);
                        $gg{$k}{$g} ++;
                        $gp{$g}{$p} ++;
                }
                my $gid=1;
                foreach my $k (keys %gg){
                        my $ggk=join(",",keys %{$gg{$k}});
                        foreach my $g (keys %{$gg{$k}}){
                                foreach my $p (keys %{$gp{$g}}){
                                        print $p,"\t$ggk\tcluster$gid\n";
                                }
                        }
                        $gid ++;
                }
        '
}

bedpe-cluster-entropy(){
usage="$FUNCTION <clustered.bedpe> <annotated.bedpe> [<annotated.bedpe>..]

##INPUT ----------------------------------------------------------------
   [ ]-------[    ]                                      
             [    ]--------[    ]   		
	                   [    ]------[    ] : cluster 1

   [ ]---a---[    ]        
             [    ]---b----[    ]             : annotated loops 1
	     

             [    ]---c----[    ]   		
	                   [    ]--d---[    ] : annotated loops 2

##OUTPUT-----------------------------------------------------------------
   [ ]---a,0-[    ]                                      
             [    ]---b,c--[    ]   		
	                   [    ]--0,d-[    ] : cluster 1
"
if [ $# -lt 1 ];then echo "$usage"; return; fi

	bedpe-join $@ |\
	perl -e 'use strict;
		my %r=();
		sub min{ my @x=@_; my @y=sort{$a<=>$b} grep {$_} @x; return $y[0];}
		sub max{ my @x=@_; my @y=sort{$b<=>$a} grep {$_} @x; return $y[0];}
		sub nrm{ my ($x,$y)=@_; return $y == 0 ? 0 : $x/$y; }
		sub ent{ my @x=@_; my $ep=0.0001; my $r=0;
			my @xx=map { $_ == 0 ? $ep : $_ } @x;
			my $s=0; map { $s+= $_; } @xx;
			map { $r+= - $_ * log($_); } map { $_/$s } @xx;
			return($r);
		}

		while(<STDIN>){chomp;my@d=split/\t/,$_;
			my @x=@d[7..$#d];
			foreach my $i (0..$#x){
				$r{$d[6]}{join("\t",@d[0..5])}[$i] += $x[$i];
			}
		}
		## todo: trans chrom
		foreach my $cluster (keys %r){
			my @m=(); ## row matrix
			my $j =0;
			my ($c,$s,$e);
			foreach my $x (keys %{$r{$cluster}}){
				my @d=split/\t/,$x;
				($c,$s,$e)=($d[0],min($s,$d[1],$d[4]),max($e,$d[2],$d[5]));

				my @y=@{$r{$cluster}{$x}};
				foreach my $i (0..$#y){
					$m[$i][$j] = $y[$i];
				}
				$j ++;
			}
			my ($nc,$ci)=split/,/,$cluster;
			#if($cluster=~/cluster800$/){
			print join("\t",$c,$s,$e,$ci,"nloops=$nc",map{ ent(@{$_}) } @m),"\n";
				#map { print join("\t",@{$_}),"\n";} @m;
			#}
		}
		
	' 
}


clusteredbedpe2bed(){
        cat $1 | groupBy -i stdin -g 8 -c 1,2,6,7,8 -o distinct,min,max,distinct,distinct | cut -f 2- |\
        cluster-short-long-density.sh - $sample $binsize $mindist 
}



bedpe-untouch-bed(){
usage="
Usage: $FUNCNAME <bedpe> <bed> 
"
	local tmpd=$( mktemp -d );
	cat $1 > $tmpd/a
	bedpe-touch-bed $tmpd/a $2 > $tmpd/b
	cat $tmpd/a | perl -e 'use strict;
		open(my $fh,"<","'$tmpd/b'") or die $!;
		my %r=();
		while(<$fh>){chomp; $r{$_}++; }
		close($fh);
		while(<STDIN>){chomp;
			next if defined $r{$_};
			print $_,"\n";
		}
	'
}

bedpe-in-bed(){
usage="
        Usage: $FUNCNAME <bedpe> <bed> [intersectBed options]
        [intersectBed options ] : default -f 1 -wa -u
"
if [ $# -lt 2 ];then echo "$usage"; return; fi

local btop=" -f 1 -wa -u "
if [ $# -gt 2 ];then btop=" ${@:3}";fi

        cat $1 | perl -ne 'chomp;my @d=split/\t/,$_;
                my ($s,$e)=($d[1],$d[5]);
                if($s > $e){ my $tmp=$e;$e=$s;$s=$tmp;}
                print $d[0],"\t",$s,"\t",$e,"\t",join("@",@d),"\n";
        ' | intersectBed -a stdin -b $2 $btop  |\
        cut -f 4 | tr "@" "\t"
}








#!/bin/bash -l

usage="
	Usage: $FUNCNAME <annotate.bedpe> [<output>.pdf]
	Draw distance and score QC plots
	output : output prefix (default=qc)
"
if [ $# -lt 1 ];then echo "$usage";exit 1;fi
input=$1;
output=${2:-qc}


## read a bedpe file obtained from annotate_loops.sh 
## draw piechart of zero and non-zero entries
## draw (log10) distance distributions for zero and non-zero intervals
cat<<-'EOF' | R --vanilla --args $input $output
args=commandArgs(trailingOnly=T);
if( args[1] == "test" ){
	library(circlize)
	d=cbind(generateRandomBed(nr=10,nc=0),generateRandomBed(nr=10,nc=3))
	d[,7]=abs(d[,7]);
	d[sample(1:nrow(d),size=round(nrow(d))/2),7]=0
}else{
	d=read.table(args[1],header=F,skip=1);
}
colnames(d)=c("c1","s1","e1","c2","s2","e2","v1","v2","col")[1:ncol(d)];
library(dplyr)
library(ggplot2)
#install.packages("cowplot", repos="http://cran.case.edu")
library(cowplot)
d=d%>% mutate(dist= 0.5* abs(s1+e1 - s2-e2 ), group= ifelse( abs(v1) > 0, "nonzeroCPM","zeroCPM") ) 



x=as.data.frame(t(table(d$group)))
colnames(x)=c("A","group","Freq")
x=x%>%arrange(desc(group))
x$ypos=cumsum(x$Freq) - 0.5*x$Freq;

## counts are added
p1=ggplot(x,aes(x = "", y = Freq, fill = group)) +
        geom_bar(stat="Identity") +
        coord_polar(theta="y",start=0,direction=1) +
        theme_void() + ylab("AQuA CPM") +
        geom_text(aes(y=ypos,label=Freq)) +
        ggtitle("AQuA CPM")


myhist=function(d,title){
	ggplot(d,aes(x=v1)) + 
	geom_histogram() + 
	xlab("CPM") + ylab("Frequency") + ggtitle(title) +
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

}
myhistlg10=function(d,title){
	ggplot(d,aes(x=log10(v1+1))) + 
	geom_histogram() + 
	xlab("log10 CPM + 1") + ylab("Frequency") + ggtitle(title) +
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

}

p2=myhist(d, "ALL CPM")
p3=myhist(d[d$group=="nonzeroCPM",], "Non-zero CPM")
p4=ggplot(d[d$group=="nonzeroCPM",], aes(x=dist,y=v1)) + geom_point() + xlab("distance") + ylab("Non-zero CPM") +
	theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p5=myhistlg10(d, "ALL CPM")
p6=myhistlg10(d[d$group=="nonzeroCPM",], "Non-zero CPM")
	
#pdf(paste0(args[2],".pdf"))
png(paste0(args[2],".png"))
plot_grid(p1,p2,p3,p4,p5,p6, align="v")
dev.off()
EOF



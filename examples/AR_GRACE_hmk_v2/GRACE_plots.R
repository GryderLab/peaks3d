
### GRACE plots #########
# Genic Rank of Active  #
# Clustered Enhancers   #
# Berkley Gryder, 2022  #
# Gryderlab.com         #
#########################
## 1. Stats on annotated and filtered BEDPE
#loops.name = "LNCaPXIP_CnR_K27ac_p-7"
#peaks3Dname = "LNCaPXIP_CnR_K27ac_p-7_LNCaPXIP_DMSO_H3K27ac_10k_4m"

peaks3Dname = "LNCaPXIP_DMSO_H3K27ac_122121_SEK_NovoG_MACS_p-7.nobl_LNCaPXIP_DMSO_H3K27ac_25k_3m_0p05";
loops.name="LNCaPXIP_DMSO_H3K27ac_122121_SEK_NovoG_MACS_p-7.nobl";
output="AR_GRACE_hmk_v2";
setwd(output);
odir="GRACE_plots/";

loops = read.table(paste0(peaks3Dname,".filtered.bedpe")) #*.filtered.bedpe
  colnames(loops)=c("left_chr","left_start","left_end","right_chr","right_start","right_end","AQuA_CPM")
  loops$distance = abs(loops$right_start - loops$left_start)
  loops$distance_log10kb = log10(loops$distance/1000)
library(ggplot2);library(viridis);library(MASS)
# Get density of points in 2 dimensions. @param x A numeric vector. @param y A numeric vector.
# @param n Create a square n by n grid to compute density. @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
loops$density <- get_density(loops$distance_log10kb, loops$AQuA_CPM, n = 50)
loopsubtitle = paste(loops.name,"\n","n = ",nrow(loops)," loops")
plot.distance_CPM = ggplot(loops) + geom_point(aes(distance_log10kb, AQuA_CPM, color = density))+
  scale_color_viridis(option="rocket")+labs(subtitle = loopsubtitle) 
  
## 2. Clustering propensity among loops
clusters = read.table(paste0(peaks3Dname,".clustered.bedGraph")) #*.clustered.bedGraph
  colnames(clusters)=c("chr","start","end","clusters")
  clusters$cluster.degree = clusters$clusters; clusters$cluster.degree=as.character( sprintf("%02d",clusters$cluster.degree));clusters$cluster.degree[clusters$clusters>9] = "10+"
plot.clusterstats = ggplot(clusters,aes(cluster.degree))+geom_bar(stat="count")+
  stat_count(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5)+labs(subtitle = paste("Clustered loops in","\n",loops.name)) 
##Plot multiple attributes
library(grid);library(gridExtra);library(cowplot);library(ggplot2)
p.save = plot_grid(plot.distance_CPM, plot.clusterstats, align = "v", nrow = 1, rel_heights = c(2/4, 2/4),  labels = "auto")
#save_plot(paste(odir,loops.name,".summary.pdf",sep=""), p.save, base_height = 5, base_aspect_ratio = 1, base_width =7)
save_plot(paste(odir,peaks3Dname,".summary.pdf",sep=""), p.save, base_height = 5, base_aspect_ratio = 1, base_width =7)
## 3. Rank order of Cluster Enhancers
density.clusters = read.table(paste0(peaks3Dname,".clustered.density.bed")) # *.clustered.density.bed
  #colnames(density.clusters)=c("chr","start","end","loops","region.name","shortCPMs","longCPMs","shortAQuA","longAQuA")
  colnames(density.clusters)=c("chr","start","end","loops","region.name","shortAQuA","longAQuA")
  density.clusters = density.clusters[,c("chr","start","end","shortAQuA","longAQuA","region.name","loops")]
  density.clusters$group = paste0("merged clusters =", nrow(density.clusters))
  
density.loops = read.table(paste0(peaks3Dname,".notclustered.density.bed")) # *.notclustered.density.bed
  colnames(density.loops)   =c("chr","start","end","shortAQuA","longAQuA")
  density.loops$region.name = paste0("loop",rownames(density.loops))
  density.loops$loops = 1
  density.loops$group = paste0("loops =", nrow(density.loops))
  
density.peaks = read.table(paste0(peaks3Dname,".unloop.density.bed")) # *.unloop.density.bed
  colnames(density.peaks)   =c("chr","start","end","shortAQuA","longAQuA")
  density.peaks$region.name = paste0("peak",rownames(density.peaks))
  density.peaks$loops = 0
  density.peaks$group = paste0("isolated peaks =",nrow(density.peaks))
  
density.all = rbind(density.clusters,density.loops,density.peaks)
my3cols <- c("#E7B800", "#2E9FDF", "#FC4E07")
cluster = ggplot(density.all, aes(x=loops,y=log2(longAQuA+1),color=group))+geom_point()+scale_color_manual(values = my3cols)+theme_bw()
scatter = ggplot(density.all, aes(x=log2(longAQuA+1),y=log2(shortAQuA+1),color=group))+geom_point(alpha=0.5)+theme_bw()+scale_color_manual(values = my3cols)+theme(legend.position = "none")
x_boxes = ggplot(density.all, aes(x=log2(longAQuA+1),y=1,color=group))+geom_boxplot()+theme_bw()+scale_color_manual(values = my3cols)+theme(legend.position = "none")
y_boxes = ggplot(density.all, aes(x=1,y=log2(shortAQuA+1),color=group))+geom_boxplot()+theme_bw()+scale_color_manual(values = my3cols)+theme(legend.position = "none")
##Print density comparisons
p.save = plot_grid(x_boxes, cluster, scatter, y_boxes,  nrow = 2, ncol = 2, labels = "auto")
#save_plot(paste(odir,loops.name,".density.pdf",sep=""), plot = p.save, base_height = 5, base_width =7 )
save_plot(paste(odir,peaks3Dname,".density.pdf",sep=""), plot = p.save, base_height = 5, base_width =7 )
#Rank
density.all$rank.CPM = rank(x = -(density.all$shortAQuA+density.all$longAQuA), ties.method = "random")
library(tidyr)
density.gather = gather(density.all, dist, AQuA, shortAQuA:longAQuA, factor_key=TRUE)
library(dplyr)
density.stats = density.gather %>% group_by(dist,group) %>% summarize(mean = mean(AQuA),sum = sum(AQuA))
p.stats = 
  ggplot(density.stats, aes(x=dist,y=sum,fill=group))+ geom_bar(position = "dodge",stat = "identity")+scale_fill_manual(values = my3cols)+theme_bw()
p.rank.CPM = 
  ggplot(density.all, aes(x=rank.CPM,y=(shortAQuA+longAQuA),color=group))+geom_point(alpha=0.3)+scale_color_manual(values = my3cols)+theme_bw()+theme(legend.position = "none")+ scale_x_reverse()
p.rank.box =
  ggplot(density.all, aes(y=(shortAQuA+longAQuA),color=group))+ geom_boxplot()+scale_color_manual(values = my3cols)+theme_bw()#+facet_wrap(~group)
p.rank.dens = 
  ggplot(density.all, aes(x=rank.CPM,fill=group))+geom_density(alpha=0.8)+scale_fill_manual(values = my3cols)+theme_bw()+theme(legend.position = "none")+ scale_x_reverse()
##Print density comparisons
p.save = plot_grid(p.rank.dens, p.stats, p.rank.CPM, p.rank.box, nrow = 2, labels = "auto")
#save_plot(paste(odir,loops.name,".ranked.pdf",sep=""), plot = p.save, base_height = 5, base_width =7 )
save_plot(paste(odir,peaks3Dname,".ranked.pdf",sep=""), plot = p.save, base_height = 5, base_width =7 )
## Write out BED file for EDEN purposes
colnames(density.all)[1] = "#chr"
#write.table(density.all,file=paste0(odir,loops.name,".ranked.bed"),append = F,row.names = F, col.names = T,quote = F, sep="\t")
write.table(density.all,file=paste0(odir,peaks3Dname,".ranked.bed"),append = F,row.names = F, col.names = T,quote = F, sep="\t")


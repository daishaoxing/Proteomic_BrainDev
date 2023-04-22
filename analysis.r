#! /usr/bin/env Rscript
# library -----------------------------------------------------------------
library(ggplot2)
library(scater)
library(Matrix)
library(cowplot)
library(scater)
library(destiny)
library(umap)
library(ggthemes)
library(dplyr)
library(patchwork)
library(circlize)
library(ggplot2)
library(Seurat)
library(scater)
library(Matrix)
library(cowplot)
library(scater)
library(destiny)
library(umap)
library(ggthemes)
library(dplyr)
library(patchwork)
library(ggsci)
library(scales)
library(tidyverse)
library(magick)
library(magrittr)
library(yyplot)
library(devtools)
library(ComplexHeatmap)
library(ggplotify)
library(circlize)
library(pheatmap)
library(hrbrthemes)
library(ggrepel)
library(ggVennDiagram)
library(RVenn)
library(devtools)
library(dplyr)
library(ggplot2)
library(loomR)
library(hdf5r)
library(dyno)
library(tidyverse)
library(ggplot2)
library(scater)
library(Matrix)
library(cowplot)
library(scater)
library(destiny)
library(umap)
library(ggthemes)
library(dplyr)
library(patchwork)
library(ggsci)
library(scales)
library(tidyverse)
library(magick)
library(magrittr)
library(yyplot)
library(devtools)
library(ggplotify)
library(PCAtools)
library(corrplot)
library(reshape2)
library(ComplexHeatmap)
library(UpSetR)
library(future)
library(paletteer)
library(data.table)
library(tidydr)
library(clusterProfiler)
######################## /========= Block1 =========/ #########################
############################## <Main Figures> #################################
##############################################################################-
# — Data preprocess --------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())

tmp_meta <- read.delim("data/metadata.txt") %>% 
  {.[.$time%in%"A",]$time<-"F50";.} %>% 
  {.[.$time%in%"B",]$time<-"F90";.} %>% 
  {.[.$time%in%"C",]$time<-"F120";.} %>% 
  {.[.$time%in%"D",]$time<-"P3";.} %>% 
  {.[.$region2%in%"1",]$region2<-"FL";.} %>% 
  {.[.$region2%in%"2",]$region2<-"TL";.} %>% 
  {.[.$region2%in%"3",]$region2<-"PL";.} %>% 
  {.[.$region2%in%"4",]$region2<-"V1";.} %>% 
  {.[.$region2%in%"5",]$region2<-"CB";.} %>% 
  {.[.$region2%in%"6",]$region2<-"STr";.} %>% 
  {.[.$region2%in%"7",]$region2<-"Hipp";.} %>% 
  {.[.$region2%in%"8",]$region2<-"MD";.} %>% 
  {.[.$region2%in%"9",]$region2<-"Amy";.}

tmp_expr <- read.delim("data/protein_Counts.Table_S2.txt") %>% 
  .[.$Gene.name!="---",] %>% 
  {.$ID<-paste(.$Gene.name,"|",.$Protein.accession,sep = "");.} %>%
  {rownames(.)<-.$ID;.} %>%
  .[,tmp_meta$sample] %>% 
  .[rowSums(.)!=0,]

save(tmp_expr,tmp_meta,file = "result/Proteomics.tmp_expr.tmp_meta.RData")

# — Figure1 ----------------------------------------------------------
# —— Fig1B Brain Region Protein Abundance PointPlot ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))

tmp_meta %<>% {.$Minor_Region<-.$region;.} %<>%
  {.[.$region%in%c("aPFC","sPFG"),]$Minor_Region<-"aPFC/sPFG";.} %<>%
  {.[.$region%in%c("PFC","pPFC","mPFG"),]$Minor_Region<-"PFC/pPFC/mPFG";.} %<>% 
  {.[.$region%in%c("TL","MT"),]$Minor_Region<-"TL/mTG";.} %<>%
  {.[.$region%in%c("PL","sPL"),]$Minor_Region<-"PL/sPL";.}


tmp_plotdf<-tmp_expr %>% rownames_to_column() %>% gather(Sample,Expr,-rowname) %>% 
  {.$Gene<-gsub("\\|.*","",.$rowname);.} %>% 
  aggregate(Expr ~ Sample, data = ., sum) %>% 
  left_join(.,tmp_meta,by=c("Sample"="sample")) %>% 
  {.$Expr<-log2(.$Expr+1);.} %>% 
  group_by(time,Minor_Region) %>% 
  dplyr::summarise(Mean_Expr = mean(Expr),SD_Expr=sd(Expr),
                   time=time,Minor_Region=Minor_Region,Expr=Expr)

Time_sort<-c("aPFC/sPFG","PFC/pPFC/mPFG","iPFG","OFC","IC","M1","aCG","sTG","TL/mTG","iTG",
             "S1","PL/sPL","V1","CB","STr","Hipp","MD","Amy")
Time_color<-c("#a5d5a7","#8fc69b","#7bbe6e","#6dbf53","#47ae51","#32994d",
              "#367941","#29a5a5","#258a8d","#277576","#4b80c2","#3470b8",
              "#1c6eb7","#dece62","#f58737","#f15c44","#b92645","#842446")

tmp_plotdf$Minor_Region<-factor(tmp_plotdf$Minor_Region,levels = Time_sort)
tmp_plotdf$time<-factor(tmp_plotdf$time,levels = c("F50","F90","F120","P3"))
(p1<-ggplot(tmp_plotdf,aes(x = Minor_Region,y=Expr,colour=Minor_Region))+
    geom_jitter(position=position_jitter(w=0.1),size=1) +
    geom_point(aes(x = Minor_Region, y = Mean_Expr), size=2, color="black", shape="-", alpha=.3) +
    geom_errorbar(mapping = aes(x = Minor_Region, y = Mean_Expr,
                                ymin = Mean_Expr - SD_Expr, ymax = Mean_Expr + SD_Expr),
                  size=.5, color="black", width=.8 ,alpha=.3)+
    # geom_boxplot(alpha=.6,outlier.size = 1)+
    scale_fill_manual(values = Time_color)+
    scale_color_manual(values = Time_color)+
    # geom_point()+
    guides(colour = guide_legend(title = "Region"))+
    guides(fill = guide_legend(title = "Region"))+
    facet_wrap(~time,ncol = 4)+
    ylab(label = "Log2 (protein abundance)")+
    theme_bw()+
    theme(axis.title = element_text(size=18,colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          # legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="black",fill = "white"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "black")))
ggsave(filename = "figure/Fig1B.Brain_Region.Protein_Abundance.PointPlot.pdf",p1,width = 7.5,height = 5)

# —— Fig1C Protein number PointPlot ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())

tmp_meta <- read.delim("data/metadata.txt") %>% 
  {.[.$time%in%"A",]$time<-"F50";.} %>% 
  {.[.$time%in%"B",]$time<-"F90";.} %>% 
  {.[.$time%in%"C",]$time<-"F120";.} %>% 
  {.[.$time%in%"D",]$time<-"P3";.}

tmp_meta %<>% {.$Minor_Region<-.$region;.} %<>%
  {.[.$region%in%c("aPFC","sPFG"),]$Minor_Region<-"aPFC/sPFG";.} %<>%
  {.[.$region%in%c("PFC","pPFC","mPFG"),]$Minor_Region<-"PFC/pPFC/mPFG";.} %<>% 
  {.[.$region%in%c("TL","MT"),]$Minor_Region<-"TL/mTG";.} %<>%
  {.[.$region%in%c("PL","sPL"),]$Minor_Region<-"PL/sPL";.}

tmp_expr <- read.delim("data/protein_Counts.Table_S2.txt") %>% 
  {.$ID<-paste(.$Gene.name,"|",.$Protein.accession,sep = "");.} %>%
  {rownames(.)<-.$ID;.} %>%
  .[,tmp_meta$sample] %>% 
  .[rowSums(.)!=0,]

tmp_plotdf<-tmp_expr %>% {.[.>0]<-1;.} %>% colSums() %>% as.data.frame() %>% 
  {colnames(.)<-c("nProtein");.} %>% rownames_to_column() %>% 
  left_join(.,tmp_meta,by=c("rowname"="sample")) %>% 
  {.$nProtein<-.$nProtein/1000;.} %>% 
  group_by(time,Minor_Region) %>% 
  dplyr::summarise(Mean_Expr = mean(nProtein),SD_Expr=sd(nProtein),
                   time=time,Minor_Region=Minor_Region,Expr=nProtein)

tmp_plotdf %>% group_by(time) %>% dplyr::summarise(mean_Expr = mean(Expr))
# # A tibble: 4 × 2
# time  mean_Expr
# <chr>     <dbl>
# 1 F120      6017.
# 2 F50       6191.
# 3 F90       5918 
# 4 P3        5885.

Time_sort<-c("aPFC/sPFG","PFC/pPFC/mPFG","iPFG","OFC","IC","M1","aCG","sTG","TL/mTG","iTG",
             "S1","PL/sPL","V1","CB","STr","Hipp","MD","Amy")
Time_color<-c("#a5d5a7","#8fc69b","#7bbe6e","#6dbf53","#47ae51","#32994d",
              "#367941","#29a5a5","#258a8d","#277576","#4b80c2","#3470b8",
              "#1c6eb7","#dece62","#f58737","#f15c44","#b92645","#842446")

tmp_plotdf$Minor_Region<-factor(tmp_plotdf$Minor_Region,levels =Time_sort)
tmp_plotdf$time<-factor(tmp_plotdf$time,levels = c("F50","F90","F120","P3"))
(p1<-ggplot(tmp_plotdf,aes(x = Minor_Region,y=Expr,colour=Minor_Region))+
    geom_jitter(position=position_jitter(w=0.1),size=1) +
    geom_point(aes(x = Minor_Region, y = Mean_Expr), size=2, color="black", shape="-", alpha=.3) +
    geom_errorbar(mapping = aes(x = Minor_Region, y = Mean_Expr,
                                ymin = Mean_Expr - SD_Expr, ymax = Mean_Expr + SD_Expr),
                  size=.5, color="black", width=.8 ,alpha=.3)+
    # geom_boxplot(alpha=.6,outlier.size = 1)+
    scale_fill_manual(values = Time_color)+
    scale_color_manual(values = Time_color)+
    # geom_point()+
    guides(colour = guide_legend(title = "Region"))+
    guides(fill = guide_legend(title = "Region"))+
    facet_wrap(~time,ncol = 4)+
    ylab(label = "(Identified proteins)/1000")+
    theme_bw()+
    theme(axis.title = element_text(size=18,colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          # legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="black",fill = "white"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "black")))
ggsave(filename = "figure/Fig1C.Brain_Region.nProteins.PointPlot.pdf",p1,width = 7.5,height = 5)

# —— Fig1D DEP number (t.test) Bubble ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))

tmp_expr<-log2(tmp_expr+1)

tmp_out<-data.frame()
for (tmp_time in unique(tmp_meta$time)) {
  # tmp_time<-"F120"
  tmp_region<-unique(tmp_meta[tmp_meta$time%in%tmp_time,]$region2)
  for (n in tmp_region) {
    # n<-"CB"
    pvalue<-c();tstat<-c();meansControl<-c();meansCase<-c();FC<-c();log2FC<-c()
    
    Case_list<-tmp_meta[tmp_meta$time%in%tmp_time&tmp_meta$region2%in%n,]$sample
    Control_list<-tmp_meta[tmp_meta$time%in%tmp_time&tmp_meta$region2%in%setdiff(tmp_region,n),]$sample
    
    edata<-tmp_expr[,c(Case_list,Control_list)] %>% .[rowSums(.)!=0,] %>% {. + 0.01}
    for (i in 1:nrow(edata)){
      #t.test
      Control_value<-edata[i,Control_list]
      Case_value<-edata[i,Case_list]
      result<-t.test(as.numeric(Case_value), as.numeric(Control_value), paired=FALSE);
      pvalue[i]<-result[["p.value"]]
      tstat[i]<-result[["statistic"]][["t"]]
      meansControl[i]<-mean(as.numeric(Control_value))
      meansCase[i]<-mean(as.numeric(Case_value))
      FC[i]<-mean(as.numeric(Case_value))/mean(as.numeric(Control_value))
      log2FC[i]<-log2(FC[i])
    }
    p_adj<-p.adjust(pvalue, method = "fdr")
    diff_df<-data.frame(rownames(edata),pvalue,tstat,meansCase,meansControl,FC,log2FC,p_adj,stringsAsFactors=FALSE) %>% 
      {names(.)[names(.) == 'rownames.edata.'] <- 'Proteomics';.}
    tmp_out<-data.frame(Versus=paste(n,"_vs_others",sep = ""),
                        Time=tmp_time,
                        nDEP_pvalue=nrow(diff_df[diff_df$pvalue<0.05,]),
                        nDEP_padj=nrow(diff_df[diff_df$p_adj<0.05,])) %>% rbind(.,tmp_out)
  }
}

tmp_out$Region<-gsub("_vs_others","",tmp_out$Versus)
tmp_out$Time<-factor(tmp_out$Time,levels = c("F50","F90","F120","P3"))
tmp_out$Region<-factor(tmp_out$Region,levels = c("Amy","MD","Hipp","STr","CB","V1","PL","TL","FL"))
write.table(tmp_out,file = "result/nDEP.summary.xls",quote = FALSE,sep = "\t",row.names=F, col.names=T)


tmp_out$Region<-factor(tmp_out$Region,levels = c("MD","Amy","Hipp","STr","CB","V1","PL","TL","FL"))
color<-c("#35703e","#4aa751","#287374","#1c6eb7","#dece62","#f58737","#e05a44","#b92645","#7e2243") %>% rev()
(p1<-ggplot(tmp_out,aes(x = Time,y=Region,size=nDEP_padj,colour=Region,ylab=''))+
    geom_point()+
    geom_text(data=tmp_out,
              aes(x=Time, y=Region, label=nDEP_padj),
              color="white",alpha=1, size=4)+
    scale_size_continuous(range=c(10,14))+
    # scale_color_gradient2(low="black",mid="#007965",high ="#ffcc29",midpoint = 0)+
    scale_color_manual(values=color)+
    theme_classic()+
    theme(axis.text=element_text(size=15, color="black"),
          axis.title = element_blank()) + NoLegend())

p2<-tmp_out %>% group_by(Region) %>% dplyr::summarise(All_DEP = sum(nDEP_padj)) %>% 
  ggplot(.,aes(x = All_DEP,y=Region))+
  geom_bar(stat = "identity")+
  theme_classic()+
  theme(axis.text=element_text(size=15, color="black"),
        axis.title = element_blank()) + NoLegend()

p3<-tmp_out %>% group_by(Time) %>% dplyr::summarise(All_DEP = sum(nDEP_padj)) %>% 
  ggplot(.,aes(x = Time,y= All_DEP))+
  geom_bar(stat = "identity")+
  theme_classic()+
  theme(axis.text=element_text(size=15, color="black"),
        axis.title = element_blank()) + NoLegend()


p<-((p3|plot_spacer())/(p1|p2))+plot_layout(heights = c(0.5,2))
ggsave(filename = "figure/Fig1D.DEP_number.bubble.pdf",width = 7,height = 6)


# —— Fig1E Subcellular location percent line ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())

protein_family<-read.delim("data/subcellular_location.tab", header=FALSE) %>% .[,c("V2","V3")] %>% 
  {colnames(.)<-c("Gene_Name","Location");.} %>% {colnames(.)<-c("Gene_Name","Family");.} %>% na.omit()

(load("result/Proteomics.tmp_expr.tmp_meta.RData"))
tmp_out_split<-tmp_expr %>% rownames_to_column() %>% gather(Sample,Expr,-rowname) %>% 
  {.$Gene<-gsub("\\|.*","",.$rowname);.} %>% 
  left_join(.,tmp_meta,by=c("Sample"="sample")) %>% 
  left_join(.,protein_family,by=c("Gene"="Gene_Name")) %>% 
  na.omit() %>% .[.$Expr!=0,] %>% .[,c("Sample","time","region2","Family","Expr")] %>% 
  aggregate(Expr ~ Sample + time + region2 + Family, data = ., sum) %>% 
  group_by(time,region2,Sample) %>% dplyr::summarise(Sum = sum(Expr),Percent=Expr/Sum*100,Family=Family) %>% 
  as.data.frame() %>% 
  {.$time<-factor(.$time,levels = c("P3","F120","F90","F50"));.} %>% 
  .[order(.$time),]

tmp_out_split$P3_vs_F120<-NA
tmp_out_split$F120_vs_F90<-NA
tmp_out_split$F90_vs_F50<-NA

for (tmp_Region in unique(tmp_out_split$region2)) {
  # tmp_Region<-"Amy"
  tmp1<-tmp_out_split[tmp_out_split$region2%in%tmp_Region,]
  for (tmp_family in unique(tmp1$Family)) {
    # tmp_family<-"Tubulin family"
    tmp2<-tmp1[tmp1$Family%in%tmp_family,]
    if (length(unique(as.character(tmp2$time)))>1) {
      tmp_list<-unique(as.character(tmp2$time)) %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
        {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
        .[.$Versus%in%intersect(c("P3_vs_F120","F120_vs_F90","F90_vs_F50"),.$Versus),]
      for (i in 1:nrow(tmp_list)) {
        # i<-2
        if (length(tmp2[tmp2$time%in%tmp_list[i,"V1"],]$Percent)>1&length(tmp2[tmp2$time%in%tmp_list[i,"V2"],]$Percent)>1) {
          tmp_out_split[tmp_out_split$region2%in%tmp_Region&tmp_out_split$Family%in%tmp_family,] %<>% 
            {.[[paste(tmp_list[i,"V1"],"_vs_",tmp_list[i,"V2"],sep = "")]]<-t.test(as.numeric(tmp2[tmp2$time%in%tmp_list[i,"V1"],]$Percent),
                                                                                   as.numeric(tmp2[tmp2$time%in%tmp_list[i,"V2"],]$Percent))[["p.value"]];.}
        }
      }
    }
  }
}
tmp_out_split[is.na(tmp_out_split)]<-1
# tmp_out_split[tmp_out_split==1]<-NA
tmp_out_split$Signal<-NA

tmp_out_split %<>% {.$Signal<-ifelse(.$P3_vs_F120<0.05|.$F120_vs_F90<0.05|.$F90_vs_F50<0.05,"Yes","No");.}
tmp_out_split[tmp_out_split==1]<-NA
tmp_out_split$anchor<-paste(tmp_out_split$region2,tmp_out_split$Family,sep = "_")

tmp_out_merge<-tmp_expr %>% rownames_to_column() %>% gather(Sample,Expr,-rowname) %>% 
  {.$Gene<-gsub("\\|.*","",.$rowname);.} %>% 
  left_join(.,tmp_meta,by=c("Sample"="sample")) %>% 
  left_join(.,protein_family,by=c("Gene"="Gene_Name")) %>% 
  na.omit() %>% .[.$Expr!=0,] %>% .[,c("Sample","time","region2","Family","Expr")] %>% 
  aggregate(Expr ~ Sample + time + region2 + Family, data = ., sum) %>% 
  group_by(time,region2,Sample) %>% dplyr::summarise(Sum = sum(Expr),Percent=Expr/Sum*100,Family=Family) %>% 
  as.data.frame() %>% 
  group_by(time,region2,Family) %>% dplyr::summarise(Mean_Percent = mean(Percent)) %>% 
  .[order(.$time,.$region2,.$Mean_Percent,decreasing = T),] %>% 
  {.$anchor<-paste(.$region2,.$Family,sep = "_");.}

tmp_head<-tmp_out_merge %>% {.[order(.$Mean_Percent,decreasing = T),]$Family} %>% unique() %>% head(.,n=20)

tmp_top<-tmp_out_merge[tmp_out_merge$Family%in%tmp_head,] %>% 
  group_by(region2,Family) %>% dplyr::summarise(SD = sd(Mean_Percent)) %>% 
  {.$anchor<-paste(.$region2,.$Family,sep = "_");.} %>% 
  .[order(.$region2,.$SD,decreasing = T),] %>% 
  {.$region2<-factor(.$region2,levels = c("FL","TL","PL","V1","CB","STr","Hipp","MD","Amy"));.} %>% 
  group_by(region2) %>% do(head(.,n=10))

tmp_out<-left_join(tmp_out_merge,unique(tmp_out_split[,c("anchor","P3_vs_F120","F120_vs_F90","F90_vs_F50","Signal")]),by=c("anchor"="anchor"))
save(tmp_out,tmp_out_split,tmp_out_merge,file = "result/Fig1E.Cellular_location.RData")


color16<-c("#4DBBD5B2","#fad390","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#6ab04c","#FF7256","#4682B4",
           "#8B658B","#f8c291","#00F5FF","#000000","#7E6148B2","#7FFF00","#78e08f","#6a89cc","#82ccdd","#b8e994")

tmp_color<-data.frame(Family=tmp_head,Color=color16,row.names = tmp_head)
tmp_Plist<-list()
tmp_width<-c()
for (i in 1:length(unique(tmp_top$region2))) {
  # i<-1
  if (i==1) {
    tmp_plotdf<-tmp_out %>% 
      .[.$anchor%in%tmp_top$anchor&.$region2%in%unique(tmp_top$region2)[i],] %>% 
      {.$time<-factor(.$time,levels = c("F50","F90","F120","P3"));.} %>% 
      {.$Family<-factor(.$Family,levels = unique(tmp_top[tmp_top$region2%in%unique(tmp_top$region2)[i],]$Family));.} %>% 
      {.$Signal<-factor(.$Signal,levels = c("Yes","No"));.}
    (p<-ggplot(data=tmp_plotdf, aes(x=time, y=Mean_Percent, group=Family, color=Family, )) + 
        geom_line(aes(linetype = Signal),size = 1) +
        # geom_text_repel(data=tmp_plotdf[tmp_plotdf$time%in%"P3",],aes(label = Family))+
        ylim(0,max(tmp_out[tmp_out$anchor%in%tmp_top$anchor,]$Mean_Percent))+
        facet_wrap(~region2,nrow = 1)+
        scale_color_manual(values = tmp_color[unique(tmp_top[tmp_top$region2%in%unique(tmp_top$region2)[i],]$Family),]$Color)+
        theme_bw()+
        theme(axis.title = element_text(size=18,colour = "black"),
              axis.title.x = element_blank(),
              axis.text = element_text(size=15,colour = "black"),
              legend.position = "none",
              legend.text= element_text(size=13),
              legend.title = element_text(size=13),
              strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
              panel.grid = element_blank(),
              strip.placement = "outside",
              strip.switch.pad.wrap = unit(0.3, "mm"),
              plot.title = element_blank(),
              strip.text = element_text(size =15,face = "bold",colour = "white")))
    tmp_Plist[[i]]<-p
    tmp_width<-c(tmp_width,length(unique(tmp_out[tmp_out$anchor%in%tmp_top$anchor&tmp_out$region2%in%unique(tmp_top$region2)[i],]$time)))
  }else{
    tmp_plotdf<-tmp_out %>% 
      .[.$anchor%in%tmp_top$anchor&.$region2%in%unique(tmp_top$region2)[i],] %>% 
      {.$time<-factor(.$time,levels = c("F50","F90","F120","P3"));.} %>% 
      {.$Family<-factor(.$Family,levels = intersect(tmp_head, unique(tmp_top[tmp_top$region2%in%unique(tmp_top$region2)[i],]$Family)));.} %>% 
      {.$Signal<-factor(.$Signal,levels = c("Yes","No"));.}
    (p<-ggplot(data=tmp_plotdf, aes(x=time, y=Mean_Percent, group=Family, color=Family, )) + 
        geom_line(aes(linetype = Signal),size = 1) +
        # geom_text_repel(data=tmp_plotdf[tmp_plotdf$time%in%"P3",],aes(label = Family))+
        ylim(0,max(tmp_out[tmp_out$anchor%in%tmp_top$anchor,]$Mean_Percent))+
        facet_wrap(~region2,nrow = 1)+
        scale_color_manual(values = tmp_color[intersect(tmp_head, unique(tmp_top[tmp_top$region2%in%unique(tmp_top$region2)[i],]$Family)),]$Color)+
        theme_bw()+
        theme(axis.title = element_blank(),
              axis.text = element_text(size=15,colour = "black"),
              axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
              legend.position = "none",
              legend.text= element_text(size=13),
              legend.title = element_text(size=13),
              strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
              panel.grid = element_blank(),
              strip.placement = "outside",
              strip.switch.pad.wrap = unit(0.3, "mm"),
              plot.title = element_blank(),
              strip.text = element_text(size =15,face = "bold",colour = "white")))
    tmp_Plist[[i]]<-p
    tmp_width<-c(tmp_width,length(unique(tmp_out[tmp_out$anchor%in%tmp_top$anchor&tmp_out$region2%in%unique(tmp_top$region2)[i],]$time)))
  }
}
pp<-Reduce("+",tmp_Plist)+plot_layout(ncol = length(tmp_Plist),widths = tmp_width,guides = "collect")
ggsave(filename = "figure/Fig1E.Subcellular-location_percent.line.pdf",pp,width = 20,height = 8)

tmp_out<-left_join(tmp_out_merge,unique(tmp_out_split[,c("anchor","P3_vs_F120","F120_vs_F90","F90_vs_F50","Signal")]),by=c("anchor"="anchor"))
tmp_plotdf<-tmp_out[tmp_out$anchor%in%tmp_top$anchor,] %>% 
  {.$time<-factor(.$time,levels = c("F50","F90","F120","P3"));.} %>% 
  {.$Family<-factor(.$Family,levels = intersect(tmp_head,unique(tmp_top$Family)));.} %>% 
  {.$Signal<-factor(.$Signal,levels = c("Yes","No"));.} %>% 
  {.$region2<-factor(.$region2,levels = c("FL","TL","PL","V1","CB","STr","Hipp","MD","Amy"));.}

(p<-ggplot(data=tmp_plotdf, aes(x=time, y=Mean_Percent, group=Family, color=Family, )) + 
    geom_line(aes(linetype = Signal),size = 1) +
    # geom_text_repel(data=tmp_plotdf[tmp_plotdf$time%in%"P3",],aes(label = Family))+
    ylim(0,max(tmp_out[tmp_out$anchor%in%tmp_top$anchor,]$Mean_Percent))+
    facet_wrap(~region2,nrow = 1)+
    scale_color_manual(values = tmp_color[intersect(tmp_head,unique(tmp_top$Family)),]$Color)+
    guides(linetype = guide_legend(ncol = 1,title = "Significance"))+
    guides(color = guide_legend(ncol = 1,title = "Cellular location"))+
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          # axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          # legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white")))
ggsave(filename = "figure/Fig1E.Subcellular-location_percent.line.legend.pdf",p,width = 20,height = 8)


setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Fig1E.Cellular_location.RData"))

tmp_head<-tmp_out_merge %>% {.[order(.$Mean_Percent,decreasing = T),]$Family} %>% unique() %>% head(.,n=20)
tmp_top<-tmp_out_merge[tmp_out_merge$Family%in%tmp_head,] %>% 
  group_by(region2,Family) %>% dplyr::summarise(SD = sd(Mean_Percent)) %>% 
  {.$anchor<-paste(.$region2,.$Family,sep = "_");.} %>% 
  .[order(.$region2,.$SD,decreasing = T),] %>% 
  {.$region2<-factor(.$region2,levels = c("FL","TL","PL","V1","CB","STr","Hipp","MD","Amy"));.} %>% 
  group_by(region2) %>% do(head(.,n=10))

tmp1<-tmp_out_merge %>% 
  group_by(region2,Family) %>% dplyr::summarise(SD = sd(Mean_Percent)) %>% 
  {.$anchor<-paste(.$region2,.$Family,sep = "_");.} %>% 
  .[order(.$region2,.$SD,decreasing = T),]

tmp_df<-left_join(tmp_out,tmp1[,c("anchor","SD")],by=c("anchor"="anchor")) %>% 
  subset(.,select=-anchor)
write.table(tmp_df,file = "result/Fig1E.Subcellular-location_percent.xls",quote = FALSE,sep = "\t",row.names=F, col.names=T)

# —— Fig1F Protein Family percent line ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())

protein_family<-read.delim("data/uniprot_human_family.tab", header=T,stringsAsFactors=FALSE,na.strings = "") %>% 
  .[,c("Gene.Names..primary.","Protein.families")] %>% {colnames(.)<-c("Gene_Name","Family");.} %>% na.omit()

(load("result/Proteomics.tmp_expr.tmp_meta.RData"))
tmp_out_split<-tmp_expr %>% rownames_to_column() %>% gather(Sample,Expr,-rowname) %>% 
  {.$Gene<-gsub("\\|.*","",.$rowname);.} %>% 
  left_join(.,tmp_meta,by=c("Sample"="sample")) %>% 
  left_join(.,protein_family,by=c("Gene"="Gene_Name")) %>% 
  na.omit() %>% .[.$Expr!=0,] %>% .[,c("Sample","time","region2","Family","Expr")] %>% 
  aggregate(Expr ~ Sample + time + region2 + Family, data = ., sum) %>% 
  group_by(time,region2,Sample) %>% dplyr::summarise(Sum = sum(Expr),Percent=Expr/Sum*100,Family=Family) %>% 
  as.data.frame() %>% 
  {.$time<-factor(.$time,levels = c("P3","F120","F90","F50"));.} %>% 
  .[order(.$time),]

tmp_out_split$P3_vs_F120<-NA
tmp_out_split$F120_vs_F90<-NA
tmp_out_split$F90_vs_F50<-NA

for (tmp_Region in unique(tmp_out_split$region2)) {
  # tmp_Region<-"Amy"
  tmp1<-tmp_out_split[tmp_out_split$region2%in%tmp_Region,]
  for (tmp_family in unique(tmp1$Family)) {
    # tmp_family<-"Tubulin family"
    tmp2<-tmp1[tmp1$Family%in%tmp_family,]
    if (length(unique(as.character(tmp2$time)))>1) {
      tmp_list<-unique(as.character(tmp2$time)) %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
        {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
        .[.$Versus%in%intersect(c("P3_vs_F120","F120_vs_F90","F90_vs_F50"),.$Versus),]
      for (i in 1:nrow(tmp_list)) {
        # i<-2
        if (length(tmp2[tmp2$time%in%tmp_list[i,"V1"],]$Percent)>1&length(tmp2[tmp2$time%in%tmp_list[i,"V2"],]$Percent)>1) {
          tmp_out_split[tmp_out_split$region2%in%tmp_Region&tmp_out_split$Family%in%tmp_family,] %<>% 
            {.[[paste(tmp_list[i,"V1"],"_vs_",tmp_list[i,"V2"],sep = "")]]<-t.test(as.numeric(tmp2[tmp2$time%in%tmp_list[i,"V1"],]$Percent),
                                                                                   as.numeric(tmp2[tmp2$time%in%tmp_list[i,"V2"],]$Percent))[["p.value"]];.}
        }
      }
    }
  }
}
tmp_out_split[is.na(tmp_out_split)]<-1
# tmp_out_split[tmp_out_split==1]<-NA
tmp_out_split$Signal<-NA

tmp_out_split %<>% {.$Signal<-ifelse(.$P3_vs_F120<0.05|.$F120_vs_F90<0.05|.$F90_vs_F50<0.05,"Yes","No");.}
tmp_out_split[tmp_out_split==1]<-NA
tmp_out_split$anchor<-paste(tmp_out_split$region2,tmp_out_split$Family,sep = "_")

tmp_out_merge<-tmp_expr %>% rownames_to_column() %>% gather(Sample,Expr,-rowname) %>% 
  {.$Gene<-gsub("\\|.*","",.$rowname);.} %>% 
  left_join(.,tmp_meta,by=c("Sample"="sample")) %>% 
  left_join(.,protein_family,by=c("Gene"="Gene_Name")) %>% 
  na.omit() %>% .[.$Expr!=0,] %>% .[,c("Sample","time","region2","Family","Expr")] %>% 
  aggregate(Expr ~ Sample + time + region2 + Family, data = ., sum) %>% 
  group_by(time,region2,Sample) %>% dplyr::summarise(Sum = sum(Expr),Percent=Expr/Sum*100,Family=Family) %>% 
  as.data.frame() %>% 
  group_by(time,region2,Family) %>% dplyr::summarise(Mean_Percent = mean(Percent)) %>% 
  .[order(.$time,.$region2,.$Mean_Percent,decreasing = T),] %>% 
  {.$anchor<-paste(.$region2,.$Family,sep = "_");.}

tmp_head<-tmp_out_merge %>% {.[order(.$Mean_Percent,decreasing = T),]$Family} %>% unique() %>% head(.,n=20)

tmp_top<-tmp_out_merge[tmp_out_merge$Family%in%tmp_head,] %>% 
  group_by(region2,Family) %>% dplyr::summarise(SD = sd(Mean_Percent)) %>% 
  {.$anchor<-paste(.$region2,.$Family,sep = "_");.} %>% 
  .[order(.$region2,.$SD,decreasing = T),] %>% 
  {.$region2<-factor(.$region2,levels = c("FL","TL","PL","V1","CB","STr","Hipp","MD","Amy"));.} %>% 
  group_by(region2) %>% do(head(.,n=10))

tmp_out<-left_join(tmp_out_merge,unique(tmp_out_split[,c("anchor","P3_vs_F120","F120_vs_F90","F90_vs_F50","Signal")]),by=c("anchor"="anchor"))
save(tmp_out,tmp_out_split,tmp_out_merge,file = "result/Fig1F.Protein-Family_percent.RData")

color16<-c("#4DBBD5B2","#fad390","#9196f7","#DC0000B2","#00A087B2",
           "#40407a","#FF82AB","#6ab04c","#FF7256","#4682B4",
           "#8B658B","#f8c291","#00F5FF","#000000","#7E6148B2","#7FFF00","#78e08f","#6a89cc","#82ccdd","#b8e994")

tmp_color<-data.frame(Family=tmp_head,Color=color16,row.names = tmp_head)
tmp_Plist<-list()
tmp_width<-c()
for (i in 1:length(unique(tmp_top$region2))) {
  # i<-1
  if (i==1) {
    tmp_plotdf<-tmp_out %>% 
      .[.$anchor%in%tmp_top$anchor&.$region2%in%unique(tmp_top$region2)[i],] %>% 
      {.$time<-factor(.$time,levels = c("F50","F90","F120","P3"));.} %>% 
      {.$Family<-factor(.$Family,levels = unique(tmp_top[tmp_top$region2%in%unique(tmp_top$region2)[i],]$Family));.} %>% 
      {.$Signal<-factor(.$Signal,levels = c("Yes","No"));.}
    (p<-ggplot(data=tmp_plotdf, aes(x=time, y=Mean_Percent, group=Family, color=Family, )) + 
        geom_line(aes(linetype = Signal),size = 1) +
        # geom_text_repel(data=tmp_plotdf[tmp_plotdf$time%in%"P3",],aes(label = Family))+
        ylim(0,max(tmp_out[tmp_out$anchor%in%tmp_top$anchor,]$Mean_Percent))+
        facet_wrap(~region2,nrow = 1)+
        scale_color_manual(values = tmp_color[unique(tmp_top[tmp_top$region2%in%unique(tmp_top$region2)[i],]$Family),]$Color)+
        theme_bw()+
        theme(axis.title = element_text(size=18,colour = "black"),
              axis.title.x = element_blank(),
              axis.text = element_text(size=15,colour = "black"),
              legend.position = "none",
              legend.text= element_text(size=13),
              legend.title = element_text(size=13),
              strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
              panel.grid = element_blank(),
              strip.placement = "outside",
              strip.switch.pad.wrap = unit(0.3, "mm"),
              plot.title = element_blank(),
              strip.text = element_text(size =15,face = "bold",colour = "white")))
    tmp_Plist[[i]]<-p
    tmp_width<-c(tmp_width,length(unique(tmp_out[tmp_out$anchor%in%tmp_top$anchor&tmp_out$region2%in%unique(tmp_top$region2)[i],]$time)))
  }else{
    tmp_plotdf<-tmp_out %>% 
      .[.$anchor%in%tmp_top$anchor&.$region2%in%unique(tmp_top$region2)[i],] %>% 
      {.$time<-factor(.$time,levels = c("F50","F90","F120","P3"));.} %>% 
      {.$Family<-factor(.$Family,levels = intersect(tmp_head, unique(tmp_top[tmp_top$region2%in%unique(tmp_top$region2)[i],]$Family)));.} %>% 
      {.$Signal<-factor(.$Signal,levels = c("Yes","No"));.}
    (p<-ggplot(data=tmp_plotdf, aes(x=time, y=Mean_Percent, group=Family, color=Family, )) + 
        geom_line(aes(linetype = Signal),size = 1) +
        # geom_text_repel(data=tmp_plotdf[tmp_plotdf$time%in%"P3",],aes(label = Family))+
        ylim(0,max(tmp_out[tmp_out$anchor%in%tmp_top$anchor,]$Mean_Percent))+
        facet_wrap(~region2,nrow = 1)+
        scale_color_manual(values = tmp_color[intersect(tmp_head, unique(tmp_top[tmp_top$region2%in%unique(tmp_top$region2)[i],]$Family)),]$Color)+
        theme_bw()+
        theme(axis.title = element_blank(),
              axis.text = element_text(size=15,colour = "black"),
              axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
              legend.position = "none",
              legend.text= element_text(size=13),
              legend.title = element_text(size=13),
              strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
              panel.grid = element_blank(),
              strip.placement = "outside",
              strip.switch.pad.wrap = unit(0.3, "mm"),
              plot.title = element_blank(),
              strip.text = element_text(size =15,face = "bold",colour = "white")))
    tmp_Plist[[i]]<-p
    tmp_width<-c(tmp_width,length(unique(tmp_out[tmp_out$anchor%in%tmp_top$anchor&tmp_out$region2%in%unique(tmp_top$region2)[i],]$time)))
  }
}
pp<-Reduce("+",tmp_Plist)+plot_layout(ncol = length(tmp_Plist),widths = tmp_width,guides = "collect")
ggsave(filename = "figure/Fig1F.Protein-Family_percent.line.pdf",pp,width = 20,height = 8)

tmp_out<-left_join(tmp_out_merge,unique(tmp_out_split[,c("anchor","P3_vs_F120","F120_vs_F90","F90_vs_F50","Signal")]),by=c("anchor"="anchor"))
tmp_plotdf<-tmp_out[tmp_out$anchor%in%tmp_top$anchor,] %>% 
  {.$time<-factor(.$time,levels = c("F50","F90","F120","P3"));.} %>% 
  {.$Family<-factor(.$Family,levels = intersect(tmp_head,unique(tmp_top$Family)));.} %>% 
  {.$Signal<-factor(.$Signal,levels = c("Yes","No"));.} %>% 
  {.$region2<-factor(.$region2,levels = c("FL","TL","PL","V1","CB","STr","Hipp","MD","Amy"));.}

(p<-ggplot(data=tmp_plotdf, aes(x=time, y=Mean_Percent, group=Family, color=Family, )) + 
    geom_line(aes(linetype = Signal),size = 1) +
    # geom_text_repel(data=tmp_plotdf[tmp_plotdf$time%in%"P3",],aes(label = Family))+
    ylim(0,max(tmp_out[tmp_out$anchor%in%tmp_top$anchor,]$Mean_Percent))+
    facet_wrap(~region2,nrow = 1)+
    scale_color_manual(values = tmp_color[intersect(tmp_head,unique(tmp_top$Family)),]$Color)+
    guides(linetype = guide_legend(ncol = 1,title = "Significance"))+
    guides(color = guide_legend(ncol = 1,title = "Protein Family"))+
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          # axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          # legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white")))
ggsave(filename = "figure/Fig1F.Protein-Family_percent.line.legend.pdf",p,width = 20,height = 8)


setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Fig1F.Protein-Family_percent.RData"))

tmp_head<-tmp_out_merge %>% {.[order(.$Mean_Percent,decreasing = T),]$Family} %>% unique() %>% head(.,n=20)
tmp_top<-tmp_out_merge[tmp_out_merge$Family%in%tmp_head,] %>% 
  group_by(region2,Family) %>% dplyr::summarise(SD = sd(Mean_Percent)) %>% 
  {.$anchor<-paste(.$region2,.$Family,sep = "_");.} %>% 
  .[order(.$region2,.$SD,decreasing = T),] %>% 
  {.$region2<-factor(.$region2,levels = c("FL","TL","PL","V1","CB","STr","Hipp","MD","Amy"));.} %>% 
  group_by(region2) %>% do(head(.,n=10))

tmp1<-tmp_out_merge %>% 
  group_by(region2,Family) %>% dplyr::summarise(SD = sd(Mean_Percent)) %>% 
  {.$anchor<-paste(.$region2,.$Family,sep = "_");.} %>% 
  .[order(.$region2,.$SD,decreasing = T),]

tmp_df<-left_join(tmp_out,tmp1[,c("anchor","SD")],by=c("anchor"="anchor")) %>% 
  subset(.,select=-anchor)
write.table(tmp_df,file = "result/Fig1F.Protein-Family_percent.line.xls",quote = FALSE,sep = "\t",row.names=F, col.names=T)








# — Figure2 ----------------------------------------------------------
# —— Fig2B Time Corr heatmap ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))

tmp_cor<-tmp_expr %>% rownames_to_column() %>% gather(Sample,Expr,-rowname) %>% 
  left_join(.,tmp_meta,by=c("Sample"="sample")) %>% 
  aggregate(Expr ~ rowname + time, data = ., sum) %>% 
  spread(.,key = "time",value = "Expr") %>% {rownames(.)<-.$rowname;.} %>% 
  subset(.,select=-rowname) %>% .[,c("F50","F90","F120","P3")] %>% 
  .[rowSums(.)!=0,] %>% 
  {log2(.+1)} %>%
  cor()

tmp_corP<-tmp_expr %>% rownames_to_column() %>% gather(Sample,Expr,-rowname) %>% 
  left_join(.,tmp_meta,by=c("Sample"="sample")) %>% 
  aggregate(Expr ~ rowname + time, data = ., sum) %>% 
  spread(.,key = "time",value = "Expr") %>% {rownames(.)<-.$rowname;.} %>% 
  subset(.,select=-rowname) %>% .[,c("F50","F90","F120","P3")] %>% 
  .[rowSums(.)!=0,] %>% 
  {log2(.+1)} %>%
  cor.mtest()

col2 = colorRampPalette(rev(c('#67001F', '#B2182B', '#D6604D', '#F4A582',
                              '#FDDBC7', '#FFFFFF', '#D1E5F0', '#92C5DE',
                              '#4393C3', '#2166AC', '#053061')))

pdf("figure/Fig2B.cor_heatmap.pdf",width = 10,height = 8)
corrplot(tmp_cor,
         p.mat = tmp_corP$p,
         method ="color",
         sig.level = c(0.001,0.01,0.05),insig = "label_sig",pch.cex = 1.5,
         addCoef.col = "white",col.lim = c(0,1),
         col=col2(50))
dev.off()


# —— Fig2D Inter-regional differences line ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))
# tmp_expr<-log10(tmp_expr+1)

# CB
tmp_title="CB"

tmp_out<-data.frame()
for (tmp_time in unique(tmp_meta[tmp_meta$region2%in%"CB",]$time)) {
  # tmp_time<-"F50"
  tmp_meta_sub<-tmp_meta[tmp_meta$region2%in%"CB",] %>% .[.$time%in%tmp_time,]
  tmp_expr_sub<-tmp_expr[,tmp_meta_sub$sample] %>% .[rowSums(.)!=0,] %>% {.<-./colSums(.)*10000;.}
  # tmp_expr_sub<-tmp_expr[,tmp_meta_sub$sample] %>% .[rowSums(.)!=0,] %>%
  #   rownames_to_column() %>% gather(Var2,Freq,-rowname) %>% 
  #   group_by(Var2) %>% dplyr::summarise(Sum = sum(Freq),rowname=rowname,Freq=Freq) %>% 
  #   {.$Percent<-.$Freq/.$Sum*10000;.} %>% .[,c("rowname","Var2","Percent")] %>% 
  #   spread(.,key = "Var2",value = "Percent") %>% as.data.frame() %>% 
  #   {rownames(.)<-.$rowname;.} %>% subset(.,select=-rowname)
  tmp_list<-unique(tmp_meta_sub$sample) %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
    {colnames(.)<-c("Case","Control");.} %>% {.$Versus<-paste(.$Case,"_vs_",.$Control,sep = "");.}
  for (n in 1:nrow(tmp_list)) {
    # n<-1
    tmp_expr_sub[[tmp_list[n,"Versus"]]]<-abs(tmp_expr_sub[[tmp_list[n,"Case"]]]-tmp_expr_sub[[tmp_list[n,"Control"]]])
  }
  
  tmp_Mean<-tmp_expr_sub[,grep("*_vs_*",colnames(tmp_expr_sub),value = T)] %>% colMeans() %>% mean()
  tmp_split<-tmp_expr_sub[,grep("*_vs_*",colnames(tmp_expr_sub),value = T)] %>% colMeans() %>% 
    as.data.frame() %>% .[,1]
  tmp_out<-data.frame(Time=tmp_time,Region=tmp_title,Mean_diff=tmp_Mean,Each_sample=tmp_split) %>% rbind(tmp_out,.)
}

tmp_out$P3_vs_F120<-NA
tmp_out$F120_vs_F90<-NA
tmp_out$F90_vs_F50<-NA
if (length(unique(as.character(tmp_out$Time)))>1) {
  tmp_list<-unique(as.character(tmp_out$Time)) %>% rev() %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
    {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
    .[.$Versus%in%intersect(c("P3_vs_F120","F120_vs_F90","F90_vs_F50"),.$Versus),]
  for (i in 1:nrow(tmp_list)) {
    # i<-1
    if (length(tmp_out[tmp_out$Time%in%tmp_list[i,"V1"],]$Each_sample)>1&length(tmp_out[tmp_out$Time%in%tmp_list[i,"V2"],]$Each_sample)>1) {
      tmp_out %<>% 
        {.[[paste(tmp_list[i,"V1"],"_vs_",tmp_list[i,"V2"],sep = "")]]<-t.test(as.numeric(tmp_out[tmp_out$Time%in%tmp_list[i,"V1"],]$Each_sample),
                                                                               as.numeric(tmp_out[tmp_out$Time%in%tmp_list[i,"V2"],]$Each_sample))[["p.value"]];.}
    }
  }
}

tmp_part1<-tmp_out %>% {.$Signal<-ifelse(.$P3_vs_F120<0.05|.$F120_vs_F90<0.05|.$F90_vs_F50<0.05,"Yes","No");.}

# CTX
rm(list=setdiff(ls(),"tmp_part1"))
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))
# tmp_expr<-log10(tmp_expr+1)

tmp_title="CTX"

tmp_out<-data.frame()
for (tmp_time in unique(tmp_meta[tmp_meta$region2%in%c("FL","TL","PL","V1"),]$time)) {
  # tmp_time<-"F50"
  tmp_meta_sub<-tmp_meta[tmp_meta$region2%in%c("FL","TL","PL","V1"),] %>% .[.$time%in%tmp_time,]
  tmp_expr_sub<-tmp_expr[,tmp_meta_sub$sample] %>% .[rowSums(.)!=0,] %>% {.<-./colSums(.)*10000;.}
  
  # tmp_expr_sub<-tmp_expr[,tmp_meta_sub$sample] %>% .[rowSums(.)!=0,] %>%
  #   rownames_to_column() %>% gather(Var2,Freq,-rowname) %>% 
  #   group_by(Var2) %>% dplyr::summarise(Sum = sum(Freq),rowname=rowname,Freq=Freq) %>% 
  #   {.$Percent<-.$Freq/.$Sum*10000;.} %>% .[,c("rowname","Var2","Percent")] %>% 
  #   spread(.,key = "Var2",value = "Percent") %>% as.data.frame() %>% 
  #   {rownames(.)<-.$rowname;.} %>% subset(.,select=-rowname)
  
  tmp_list<-unique(tmp_meta_sub$sample) %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
    {colnames(.)<-c("Case","Control");.} %>% {.$Versus<-paste(.$Case,"_vs_",.$Control,sep = "");.} %>% 
    left_join(.,tmp_meta_sub[,c("sample","region2")],by=c("Case"="sample")) %>% {names(.)[names(.) == 'region2'] <- 'Case_Region';.} %>% 
    left_join(.,tmp_meta_sub[,c("sample","region2")],by=c("Control"="sample")) %>% {names(.)[names(.) == 'region2'] <- 'Control_Region';.} %>% 
    .[.$Case_Region!=.$Control_Region,]
  for (n in 1:nrow(tmp_list)) {
    # n<-1
    tmp_expr_sub[[tmp_list[n,"Versus"]]]<-abs(tmp_expr_sub[[tmp_list[n,"Case"]]]-tmp_expr_sub[[tmp_list[n,"Control"]]])
  }
  tmp_Mean<-tmp_expr_sub[,grep("*_vs_*",colnames(tmp_expr_sub),value = T)] %>% colMeans() %>% mean()
  tmp_split<-tmp_expr_sub[,grep("*_vs_*",colnames(tmp_expr_sub),value = T)] %>% colMeans() %>% 
    as.data.frame() %>% .[,1]
  tmp_out<-data.frame(Time=tmp_time,Region=tmp_title,Mean_diff=tmp_Mean,Each_sample=tmp_split) %>% rbind(tmp_out,.)
}

tmp_out$P3_vs_F120<-NA
tmp_out$F120_vs_F90<-NA
tmp_out$F90_vs_F50<-NA
if (length(unique(as.character(tmp_out$Time)))>1) {
  tmp_list<-unique(as.character(tmp_out$Time)) %>% rev() %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
    {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
    .[.$Versus%in%intersect(c("P3_vs_F120","F120_vs_F90","F90_vs_F50"),.$Versus),]
  for (i in 1:nrow(tmp_list)) {
    # i<-1
    if (length(tmp_out[tmp_out$Time%in%tmp_list[i,"V1"],]$Each_sample)>1&length(tmp_out[tmp_out$Time%in%tmp_list[i,"V2"],]$Each_sample)>1) {
      tmp_out %<>% 
        {.[[paste(tmp_list[i,"V1"],"_vs_",tmp_list[i,"V2"],sep = "")]]<-t.test(as.numeric(tmp_out[tmp_out$Time%in%tmp_list[i,"V1"],]$Each_sample),
                                                                               as.numeric(tmp_out[tmp_out$Time%in%tmp_list[i,"V2"],]$Each_sample))[["p.value"]];.}
    }
  }
}

tmp_part2<-tmp_out %>% {.$Signal<-ifelse(.$P3_vs_F120<0.05|.$F120_vs_F90<0.05|.$F90_vs_F50<0.05,"Yes","No");.}

# sCTX
rm(list=setdiff(ls(),c("tmp_part1","tmp_part2")))
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))
# tmp_expr<-log10(tmp_expr+1)

tmp_title="sCTX"

tmp_out<-data.frame()
for (tmp_time in unique(tmp_meta[tmp_meta$region2%in%c("STr","Hipp","MD","Amy"),]$time)) {
  # tmp_time<-"F90"
  tmp_meta_sub<-tmp_meta[tmp_meta$region2%in%c("STr","Hipp","MD","Amy"),] %>% .[.$time%in%tmp_time,]
  tmp_expr_sub<-tmp_expr[,tmp_meta_sub$sample] %>% .[rowSums(.)!=0,] %>% {.<-./colSums(.)*10000;.}
  
  # tmp_expr_sub<-tmp_expr[,tmp_meta_sub$sample] %>% .[rowSums(.)!=0,] %>%
  #   rownames_to_column() %>% gather(Var2,Freq,-rowname) %>% 
  #   group_by(Var2) %>% dplyr::summarise(Sum = sum(Freq),rowname=rowname,Freq=Freq) %>% 
  #   {.$Percent<-.$Freq/.$Sum*10000;.} %>% .[,c("rowname","Var2","Percent")] %>% 
  #   spread(.,key = "Var2",value = "Percent") %>% as.data.frame() %>% 
  #   {rownames(.)<-.$rowname;.} %>% subset(.,select=-rowname)
  
  tmp_list<-unique(tmp_meta_sub$sample) %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
    {colnames(.)<-c("Case","Control");.} %>% {.$Versus<-paste(.$Case,"_vs_",.$Control,sep = "");.} %>% 
    left_join(.,tmp_meta_sub[,c("sample","region2")],by=c("Case"="sample")) %>% {names(.)[names(.) == 'region2'] <- 'Case_Region';.} %>% 
    left_join(.,tmp_meta_sub[,c("sample","region2")],by=c("Control"="sample")) %>% {names(.)[names(.) == 'region2'] <- 'Control_Region';.} %>% 
    .[.$Case_Region!=.$Control_Region,]
  for (n in 1:nrow(tmp_list)) {
    # n<-1
    tmp_expr_sub[[tmp_list[n,"Versus"]]]<-abs(tmp_expr_sub[[tmp_list[n,"Case"]]]-tmp_expr_sub[[tmp_list[n,"Control"]]])
  }
  tmp_Mean<-tmp_expr_sub[,grep("*_vs_*",colnames(tmp_expr_sub),value = T)] %>% colMeans() %>% mean()
  tmp_split<-tmp_expr_sub[,grep("*_vs_*",colnames(tmp_expr_sub),value = T)] %>% colMeans() %>% 
    as.data.frame() %>% .[,1]
  tmp_out<-data.frame(Time=tmp_time,Region=tmp_title,Mean_diff=tmp_Mean,Each_sample=tmp_split) %>% rbind(tmp_out,.)
}

tmp_out$P3_vs_F120<-NA
tmp_out$F120_vs_F90<-NA
tmp_out$F90_vs_F50<-NA
if (length(unique(as.character(tmp_out$Time)))>1) {
  tmp_list<-unique(as.character(tmp_out$Time)) %>% rev() %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
    {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
    .[.$Versus%in%intersect(c("P3_vs_F120","F120_vs_F90","F90_vs_F50"),.$Versus),]
  for (i in 1:nrow(tmp_list)) {
    # i<-1
    if (length(tmp_out[tmp_out$Time%in%tmp_list[i,"V1"],]$Each_sample)>1&length(tmp_out[tmp_out$Time%in%tmp_list[i,"V2"],]$Each_sample)>1) {
      tmp_out %<>% 
        {.[[paste(tmp_list[i,"V1"],"_vs_",tmp_list[i,"V2"],sep = "")]]<-t.test(as.numeric(tmp_out[tmp_out$Time%in%tmp_list[i,"V1"],]$Each_sample),
                                                                               as.numeric(tmp_out[tmp_out$Time%in%tmp_list[i,"V2"],]$Each_sample))[["p.value"]];.}
    }
  }
}

tmp_part3<-tmp_out %>% {.$Signal<-ifelse(.$P3_vs_F120<0.05|.$F120_vs_F90<0.05|.$F90_vs_F50<0.05,"Yes","No");.}


tmp_merge<-rbind(tmp_part1,tmp_part2,tmp_part3)
save(tmp_merge,file = "result/Fig2D.Inter-regional.differences.RData")

load(file = "result/Fig2D.Inter-regional.differences.RData")

(p<-tmp_merge[,c("Time","Region","Mean_diff","Signal")] %>% unique() %>% 
    {.$Time<-factor(.$Time,levels = c("F50","F90","F120","P3"));.} %>% 
    {.$Signal<-factor(.$Signal,levels = c("Yes","No"));.} %>%
    ggplot(data=., aes(x=Time, y=Mean_diff, group=Region, color=Region, )) + 
    geom_line(aes(linetype = Signal),size = 1) +
    geom_point()+
    # geom_text_repel(data=tmp_plotdf[tmp_plotdf$time%in%"P3",],aes(label = Family))+
    # ylim(0,max(tmp_out[tmp_out$anchor%in%tmp_top$anchor,]$Mean_Percent))+
    # facet_wrap(~region2,nrow = 1)+
    guides(linetype = guide_legend(ncol = 1,title = "Significance"))+
    scale_color_manual(values = c("#f1c40f","#2980b9","#c0392b"))+
    ylab(label = "Mean difference")+
    theme_bw()+
    theme(axis.title = element_text(size=18,colour = "black"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          # legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white")))
ggsave(filename = "figure/Fig2D.Inter-regional_differences.line.pdf",p,width = 8,height = 6)


# — Figure3 ----------------------------------------------------------
# —— Fig3B Proteomics Distance by Time Heatmap ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))

tmp_own_meta <- read.delim("data/metadata.txt") %>% 
  {.[.$time%in%"A",]$time<-"F50";.} %>% 
  {.[.$time%in%"B",]$time<-"F90";.} %>% 
  {.[.$time%in%"C",]$time<-"F120";.} %>% 
  {.[.$time%in%"D",]$time<-"P3";.} %>% 
  {.[.$region2%in%"1",]$region2<-"FL";.} %>% 
  {.[.$region2%in%"2",]$region2<-"TL";.} %>% 
  {.[.$region2%in%"3",]$region2<-"PL";.} %>% 
  {.[.$region2%in%"4",]$region2<-"V1";.} %>% 
  {.[.$region2%in%"5",]$region2<-"CB";.} %>% 
  {.[.$region2%in%"6",]$region2<-"STr";.} %>% 
  {.[.$region2%in%"7",]$region2<-"Hipp";.} %>% 
  {.[.$region2%in%"8",]$region2<-"MD";.} %>% 
  {.[.$region2%in%"9",]$region2<-"Amy";.}


Time_sort<-c("FL","TL","PL","V1","CB","STr","Hipp","MD","Amy")

tmp_own_meta$Minor_Region<-factor(tmp_own_meta$region2,levels = Time_sort)

tmp_input<-data.frame()
for (tmp_time in unique(tmp_meta$time)) {
  # tmp_time<-"F50"
  tmp_samples<-tmp_meta[tmp_meta$time%in%tmp_time,]$sample
  tmp_input<-tmp_expr[,tmp_samples] %>% 
    .[rowSums(.)!=0,] %>% 
    {log2(.+1)} %>% 
    t() %>% dist() %>% as.matrix() %>% as.data.frame() %>% 
    rownames_to_column() %>% gather(rowname2,Corr,-rowname) %>% 
    left_join(.,tmp_own_meta[tmp_own_meta$sample%in%tmp_samples,c("sample","Minor_Region")],by=c("rowname"="sample")) %>% 
    {names(.)[names(.) == 'Minor_Region'] <- 'Sample1';.} %>% 
    left_join(.,tmp_own_meta[tmp_own_meta$sample%in%tmp_samples,c("sample","Minor_Region")],by=c("rowname2"="sample")) %>% 
    {names(.)[names(.) == 'Minor_Region'] <- 'Sample2';.} %>% 
    .[.$Sample1%in%"CB",] %>% {.$Versus<-paste(.$Sample1,"_",.$Sample2,sep = "");.} %>% 
    {.$Time<-tmp_time;.} %>% 
    group_by(Sample1,Sample2) %>% dplyr::summarise(Corr_mean = mean(Corr),Corr=Corr,
                                                   Versus=Versus,Time=Time,) %>% as.data.frame() %>% 
    rbind(tmp_input,.)
}

tmp_input$P3_vs_F120<-NA
tmp_input$F120_vs_F90<-NA
tmp_input$F90_vs_F50<-NA
for (tmp_Versus in unique(tmp_input$Versus)) {
  # tmp_Versus<-"CB_FL"
  tmp1<-tmp_input[tmp_input$Versus%in%tmp_Versus,]
  if (length(unique(as.character(tmp1$Time)))>1) {
    tmp_list<-unique(as.character(tmp1$Time)) %>% rev() %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
      {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
      .[.$Versus%in%intersect(c("P3_vs_F120","F120_vs_F90","F90_vs_F50"),.$Versus),]
    for (i in 1:nrow(tmp_list)) {
      # i<-1
      if (length(tmp1[tmp1$Time%in%tmp_list[i,"V1"],]$Corr)>1&length(tmp1[tmp1$Time%in%tmp_list[i,"V2"],]$Corr)>1) {
        tmp_input[tmp_input$Versus%in%tmp_Versus,] %<>% 
          {.[[tmp_list[i,"Versus"]]]<-t.test(as.numeric(tmp1[tmp1$Time%in%tmp_list[i,"V1"],]$Corr),
                                             as.numeric(tmp1[tmp1$Time%in%tmp_list[i,"V2"],]$Corr))[["p.value"]];.}
      }
    }
  }
}

tmp_input %<>% {.$Signal<-ifelse(.$P3_vs_F120<0.05|.$F120_vs_F90<0.05|.$F90_vs_F50<0.05,"Yes","No");.}
save(tmp_input,file = "result/Fig3B.Proteomics_Distance.by_Time.line.RData")


(p<-tmp_input[,c("Sample2","Time","Corr_mean")] %>% unique() %>% spread(.,key = "Time",value = "Corr_mean") %>% 
    {rownames(.)<-.$Sample2;.} %>% subset(.,select=-Sample2) %>% 
    .[c("FL","TL","PL","V1","STr","Hipp","Amy","MD"),c("F50","F90","F120","P3")] %>% 
    pheatmap::pheatmap(.,cluster_cols = F,cluster_rows = F,display_numbers = T,number_color = "black",fontsize = 15,
                       color = colorRampPalette(rev(c('#B2182B', '#D6604D', '#F4A582',
                                                      '#FDDBC7')))(50)) %>% as.ggplot())

ggsave(filename = "figure/Fig3B.Proteomics_Distance.by_Time.Heatmap.pdf",p,width = 6,height = 6)






# —— Fig3C CB-specific Marker heatmap ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))

tmp_expr<-tmp_expr %>% rownames_to_column() %>% {.$Gene_name<-gsub("\\|.*","",.$rowname);.} %>% 
  .[,c("Gene_name",tmp_meta$sample)] %>% as.data.table() %>% 
  .[,lapply(.SD,mean),"Gene_name"] %>% as.data.frame() %>% {rownames(.)<-.$Gene_name;.} %>% subset(.,select=-Gene_name)


myobj<- CreateSeuratObject(counts = tmp_expr)
myobj<- tmp_meta %>% {rownames(.)<-.$sample;.} %>% AddMetaData(myobj, .)

seurat.combined <- NormalizeData(myobj)
seurat.combined <- FindVariableFeatures(seurat.combined, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seurat.combined)
seurat.combined <- ScaleData(seurat.combined, features = all.genes)
seurat.combined <- RunPCA(seurat.combined, features = VariableFeatures(object = seurat.combined), npcs = 30, verbose = T)
seurat.combined <- FindNeighbors(seurat.combined, dims = 1:30)
seurat.combined <- FindClusters(seurat.combined,resolution = seq(0.1,1,0.1))
seurat.combined <- RunUMAP(seurat.combined, seed.use = 42, reduction = "pca", dims = 1:30)
seurat.combined <- RunTSNE(seurat.combined, seed.use = 2, reduction = "pca", dims = 1:30)

tmp_Marker<-data.frame()
for (tmp_time in unique(myobj$time)) {
  # tmp_time<-"F50"
  subObj<-subset(seurat.combined,time%in%tmp_time) %>% 
    {.$region2<-factor(.$region2,levels = intersect(c("CB","FL","TL","PL","V1","STr","Hipp","MD","Amy"),unique(.$region2)));.}
  
  Idents(subObj)<-"region2"
  all.markers <- FindAllMarkers(subObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>% 
    .[.$p_val<0.05&.$cluster%in%"CB",]
  tmp_Marker<- all.markers %>% {.$Time<-tmp_time;.} %>% rbind(tmp_Marker,.)
  (p<-DoHeatmap(subObj, features = unique(head(all.markers$gene,n=20)),group.by = "region2",label = 10,
                raster = F,draw.lines=F,group.colors = c("#dece62","#4aa751","#35703e","#287374","#1c6eb7","#f58737","#e05a44","#b92645","#7e2243"))+
      scale_fill_gradientn(colors = colorRampPalette(c("#1a2a6c", "white", "#c21e20"))(20),na.value = "white"))
  
  
  ggsave(paste("figure/Fig3C.",tmp_time,".CB-specific_Marker.heatmap.pdf",sep = ""),p,width = 8,height = 6)
  
}
write.table(tmp_Marker,file = "result/Fig3C.CB-specific.Marker.xls",quote = FALSE,sep = "\t",row.names = F)


# —— Fig3DE Granule&Purkinje Change by Time in cortical point add Line ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
ortholog <- read.delim("data/Homo_MFA_MMU_Pig_Mouse.Gid_Gname.xls", header=T,stringsAsFactors=FALSE,na.strings = "")

(load("result/Proteomics.tmp_expr.tmp_meta.RData"))
tmp_meta<-tmp_meta[tmp_meta$region2%in%"CB",]

tmp_expr<-tmp_expr %>% rownames_to_column() %>% {.$Gene_name<-gsub("\\|.*","",.$rowname);.} %>% 
  .[,c("Gene_name",tmp_meta$sample)] %>% as.data.table() %>% 
  .[,lapply(.SD,mean),"Gene_name"] %>% as.data.frame() %>% 
  {rownames(.)<-.$Gene_name;.} %>% subset(.,select=-Gene_name) %>% {log10(.+1)}

# disease_Gene<-data.frame(Gene=c("APC", "BRCA2", "PALB2", "PTCH1", "SUFU", "TP53",
#                                 "CDH15", "CALB2", "RBFOX3", "RELN", "ZIC1", "GABRA6",
#                                 "PCP2", "PVALB", "CABLES1", "ITPR1", "GRID2IP"),
#                          Type=c(rep("Medulloblastomas",6),
#                                 rep("Granule",6),
#                                 rep("Purkinje",5)))

disease_Gene<-data.frame(Gene=c("CDH15", "CALB2", "RBFOX3", "RELN", "ZIC1", "GABRA6",
                                "PCP2", "PVALB", "CABLES1", "ITPR1", "GRID2IP"),
                         Type=c(rep("Granule",6),
                                rep("Purkinje",5)))


setdiff(disease_Gene$Gene,rownames(tmp_expr))
# "BRCA2"   "PALB2"   "PTCH1"   "SUFU"    "TP53"    "RBFOX3"  "CABLES1"
disease_Gene<-disease_Gene[disease_Gene$Gene%in%intersect(disease_Gene$Gene,rownames(tmp_expr)),]


tmp_out<-data.frame()
for (tmp_Region in c("CB")) {
  # tmp_Region<-"CB"
  for (tmp_Time in unique(tmp_meta$time)) {
    # tmp_Time<-"F50"
    tmp_Sample<-tmp_meta[tmp_meta$region2%in%tmp_Region&tmp_meta$time%in%tmp_Time,]
    tmp_out<-tmp_expr[unique(disease_Gene$Gene),tmp_Sample$sample] %>% rownames_to_column() %>% gather(Samples,Expr,-rowname) %>% 
      left_join(.,disease_Gene,by=c("rowname"="Gene")) %>% 
      group_by(rowname) %>% dplyr::summarise(mean_Expr = mean(Expr),
                                             rowname=rowname,Expr=Expr,
                                             Samples=Samples,Type=Type) %>% as.data.frame() %>%
      left_join(.,tmp_meta,by=c("Samples"="sample")) %>% rbind(tmp_out,.)
  }
}

tmp_out$P3_vs_F120<-NA
tmp_out$F120_vs_F90<-NA
tmp_out$F90_vs_F50<-NA
for (tmp_Region in c("CB")) {
  # tmp_Region<-"FL"
  for (tmp_Type in unique(tmp_out$Type)) {
    # tmp_Type<-"AD"
    if (length(unique(as.character(tmp_out[tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$time)))>1) {
      tmp_list<-unique(as.character(tmp_out[tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$time)) %>% rev() %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
        {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
        .[.$Versus%in%intersect(c("P3_vs_F120","F120_vs_F90","F90_vs_F50"),.$Versus),]
      for (i in 1:nrow(tmp_list)) {
        # i<-1
        if (length(tmp_out[tmp_out$time%in%tmp_list[i,"V1"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$Expr)>1&
            length(tmp_out[tmp_out$time%in%tmp_list[i,"V2"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$Expr)>1) {
          tmp_out[tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,] %<>% 
            {.[[tmp_list[i,"Versus"]]]<-t.test(as.numeric(tmp_out[tmp_out$time%in%tmp_list[i,"V1"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$Expr),
                                               as.numeric(tmp_out[tmp_out$time%in%tmp_list[i,"V2"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$Expr))[["p.value"]];.}
        }
      }
    }
  }
}

tmp_plotdf<-tmp_out %>% {.$Signal<-ifelse(.$P3_vs_F120<0.05|.$F120_vs_F90<0.05|.$F90_vs_F50<0.05,"Yes","No");.}
# save(tmp_plotdf,file = "result/Fig5D.Disease_Change_by_Time_in_cortical.v7.RData")
# write.table(tmp_plotdf[order(tmp_plotdf$region2,tmp_plotdf$Type,tmp_plotdf$time),],
#             file = "result/Fig3DE.Granule_Purkinje.Change_by_Time_in_CB.v1.xls",
#             quote = FALSE,sep = "\t",row.names=F, col.names=T)

pd <- position_dodge(0.5)

# load("result/Fig5D.Disease_Change_by_Time_in_cortical.v7.RData")
tmp_plotdf$time<-factor(tmp_plotdf$time,levels = c("F50","F90","F120","P3"))
tmp_plotdf$Type<-factor(tmp_plotdf$Type,levels = c("Medulloblastomas","Granule","Purkinje"))
tmp_plotdf$Signal<-factor(tmp_plotdf$Signal,levels = c("Yes","No"))

tmp_plotdf<-tmp_plotdf %>% group_by(time,rowname) %>% 
  dplyr::summarise(All_mean = mean(Expr),All_sd = sd(Expr),
                   Expr=Expr,time=time,Type=Type,region2=region2,Signal=Signal)

tmp_plotdf$rowname<-factor(tmp_plotdf$rowname,levels = disease_Gene$Gene)

color14<-c("#633221","#B8473D","#D59681","#BF7D40","#D1B175",
           "#3D81B8","#9FDFD6","#B6CFE7","#BDBAE8","#1C1D54",
           "#174527","#4B6B24","#47541C","#53C671")
for (i in unique(tmp_plotdf$Type)) {
  # i<-"Granule"
  tmp<-tmp_plotdf[tmp_plotdf$Type%in%i,]
  (p<-ggplot(tmp, aes(x=time, y=All_mean, colour=rowname, group=rowname)) + 
      # geom_errorbar(aes(ymin=All_mean-All_sd, ymax=All_mean+All_sd, colour=rowname), width=.4, position=pd) +
      geom_line(aes(group=rowname,linetype = Signal),position=pd) +
      geom_point(aes(x=time, y=Expr, colour=rowname, group=rowname),position=pd,size=1) +
      geom_point(aes(x=time, y=All_mean, colour=rowname, group=rowname),position=pd, size=3, shape="-", fill="white")+
      geom_text_repel(data=unique(tmp[tmp$time%in%"P3",c("rowname","time","All_mean","Type")]),
                      aes(label = rowname),max.overlaps=Inf)+
      # stat_summary(fun=mean, geom="point")+
      # stat_summary(fun=median, geom="line", aes(group=Type,linetype = Signal))+
      scale_fill_manual(values = c("#3578ad","#4ca64a","#8ecfc9","#d5211e","#ee7b1b"))+
      scale_color_manual(values = c("#3578ad","#4ca64a","#8ecfc9","#d5211e","#ee7b1b"))+
      guides(linetype = guide_legend(ncol = 1,title = "Significance"))+
      guides(color = guide_legend(ncol = 1,title = "Protein"))+
      # geom_point()+
      # facet_wrap(~Type)+
      ylab(label = "Mean Log10 (protein abundance)")+
      theme_bw()+
      theme(axis.title = element_text(size=18,colour = "black"),
            axis.title.x = element_blank(),
            axis.text = element_text(size=15,colour = "black"),
            # legend.position = "none",
            legend.text= element_text(size=13),
            legend.title = element_text(size=13),
            strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
            panel.grid = element_blank(),
            strip.placement = "outside",
            strip.switch.pad.wrap = unit(0.3, "mm"),
            plot.title = element_blank(),
            strip.text = element_text(size =15,face = "bold",colour = "white")))
  ggsave(filename = paste("figure/Fig3DE.",i,".Change_by_Time_in_CB.Pointplot_add_line.pdf",sep = ""),p,width = 7,height = 5)
  
}
# — Figure4 ----------------------------------------------------------
# —— Fig4E Gene sets Change by Time in cortical & subcortical point add Line  ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))
tmp_meta <- tmp_meta %>% 
  mutate(tmp_region=region2,region2=case_when(
    tmp_region %in% c("FL","TL","PL","V1") ~ "Cortical",
    tmp_region %in% c("STr","Hipp","Amy","MD") ~ "Subcortical",
    tmp_region %in% c("CB") ~ "CB"
  ))

tmp_expr<-tmp_expr %>% rownames_to_column() %>% {.$Gene_name<-gsub("\\|.*","",.$rowname);.} %>% 
  .[,c("Gene_name",tmp_meta$sample)] %>% as.data.table() %>% 
  .[,lapply(.SD,mean),"Gene_name"] %>% as.data.frame() %>% 
  {rownames(.)<-.$Gene_name;.} %>% subset(.,select=-Gene_name) %>% {log10(.+1)}

disease <- read.table("data/4type.Geneset.txt", quote="\"", comment.char="")
disease_Gene<-disease$V1 %>% intersect(.,rownames(tmp_expr)) %>%{ disease[disease$V1%in%.,]} %>% 
  {colnames(.)<-c("Gene","Type");.}

tmp_out<-data.frame()
for (tmp_Region in c("Cortical","Subcortical")) {
  # tmp_Region<-"Cortical"
  for (tmp_Time in unique(tmp_meta$time)) {
    # tmp_Time<-"F50"
    if (nrow(tmp_meta[tmp_meta$region2%in%tmp_Region&tmp_meta$time%in%tmp_Time,])!=0) {
      tmp_Sample<-tmp_meta[tmp_meta$region2%in%tmp_Region&tmp_meta$time%in%tmp_Time,]
      tmp_out<-tmp_expr[unique(disease_Gene$Gene),tmp_Sample$sample] %>% rownames_to_column() %>% gather(Samples,Expr,-rowname) %>% 
        left_join(.,disease_Gene,by=c("rowname"="Gene")) %>% 
        group_by(Samples,Type) %>% dplyr::summarise(mean_Expr = mean(Expr)) %>% as.data.frame() %>%
        left_join(.,tmp_meta,by=c("Samples"="sample")) %>% rbind(tmp_out,.)
    }
  }
}

tmp_out$P3_vs_F120<-NA
tmp_out$F120_vs_F90<-NA
tmp_out$F90_vs_F50<-NA
for (tmp_Region in c("Cortical","Subcortical")) {
  # tmp_Region<-"Cortical"
  for (tmp_Type in unique(tmp_out$Type)) {
    # tmp_Type<-"Dendrite_development"
    if (length(unique(as.character(tmp_out[tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$time)))>1) {
      tmp_list<-unique(as.character(tmp_out[tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$time)) %>% rev() %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
        {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
        .[.$Versus%in%intersect(c("P3_vs_F120","F120_vs_F90","F90_vs_F50"),.$Versus),]
      for (i in 1:nrow(tmp_list)) {
        # i<-1
        if (length(tmp_out[tmp_out$time%in%tmp_list[i,"V1"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$mean_Expr)>1&
            length(tmp_out[tmp_out$time%in%tmp_list[i,"V2"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$mean_Expr)>1) {
          tmp_out[tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,] %<>% 
            {.[[tmp_list[i,"Versus"]]]<-t.test(as.numeric(tmp_out[tmp_out$time%in%tmp_list[i,"V1"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$mean_Expr),
                                               as.numeric(tmp_out[tmp_out$time%in%tmp_list[i,"V2"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$mean_Expr))[["p.value"]];.}
        }
      }
    }
  }
}

tmp_plotdf<-tmp_out %>% {.$Signal<-ifelse(.$P3_vs_F120<0.05|.$F120_vs_F90<0.05|.$F90_vs_F50<0.05,"Yes","No");.}
# save(tmp_plotdf,file = "result/Fig5D.Disease_Change_by_Time_in_cortical.v7.RData")
# write.table(tmp_plotdf[order(tmp_plotdf$region2,tmp_plotdf$Type,tmp_plotdf$time),],
#             file = "result/Fig5D.Disease_Change_by_Time_in_cortical.v7.xls",
#             quote = FALSE,sep = "\t",row.names=F, col.names=T)
# 
# load("result/Fig5D.Disease_Change_by_Time_in_cortical.v7.RData")
tmp_plotdf$time<-factor(tmp_plotdf$time,levels = c("F50","F90","F120","P3"))
tmp_plotdf$Type<-factor(tmp_plotdf$Type,levels = c("Dendrite_development","Neuron_differentiation","Neuronal_activity","Synapse_assembly"))
tmp_plotdf$Signal<-factor(tmp_plotdf$Signal,levels = c("Yes","No"))

tmp_plotdf<-tmp_plotdf %>% group_by(time,Type,region2) %>% 
  dplyr::summarise(All_mean = mean(mean_Expr),All_sd = sd(mean_Expr),
                   mean_Expr=mean_Expr,time=time,Type=Type,region2=region2,Signal=Signal)

pd <- position_dodge(0.5)
(p1<-ggplot(tmp_plotdf[tmp_plotdf$region2%in%"Cortical",], aes(x=time, y=All_mean, colour=Type, group=Type)) + 
    geom_errorbar(aes(ymin=All_mean-All_sd, ymax=All_mean+All_sd, colour=Type), width=.4, position=pd) +
    geom_line(aes(group=Type,linetype = Signal),position=pd) +
    geom_point(position=pd, size=3, shape=21, fill="white",alpha=.3)+
    # stat_summary(fun=mean, geom="point")+
    # stat_summary(fun=median, geom="line", aes(group=Type,linetype = Signal))+
    scale_fill_manual(values = c("#3578ad","#4ca64a","#8ecfc9","#d5211e","#ee7b1b"))+
    scale_color_manual(values = c("#3578ad","#4ca64a","#8ecfc9","#d5211e","#ee7b1b"))+
    guides(linetype = guide_legend(ncol = 1,title = "Significance"))+
    # geom_point()+
    facet_wrap(~region2)+
    ylab(label = "Mean Log10 (protein abundance)")+
    ylim(1.3,5.8)+
    theme_bw()+
    theme(axis.title = element_text(size=18,colour = "black"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          # legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white")))

(p2<-ggplot(tmp_plotdf[tmp_plotdf$region2%in%"Subcortical",], aes(x=time, y=All_mean, colour=Type, group=Type)) + 
    geom_errorbar(aes(ymin=All_mean-All_sd, ymax=All_mean+All_sd, colour=Type), width=.4, position=pd) +
    geom_line(aes(group=Type,linetype = Signal),position=pd) +
    geom_point(position=pd, size=3, shape=21, fill="white",alpha=.3)+
    # stat_summary(fun=mean, geom="point")+
    # stat_summary(fun=median, geom="line", aes(group=Type,linetype = Signal))+
    scale_fill_manual(values = c("#3578ad","#4ca64a","#8ecfc9","#d5211e","#ee7b1b"))+
    scale_color_manual(values = c("#3578ad","#4ca64a","#8ecfc9","#d5211e","#ee7b1b"))+
    guides(linetype = guide_legend(ncol = 1,title = "Significance"))+
    # geom_point()+
    facet_wrap(~region2)+
    ylab(label = "Mean Log10 (protein abundance)")+
    ylim(1.3,5.8)+
    theme_bw()+
    theme(axis.title = element_text(size=18,colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          # legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white")))

p<-(p1|p2)+plot_layout(guides = "collect",widths = c(4,3))
ggsave(filename = "figure/Fig4E.4Type_Genesets_Change_by_Time.Pointplot_add_line.pdf",p,width = 10,height = 5)


# — Figure5 ----------------------------------------------------------
# —— Fig5B Protein Change by Time in cortical BarPlot ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))

tmp_expr_Counts<-tmp_expr %>% {log2(.+1)}
tmp_expr[tmp_expr>0]<-1

(load("data/human_enrichment.RData"))

tmp_plotdf<-data.frame()
for (tmp_Region in c("FL","TL","PL","V1")) {
  # tmp_Region<-"FL"
  tmp_list<-unique(as.character(tmp_meta$time)) %>% rev() %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
    {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
    .[.$Versus%in%intersect(c("P3_vs_F120","F120_vs_F90","F90_vs_F50"),.$Versus),]
  for (n in 1:nrow(tmp_list)) {
    # n=1
    pvalue<-c();tstat<-c();meansControl<-c();meansCase<-c();FC<-c();log2FC<-c()
    
    Case_list<-tmp_meta[tmp_meta$time%in%tmp_list[n,"V2"],]$sample
    Control_list<-tmp_meta[tmp_meta$time%in%tmp_list[n,"V1"],]$sample
    
    edata<-tmp_expr_Counts[,c(Case_list,Control_list)] %>% .[rowSums(.)!=0,] %>% {. + 0.01}
    for (tmp_row in 1:nrow(edata)){
      #t.test
      Control_value<-edata[tmp_row,Control_list]
      Case_value<-edata[tmp_row,Case_list]
      result<-t.test(as.numeric(Case_value), as.numeric(Control_value), paired=FALSE);
      pvalue[tmp_row]<-result[["p.value"]]
      tstat[tmp_row]<-result[["statistic"]][["t"]]
      meansControl[tmp_row]<-mean(as.numeric(Control_value))
      meansCase[tmp_row]<-mean(as.numeric(Case_value))
      FC[tmp_row]<-mean(as.numeric(Case_value))/mean(as.numeric(Control_value))
      log2FC[tmp_row]<-log2(FC[tmp_row])
    }
    p_adj<-p.adjust(pvalue, method = "fdr")
    diff_df<-data.frame(rownames(edata),pvalue,tstat,meansCase,meansControl,FC,log2FC,p_adj,stringsAsFactors=FALSE) %>% 
      .[.$p_adj<0.05,]
    
    tmp_split<-tmp_meta[tmp_meta$time%in%c(tmp_list[n,"V1"],tmp_list[n,"V2"])&tmp_meta$region2%in%tmp_Region,] %>% 
      {unique(.$sample)} %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
      {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
      left_join(.,tmp_meta[,c("sample","time")],by=c("V1"="sample")) %>% {names(.)[names(.) == 'time'] <- 'Case_time';.} %>% 
      left_join(.,tmp_meta[,c("sample","time")],by=c("V2"="sample")) %>% {names(.)[names(.) == 'time'] <- 'Control_time';.} %>% 
      .[.$Case_time!=.$Control_time,] %>% {.$Time<-paste(.$Control_time,"_vs_",.$Case_time,sep = "");.}
    for (i in 1:nrow(tmp_split)) {
      # i<-1
      tmp_New<-{tmp_expr[,tmp_split[i,"V2"],drop=F]-tmp_expr[,tmp_split[i,"V1"],drop=F]} %>% 
        {colnames(.)<-"New";.} %>% .[.$New==1,,drop=F] %>% 
        {.$Region<-tmp_Region;.$Time<-tmp_split[i,"Time"];.$Versus<-tmp_split[i,"Versus"];.} %>% 
        {.$Gene<-rownames(.);.} %>% subset(.,select=-New) %>% 
        .[intersect(rownames(.),diff_df$`rownames.edata.`),]
      
      tmp_Gene<-tmp_New$Gene %>% gsub("\\|.*","",.)
      term2gene <- gobp[,c('GO','SYMBOL')]
      term2gene<-term2gene[!is.na(term2gene$SYMBOL),]
      term2name <- gobp[,c('GO','Name')]
      tmp_enricher<-enricher(gene = tmp_Gene,
                             pvalueCutoff = 0.05,pAdjustMethod = "fdr",
                             minGSSize = 10,maxGSSize = Inf,
                             qvalueCutoff = 0.2,
                             TERM2GENE = term2gene,TERM2NAME = term2name) %>% 
        .@result %>% as.data.frame() %>% .[.$pvalue<0.05,] %>% 
        {.$Region<-tmp_Region;.$Time<-tmp_split[i,"Time"];.$Versus<-tmp_split[i,"Versus"];.}
      
      
      tmp_plotdf<-data.frame(Region=tmp_Region,Time=tmp_split[i,"Time"],Versus=tmp_split[i,"Versus"],
                             nGene=nrow(tmp_New),nTerm=nrow(tmp_enricher)) %>% rbind(tmp_plotdf,.)
      
    }
  }
}

save(tmp_plotdf,file = "result/Fig5B.Protein_Change.by_Time.in_cortical.RData")

load("result/Fig5B.Protein_Change.by_Time.in_cortical.RData")

PFC<-c("S50_1_1L","S50_2_1L","S50_3_2L",
       "S90_1_2L","S90_2_6L","S90_2_6R",
       "S120_3_6L","S120_4_6L","S120_4_6R",
       "SP3_1_6L","SP3_2_6L","SP3_2_6R")

tmp_plotdf %<>% {.$tmp1<-gsub("_vs_.*","",.$Versus);.} %<>%
  {.$tmp2<-gsub(".*_vs_","",.$Versus);.} %<>%
  .[.$Region%in%c("TL","PL","V1")|(.$tmp1%in%PFC&.$tmp2%in%PFC),]

tmp_plotdf$Time<-factor(tmp_plotdf$Time,levels = rev(c("P3_vs_F120","F120_vs_F90","F90_vs_F50")))
levels(tmp_plotdf$Time) %<>% {.[. %in% c("P3_vs_F120")] <- "F120_vs_P3";.}
levels(tmp_plotdf$Time) %<>% {.[. %in% c("F120_vs_F90")] <- "F90_vs_F120";.}
levels(tmp_plotdf$Time) %<>% {.[. %in% c("F90_vs_F50")] <- "F50_vs_F90";.}

tmp_plotdf<-tmp_plotdf %>% group_by(Time,Region) %>% 
  dplyr::summarise(Mean_nGene = mean(nGene),SD_nGene=sd(nGene),
                   Mean_nTerm = mean(nTerm),SD_nTerm=sd(nTerm),
                   Time=Time,Region=Region,nGene=nGene,nTerm=nTerm) %>% 
  {.[.$Region%in%"FL",]$Region<-"PFC";.}

tmp_plotdf$Region<-factor(tmp_plotdf$Region,levels = c("PFC","TL","PL","V1"))

my_comparisons=list(c("PFC","PL"),c("PFC","TL"),c("PFC","V1"))
(p1<-ggplot(tmp_plotdf, aes(x=Region, y=nGene, fill=Region)) + 
  geom_bar(aes(x=Region, y=Mean_nGene, fill=Region),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Mean_nGene - SD_nGene, ymax=Mean_nGene + SD_nGene),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  facet_wrap(~Time,ncol = 4)+
  ggpubr::stat_compare_means(aes(group=Region),comparisons = my_comparisons,
                             symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                              symbols = c("****", "***", "**", "*", "ns")),
                             method = "t.test",
                             label = "p.signif",
                             size = 5)+
  scale_fill_manual(values = c("#6cb943","#028589","#376ab2","#1866b1"))+
  scale_color_manual(values = c("#6cb943","#028589","#376ab2","#1866b1"))+
  # geom_point()+
  ylab(label = "New proteins")+
  theme_bw()+
  theme(axis.title = element_text(size=18,colour = "black"),
        axis.title.x = element_blank(),
        axis.text = element_text(size=15,colour = "black"),
        axis.text.x = element_text(angle = 30,hjust = 1),
        legend.position = "none",
        legend.text= element_text(size=13),
        legend.title = element_text(size=13),
        strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
        panel.grid = element_blank(),
        strip.placement = "outside",
        strip.switch.pad.wrap = unit(0.3, "mm"),
        plot.title = element_blank(),
        strip.text = element_text(size =15,face = "bold",colour = "white")))

my_comparisons=list(c("PFC","PL"),c("PFC","TL"),c("PFC","V1"))
(p2<-ggplot(tmp_plotdf, aes(x=Region, y=nTerm, fill=Region)) + 
    geom_bar(aes(x=Region, y=Mean_nTerm, fill=Region),position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Mean_nTerm - SD_nTerm, ymax=Mean_nTerm + SD_nTerm),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    facet_wrap(~Time,ncol = 4)+
    ggpubr::stat_compare_means(aes(group=Region),comparisons = my_comparisons,
                               symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                symbols = c("****", "***", "**", "*", "ns")),
                               method = "t.test",
                               label = "p.signif",
                               size = 5)+
    scale_fill_manual(values = c("#6cb943","#028589","#376ab2","#1866b1"))+
    scale_color_manual(values = c("#6cb943","#028589","#376ab2","#1866b1"))+
    # geom_point()+
    ylab(label = "New biological processes")+
    theme_bw()+
    theme(axis.title = element_text(size=18,colour = "black"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          axis.text.x = element_text(angle = 30,hjust = 1),
          legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white")))

(p<-p1/p2)
ggsave(filename = "figure/Fig5B.Protein_Change_by_Time.in_cortical.BarPlot.pdf",p,width = 7,height = 9)

# —— Fig5C Temporal Change by Time in cortical BarPlot ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))

tmp_expr<-tmp_expr %>% rownames_to_column() %>% {.$Gene_name<-gsub("\\|.*","",.$rowname);.} %>% 
  .[,c("Gene_name",tmp_meta$sample)] %>% as.data.table() %>% 
  .[,lapply(.SD,mean),"Gene_name"] %>% as.data.frame() %>% 
  {rownames(.)<-.$Gene_name;.} %>% subset(.,select=-Gene_name) %>% {log10(.+1)}

synapse <- read.table("data/synapse_assembly.txt", quote="\"", comment.char="") %>% .$V1 %>% 
  intersect(.,gsub("\\|.*","",rownames(tmp_expr)))

tmp_plotdf<-data.frame()
for (tmp_Region in c("FL","TL","PL","V1")) {
  # tmp_Region<-"FL"
  tmp_list<-unique(as.character(tmp_meta$time)) %>% rev() %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
    {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
    .[.$Versus%in%intersect(c("P3_vs_F120","F120_vs_F90","F90_vs_F50"),.$Versus),] %>% 
    {.$Days<-c(33,30,40);.}
  for (n in 1:nrow(tmp_list)) {
    # n=1
    tmp_split<-tmp_meta[tmp_meta$time%in%c(tmp_list[n,"V1"],tmp_list[n,"V2"])&tmp_meta$region2%in%tmp_Region,] %>% 
      {unique(.$sample)} %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
      {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
      left_join(.,tmp_meta[,c("sample","time")],by=c("V1"="sample")) %>% {names(.)[names(.) == 'time'] <- 'Case_time';.} %>% 
      left_join(.,tmp_meta[,c("sample","time")],by=c("V2"="sample")) %>% {names(.)[names(.) == 'time'] <- 'Control_time';.} %>% 
      .[.$Case_time!=.$Control_time,] %>% {.$Time<-paste(.$Control_time,"_vs_",.$Case_time,sep = "");.} %>% 
      {.$Days<-tmp_list[n,"Days"];.}
    for (i in 1:nrow(tmp_split)) {
      # i<-1
      tmp_Total<-{(colSums(tmp_expr[,tmp_split[i,"V2"],drop=F])-colSums(tmp_expr[,tmp_split[i,"V1"],drop=F]))/log2(tmp_split[i,"Days"])} %>% 
        as.numeric()
      
      tmp_Synaptic<-{(colSums(tmp_expr[synapse,tmp_split[i,"V2"],drop=F])-colSums(tmp_expr[synapse,tmp_split[i,"V1"],drop=F]))/log2(tmp_split[i,"Days"])} %>% 
        as.numeric()
      
      tmp_plotdf<-data.frame(Region=tmp_Region,Time=tmp_split[i,"Time"],Versus=tmp_split[i,"Versus"],
                             Total=tmp_Total,Synaptic=tmp_Synaptic) %>% rbind(tmp_plotdf,.)
      
    }
  }
}

save(tmp_plotdf,file = "result/Fig5C.Temporal_Change.by_Time.in_cortical.RData")
load("result/Fig5C.Temporal_Change.by_Time.in_cortical.RData")

PFC<-c("S50_1_1L","S50_2_1L","S50_3_2L",
       "S90_1_2L","S90_2_6L","S90_2_6R",
       "S120_3_6L","S120_4_6L","S120_4_6R",
       "SP3_1_6L","SP3_2_6L","SP3_2_6R")

tmp_plotdf %<>% {.$tmp1<-gsub("_vs_.*","",.$Versus);.} %<>%
  {.$tmp2<-gsub(".*_vs_","",.$Versus);.} %<>%
  .[.$Region%in%c("TL","PL","V1")|(.$tmp1%in%PFC&.$tmp2%in%PFC),]

tmp_plotdf$Time<-factor(tmp_plotdf$Time,levels = rev(c("P3_vs_F120","F120_vs_F90","F90_vs_F50")))

levels(tmp_plotdf$Time) %<>% {.[. %in% c("P3_vs_F120")] <- "F120_vs_P3";.}
levels(tmp_plotdf$Time) %<>% {.[. %in% c("F120_vs_F90")] <- "F90_vs_F120";.}
levels(tmp_plotdf$Time) %<>% {.[. %in% c("F90_vs_F50")] <- "F50_vs_F90";.}

tmp_plotdf<-tmp_plotdf %>% group_by(Time,Region) %>% 
  dplyr::summarise(Mean_Total = mean(Total),SD_Total=sd(Total),
                   Mean_Synaptic = mean(Synaptic),SD_Synaptic=sd(Synaptic),
                   Time=Time,Region=Region,Total=Total,Synaptic=Synaptic) %>% 
  {.[.$Region%in%"FL",]$Region<-"PFC";.}

tmp_plotdf$Region<-factor(tmp_plotdf$Region,levels = c("PFC","TL","PL","V1"))

my_comparisons=list(c("PFC","PL"),c("PFC","TL"),c("PFC","V1"))
(p1<-ggplot(tmp_plotdf, aes(x=Region, y=Total, fill=Region)) + 
    geom_bar(aes(x=Region, y=Mean_Total, fill=Region),position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Mean_Total - SD_Total, ymax=Mean_Total + SD_Total),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    facet_wrap(~Time,ncol = 4)+
    ggpubr::stat_compare_means(aes(group=Region),comparisons = my_comparisons,
                               symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                symbols = c("****", "***", "**", "*", "ns")),
                               method = "t.test",
                               label = "p.signif",
                               size = 5)+
    scale_fill_manual(values = c("#6cb943","#028589","#376ab2","#1866b1"))+
    scale_color_manual(values = c("#6cb943","#028589","#376ab2","#1866b1"))+
    scale_y_continuous(breaks=c(-400,-200,0,200,400))+
    # geom_point()+
    ylab(label = "Total proteins/log2(days)")+
    theme_bw()+
    theme(axis.title = element_text(size=18,colour = "black"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          axis.text.x = element_text(angle = 30,hjust = 1),
          legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white")))

my_comparisons=list(c("PFC","PL"),c("PFC","TL"),c("PFC","V1"))
(p2<-ggplot(tmp_plotdf, aes(x=Region, y=Synaptic, fill=Region)) + 
    geom_bar(aes(x=Region, y=Mean_Synaptic, fill=Region),position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=Mean_Synaptic - SD_Synaptic, ymax=Mean_Synaptic + SD_Synaptic),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
    facet_wrap(~Time,ncol = 4)+
    ggpubr::stat_compare_means(aes(group=Region),comparisons = my_comparisons,
                               symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
                                                symbols = c("****", "***", "**", "*", "ns")),
                               method = "t.test",
                               label = "p.signif",
                               size = 5)+
    scale_fill_manual(values = c("#6cb943","#028589","#376ab2","#1866b1"))+
    scale_color_manual(values = c("#6cb943","#028589","#376ab2","#1866b1"))+
    ylab(label = "Synaptic proteins/log2(days)")+
    theme_bw()+
    theme(axis.title = element_text(size=18,colour = "black"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          axis.text.x = element_text(angle = 30,hjust = 1),
          legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white")))


(p<-p1/p2)
ggsave(filename = "figure/Fig5C.Temporal_Change_by_Time.in_cortical.BarPlot.pdf",p,width = 7,height = 9)


# —— Fig5D Disease Change by Time in cortical point add Line ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())
(load("result/Proteomics.tmp_expr.tmp_meta.RData"))
tmp_meta<-tmp_meta[tmp_meta$region2%in%"V1"|tmp_meta$sample%in%c("S50_1_1L","S50_2_1L","S50_3_2L","S90_1_2L","S90_2_6L","S90_2_6R","S120_3_6L","S120_4_6L","S120_4_6R","SP3_1_6L","SP3_2_6L","SP3_2_6R"),]

tmp_expr<-tmp_expr %>% rownames_to_column() %>% {.$Gene_name<-gsub("\\|.*","",.$rowname);.} %>% 
  .[,c("Gene_name",tmp_meta$sample)] %>% as.data.table() %>% 
  .[,lapply(.SD,mean),"Gene_name"] %>% as.data.frame() %>% 
  {rownames(.)<-.$Gene_name;.} %>% subset(.,select=-Gene_name) %>% {log10(.+1)}

disease <- read.table("data/disease_gene.v2.txt", quote="\"", comment.char="") %>% 
  .[.$V2%in%c("AD","ASD","MDD","PD","SCZ"),]
disease_Gene<-disease$V1 %>% intersect(.,rownames(tmp_expr)) %>%{ disease[disease$V1%in%.,]} %>% 
  {colnames(.)<-c("Gene","Type");.}

tmp_out<-data.frame()
for (tmp_Region in c("FL","V1")) {
  # tmp_Region<-"FL"
  for (tmp_Time in unique(tmp_meta$time)) {
    # tmp_Time<-"F50"
    tmp_Sample<-tmp_meta[tmp_meta$region2%in%tmp_Region&tmp_meta$time%in%tmp_Time,]
    tmp_out<-tmp_expr[unique(disease_Gene$Gene),tmp_Sample$sample] %>% rownames_to_column() %>% gather(Samples,Expr,-rowname) %>% 
      left_join(.,disease_Gene,by=c("rowname"="Gene")) %>% 
      group_by(Samples,Type) %>% dplyr::summarise(mean_Expr = mean(Expr)) %>% as.data.frame() %>%
      left_join(.,tmp_meta,by=c("Samples"="sample")) %>% rbind(tmp_out,.)
  }
}

tmp_out$P3_vs_F120<-NA
tmp_out$F120_vs_F90<-NA
tmp_out$F90_vs_F50<-NA
for (tmp_Region in c("FL","V1")) {
  # tmp_Region<-"FL"
  for (tmp_Type in unique(tmp_out$Type)) {
    # tmp_Type<-"AD"
    if (length(unique(as.character(tmp_out[tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$time)))>1) {
      tmp_list<-unique(as.character(tmp_out[tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$time)) %>% rev() %>% combn(.,2) %>% t() %>% as.data.frame() %>% 
        {.$Versus<-paste(.$V1,"_vs_",.$V2,sep = "");.} %>% 
        .[.$Versus%in%intersect(c("P3_vs_F120","F120_vs_F90","F90_vs_F50"),.$Versus),]
      for (i in 1:nrow(tmp_list)) {
        # i<-1
        if (length(tmp_out[tmp_out$time%in%tmp_list[i,"V1"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$mean_Expr)>1&
            length(tmp_out[tmp_out$time%in%tmp_list[i,"V2"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$mean_Expr)>1) {
          tmp_out[tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,] %<>% 
            {.[[tmp_list[i,"Versus"]]]<-t.test(as.numeric(tmp_out[tmp_out$time%in%tmp_list[i,"V1"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$mean_Expr),
                                               as.numeric(tmp_out[tmp_out$time%in%tmp_list[i,"V2"]&tmp_out$Type%in%tmp_Type&tmp_out$region2%in%tmp_Region,]$mean_Expr))[["p.value"]];.}
        }
      }
    }
  }
}

tmp_plotdf<-tmp_out %>% {.$Signal<-ifelse(.$P3_vs_F120<0.05|.$F120_vs_F90<0.05|.$F90_vs_F50<0.05,"Yes","No");.}
save(tmp_plotdf,file = "result/Fig5D.Disease_Change_by_Time_in_cortical.RData")
write.table(tmp_plotdf[order(tmp_plotdf$region2,tmp_plotdf$Type,tmp_plotdf$time),],
            file = "result/Fig5D.Disease_Change_by_Time_in_cortical.xls",
            quote = FALSE,sep = "\t",row.names=F, col.names=T)

load("result/Fig5D.Disease_Change_by_Time_in_cortical.RData")
tmp_plotdf$time<-factor(tmp_plotdf$time,levels = c("F50","F90","F120","P3"))
tmp_plotdf$Type<-factor(tmp_plotdf$Type,levels = c("ASD","MDD","SCZ","AD","PD"))
tmp_plotdf$Signal<-factor(tmp_plotdf$Signal,levels = c("Yes","No"))

tmp_plotdf<-tmp_plotdf %>% group_by(time,Type,region2) %>% 
  dplyr::summarise(All_mean = mean(mean_Expr),All_sd = sd(mean_Expr),
                   mean_Expr=mean_Expr,time=time,Type=Type,region2=region2,Signal=Signal) %>% 
  {.[.$region2%in%"FL",]$region2<-"PFC";.}


pd <- position_dodge(0.5)
(p<-ggplot(tmp_plotdf, aes(x=time, y=All_mean, colour=Type, group=Type)) + 
    geom_errorbar(aes(ymin=All_mean-All_sd, ymax=All_mean+All_sd, colour=Type), width=.4, position=pd) +
    geom_line(aes(group=Type,linetype = Signal),position=pd) +
    geom_point(position=pd, size=3, shape=21, fill="white",alpha=.3)+
    # stat_summary(fun=mean, geom="point")+
    # stat_summary(fun=median, geom="line", aes(group=Type,linetype = Signal))+
    scale_fill_manual(values = c("#3578ad","#4ca64a","#8ecfc9","#d5211e","#ee7b1b"))+
    scale_color_manual(values = c("#3578ad","#4ca64a","#8ecfc9","#d5211e","#ee7b1b"))+
    guides(linetype = guide_legend(ncol = 1,title = "Significance"))+
    # geom_point()+
    facet_wrap(~region2)+
    ylab(label = "Mean Log10 (protein abundance)")+
    theme_bw()+
    theme(axis.title = element_text(size=18,colour = "black"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          # legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white")))
ggsave(filename = "figure/Fig5D.Disease_Change_by_Time_in_cortical.Pointplot_add_line.pdf",p,width = 13,height = 5)









# — Figure6 ----------------------------------------------------------
# —— Fig6C human_monkey_mouse Comparison enricher ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())

Corr <- read.delim("data/Protein_cor.tab") %>% 
  .[,c("gene","Mouse_vs_Human","Monkey_vs_Mouse","MonkeyD_vs_Human")]

tmp1<-Corr[Corr$MonkeyD_vs_Human>0.75&Corr$MonkeyD_vs_Human<=1,]$gene
tmp2<-Corr[Corr$Monkey_vs_Mouse>0.75&Corr$Monkey_vs_Mouse<=1,]$gene
tmp3<-Corr[Corr$Mouse_vs_Human>0.75&Corr$Mouse_vs_Human<=1,]$gene

primate_specific<-setdiff(tmp1,tmp2) %>% setdiff(.,tmp3)

tmp4<-Corr[Corr$MonkeyD_vs_Human> -0.75&Corr$MonkeyD_vs_Human<= -0.5,]$gene
tmp5<-Corr[Corr$Mouse_vs_Human> -0.75&Corr$Mouse_vs_Human<= -0.5,]$gene

human_specific<-c(tmp4,tmp5) %>% unique()

tmp6<-Corr[Corr$MonkeyD_vs_Human>0.75&Corr$MonkeyD_vs_Human<=1,]$gene
tmp7<-Corr[Corr$Monkey_vs_Mouse>0.75&Corr$Monkey_vs_Mouse<=1,]$gene

primate_Common<-intersect(tmp6,tmp7) %>% unique()


tmp_input<-data.frame(Gene=c(primate_specific,human_specific,primate_Common),
                      Group=c(rep("primate_specific",length(primate_specific)),
                              rep("human_specific",length(human_specific)),
                              rep("primate_Common",length(primate_Common))))



(load("data/human_enrichment.RData"))

tmp_enricher<-data.frame()
for (tmp_time in unique(tmp_input$Group)) {
  # tmp_time<-"primate_specific"
  tmp_Gene<-tmp_input[tmp_input$Group%in%tmp_time,]$Gene %>% unique()
  term2gene <- gobp[,c('GO','SYMBOL')]
  term2gene<-term2gene[!is.na(term2gene$SYMBOL),]
  term2name <- gobp[,c('GO','Name')]
  
  tmp_enricher<-enricher(gene = tmp_Gene,
                         pvalueCutoff = 0.05,pAdjustMethod = "fdr",
                         minGSSize = 10,maxGSSize = Inf,
                         qvalueCutoff = 0.2,
                         TERM2GENE = term2gene,TERM2NAME = term2name) %>% 
    .@result %>% as.data.frame() %>% .[.$pvalue<0.05,] %>% 
    {.$Group<-tmp_time;.} %>% rbind(tmp_enricher,.)
}
write.table(tmp_enricher,
            file = "result/Fig6C.human_monkey_mouse.Comparison.enricher.xls",
            quote = FALSE,sep = "\t",row.names=F, col.names=T)

(p1<-tmp_enricher %>% {.$Group<-factor(.$Group,levels = c("human_specific","primate_specific","primate_Common"));.} %>% 
    .[order(.$Group),] %>% 
    group_by(Group) %>% do(head(.,n=15)) %>% 
    {.[.$pvalue<1e-20,]$pvalue<-1e-20;.} %>% 
    {.$Description<-factor(.$Description,levels = unique(.$Description));.} %>% 
    # group_by(Corr,Group) %>% do(head(.,n=10)) %>% 
    ggplot(.,aes(x = Group,y=Description,size=Count,colour=-log10(pvalue),ylab=''))+
    geom_point()+
    # facet_wrap(~Corr)+
    scale_size_continuous(range=c(4,7))+
    scale_color_gradientn(colors = rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))+
    # labs(title = "Protein family")+
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text = element_text(size=15,colour = "black"),
          # axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),
          # legend.position = "none",
          legend.text= element_text(size=13),
          legend.title = element_text(size=13),
          strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
          panel.grid = element_blank(),
          strip.placement = "outside",
          strip.switch.pad.wrap = unit(0.3, "mm"),
          plot.title = element_blank(),
          strip.text = element_text(size =15,face = "bold",colour = "white"))+
    RotatedAxis())
ggsave(filename = "figure/Fig6C.human_monkey_mouse.Comparison.enricher.Bubble.pdf",p1,width = 11,height = 10)


# —— Fig6D RNA_Proteomics Period 8/Period 7 Log2FC Bubble -------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())

library(data.table)
tmp_own_meta <- read.delim("data/metadata.txt") %>% 
  {.[.$time%in%"A",]$time<-"F50";.} %>% 
  {.[.$time%in%"B",]$time<-"F90";.} %>% 
  {.[.$time%in%"C",]$time<-"F120";.} %>% 
  {.[.$time%in%"D",]$time<-"P3";.} %>% 
  {.[.$region2%in%"1",]$region2<-"FL";.} %>% 
  {.[.$region2%in%"2",]$region2<-"TL";.} %>% 
  {.[.$region2%in%"3",]$region2<-"PL";.} %>% 
  {.[.$region2%in%"4",]$region2<-"V1";.} %>% 
  {.[.$region2%in%"5",]$region2<-"CB";.} %>% 
  {.[.$region2%in%"6",]$region2<-"STr";.} %>% 
  {.[.$region2%in%"7",]$region2<-"Hipp";.} %>% 
  {.[.$region2%in%"8",]$region2<-"MD";.} %>% 
  {.[.$region2%in%"9",]$region2<-"Amy";.}

tmp_own_expr <- read.delim("data/protein_Counts.Table_S2.txt") %>% 
  .[.$Gene.name!="---",] %>% 
  .[,c("Gene.name",tmp_own_meta$sample)] %>% as.data.table() %>% 
  .[,lapply(.SD,mean),"Gene.name"] %>% as.data.frame() %>% 
  {rownames(.)<-.$Gene.name;.} %>% .[,tmp_own_meta$sample] %>% 
  .[rowSums(.)!=0,]

tmp_ref_meta <- read.delim("data/ref.meta.txt")
tmp_ref_RPKM <- read.delim("data/nhp_development_RPKM_rmTechRep.digits2.txt",row.names = 1) %>% 
  rownames_to_column() %>% 
  {.$Gene_name<-gsub(".*\\|","",.$rowname);.} %>% 
  .[,c("Gene_name",tmp_ref_meta$Sample)] %>% as.data.table() %>% 
  .[,lapply(.SD,mean),"Gene_name"] %>% as.data.frame() %>% 
  {rownames(.)<-.$Gene_name;.} %>% .[,tmp_ref_meta$Sample]

tmp_Gene<-intersect(rownames(tmp_own_expr),rownames(tmp_ref_RPKM))
tmp_Region<-c("mPFG","iPFG","OFC","M1","aCG","sTG",
              "iTG","S1","sPL","V1","CB","STr",
              "Hipp","MD","Amy")
tmp_ref_sample<-tmp_ref_meta[tmp_ref_meta$Predicted.period%in%c(7,8),] %>% 
  {.$Region_merge<-.$Region;.} %>% 
  {.[.$Region_merge%in%"AMY",]$Region_merge<-"Amy";.} %>% 
  {.[.$Region_merge%in%"CBC",]$Region_merge<-"CB";.} %>% 
  {.[.$Region_merge%in%"HIP",]$Region_merge<-"Hipp";.} %>% 
  {.[.$Region_merge%in%"MD",]$Region_merge<-"MD";.} %>% 
  {.[.$Region_merge%in%"A1C",]$Region_merge<-"None";.} %>% 
  {.[.$Region_merge%in%"DFC",]$Region_merge<-"mPFG";.} %>% 
  {.[.$Region_merge%in%"IPC",]$Region_merge<-"sPL";.} %>% 
  {.[.$Region_merge%in%"ITC",]$Region_merge<-"iTG";.} %>% 
  {.[.$Region_merge%in%"M1C",]$Region_merge<-"M1";.} %>% 
  {.[.$Region_merge%in%"MFC",]$Region_merge<-"aCG";.} %>% 
  {.[.$Region_merge%in%"OFC",]$Region_merge<-"OFC";.} %>% 
  {.[.$Region_merge%in%"S1C",]$Region_merge<-"S1";.} %>% 
  {.[.$Region_merge%in%"STC",]$Region_merge<-"sTG";.} %>% 
  {.[.$Region_merge%in%"V1C",]$Region_merge<-"V1";.} %>% 
  {.[.$Region_merge%in%"VFC",]$Region_merge<-"iPFG";.} %>% 
  {.[.$Region_merge%in%"STR",]$Region_merge<-"STr";.} %>% 
  .[.$Region_merge%in%tmp_Region,] %>% 
  {.$Region_merge<-factor(.$Region_merge,levels = tmp_Region);.}

tmp_own_sample<-tmp_own_meta[tmp_own_meta$time%in%c("F120","P3"),] %>% 
  .[.$region%in%tmp_Region,] %>% 
  {.$Region_merge<-factor(.$region,levels = tmp_Region);.}

tmp_meta1<-tmp_own_sample[,c("sample","time","Region_merge","region2")] %>% 
  {.$Species<-"Macaca fascicularis";.} %>% 
  {.$Type<-"Proteomics";.} %>% 
  {colnames(.)<-c("Samples","Time","Minor_Region","Major_Region","Species","Type");.} %>% 
  .[,c("Samples","Species","Time","Minor_Region","Major_Region","Type")]

tmp_meta2<-tmp_ref_sample[,c("Sample","Species","Age","Region_merge")] %>% 
  left_join(.,unique(tmp_meta1[,c("Minor_Region","Major_Region")]),by=c("Region_merge"="Minor_Region")) %>% 
  {.$Type<-"RNA-seq";.} %>% 
  {colnames(.)<-c("Samples","Species","Time","Minor_Region","Major_Region","Type");.} %>% 
  .[,c("Samples","Species","Time","Minor_Region","Major_Region","Type")] %>% 
  {.[.$Species%in%"Human",]$Species<-"Homo sapiens";.} %>% 
  {.[.$Species%in%"Macaque",]$Species<-"Macaca mulatta";.}

tmp_meta_final<-rbind(tmp_meta1,tmp_meta2)
# tmp_ref_expr_Counts<-tmp_ref_Counts[tmp_Gene,tmp_meta_final[tmp_meta_final$Type%in%"RNA-seq",]$Samples]
tmp_ref_expr_RPKM<-tmp_ref_RPKM[tmp_Gene,tmp_meta_final[tmp_meta_final$Type%in%"RNA-seq",]$Samples]

tmp_own_expr_Proteomics<-tmp_own_expr[tmp_Gene,tmp_meta_final[tmp_meta_final$Type%in%"Proteomics",]$Samples]


tmp_title<-"Macaca mulatta"

tmp1<-tmp_own_expr_Proteomics %>% 
  log1p() %>% 
  {(.-min(.))/(max(.)-min(.))} %>% 
  {.+0.01} %>% 
  rownames_to_column() %>% gather(Samples,Expr,-rowname) %>% 
  left_join(.,tmp_meta_final[,c("Time","Samples","Major_Region","Species")],by=c("Samples"="Samples")) %>% 
  group_by(rowname,Time,Major_Region,Species) %>% dplyr::summarise(Proteomics = mean(Expr)) %>% as.data.frame() %>% 
  {.$anchor<-paste(.$rowname,.$Major_Region,sep = "_");.}

tmp_gourp1<-tmp1[tmp1$Time%in%unique(tmp1$Time)[1],c("anchor","Proteomics")] %>%
  {colnames(.)<-c("anchor",unique(tmp1$Time)[1]);.} %>% 
  left_join(.,tmp1[tmp1$Time%in%unique(tmp1$Time)[2],c("anchor","rowname","Major_Region","Proteomics")],by=c("anchor"="anchor")) %>% 
  {names(.)[names(.) == "Proteomics"] <- unique(tmp1$Time)[2];.} %>% 
  {.$Log2FC<-log2(.[[unique(tmp1$Time)[2]]]/.[[unique(tmp1$Time)[1]]]);.} %>% 
  {.$Versus<-paste(unique(tmp1$Time)[2],"/",unique(tmp1$Time)[1],sep = "");.}

tmp_meta<-tmp_meta_final[tmp_meta_final$Species%in%tmp_title,]

tmp2<-tmp_ref_expr_RPKM[,tmp_meta$Samples] %>% 
  {(.-min(.))/(max(.)-min(.))} %>% 
  {.+0.01} %>% 
  rownames_to_column() %>% gather(Samples,Expr,-rowname) %>% 
  left_join(.,tmp_meta_final[,c("Time","Samples","Major_Region","Species")],by=c("Samples"="Samples")) %>% 
  group_by(rowname,Time,Major_Region,Species) %>% dplyr::summarise(RNA = mean(Expr)) %>% as.data.frame() %>% 
  {.$anchor<-paste(.$rowname,.$Major_Region,sep = "_");.}

tmp_gourp2<-tmp2[tmp2$Time%in%unique(tmp2$Time)[1],c("anchor","RNA")] %>%
  {colnames(.)<-c("anchor",unique(tmp2$Time)[1]);.} %>% 
  left_join(.,tmp2[tmp2$Time%in%unique(tmp2$Time)[2],c("anchor","rowname","Major_Region","RNA")],by=c("anchor"="anchor")) %>% 
  {names(.)[names(.) == "RNA"] <- unique(tmp2$Time)[2];.} %>% 
  {.$Log2FC<-log2(.[[unique(tmp2$Time)[2]]]/.[[unique(tmp2$Time)[1]]]);.} %>% 
  {.$Versus<-paste(unique(tmp2$Time)[2],"/",unique(tmp2$Time)[1],sep = "");.}

color<-c("#35703e","#4aa751","#287374","#1c6eb7","#dece62","#f58737","#e05a44","#b92645","#7e2243") %>% rev()
plotdf<-tmp_gourp1[,c("anchor","rowname","Major_Region","Log2FC")] %>% 
  {colnames(.)<-c("anchor","Gene","Major_Region","Proteomics");.} %>% 
  left_join(.,tmp_gourp2[,c("anchor","Log2FC")],by=c("anchor"="anchor")) %>% 
  {names(.)[names(.) == "Log2FC"] <- "RNA";.} %>% 
  {.$Major_Region<-factor(.$Major_Region,levels = c("FL","TL","PL","V1","CB","STr","Hipp","MD","Amy"));.} %>% 
  {.$Label<-"None";.} %>% 
  {.[.$RNA>=.$Proteomics-1&.$RNA<=.$Proteomics+1&.$RNA>=-1&.$RNA<=1&.$Proteomics>=-1&.$Proteomics<=1,]$Label<-"Type1";.} %>% 
  {.[.$RNA>=.$Proteomics-1&.$RNA<=.$Proteomics+1&.$Label!="Type1",]$Label<-"Type2";.} %>% 
  {.[.$RNA>=-.$Proteomics-1&.$RNA<=-.$Proteomics+1&.$Label!="Type1",]$Label<-"Type3";.} %>% 
  {.[.$RNA>=1&.$Proteomics>=1&.$Label%in%"None",]$Label<-"Type4";.} %>% 
  {.[.$RNA<=-1&.$Proteomics<=-1&.$Label%in%"None",]$Label<-"Type4";.} %>% 
  {.[.$RNA>-1&.$RNA<1&.$Label%in%"None",]$Label<-"Type5";.} %>% 
  {.[.$Proteomics>-1&.$Proteomics<1&.$Label%in%"None",]$Label<-"Type6";.} %>% 
  {.[.$Label%in%"None",]$Label<-"Type3";.}

plotdf$Label2<-as.character(plotdf$Label)
plotdf[plotdf$Label%in%"Type1",]$Label2<-"Type1"
plotdf[plotdf$Label%in%"Type5",]$Label2<-"Type2"
plotdf[plotdf$Label%in%"Type6",]$Label2<-"Type3"
plotdf[plotdf$Label%in%"Type2",]$Label2<-"Type4"
plotdf[plotdf$Label%in%"Type3",]$Label2<-"Type5"
plotdf[plotdf$Label%in%"Type4",]$Label2<-"Type6"
plotdf$Label2<-factor(plotdf$Label2,levels = c("Type1","Type2","Type3","Type4","Type5","Type6"))

save(plotdf,file = paste("result/Fig6D.Macaca_fascicularis.",gsub(" ","_",tmp_title),".Period8_vs_Period7.Gene_type.RData",sep = ""))

tmp_Plist<-list()
tmp_summary<-data.frame()
color6<-c("grey80","#6ab04c","#FFA500","#DC0000B2","#40407a","#7E6148B2")
for (i in levels(plotdf$Major_Region)) {
  # i<-"FL"
  tmp_plotdf<- plotdf[plotdf$Major_Region%in%i,]
  
  (p1<-ggplot(data=tmp_plotdf,aes(x=Proteomics,y=RNA,color=Label2))+
      geom_point(size=1,alpha=0.9)+
      xlim(min(tmp_plotdf$Proteomics,tmp_plotdf$RNA),max(tmp_plotdf$Proteomics,tmp_plotdf$RNA))+
      ylim(min(tmp_plotdf$Proteomics,tmp_plotdf$RNA),max(tmp_plotdf$Proteomics,tmp_plotdf$RNA))+
      geom_hline(yintercept = c(-1,1),lty=2,lwd=0.6,alpha=0.8)+
      geom_vline(xintercept = c(-1,1),lty=2,lwd=0.6,alpha=0.8)+
      geom_abline(intercept = c(-1,1),lty=2,lwd=0.6,alpha=0.8,slope = 1)+
      # geom_abline(intercept = c(-1,1),lty=2,lwd=0.6,alpha=0.8,slope = -1)+
      # labs(title = paste("Macaca_fascicularis & ",tmp_title,sep = ""))+
      labs(x="log2 fold-change by proteomics",y="log2 fold-change by RNA-seq")+
      scale_colour_manual(values = color6)+
      facet_wrap(~Major_Region,ncol = 5)+
      guides(color = guide_legend(override.aes = list(size = 5)))+
      theme_bw()+
      theme(axis.title = element_text(size=15),
            # axis.title.x = element_blank(),
            axis.text = element_text(size=15,colour = "black"),
            legend.text= element_text(size=13),
            legend.title = element_text(size=13),
            strip.background = element_rect(color="#7F8487",fill = "#7F8487"),
            panel.grid = element_blank(),
            strip.placement = "outside",
            strip.switch.pad.wrap = unit(0.3, "mm"),
            plot.title = element_text(size=18),
            strip.text = element_text(size =15,face = "bold",colour = "white")))
  
  (p2<-table(tmp_plotdf$Label2) %>% as.data.frame() %>% 
      {.$Var1<-factor(.$Var1,levels = c("Type1","Type2","Type3","Type4","Type5","Type6"));.} %>% 
      # .[.$Var1!="Type1",] %>% 
      ggplot(., aes(x="", y=Freq, fill=Var1)) +
      geom_bar(stat="identity", width=1, color="white") +
      scale_fill_manual(values = color6)+
      coord_polar("y", start=0) +
      theme_void()+
      NoLegend())
  
  p<-p1+inset_element(p2, left = 1, bottom = 0.08, right = 0.7, top = 0.2)
  tmp_Plist[[i]]<-p
  
  tmp_summary<-table(tmp_plotdf$Label2) %>% as.data.frame() %>% 
    {.$Sum<-sum(.$Freq);.} %>% {.$Percent<-(.$Freq/.$Sum)*100;.} %>% .[,c("Var1","Percent")] %>% 
    {.$Region<-i;.} %>% rbind(tmp_summary,.)
  
}

tmp_out<-tmp_summary %>% spread(.,key = "Region",value = "Percent") %>% {rownames(.)<-.$Var1;.} %>% subset(.,select=-Var1)
write.table(tmp_out[,c("FL","TL","PL","V1","STr","Hipp","MD","Amy","CB")],file = "result/Fig6D.Period8_vs_Period7.percent.summary.xls",quote = FALSE,sep = "\t",row.names=TRUE, col.names=NA)

p<-(tmp_Plist$FL+tmp_Plist$TL+tmp_Plist$PL+tmp_Plist$V1)+plot_layout(ncol = 2,guides = "collect")
ggsave(filename = paste("figure/Fig6D.Cortical.Period8_vs_Period7.RNA_Proteomics.v5.pdf",sep = ""),
       p,width = 9,height = 8)

p<-(tmp_Plist$STr+tmp_Plist$Hipp+tmp_Plist$MD+tmp_Plist$Amy)+plot_layout(ncol = 2,guides = "collect")
ggsave(filename = paste("figure/Fig6D.SubCortical.Period8_vs_Period7.RNA_Proteomics.v5.pdf",sep = ""),
       p,width = 9,height = 8)

p<-tmp_Plist$CB
ggsave(filename = paste("figure/Fig6D.CB.Period8_vs_Period7.RNA_Proteomics.v5.pdf",sep = ""),
       p,width = 5,height = 4)

######################## /========= Block2 =========/ #########################
############################# <2023.4.3 review> ###############################
##############################################################################-
# — Figure6 ----------------------------------------------------------
# —— Fig6A P3 top protein family cor with monkey mouse heatmap ----------------------------------------------------------
setwd("/data1/zhur/proj/Siwei/Brain/R_analysis_Final")
rm(list=ls())

protein_family<-read.delim("data/uniprot_human_family.tab", header=T,stringsAsFactors=FALSE,na.strings = "") %>% 
  .[,c("Gene.Names..primary.","Protein.families")] %>% {colnames(.)<-c("Gene_Name","Family");.} %>% na.omit()

(load("result/Fig1F.Protein-Family_percent.RData"))
tmp_head<-tmp_out_merge %>% {.[order(.$Mean_Percent,decreasing = T),]$Family} %>% unique() %>% head(.,n=11)

expr_all <- read.delim("data/all_expr.tab", row.names=1) %>%
  {.[,grep("^Monkey_D|^Human|^Mouse",colnames(.),value = T)]} %>%
  {colnames(.)<-gsub("Monkey_D","Monkey_",colnames(.));.} %>% 
  .[intersect(rownames(.),protein_family[protein_family$Family%in%c(tmp_head),]$Gene_Name),]

tmp_ph<-data.frame()
for (i in rownames(expr_all)) {
  # i<-"AAK1"
  tmp_cor<-expr_all[i,] %>% t() %>% as.data.frame() %>% rownames_to_column() %>% 
    separate(., col = rowname, into = c("Species", "Region"), sep = "_") %>% 
    {colnames(.)<-c("Species", "Region", "Expr");.} %>% 
    spread(.,key = "Species",value = "Expr") %>% 
    {rownames(.)<-.[,1];.} %>% .[,-1] %>% 
    cor(.,method = "spearman")
  
  tmp_ph<-data.frame(Gene=i,
                     Human_vs_Monkey=tmp_cor[grep("Human",rownames(tmp_cor),value = T),grep("Monkey",rownames(tmp_cor),value = T)],
                     Human_vs_Mouse=tmp_cor[grep("Human",rownames(tmp_cor),value = T),grep("Mouse",rownames(tmp_cor),value = T)],
                     Monkey_vs_Mouse=tmp_cor[grep("Monkey",rownames(tmp_cor),value = T),grep("Mouse",rownames(tmp_cor),value = T)]) %>% 
    rbind(.,tmp_ph)
}

tmp1<-tmp_ph %>% 
  gather(Versus,Corr,-Gene) %>% 
  left_join(.,protein_family[protein_family$Family%in%c(tmp_head),],by=c("Gene"="Gene_Name")) %>% na.omit() %>% 
  {.$Family<-factor(.$Family,levels = intersect(tmp_head,unique(.$Family)));.}

ann_colors<-list(
  `Versus`=c(`Human_vs_Monkey`="#3381ae",
             `Human_vs_Mouse`="#f7a369",
             `Monkey_vs_Mouse`="#cb2a2d"))


i<-c("Tubulin family","Metallo-dependent hydrolases superfamily","ALB/AFP/VDB family","Intermediate filament family")
annotation_row<-tmp1[tmp1$Family%in%i&tmp1$Versus%in%c("Human_vs_Monkey"),] %>% 
  .[order(.$Family,-.$Corr),] %>% .[,c("Gene","Family")] %>% unique() %>% 
  {rownames(.)<-.$Gene;.}
(p1<-tmp1[tmp1$Family%in%i,] %>% .[,c("Gene","Versus","Corr")] %>% 
    spread(.,key = "Versus",value = "Corr") %>% {rownames(.)<-.[,1];.} %>% .[,-1] %>% 
    .[annotation_row$Gene,c("Human_vs_Monkey","Human_vs_Mouse","Monkey_vs_Mouse")] %>% 
    pheatmap(.,cluster_cols = F,cluster_rows = F,display_numbers = F,number_color = "black",fontsize = 12,
             show_colnames = F,annotation_legend = T,border_color=NA,
             # breaks = c(seq(-1,0,by=0.25),seq(0.25,1,by=0.25)),
             annotation_col = data.frame(row.names = c("Human_vs_Monkey","Human_vs_Mouse","Monkey_vs_Mouse"),
                                         Versus=c("Human_vs_Monkey","Human_vs_Mouse","Monkey_vs_Mouse")),
             annotation_row = annotation_row[,"Family",drop=F],
             annotation_colors = ann_colors,
             color = colorRampPalette(rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))(50)) %>% as.ggplot())
write.table(annotation_row,file = "result/Fig6A.top_protein_family.cor.heatmap.part1.Gene_list.xls",quote = FALSE,sep = "\t",row.names=F, col.names=T)

i<-c("Protein kinase superfamily")
annotation_row<-tmp1[tmp1$Family%in%i&tmp1$Versus%in%c("Human_vs_Monkey"),] %>% 
  .[order(.$Family,-.$Corr),] %>% .[,c("Gene","Family")] %>% unique() %>% 
  {rownames(.)<-.$Gene;.}
(p2<-tmp1[tmp1$Family%in%i,] %>% .[,c("Gene","Versus","Corr")] %>% 
    spread(.,key = "Versus",value = "Corr") %>% {rownames(.)<-.[,1];.} %>% .[,-1] %>% 
    .[annotation_row$Gene,c("Human_vs_Monkey","Human_vs_Mouse","Monkey_vs_Mouse")] %>% 
    pheatmap(.,cluster_cols = F,cluster_rows = F,display_numbers = F,number_color = "black",fontsize = 12,
             show_colnames = F,show_rownames = F,annotation_legend = T,border_color=NA,
             # breaks = c(seq(-1,0,by=0.25),seq(0.25,1,by=0.25)),
             annotation_col = data.frame(row.names = c("Human_vs_Monkey","Human_vs_Mouse","Monkey_vs_Mouse"),
                                         Versus=c("Human_vs_Monkey","Human_vs_Mouse","Monkey_vs_Mouse")),
             annotation_row = annotation_row[,"Family",drop=F],
             annotation_colors = ann_colors,
             color = colorRampPalette(rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))(50)) %>% as.ggplot())

write.table(annotation_row,file = "result/Fig6A.top_protein_family.cor.heatmap.part2.Gene_list.xls",quote = FALSE,sep = "\t",row.names=F, col.names=T)

tmp_out<-tmp1[tmp1$Family%in%i,] %>% .[,c("Gene","Versus","Corr")] %>% 
  spread(.,key = "Versus",value = "Corr") %>% {rownames(.)<-.[,1];.} %>% .[,-1] %>% 
  .[annotation_row$Gene,c("Human_vs_Monkey","Human_vs_Mouse","Monkey_vs_Mouse")]
write.table(tmp_out,file = "result/Fig6A.top_protein_family.cor.heatmap.kinase.Gene_corr_list.xls",quote = FALSE,sep = "\t",row.names=T, col.names=NA)


i<-c("Spectrin family","MAP1 family","Heat shock protein 70 family","Immunoglobulin superfamily","Cation transport ATPase (P-type) (TC 3.A.3) family")
annotation_row<-tmp1[tmp1$Family%in%i&tmp1$Versus%in%c("Human_vs_Monkey"),] %>% 
  .[order(.$Family,-.$Corr),] %>% .[,c("Gene","Family")] %>% unique() %>% 
  {rownames(.)<-.$Gene;.}
(p3<-tmp1[tmp1$Family%in%i,] %>% .[,c("Gene","Versus","Corr")] %>% 
    spread(.,key = "Versus",value = "Corr") %>% {rownames(.)<-.[,1];.} %>% .[,-1] %>% 
    .[annotation_row$Gene,c("Human_vs_Monkey","Human_vs_Mouse","Monkey_vs_Mouse")] %>% 
    pheatmap(.,cluster_cols = F,cluster_rows = F,display_numbers = F,number_color = "black",fontsize = 12,
             show_colnames = F,show_rownames = T,annotation_legend = T,border_color=NA,
             # breaks = c(seq(-1,0,by=0.25),seq(0.25,1,by=0.25)),
             annotation_col = data.frame(row.names = c("Human_vs_Monkey","Human_vs_Mouse","Monkey_vs_Mouse"),
                                         Versus=c("Human_vs_Monkey","Human_vs_Mouse","Monkey_vs_Mouse")),
             annotation_row = annotation_row[,"Family",drop=F],
             annotation_colors = ann_colors,
             color = colorRampPalette(rev(c("#810000","#CE1212","#F05454","#ffd06f","#ffe6b7","#aadce0","#72bcd5","#528fad","#376795","#1e466e")))(50)) %>% as.ggplot())
write.table(annotation_row,file = "result/Fig6A.top_protein_family.cor.heatmap.part3.Gene_list.xls",quote = FALSE,sep = "\t",row.names=F, col.names=T)

p<-(p1|p2|p3)+plot_layout(widths = c(2.6,1.95,2.85))
ggsave("figure/Fig6A.top_protein_family.cor.heatmap.pdf",p,width = 22,height = 9)
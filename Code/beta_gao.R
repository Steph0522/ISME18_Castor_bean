#beta
library(readxl)
library(tidyverse)
library(betapart)
library(hilldiv)
library(qiime2R)
log_norm <- function(otu) {log(otu + 1)}
relabunda<- function(x){(as.data.frame(t(t(x)/colSums(x)))*100)}

meta<-read.delim(file = "Data//FINALMAP", check.names = F) %>% filter(Month==2) %>% filter(!Type_of_soil=="Uncultivated")

table_bac<- read_qza("Data/dada2.200.trimmed.single.filtered.cleaned.qza")$data %>% as.data.frame() %>% t() %>% as.data.frame() %>% 
  mutate_if(is.numeric,  ~ifelse(.x==0, 0, 1)) 
table_fung<- read_qza("Data/merge_table_240_noplant_filtered_nous.qza")$data %>% 
  as.data.frame()%>%  t() %>% as.data.frame() %>% 
  mutate_if(is.numeric,  ~ifelse(.x==0, 0, 1)) 


#tablas

fung<-table_fung%>% select_if(is.numeric) %>% 
  t() %>% relabunda() %>% log_norm() %>% t() %>% as.data.frame(
      ) %>% rownames_to_column(var="#SampleID") %>% inner_join(meta) %>% dplyr::select(
        -BarcodeSequence:-Description) %>% column_to_rownames(var = "#SampleID")
env<-table_fung%>% select_if(is.numeric) %>% 
  t() %>% relabunda() %>% log_norm() %>% t() %>% as.data.frame(
  ) %>% rownames_to_column(var="#SampleID") %>% inner_join(meta) %>% dplyr::select(
    `#SampleID`, BarcodeSequence:Description) %>% column_to_rownames(var = "#SampleID")

bact<-table_bac%>% select_if(is.numeric) %>% 
  t() %>% relabunda() %>% log_norm() %>% t() %>% as.data.frame(
  ) %>% rownames_to_column(var="#SampleID") %>% inner_join(meta) %>% dplyr::select(
    -BarcodeSequence:-Description) %>% column_to_rownames(var = "#SampleID")
envb<-table_bac%>% select_if(is.numeric) %>% 
  t() %>% relabunda() %>% log_norm() %>% t() %>% as.data.frame(
  ) %>% rownames_to_column(var="#SampleID") %>% inner_join(meta) %>% dplyr::select(
    `#SampleID`, BarcodeSequence:Description) %>% column_to_rownames(var = "#SampleID")


 #betapart
library(betapart)
library(vegan)
library(colorRamps)
fungb<- fung
bacb<- bact

#conver to incidence
fungb[fungb>0]=1 
bacb[bacb>0]=1 

#betapart function 
library(betapart)
library(vegan)
fd<-beta.pair(fungb, index.family = "jaccard")
fd1<-beta.pair(bacb, index.family = "jaccard")


fung.dist.jac<-fd$beta.jac
fung.dist.jac1<-fd1$beta.jac

fung.dist.jtu<-fd$beta.jtu
fung.dist.jtu1<-fd1$beta.jtu

fung.dist.jne<-fd$beta.jne
fung.dist.jne1<-fd1$beta.jne

jd<-betadisper(fung.dist.jac,factor(env$Type_of_soil))
jd1<-betadisper(fung.dist.jac1,factor(envb$Type_of_soil))

bd<-betadisper(fung.dist.jtu,factor(env$Type_of_soil))
bd1<-betadisper(fung.dist.jtu1,factor(envb$Type_of_soil))

nd<-betadisper(fung.dist.jne,factor(env$Type_of_soil))
nd1<-betadisper(fung.dist.jne1,factor(envb$Type_of_soil))

adj<- adonis2(fung.dist.jac~env$Type_of_soil)
adj1<- adonis2(fung.dist.jac1~envb$Type_of_soil)

adt<-adonis2(fung.dist.jtu~env$Type_of_soil)
adt1<-adonis2(fung.dist.jtu1~envb$Type_of_soil)

adn<-adonis2(fung.dist.jne~env$Type_of_soil)
adn1<-adonis2(fung.dist.jne1~envb$Type_of_soil)

 pdf(file="beta_gao_explo_bact_fung.pdf", width =6, height =6.5)
 library(viridis)
 pal<- viridis(6, option = "D") 
pal = c("#800000", "#808000", "#008000", "#D35400", "#2E4053")
 #plots
 par(mfrow=c(3,2),mar=c(2, 2, 0.5, 0.5))
 color=rgb(0,0,0,alpha=0.5) 
 
 plot(jd, col=pal, 
      hull = FALSE,cex = 1,label.cex = 0.8,
      seg.col=pal, main="", xlab="PC1", ylab="PC2", sub.caption="")
 text(-0.3, -0.3, sprintf("F=%.2f, P=%.3f", adj$F[1], adj$`Pr(>F)`[1]), col="Gray21", cex =1, font =2)
 title("Jaccard fungi", adj = 0.031, line = -1)
 
 plot(jd1, col=pal, 
      hull = FALSE,cex = 1,label.cex = 0.8,
      seg.col=pal, main="", xlab="PC1", ylab="PC2", sub.caption=" ")
 text(-.2, -0.25, sprintf("F=%.2f, P=%.3f", adj1$F[1], adj1$`Pr(>F)`[1]), col="Gray21", cex =1, font =2)
 title("Jaccard bacteria", adj = 0.031, line = -1)
 
 plot(bd, col=pal, 
      hull = FALSE,cex = 1,label.cex = 0.8,
      seg.col=pal, main="", xlab="PC1", ylab="PC2", sub.caption=" ")
 text(-.4, -0.3, sprintf("F=%.2f, P=%.3f", adt$F[1], adt$`Pr(>F)`[1]), col="Gray21", cex =1, font =2)
 title("Turnover fungi", adj = 0.031, line = -1)
 
 plot(bd1, col=pal, 
      hull = FALSE,cex = 1,label.cex = 0.8,
      seg.col=pal, main="", xlab="PC1", ylab="PC2", sub.caption=" ")
 text(-.4, -0.3, sprintf("F=%.2f, P=%.3f", adt1$F[1], adt1$`Pr(>F)`[1]), col="Gray21", cex =1, font =2)
 title("Turnover bacteria", adj = 0.031, line = -1)
 
 plot(nd, col=pal, 
      hull = FALSE,cex = 1,label.cex = 0.8,
      seg.col=pal, main="", xlab="PC1", ylab="PC2", sub.caption=" ")
 text(-0.55, -0.2, sprintf("F=%.2f, P=%.3f", adn$F[1], adn$`Pr(>F)`[1]), col="Gray21", cex =1, font =2)
 title("Nestedness fungi", adj = 0.031, line = -1)
 
 plot(nd1, col=pal, 
      hull = FALSE,cex = 1,label.cex = 0.8,
      seg.col=pal, main="", xlab="PC1", ylab="PC2", sub.caption=" ")
 text(-0.05, -0.1, sprintf("F=%.2f, P=%.3f", adn1$F[1], adn1$`Pr(>F)`[1]), col="Gray21", cex =1, font =2)
 title("Nestedness bacteria", adj = 0.031, line = -1)
 
 
 #plots

 dev.off()
 
 
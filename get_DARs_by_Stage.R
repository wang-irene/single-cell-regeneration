#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
obj<-readRDS(paste(args[1]))

# load packages
library(Signac,quietly = T)
library(Seurat, quietly = T)
library(GenomeInfoDb, quietly = T)
library(plyr, quietly = T)
library(dplyr, quietly = T)
library(tidyverse, quietly = T)
library(GenomicRanges, quietly = T)
library(openxlsx, quietly = T)
library(BSgenome.Drerio.UCSC.danRer11, quietly = T)
library(future, quietly = T)
set.seed(1234)
plan("multicore", workers = 6)
options(future.globals.maxSize = 120000 * 1024^2)

# function
GetDARs <- function(object, assay, ident1, ident2) {
    allDARs <- FindMarkers(object, 
                assay = assay, 
                ident.1 = ident1, 
                ident.2 = ident2, 
                logfc.threshold = 0.20,
                min.pct = 0.25,
                test.use = 'LR',
                latent.vars = paste("nFeature", assay, sep = "_")) %>% filter(p_val_adj < 0.05)
    posDARs <- allDARs %>% filter(avg_log2FC > 0.2)
    negDARs <- allDARs %>% filter(avg_log2FC < 0.2)
    return(list(All_DARs = allDARs, Closed_DARs = posDARs, Open_DARs = negDARs))
}

# DARs

print("Get DARs")

sample="WNNcelltype"
assay="peaks"
DefaultAssay(obj) <- assay

#obj <- RunTFIDF(obj) #TF-IDF normalization again in the subset data.
#obj <- FindTopFeatures(obj, min.cutoff = 20)
#annotation<-readRDS("/scratch/ichen/FinRegen_10xMultiome/annotations/Danio_rerio_GRCz11.104_Ensembl_gtf_annotation_onlyChr_wEGFP.rds")
#seqlevelsStyle(annotation) <- "UCSC"
#genome(annotation) <- "danRer11"
#Annotation(obj) <- annotation

Idents(obj)<-obj$Stage
idents=unique(Idents(obj))

## for loop to compare between stages
#for (i in 1:length(idents)){
#i=4
DAR_assay="peaks"
#celltype= idents[[i]]
#celltype="Mesenchymal_1"
print(paste("Find DARs for ",obj,"...",sep=""))

### prevent writing results from last loop to current loop
dar_0vs1<-list()
dar_0vs2<-list()
dar_0vs4<-list()
dar_0vs6<-list()
dar_1vs2<-list()
dar_1vs4<-list()
dar_1vs6<-list()
dar_2vs4<-list()
dar_2vs6<-list()
dar_4vs6<-list()


print(paste("Find DARs for ",obj," 0vs1...",sep=""))
tryCatch(dar_0vs1<-GetDARs(object=obj,assay=DAR_assay,ident1="mfin",ident2="1dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",obj," 0vs2...",sep=""))
tryCatch(dar_0vs2<-GetDARs(object=obj,assay=DAR_assay,ident1="mfin",ident2="2dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",obj," 0vs4...",sep=""))
tryCatch(dar_0vs4<-GetDARs(object=obj,assay=DAR_assay,ident1="mfin",ident2="4dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",obj," 0vs6...",sep=""))
tryCatch(dar_0vs6<-GetDARs(object=obj,assay=DAR_assay,ident1="mfin",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",obj," 1vs2...",sep=""))
tryCatch(dar_1vs2<-GetDARs(object=obj,assay=DAR_assay,ident1="1dpa",ident2="2dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",obj," 1vs4...",sep=""))
tryCatch(dar_1vs4<-GetDARs(object=obj,assay=DAR_assay,ident1="1dpa",ident2="4dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",obj," 1vs6...",sep=""))
tryCatch(dar_1vs6<-GetDARs(object=obj,assay=DAR_assay,ident1="1dpa",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",obj," 2vs4...",sep=""))
tryCatch(dar_2vs4<-GetDARs(object=obj,assay=DAR_assay,ident1="2dpa",ident2="4dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",obj," 2vs6...",sep=""))
tryCatch(dar_2vs6<-GetDARs(object=obj,assay=DAR_assay,ident1="2dpa",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

print(paste("Find DARs for ",obj," 4vs6...",sep=""))
tryCatch(dar_4vs6<-GetDARs(object=obj,assay=DAR_assay,ident1="4dpa",ident2="6dpa"), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

dar.list<-list()
dar.list<-list(
    dar_0vs1=dar_0vs1$All_DARs, dar_0vs1_close=dar_0vs1$Closed_DARs, dar_0vs1_open=dar_0vs1$Open_DARs,
    dar_0vs2=dar_0vs2$All_DARs, dar_0vs2_close=dar_0vs2$Closed_DARs, dar_0vs2_open=dar_0vs2$Open_DARs,
    dar_0vs4=dar_0vs4$All_DARs, dar_0vs4_close=dar_0vs4$Closed_DARs, dar_0vs4_open=dar_0vs4$Open_DARs,
    dar_0vs6=dar_0vs6$All_DARs, dar_0vs6_close=dar_0vs6$Closed_DARs, dar_0vs6_open=dar_0vs6$Open_DARs,
    dar_1vs2=dar_1vs2$All_DARs, dar_1vs2_close=dar_1vs2$Closed_DARs, dar_1vs2_open=dar_1vs2$Open_DARs, 
    dar_1vs4=dar_1vs4$All_DARs, dar_1vs4_close=dar_1vs4$Closed_DARs, dar_1vs4_open=dar_1vs4$Open_DARs, 
    dar_1vs6=dar_1vs6$All_DARs, dar_1vs6_close=dar_1vs6$Closed_DARs, dar_1vs6_open=dar_1vs6$Open_DARs, 
    dar_2vs4=dar_2vs4$All_DARs, dar_2vs4_close=dar_2vs4$Closed_DARs, dar_2vs4_open=dar_2vs4$Open_DARs, 
    dar_2vs6=dar_2vs6$All_DARs, dar_2vs6_close=dar_2vs6$Closed_DARs, dar_2vs6_open=dar_2vs6$Open_DARs, 
    dar_4vs6=dar_4vs6$All_DARs, dar_4vs6_close=dar_4vs6$Closed_DARs, dar_4vs6_open=dar_4vs6$Open_DARs
    )

require(openxlsx)
print(paste("write DARs to xlsx ",obj,"...",sep=""))
write.xlsx(dar.list, file = paste(sample,obj,"DARs_stage_comparison.xlsx",sep="_"), sheetName = names(dar.list), rowNames = F,overwrite=T)
stat <-lapply(dar.list, function(x) {print(nrow(x))}) %>% unlist() %>% data.frame() 
write.table(stat,file=paste(sample,obj,"DARs_stage_comparison_stats.txt",sep="_"),quote=F,sep="\t")

# Initiate for loop
for (i in c(10, 11, 12, 13, 14, 15)) { ## Replace w/ the dims you want to test
 
    dim_name = paste("Dim1-",i,sep="")
    dim=1:i
    hem <- RunUMAP(hem, assay="integrated", dims = dim) ## Replace w/ correct object 
    hem <- FindNeighbors(hem, dims = dim)
    hem <- FindClusters(hem, resolution = resolution_to_use)
 
 
# UMAP plots with different resolution (per each dim), as well as one umap colored by stage
    pdf(paste(sample,dim_name,"umap_rescombos.pdf",sep="_"), width=24, height=24)
 
    p0<-DimPlot(hem,reduction="umap",group.by = 'orig.ident', label = FALSE, pt.size = 0.1)+ 
      ggplot2::ggtitle(paste(sample,dim_name,"BySample",sep="_"))
    p1<-DimPlot(hem,reduction="umap",group.by = 'integrated_snn_res.1', label = TRUE, pt.size = 0.1)+ 
      ggplot2::ggtitle(paste(sample,dim_name,"Res1_ByCluster",sep="_"))
    p2<-DimPlot(hem, reduction="umap",group.by = 'integrated_snn_res.0.9', label = TRUE, pt.size = 0.1)+ 
      ggplot2::ggtitle(paste(sample,dim_name,"Res0.9_ByCluster",sep="_"))
    p3<-DimPlot(hem, reduction="umap",group.by = 'integrated_snn_res.0.8', label = TRUE, pt.size = 0.1)+ 
      ggplot2::ggtitle(paste(sample,dim_name,"Res0.8_ByCluster",sep="_"))
    p4<-DimPlot(hem, reduction="umap",group.by = 'integrated_snn_res.0.7', label = TRUE, pt.size = 0.1)+ 
      ggplot2::ggtitle(paste(sample,dim_name,"Res0.7_ByCluster",sep="_"))
    p5<-DimPlot(hem, reduction="umap",group.by = 'integrated_snn_res.0.6', label = TRUE, pt.size = 0.1)+ 
      ggplot2::ggtitle(paste(sample,dim_name,"Res0.6_ByCluster",sep="_"))
    p6<-DimPlot(hem, reduction="umap",group.by = 'integrated_snn_res.0.5', label = TRUE, pt.size = 0.1)+ 
      ggplot2::ggtitle(paste(sample,dim_name,"Res0.5_ByCluster",sep="_"))
    p7<-DimPlot(hem, reduction="umap",group.by = 'integrated_snn_res.0.4', label = TRUE, pt.size = 0.1)+ 
      ggplot2::ggtitle(paste(sample,dim_name,"Res0.4_ByCluster",sep="_"))
    p8<-DimPlot(hem, reduction="umap",group.by = 'integrated_snn_res.0.3', label = TRUE, pt.size = 0.1)+ 
      ggplot2::ggtitle(paste(sample,dim_name,"Res0.3_ByCluster",sep="_"))
    p9<-DimPlot(hem, reduction="umap",group.by = 'integrated_snn_res.0.2', label = TRUE, pt.size = 0.1)+ 
      ggplot2::ggtitle(paste(sample,dim_name,"Res0.2_ByCluster",sep="_"))
    p10<-DimPlot(hem, reduction="umap",group.by = 'integrated_snn_res.0.1', label = TRUE, pt.size = 0.1)+ 
      ggplot2::ggtitle(paste(sample,dim_name,"Res0.1_ByCluster",sep="_"))
 
    print(
      p0+p1+p2+p3+p4+p5+p6+p7+p8+p9+p10
    )
    dev.off()
  }

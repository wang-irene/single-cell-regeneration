# Load packages
library(JASPAR2020) ## Depends if you choose to use this as pfm, see below
library(TFBSTools)
library(BSgenome.Drerio.UCSC.danRer11)
library(patchwork)
library(motifmatchr)
set.seed(1234)

# Read in file, here I use the whole sc object, subsetted into atac
integrated_atac <- readRDS("/scratch/ichen/FinRegen_10xMultiome/Signac_processed/integration/rPCA_Seurat_integration/integrated_rPCA_pctMT10/redo2/integrated_dietseurat_atac_only.rds")
DefaultAssay(integrated_atac) <- "peaks"
Idents(integrated_atac) <- integrated_atac$celltype_wnn

# Add motif matrix -- get list of motif position frequency matrices from JASPAR database, PATH 1:
pfm <- getMatrixSet(x=JASPAR2020, opts=list(collection="CORE", tax_group='vertebrates', all_versions=F)) ## Our server's lab version only supports 2020

# Add motif matrix -- personal file in order to have updated motif names, PATH 2: 
pfm <- readRDS("/scratch/chemauer/STREAM/JASPAR_matrix.rds")

# Footprinting
integrated_atac <- AddMotifs(object = integrated_atac, assay="peaks", genome=BSgenome.Drerio.UCSC.danRer11, pfm=pfm)
motif_names <- c("MA0853.1", "MA0879.2", "MA1629.1", "MA0716.1", "MA0701.2", "MA0655.1", "MA1990.1", "MA1491.2", "MA1578.1", "MA0734.3", "MA0736.1")
integrated_atac <- Footprint(object=integrated_atac, motif.name=motif_names, genome=BSgenome.Drerio.UCSC.danRer11)
 
# Create plot
PlotFootprint(integrated_atac, features=motif_names)

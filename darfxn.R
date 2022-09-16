# Function to compile all 0 vs x days, where x is a post injury upregulated DAR day. Can be manipulated for DEGs as well.
darfxn <- function(dar_file, name) {
    dar_idents<-c("dar_0vs1_open", "dar_0vs2_open", "dar_0vs4_open", "dar_0vs6_open")
    dar.df.list <- lapply(dar_idents, function(x) read.xlsx(xlsxFile=dar_file, sheet=x, colNames=TRUE, cols=1))
    names(dar.df.list) <- dar_idents
    lapply(dar_idents, function(x) write.table(dar.df.list[[x]], paste(x, name, sep="_"), quote=FALSE, sep="-", row.names=FALSE, col.names=FALSE))
}

# Usage
wnn_names <- c("All_MES", "C12", "Pigment", "ALL_HEM", "BE", "Endothelial", "ALL_MUC", "C18", "SE", "MLC", "IE")
for (i in wnn_names) {
    file_path <- paste0("/scratch/iwang/celltype_dars/logfc0.25_minpct0.05_newparams/WNNcelltype_, i, "_DARs_stage_comparison.xlsx")
    darfxn(file_path, i) ## See yesterday's work for darfxn()
}

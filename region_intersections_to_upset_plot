# Read in individual strung, merged files
se <- read.csv("/scratch/iwang/celltype_dars/logfc0.25_minpct0.05_newparams/mergestuff/strung_merged_all_SE_0vsx", header=FALSE)
ie <- read.csv("/scratch/iwang/celltype_dars/logfc0.25_minpct0.05_newparams/mergestuff/strung_merged_all_IE_0vsx", header=FALSE)
be <- read.csv("/scratch/iwang/celltype_dars/logfc0.25_minpct0.05_newparams/mergestuff/strung_merged_all_BE_0vsx", header=FALSE)
ep_muc <- read.csv("/scratch/iwang/celltype_dars/logfc0.25_minpct0.05_newparams/mergestuff/strung_merged_all_muc_0vsx", header=FALSE)
mesen <- read.csv("/scratch/iwang/celltype_dars/logfc0.25_minpct0.05_newparams/mergestuff/strung_merged_all_mes_0vsx", header=FALSE)
hema <- read.csv("/scratch/iwang/celltype_dars/logfc0.25_minpct0.05_newparams/mergestuff/strung_merged_all_hem_0vsx", header=FALSE)

# Read in files into list
to_intersect_list <- list(Superficial_Epithelial=se$V1, Intermediate_Epithelial=ie$V1, Basal_Epithelial=be$V1, Epidermal_Mucous=ep_muc$V1, Mesenchymal=mesen$V1, Hematopoietic=hema$V1)

# Create upset plot
pdf("celltype_intersections_complexupsetplot.pdf",width=22, height=10)
upset(
    fromList(to_intersect_list),
    colnames(fromList(to_intersect_list))[1:6],
    name='Cell Type',
    width_ratio=0.12,
    #min_size=10, ## Made a second plot that filtered for min number of intersections
    sort_intersections_by='degree', 
    set_sizes=(upset_set_size() + theme(axis.ticks.x=element_line(), axis.text.x=element_text(size=14, angle=20))),
    height_ratio=0.6,
    base_annotations=list('Intersection size'=intersection_size(text = list(size=4.2), bar_number_threshold=1)),
    matrix=(
      intersection_matrix(geom=geom_point(shape='circle filled', size=5.2))
      + scale_color_manual(
          values=c('Superficial_Epithelial'="#8B2DB2", 'Intermediate_Epithelial'="#CE6DBD",
          'Basal_Epithelial'="#9C9EDE",'Epidermal_Mucous'="#F28E2BFF",'Hematopoietic'="#59A14FFF",'Mesenchymal'="#B6992DFF"),
      )
    ),
    queries=list(
      upset_query(set='Superficial_Epithelial', fill="#8B2DB2"),
      upset_query(set='Intermediate_Epithelial', fill="#CE6DBD"),
      upset_query(set='Basal_Epithelial', fill="#9C9EDE"),
      upset_query(set='Epidermal_Mucous', fill="#F28E2BFF"),
      upset_query(set='Hematopoietic', fill="#59A14FFF"),
      upset_query(set='Mesenchymal', fill="#B6992DFF")
    ),
    themes=upset_default_themes(text=element_text(size=22))
  )
dev.off()

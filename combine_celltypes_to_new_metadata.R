integrated_atac[["celltype_wnn_merge"]] <- mapvalues(integrated_atac$celltype_wnn, ## "celltype_wnn_merge" is new metadata column name, "celltype_wnn" is what you are extracting from
  from = c("Superficial_Epithelial", "Intermediate_Epithelial", "Basal_Epithelial", "Epidermal_Mucous_1", "Epidermal_Mucous_2", "Hematopoietic_1", "Hematopoietic_2", "Mesenchymal_1", "Mesenchymal_2", "Metaphocyte", "Pigment", "Endothelial", "C12", "C18"), 
  to = c("Superficial_Epithelial",
    "Intermediate_Epithelial",
    "Basal_Epithelial",
    "Epidermal_Mucous", "Epidermal_Mucous",
    "Hematopoietic", "Hematopoietic",
    "Mesenchymal","Mesenchymal",
    "Metaphocyte",
    "Pigment",
    "Endothelial",
    "C12", 
    "C18"
    )) %>% factor( 
    levels = c("Superficial_Epithelial",
    "Intermediate_Epithelial",  
    "Basal_Epithelial",
    "Epidermal_Mucous",
    "Hematopoietic",
    "Mesenchymal",
    "Pigment",
    "Endothelial",
    "Metaphocyte",
    "C12","C18"
    ))

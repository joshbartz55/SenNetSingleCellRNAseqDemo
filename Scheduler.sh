#!/bin/bash

Rscript - <<'EOF'
##############################################################################
#------------------------------ Load Libraries ------------------------------#
##############################################################################

# Make sure any required libraries are loaded
# library(Seurat)
# library(yourGCLpackage)

##############################################################################
#---------------------------------- Load Data -------------------------------#
##############################################################################
load('/projects/standard/dong0265/bartz209/SenNetFigure/data/integrated_data.RDATA')
integrated_data$predicted.celltype = gsub(' ','_',integrated_data$predicted.celltype)
cluster_names = gsub(' ','_',cluster_names)
##############################################################################
#-------------------------------- Calculate GCL ----------------------------#
##############################################################################

# Loop over clusters
for(i in 1:length(cluster_names)){
  cell_type <- cluster_names[[i]]
  ct_subset <- subset(integrated_data, subset = predicted.celltype == cell_type)
  print(cell_type)
  
  sample_ids <- unique(ct_subset$sample_id)
  
  for(j in 1:length(sample_ids)){
    sample <- sample_ids[[j]]
    print(sample)
    
    sample_ct_subset <- subset(ct_subset, subset = sample_id == sample)
    
    filename <- paste0('/projects/standard/dong0265/bartz209/SenNetFigure/data/',
                       cell_type,'_',sample,'.Rdata')
    
    save(sample_ct_subset, file = filename)
    
    # Run setGCL function
    cmd <- paste("sbatch SubGCL.sh", cell_type, sample)
    system(cmd)
  }
}
EOF


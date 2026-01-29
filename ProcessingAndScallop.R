library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)
library(SeuratObject)
library(SeuratDisk)
library(purrr)
library(foreach)
library(rdist)
library(doParallel)
library(dplyr)
########################################################################################################################
#---------------------------------------------------FUNCTION DEFINITIONS-----------------------------------------------#
########################################################################################################################

setScallop = function(seurat_object){
  filename = "seurat_obj.h5seurat"
  SaveH5Seurat(seurat_object, filename = filename, overwrite = T)
  Convert(filename, dest = "h5ad", overwrite = TRUE)
  filename = gsub('seurat_obj.h5seurat','seurat_obj.h5ad',filename)
  cmd = sprintf(
    "/opt/homebrew/bin/python3.11 'Scallop_Runner.py' '%s'",
    filename
  )
  system(cmd, wait = TRUE)
  scallop_scores = read.csv(gsub('seurat_obj.h5ad','scallop_scores.csv',filename))
  seurat_object@misc$Scallop = scallop_scores
  return(seurat_object) 
}x

addUniqueName = function(name, name_list, i = 1, modded = FALSE){
  count = i + 1
  if(name %in% name_list){
    if(modded){
      name = substr(name,1,nchar(name)-2)
    }
    name_list = addUniqueName(paste(name, count), name_list, i = count, modded = TRUE)
  }
  else{
    name_list = c(name_list, name)
  }
  return(name_list)
}

assignID = function(cluster, marker_map){
  ret = cbind(marker_map, matches = rep(0,nrow(marker_map)))
  
  for(i in 1:nrow(cluster)){
    marker = cluster[i,]
    marker_name = marker$gene
    for(i in 1:nrow(ret)){
      entry = ret[i,]
      marker_list = unlist(strsplit(entry$markerEnsemblIDs, split=" "))
      if(marker_name %in% marker_list){
        ret[i,]$matches = ret[i,]$matches + 1
      }
    }
  }
  return(ret)
}


##############################################################################
#---------------------------------- Load Data -------------------------------#
##############################################################################
h5ad_files = list.files('Research/SenNetPortalProject/data', pattern = '*.h5ad', full.names = T)

all_data <- lapply(h5ad_files, function(f) {
  sce <- readH5AD(f)
  
  rowData(sce)$hugo_symbol <- NULL
  
  obj <- as.Seurat(
    sce,
    counts = "X",
    data = NULL
  )

  return(obj)
})

all_data[[1]]$sample_id = 1
all_data[[2]]$sample_id = 2
all_data[[3]]$sample_id = 3

#############################################################################
#--------------------------------- CLeaning --------------------------------#
#############################################################################
cleaned_data = lapply(X = all_data, FUN = function(x){
  x[["RNA"]] = x[["originalexp"]]   # copy assay
  DefaultAssay(x) ="RNA"       
  x[["originalexp"]] = NULL          # remove old assay
  
  feat_quant = quantile(x@meta.data[["nFeature_RNA"]], prob=0.95)
  selected_c = WhichCells(x, expression = (nFeature_RNA > 200 & nFeature_RNA < feat_quant[[1]]))
  selected_f = rownames(x)
  filtered_data = subset(x,features = selected_f, cells = selected_c)
  filtered_data[["percent.mt"]] = PercentageFeatureSet(filtered_data, pattern = "^MT-")
  filtered_data = subset(filtered_data, percent.mt < 5)
  normalizedData = NormalizeData(filtered_data)
  print(ncol(x))
  x = FindVariableFeatures(normalizedData, selection.method = "vst", nfeatures = 2000)
})

##############################################################################
#-------------------------------- Integration -------------------------------#
##############################################################################
features = SelectIntegrationFeatures(object.list = cleaned_data)
anchors = FindIntegrationAnchors(object.list = cleaned_data, anchor.features = features, k.filter = 50)
integrated_data = IntegrateData(anchorset = anchors, k.weight = 30)
DefaultAssay(integrated_data) = "integrated"
##############################################################################
#--------------------------------- Clustering -------------------------------#
##############################################################################
integrated_data = ScaleData(integrated_data, verbose = FALSE)
integrated_data = RunPCA(integrated_data, npcs = 30, verbose = FALSE)
integrated_data = RunUMAP(integrated_data, reduction = "pca", dims = 1:30)
integrated_data = FindNeighbors(integrated_data, reduction = "pca", dims = 1:30)
integrated_data = FindClusters(integrated_data, resolution = 0.5)
##############################################################################
#------------------------- Cell Type Identification -------------------------#
##############################################################################
DefaultAssay(integrated_data) = "RNA"
found_markers = FindAllMarkers(integrated_data)
found_markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
found_markers$gene = sub("\\..*$", "", found_markers$gene)

all_markers = read.csv('/Users/joshbartz/Desktop/Research/DBA/ReferenceData/cellMarkers.csv')
#Separate and format markers dictionary for specified tissue type
tissue_markers <- all_markers[all_markers$tissueType=='Lung',]
# Return top num_markers for cluster specified 'x'
gen_marker_table = function(x){
  found_markers[found_markers$cluster == x, ] %>%
    head(n=100)
}
# Create a data frame of results for clusters
max_cluster = length(unique(found_markers$cluster)) -1
top_markers <- map_dfr(0:max_cluster, gen_marker_table)
#idenyify cell types of clusters
cluster_names = c()
cluster2_names = c()
cluster3_names = c()
for(i in 0:max_cluster){
  cluster_markers = top_markers[top_markers$cluster == i,]
  id_matrix = assignID(cluster_markers, tissue_markers)
  id_matrix = id_matrix[order(id_matrix$matches, decreasing = TRUE), ]
  cluster_id  = if(id_matrix[1,]$matches != 0) id_matrix[1,]$cellName else 'Unkown'
  cluster_id2 = if(id_matrix[2,]$matches != 0) id_matrix[2,]$cellName else 'NA'
  if(nrow(id_matrix) >= 3){
    cluster_id3 = if(id_matrix[3,]$matches != 0) id_matrix[3,]$cellName else 'NA'
  }else{
    cluster_id3 = 'NA' 
  }
  cluster_names = addUniqueName(cluster_id, cluster_names)
  cluster2_names = c(cluster2_names, cluster_id2)
  cluster3_names = c(cluster3_names, cluster_id3)
}
#apply assigned cell types to clusters
names(cluster_names) <- levels(integrated_data)
integrated_data <- RenameIdents(integrated_data, cluster_names)
umap = DimPlot(integrated_data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
integrated_data$predicted.celltype = Idents(integrated_data)
save(integrated_data, cluster_names, file = 'Research/SenNetPortalProject/data/integrated_data.RDATA')
##############################################################################
#------------------------------- Calculate Scallop --------------------------#
##############################################################################
integrated_data$predicted.celltype = gsub(' ','_',integrated_data$predicted.celltype)
cluster_names = gsub(' ','_',cluster_names)

for(i in 1:length(cluster_names)){
  cell_type = cluster_names[[i]]
  ct_subset = subset(integrated_data, subset = predicted.celltype == cell_type)
  print(cell_type)
  sample_ids = unique(ct_subset$sample_id)
  for(j in 1:length(sample_ids)){
    sample = sample_ids[[j]]
    print(sample)
    sample_ct_subset = subset(ct_subset, subset = sample_id == sample)
    sample_ct_subset = setScallop(sample_ct_subset)
    write.csv(sample_ct_subset@misc$Scallop, file = paste0('Research/SenNetPortalProject/ScallopResults/',cell_type,'_',sample,'_Scallop.csv'))
  }
}





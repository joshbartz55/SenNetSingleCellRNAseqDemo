library(Seurat)
library(purrr)
library(foreach)
library(rdist)
library(doParallel)
library(dplyr)
########################################################################################################################
#---------------------------------------------------FUNCTION DEFINITIONS-----------------------------------------------#
########################################################################################################################


Vn = function(matrix1, matrix2, n){
  return((1/(n*(n-3)))*(sum(matrix1*matrix2) - (n/(n-2))*sum(diag(matrix1)*diag(matrix2))))
}

Rn = function(matrix1, matrix2, n){
  return(Vn(matrix1, matrix2,n)/sqrt(Vn(matrix1,matrix1,n)*Vn(matrix2,matrix2,n)))
}

calculateGCL = function(popData, num_divisions=50){
  popData=matrix(as.numeric(unlist(popData)), nrow = nrow(popData))
  num_genes = nrow(popData)
  num_cells = ncol(popData)
  
  num_cores = 10
  cl = makeCluster(num_cores)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("Vn", "Rn"))
  
  gcl = foreach(i=1:num_divisions, .combine = 'c', .packages = "rdist") %dopar%{
    #Split genes into two random groups
    random_permutation = sample(num_genes)
    G1 = t.default(popData[random_permutation[1:floor(num_genes/2)],])
    G2 = t.default(popData[random_permutation[(floor(num_genes/2)+1):num_genes],])
    
    #Create and Populate Matrices
    d1 = cdist(G1,G1)
    mean1 = colMeans(d1)
    elementMeans1 = mean(d1) 
    Aij = d1 - as.matrix(mean1)%*%matrix(1,1,num_cells)
    Aij = Aij - matrix(1,num_cells,1)%*%mean1
    Aij = Aij + elementMeans1
    Aij = Aij -d1/num_cells
    Aij[seq(1,length(Aij), num_cells+1)] = mean1 -elementMeans1
    Aij = (num_cells/(num_cells-1))* Aij
    
    #Repeat for second split
    d2 = cdist(G2,G2)
    mean2 = colMeans(d2)
    elementMeans2 = mean(d2) 
    Bij = d2 - as.matrix(mean2)%*%matrix(1,1,num_cells)
    Bij = Bij - matrix(1,num_cells,1)%*%mean2
    Bij = Bij + elementMeans2
    Bij = Bij - d2/num_cells
    Bij[seq(1,length(Bij), num_cells+1)] = mean1 -elementMeans1
    Bij = (num_cells/(num_cells-1))* Bij
    
    rn = Rn(Bij,Aij,num_cells)
    return(rn)
  }
  stopCluster(cl)
  return(gcl)
}

setGCL = function(seurat_object){
  dataset = as.matrix(seurat_object@assays$RNA@counts)
  gcl_values = calculateGCL(dataset)
  seurat_object@misc$GCL = gcl_values
  return(seurat_object)
}
########################################################################################################################
#-------------------------------------------------------CALCULATE GCL--------------------------------------------------#
########################################################################################################################

args = commandArgs(trailingOnly = TRUE)

cell_type = args[[1]]
sample_id = args[[2]]
filename = paste0('/projects/standard/dong0265/bartz209/SenNetFigure/data/',cell_type,'_',sample_id,'.Rdata')
print(paste0('Starting Analysis of ', filename))
load(filename)
outfile = paste0('/projects/standard/dong0265/bartz209/SenNetFigure/results/',cell_type,'_',sample_id,'.csv')
res = setGCL(sample_ct_subset)
write.csv(res@misc$GCL, outfile)





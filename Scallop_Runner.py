import sys 
print(f"Running Python {sys.version}")
sys.stdout.flush()

import os
import pandas as pd
import numpy as np
import anndata as ad
from sklearn.preprocessing import MinMaxScaler, RobustScaler
from scipy import sparse
# Single-cell
import scanpy as sc
import scanpy.external as sce
import anndata 
import scallop as sl

if len(sys.argv) < 2:
    print("Usage: python3 Scallop_Runner.py <path_to_h5ad>")
    sys.exit(1)
  
adata_path = sys.argv[1]   
outdir = adata_path.replace('seurat_obj.h5ad', 'scallop_scores.csv')
print(f"Reading AnnData from: {adata_path}") 
adata = ad.read_h5ad("/Users/joshbartz/Desktop/seurat_obj.h5ad")

scalliwag_sample = sl.Scallop(adata)
res_df = sl.tl.getScore(scalliwag_sample, res=1.2, n_trials=30, frac_cells=0.95, do_return=True)
res_df.to_csv(outdir)
print(f"Saved Scallop scores to {outdir}")
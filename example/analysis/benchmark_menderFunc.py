# IMPORTANT: requires python 3.10
import warnings
warnings.filterwarnings("ignore")
import time
# from pathlib import Path
# from operator import itemgetter
import scanpy as sc
import pandas as pd
import numpy as np
# from joblib import Parallel, delayed
from sklearn.metrics import *
import MENDER


def load_sample(file_list):
    # function to load anndata files 
    # and add them to a list
    adata_list = []
    for file in file_list:
        print(file)
        adata_list.append(sc.read_h5ad(file))
    return adata_list


def nscluster_batch(in_adata, batch):
    # function to perform Leiden clustering
    # on batch-corrected samples in anndata file
    if batch == 'True':
        sc.pp.neighbors(in_adata, use_rep = 'X_pca_harmony')
        sc.tl.umap(in_adata)
        print(in_adata)
    elif batch == 'False':
        sc.pp.neighbors(in_adata, use_rep = 'X_pca')
        sc.tl.umap(in_adata)
        print(in_adata)
    # Leiden clustering
    sc.tl.leiden(in_adata, resolution = 2, key_added = 'ct')
    in_adata.obs['ct'] = in_adata.obs['ct'].astype('category')
    return in_adata

def mender_run(file_list, scale, mode, radius, nSeed, batch, msm_res):
    # function to perform mender
    adata_list = load_sample(file_list)
    print("Data loaded.")
    # merge all anndata files
    adata_raw = adata_list[0].concatenate(adata_list[1:])
    adata_raw.obs['slice_id'] = adata_raw.obs['sampleID'].astype('category')
    # perform Leiden clustering
    adata = nscluster_batch(adata_raw.copy(), batch)
    print("Initial clustering performed.")
    # UMAP visualisation
    # sc.pl.umap(adata, color = ["domainAnnotations", "ct", "sampleID"], size = 6)
    # can also use cell annotations in place of Leiden clusters
    msm = MENDER.MENDER(adata, batch_obs = 'sampleID', ct_obs = 'ct', random_seed = int(nSeed))
    print("Mender object created.")
    # set the MENDER parameters
    msm.prepare()
    print("Mender run prepared.")
    msm.set_MENDER_para(n_scales = int(scale), nn_mode = mode, nn_para = int(radius))
    print("Mender parameters set.")
    # construct the context representation
    # the number of processings
    # msm.run_representation_mp()
    msm.run_representation_mp(mp = 1)
    print("Mender context representation constructed.")
    # set the spatial clustering parameter
    msm.run_clustering_normal(msm_res)
    print("Mender clustering performed.")
    return msm.adata_MENDER
    

# -*- coding: utf-8 -*-
# @Author: dongqing
# @Date:   2023-04-13 22:38:41
# @Last Modified by:   dongqing
# @Last Modified time: 2023-04-14 21:52:36



import os
import scipy
import gensim

import numpy as np
import pandas as pd
import argparse as ap

from gensim.models import LdaModel
from sklearn.preprocessing import StandardScaler

from STRIDE.ModelTrain import scProcess, stProcess, scLDA, ModelRetrieve
from STRIDE.Deconvolution import SpatialDeconvolve

scrna_dir = "Data/scRNA"
st_dir = "Data/ST_simulation/FigS1_simulation"
res_dir = "Result/FigS1_simulation/STRIDE"

st_count_file = os.path.join(st_dir, "BRCA_EMTAB8107_simulation_20211109_count_gene_count.h5")
sc_scale_factor = 20000
st_scale_factor = 20000
out_dir = "Simulation_20211109_20230413"
out_prefix = "Simulation_20211109"
normalize = True
model_dir = os.path.join(res_dir, "Marker_gene_overlapped_rm_stride1b/model")

print("Reading spatial count matrix...")
st_info = stProcess(st_count_file, st_scale_factor)
st_count_mat = st_info["scale_matrix"]
st_count_genes = st_info["genes"]
st_count_spots = st_info["spots"]
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

print("Retrieving the pre-trained model...")
modelretrieve_res = ModelRetrieve(model_dir, st_count_genes)
genes_dict = modelretrieve_res["genes_dict"]
print("Deconvolving spatial transcriptomics...")
ntopics_selected = 29
model_selected = "BayesNorm"
spot_celltype_array_norm_df = SpatialDeconvolve(st_count_mat, st_count_genes, st_count_spots, genes_dict, model_selected, ntopics_selected, normalize, out_dir, model_dir, out_prefix)

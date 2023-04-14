# -*- coding: utf-8 -*-
# @Author: dongqing
# @Date:   2023-04-13 22:25:13
# @Last Modified by:   dongqing
# @Last Modified time: 2023-04-14 21:38:11


import os
import scipy
import gensim

import numpy as np
import pandas as pd
import argparse as ap

from gensim.corpora.dictionary import Dictionary

from sklearn.preprocessing import StandardScaler

from STRIDE.ModelTrain import scProcess, stProcess, scLDA, ModelRetrieve
from STRIDE.Deconvolution import SpatialDeconvolve
from STRIDE.utility.IO import read_10X_h5, read_count, write_10X_h5

scrna_dir = "Data/scRNA"
st_dir = "Data/ST_simulation/Fig2_simulation"
res_dir = "Result/Fig2_simulation/STRIDE"
st_count_file = os.path.join(st_dir, "BRCA_EMTAB8107_sim_scaled_2e4_gene_count.h5")

out_dir = "Simulation_1221_20230414"
out_prefix = "Simulation_1221"
normalize = True
model_dir = os.path.join(res_dir, "Marker_gene_scaled_norm/model")
gene_use = os.path.join(scrna_dir, "BRCA_EMTAB8107_markers_top500.txt")
genes_dict_file = os.path.join(model_dir, "../Gene_dict.txt")

print("Reading spatial count matrix...")
st_count = read_10X_h5(st_count_file)
st_count_mat = st_count.matrix
st_count_genes = st_count.names.tolist()
st_count_genes = [i.decode() for i in st_count_genes]
st_count_spots = st_count.barcodes.tolist()
st_count_spots = [i.decode() for i in st_count_spots]

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

genes_dict = Dictionary.load_from_text(genes_dict_file)
ntopics_selected = 28
model_selected = "BayesNorm"
spot_celltype_array_norm_df = SpatialDeconvolve(st_count_mat, st_count_genes, st_count_spots, genes_dict, model_selected, ntopics_selected, normalize, out_dir, model_dir, out_prefix)

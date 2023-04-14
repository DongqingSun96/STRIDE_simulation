# STRIDE_simulation
This repository contains all the results and codes to generate the performance evaluation plots in main Figure 2 (Fig 2D-G) and supplementary Figure S1 (Fig S1D-G) of [STRIDE paper](https://academic.oup.com/nar/article/50/7/e42/6543547).

*To run the codes, please clone the repository and download the [Data](https://drive.google.com/drive/folders/1rQ4IsmHRmLflsmEKHjEvq5rq8t4PQwBZ?usp=share_link) folder to the `STRIDE_simulation` filder.*

## About files

### Data: three simulated ST datasets and the reference scRNA-seq
* [Data/ST_simulation](Data/ST_simulation) contains three simulated ST datasets used in STRIDE paper to compare STRIDE's deconvolution performance with other tools.
  Folder | Explanation
  --- | --- 
  `Fig2_simulation` | Random simulation, named as `ST1` (Fig 2D-F and Fig S1D).
  `FigS1_simulation` | Complex TME, named as `ST2` (Fig S1E-G).
  `Fig2_seqdepth_simulation` | Benchmark different sequencing depths, named as `ST3` (Fig2G).

  Generally, each folder contains three files.
  File | Explanation
  --- | --- 
  `*_gene_count.h5` | Spot-level gene expression file, with genes as rows and spots as columns.
  `*.data.rds` | The simulated ST data, in which `$cell_composition` stores cell compositions for each simulated spot and `$topic_profiles` stores merged gene expression for each spot.
  `*_res.rds` | Processed `Seurat` object for the simulated ST data.

* [Data/scRNA](Data/scRNA) contains the scRNA-seq reference used to simulate ST data.

  File | Explanation
  --- | --- 
  `BRCA_EMTAB8107_celltype.txt` | Cell-type annotation file all the cells included in `BRCA_EMTAB8107_count_gene_count.h5`.
  `BRCA_EMTAB8107_count_gene_count.h5` | Single-cell gene expression file, with genes as rows and cells as columns.
  `BRCA_EMTAB8107_markers_top500.txt` | Gene list used to run STRIDE deconvolution on `ST1`.
  `BRCA_EMTAB8107_markers_top500_overlapped_rm.txt` | Gene list used to run STRIDE deconvolution on `ST2` and `ST3`.
  `BRCA_EMTAB8107_res.rds` | Processed `Seurat` object for the reference scRNA-seq data.


### Result: deconvolution results by different algorithms to re-generate the performance comparison plots
* [Result](Result) contains deconvolution results by different algorithms to re-generate the performance comparison plots.
  
  Folder | Explanation
  --- | --- 
  `Fig2_simulation` | `ST1` (Fig 2D-F and Fig S1D)
  `FigS1_simulation` | `ST2` (Fig S1E-G)
  `Fig2_seqdepth_simulation` | `ST3` (Fig2G)

  The trained models to re-run STRIDE on the three simulated datasets are also included.

  Folder | Explanation
  --- | --- 
  `Fig2_simulation/STRIDE/Marker_gene_scaled_norm` | Trained STRIDE model for `ST1`.
  `FigS1_simulation/STRIDE/Marker_gene_overlapped_rm_stride1b` | Trained STRIDE model for `ST2` and `ST3`.

### Code: codes to re-run STRIDE deconvolution with provided models and to re-generate the performance comparison plots
* [Code/Evaluation](Code/Evaluation) contains the codes to re-generate the performance evaluation plots from deconvolution results by different algorithms.
  
  File | Explanation
  --- | --- 
  `Fig2_simulation.R` | Fig 2D-F and Fig S1D
  `FigS1_simulation.R` | Fig S1E-G
  `Fig2_seqdepth_simulation.R` | Fig2G

* [Code/STRIDE](Code/STRIDE) contains the codes to re-run STRIDE on the three simulated datasets.

  File | Explanation
  --- | --- 
  `Fig2_simulation_STRIDE.py` | Run STRIDE deconvolution on `ST1`.
  `FigS1_simulation_STRIDE.py` | Run STRIDE deconvolution on `ST2`.
  `Fig2_seqdepth_simulation_STRIDE.sh` | Run STRIDE deconvolution on `ST3`.

  *Note: To run STRIDE, please install `STRIDE 0.0.2b0` by `pip install stridespatial`.*


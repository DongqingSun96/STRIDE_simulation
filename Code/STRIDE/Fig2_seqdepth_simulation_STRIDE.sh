# @Author: dongqing
# @Date:   2023-04-14 21:43:47
# @Last Modified by:   dongqing
# @Last Modified time: 2023-04-14 21:57:26


samples=("1000" "2500" "5000" "10000" "15000" "20000")
for sample in ${samples[*]}; do
    echo $sample
    STRIDE deconvolve --sc-count Data/scRNA/BRCA_EMTAB8107_count_gene_count.h5 \
    --sc-celltype Data/scRNA/BRCA_EMTAB8107_celltype.txt \
    --st-count Data/ST_simulation/Fig2_seqdepth_simulation/${sample}/BRCA_EMTAB8107_simulation_20211110_seqdepth_${sample}_count_gene_count.h5 \
    --model-dir Result/FigS1_simulation/STRIDE/Marker_gene_overlapped_rm_stride1b/model \
    --gene-use Data/scRNA/BRCA_EMTAB8107_markers_top500_overlapped_rm.txt \
    --outdir Simulation_20211110_seqdepth_20230413/${sample}/ --outprefix Simulation_20211110_seqdepth_${sample} --normalize;
done

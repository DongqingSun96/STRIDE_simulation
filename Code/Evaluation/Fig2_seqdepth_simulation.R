library(ggplot2)
library(ggsci)

setwd("STRIDE_simulation")

# compare STRIDE with simulation
# cell-level correlation 
{
  data_prefix_dir = "Data/ST_simulation/Fig2_seqdepth_simulation"
  res_prefix_dir = "Result/Fig2_seqdepth_simulation"
  
  deconv_comp_df_list = list()
  for (reads in c(1000, 2500, 5000, 10000, 15000, 20000)){
    print(reads)
    sim.res = readRDS(file.path(data_prefix_dir, as.character(reads), paste0("BRCA_EMTAB8107_simulation_20211110_seqdepth_", reads,".rds")))
    sim.cell.comp = sim.res$cell_composition
    sim.cell.comp = sim.cell.comp[,sort(names(sim.cell.comp))]
    sim.cell.comp.frac = sim.cell.comp/rowSums(sim.cell.comp)
  
    NMF_deconv = read.table(file.path(res_prefix_dir, reads, "NMFreg", paste0("Simulation_20211110_", reads, "_NMFreg_deconv.txt")), 
                            header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
    colnames(NMF_deconv) = colnames(sim.cell.comp.frac)
    corr_nmfreg = cor(t(NMF_deconv), t(sim.cell.comp.frac))
    corr_nmfreg = diag(corr_nmfreg)
    median(corr_nmfreg, na.rm = TRUE)
    
    celltype.predictions.score = read.table(file.path(res_prefix_dir, as.character(reads), "CCA", paste0("BRCA_EMTAB8107_simulation_20211110_", reads,"_celltype.predictions.score.txt")), header = TRUE, row.names = 1, sep = "\t")
    corr_cca = cor(t(celltype.predictions.score), t(sim.cell.comp.frac))
    corr_cca = diag(corr_cca)
    median(corr_cca, na.rm = TRUE)
    
    spotlight_ls = readRDS(file.path(res_prefix_dir, as.character(reads), "SPOTlight",  "spotlight_res.rds"))
    corr_spotlight = cor(t(spotlight_ls[[2]][,-ncol(spotlight_ls[[2]])]), t(sim.cell.comp.frac))
    corr_spotlight = diag(corr_spotlight)
    median(corr_spotlight, na.rm = TRUE)
    
    RCTD_res = read.table(file.path(res_prefix_dir, as.character(reads), "RCTD", paste0("BRCA_EMTAB8107_simulation_20211110_", reads, "_RCTD.txt")), header = TRUE, row.names = 1, sep = "\t")
    corr_rctd = cor(t(RCTD_res), t(sim.cell.comp.frac))
    corr_rctd = diag(corr_rctd)
    median(corr_rctd, na.rm = TRUE)
    
    STRIDE_deconv = read.table(file.path(res_prefix_dir, as.character(reads), "STRIDE", paste0("Simulation_20211110_", reads, "_spot_celltype_frac.txt")), sep = "\t", header = TRUE, row.names = 1)
    corr_STRIDE = cor(t(as.matrix(STRIDE_deconv)), t(sim.cell.comp.frac))
    corr_STRIDE = diag(corr_STRIDE)
    median(corr_STRIDE, na.rm = TRUE)
    
    cell2location_deconv = read.table(file.path(res_prefix_dir, as.character(reads),"Cell2location", "W_mRNA_count.csv"),
                                      header = TRUE, row.names = 1, check.names = FALSE, sep = ",")
    
    colnames(cell2location_deconv) = gsub("mean_nUMI_factors", "", colnames(cell2location_deconv))
    cell2location_deconv = cell2location_deconv/rowSums(cell2location_deconv)
    corr_cell2location = cor(t(cell2location_deconv), t(sim.cell.comp.frac))
    corr_cell2location = diag(corr_cell2location)
    median(corr_cell2location, na.rm = TRUE)
    
    deconv_other_df = data.frame(Corr = c(corr_spotlight, corr_nmfreg, corr_cca, corr_rctd, corr_cell2location),
                                 Method = rep(c("SPOTlight", "NMFreg", "CCA", "RCTD", "Cell2location"), each = length(corr_spotlight)))
    lda_res = data.frame(Corr = corr_STRIDE,
                         Method = "STRIDE")
    deconv_comp_df = rbind(lda_res, deconv_other_df)
    deconv_comp_df$Method = factor(deconv_comp_df$Method, 
                                   levels = c("STRIDE","RCTD","Cell2location",  "NMFreg","SPOTlight", "CCA"))
    
    deconv_comp_df$Spot = paste("spot", 1:nrow(sim.cell.comp.frac), sep = "_")
    deconv_comp_df$Counts = reads
    # rownames(sim.cell.comp.frac) = paste("spot", 1:nrow(sim.cell.comp.frac), sep = "_")
    deconv_comp_df_list[[as.character(reads)]] = deconv_comp_df
  }
  deconv_comp_df_rbind = do.call(rbind, deconv_comp_df_list)
  deconv_comp_df_rbind$Counts = factor(deconv_comp_df_rbind$Counts, unique(deconv_comp_df_rbind$Counts))
  data_prefix = "Simulation_20211110_seqdepth"
  p = ggplot(deconv_comp_df_rbind, aes(x = Counts, y = Corr, fill = Method)) + 
    geom_boxplot(outlier.size = 0.1, weight = 0.5) +
    scale_fill_manual(values = pal_simpsons("springfield")(16)[c(2,1,3:4,6:7)]) +
    # scale_fill_simpsons(palette = c("springfield"), alpha = 0.8) + 
    # scale_fill_jco() + 
    facet_grid(cols = vars(Method)) + 
    labs(x = "Sequencing depth",
         y = "Pearson's correlation", 
         title = "Spot-level correlation") + 
    scale_y_continuous(limits = c(0,1)) + 
    theme_bw(base_size = 14) +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  ggsave(paste(data_prefix, "corr_compare_methods.pdf", sep = "_"), p, width = 9, height = 4.55)
}

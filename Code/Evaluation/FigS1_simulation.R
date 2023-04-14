library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggsci)

setwd("STRIDE_simulation")

{

  # comparison
  {
    # Fig S1E cell-level correlation (choose the topic number and compare differnt methods)
    {
      data_prefix_dir = "Data/ST_simulation/FigS1_simulation"
      res_prefix_dir = "Result/FigS1_simulation"
      data_prefix = "Simulation_20211109"
      
      sim.res = readRDS(file.path(data_prefix_dir, "BRCA_EMTAB8107_simulation_20211109.rds"))
      sim.cell.comp = sim.res$cell_composition
      sim.cell.comp = sim.cell.comp[,sort(names(sim.cell.comp))]
      sim.cell.comp.frac = sim.cell.comp/rowSums(sim.cell.comp)

      NMF_deconv = read.table(file.path(res_prefix_dir, "NMFreg", paste0("Simulation_20211109", "_NMFreg_deconv.txt")), 
                              header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
      colnames(NMF_deconv) = colnames(sim.cell.comp.frac)
      corr_nmfreg = cor(t(NMF_deconv), t(sim.cell.comp.frac))
      corr_nmfreg = diag(corr_nmfreg)
      median(corr_nmfreg, na.rm = TRUE)
      
      celltype.predictions.score = read.table(file.path(res_prefix_dir, "CCA",paste0("BRCA_EMTAB8107_simulation_20211109", "_celltype.predictions.score.txt")), header = TRUE, row.names = 1, sep = "\t")
      corr_cca = cor(t(celltype.predictions.score), t(sim.cell.comp.frac))
      corr_cca = diag(corr_cca)
      median(corr_cca, na.rm = TRUE)
      
      spotlight_ls = readRDS(file.path(res_prefix_dir, "SPOTlight/spotlight_res.rds"))
      corr_spotlight = cor(t(spotlight_ls[[2]][,-ncol(spotlight_ls[[2]])]), t(sim.cell.comp.frac))
      corr_spotlight = diag(corr_spotlight)
      median(corr_spotlight, na.rm = TRUE)
      
      RCTD_res = read.table(file.path(res_prefix_dir, "RCTD", paste0("BRCA_EMTAB8107_simulation_20211109","_RCTD.txt")), header = TRUE, row.names = 1, sep = "\t")
      corr_rctd = cor(t(RCTD_res), t(sim.cell.comp.frac))
      corr_rctd = diag(corr_rctd)
      median(corr_rctd, na.rm = TRUE)
      
      STRIDE_deconv = read.table(file.path(res_prefix_dir, "STRIDE", paste0(data_prefix, "_spot_celltype_frac.txt")), sep = "\t", header = TRUE, row.names = 1)
      corr_STRIDE = cor(t(as.matrix(STRIDE_deconv)), t(sim.cell.comp.frac))
      corr_STRIDE = diag(corr_STRIDE)
      median(corr_STRIDE, na.rm = TRUE)
      
      cell2location_deconv = read.table(file.path(res_prefix_dir,  "Cell2location", "W_mRNA_count.csv"), 
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
                                     levels = c("STRIDE", "RCTD","Cell2location", "NMFreg", "SPOTlight", "CCA"))
      
      deconv_comp_df$Spot = paste("spot", 1:nrow(sim.cell.comp.frac), sep = "_")
      rownames(sim.cell.comp.frac) = paste("spot", 1:nrow(sim.cell.comp.frac), sep = "_")
      
      p = ggplot(deconv_comp_df, aes(x = Method, y = Corr, fill = Method)) + 
        geom_boxplot(outlier.size = 0.1, weight = 0.5) +  
        scale_fill_manual(values = pal_simpsons("springfield")(16)[c(2,1,3:4,6:7)]) + 
        labs(y = "Pearson's correlation", 
             title = "Spot-level correlation") + 
        scale_y_continuous(limits = c(0,1)) + 
        theme_bw(base_size = 14) +
        theme(legend.position = "none", 
              axis.title.x = element_blank(),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      ggsave(paste(data_prefix, "corr_compare_methods.pdf", sep = "_"), p, width = 4.23, height = 4.55)
      
    }
    
    # Fig S1F
    {
      RootMeanSquareErrorPresence = function(x,y){
        y_above0 = which(y > 0)
        abs_error = x[y_above0] - y[y_above0]
        square_error = abs_error^2
        return(sum(square_error)/length(square_error))
      }
      
      RootMeanSquareErrorAbsence = function(x,y){
        y_0 = which(y == 0)
        abs_error = x[y_0] - y[y_0]
        square_error = abs_error^2
        return(sum(square_error)/length(square_error))
      }
      
      rmse_nmfreg_pre = sapply(1:nrow(NMF_deconv), function(x){
        return(RootMeanSquareErrorPresence(NMF_deconv[x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_nmfreg_abs = sapply(1:nrow(NMF_deconv), function(x){
        return(RootMeanSquareErrorAbsence(NMF_deconv[x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_nmfreg = c(rmse_nmfreg_pre, rmse_nmfreg_abs)
      rmse_nmfreg_presence = rep(c("Present", "Absent"), each = length(rmse_nmfreg_pre))
      
      rmse_rctd_pre = sapply(1:nrow(RCTD_res), function(x){
        return(RootMeanSquareErrorPresence(RCTD_res[x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_rctd_abs = sapply(1:nrow(RCTD_res), function(x){
        return(RootMeanSquareErrorAbsence(RCTD_res[x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_rctd = c(rmse_rctd_pre, rmse_rctd_abs)
      rmse_rctd_presence = rep(c("Present", "Absent"), each = length(rmse_rctd_pre))
      
      
      rmse_spotlight_pre = sapply(1:nrow(spotlight_ls[[2]][,-ncol(spotlight_ls[[2]])]), function(x){
        return(RootMeanSquareErrorPresence(spotlight_ls[[2]][,-ncol(spotlight_ls[[2]])][x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_spotlight_abs = sapply(1:nrow(spotlight_ls[[2]][,-ncol(spotlight_ls[[2]])]), function(x){
        return(RootMeanSquareErrorAbsence(spotlight_ls[[2]][,-ncol(spotlight_ls[[2]])][x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_spotlight = c(rmse_spotlight_pre, rmse_spotlight_abs)
      rmse_spotlight_presence = rep(c("Present", "Absent"), each = length(rmse_spotlight_pre))
      
      
      rmse_cca_pre = sapply(1:nrow(celltype.predictions.score), function(x){
        return(RootMeanSquareErrorPresence(celltype.predictions.score[x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_cca_abs = sapply(1:nrow(celltype.predictions.score), function(x){
        return(RootMeanSquareErrorAbsence(celltype.predictions.score[x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_cca = c(rmse_cca_pre, rmse_cca_abs)
      rmse_cca_presence = rep(c("Present", "Absent"), each = length(rmse_cca_pre))
      
      rmse_STRIDE_pre = sapply(1:nrow(STRIDE_deconv), function(x){
        return(RootMeanSquareErrorPresence(STRIDE_deconv[x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_STRIDE_abs = sapply(1:nrow(STRIDE_deconv), function(x){
        return(RootMeanSquareErrorAbsence(STRIDE_deconv[x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_STRIDE = c(rmse_STRIDE_pre, rmse_STRIDE_abs)
      rmse_STRIDE_presence = rep(c("Present", "Absent"), each = length(rmse_STRIDE_pre))
      
      rmse_cell2location_pre = sapply(1:nrow(cell2location_deconv), function(x){
        return(RootMeanSquareErrorPresence(cell2location_deconv[x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_cell2location_abs = sapply(1:nrow(cell2location_deconv), function(x){
        return(RootMeanSquareErrorAbsence(cell2location_deconv[x,], sim.cell.comp.frac[x, ]))
      })
      
      rmse_cell2location = c(rmse_cell2location_pre, rmse_cell2location_abs)
      rmse_cell2location_presence = rep(c("Present", "Absent"), each = length(rmse_cell2location_pre))
      
      
      deconv_comp_df_rmse = data.frame(RMSE = c(rmse_spotlight, rmse_nmfreg, rmse_cca, rmse_rctd, rmse_STRIDE, rmse_cell2location),
                                       Presence = c(rmse_spotlight_presence, rmse_nmfreg_presence, rmse_cca_presence, rmse_rctd_presence, rmse_STRIDE_presence, rmse_cell2location_presence),
                                       Method = rep(c("SPOTlight", "NMFreg", "CCA", "RCTD", "STRIDE", "Cell2location"), each = length(rmse_spotlight)))
      deconv_comp_df_rmse$Method = factor(deconv_comp_df_rmse$Method, 
                                          levels = c("STRIDE","RCTD", "Cell2location", "NMFreg", "SPOTlight", "CCA"))
      
      mycolor = pal_npg(palette = c("nrc"), alpha = 0.8)(10)[c(5,4)]
      p = ggplot(deconv_comp_df_rmse, aes(x = Method, y = RMSE, fill = Presence)) + 
        geom_boxplot(outlier.size = 0.2, weight = 0.5) + theme_bw(base_size = 14) + 
        # scale_fill_simpsons(palette = c("springfield"), alpha = 0.8) + 
        scale_fill_manual(values = mycolor) +
        theme(axis.title.x = element_blank(),
              plot.title = element_text(hjust = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
              legend.title = element_blank()) + 
        labs(y = "Root Mean Squared Error",
             title = "Spot-level RMSE") + 
        scale_y_continuous(limits = c(0,0.25))
      ggsave(paste(data_prefix, "rmse_presence_compare_methods.pdf", sep = "_"), p, width = 5.4, height = 4.55)
      
    }
    
    
    
    # Fig S1G celltype-level correlation
    {
      corr_celltype_STRIDE = cor(as.matrix(STRIDE_deconv), sim.cell.comp.frac)
      corr_celltype_STRIDE = diag(corr_celltype_STRIDE)
      median(corr_celltype_STRIDE, na.rm = TRUE)
      
      corr_nmfreg_celltype = cor(NMF_deconv, sim.cell.comp.frac)
      corr_nmfreg_celltype = diag(corr_nmfreg_celltype)
      
      corr_cca_celltype = cor(celltype.predictions.score, sim.cell.comp.frac)
      corr_cca_celltype = diag(corr_cca_celltype)
      
      corr_spotlight_celltype = cor(spotlight_ls[[2]][,-ncol(spotlight_ls[[2]])], sim.cell.comp.frac)
      corr_spotlight_celltype = diag(corr_spotlight_celltype)
      names(corr_spotlight_celltype) = colnames(spotlight_ls[[2]][,-ncol(spotlight_ls[[2]])])
      
      corr_rctd_celltype = cor(RCTD_res, sim.cell.comp.frac)
      corr_rctd_celltype = diag(corr_rctd_celltype)
      
      corr_cell2location_celltype = cor(cell2location_deconv, sim.cell.comp.frac)
      corr_cell2location_celltype = diag(corr_cell2location_celltype)
      
      
      
      deconv_celltype_df = data.frame(Corr = c(corr_spotlight_celltype, corr_nmfreg_celltype, corr_cca_celltype, corr_rctd_celltype, corr_celltype_STRIDE, corr_cell2location_celltype),
                                      Celltype = rep(names(corr_spotlight_celltype), 6),
                                      Method = rep(c("SPOTlight", "NMFreg", "CCA", "RCTD", "STRIDE", "Cell2location"), each = length(corr_spotlight_celltype)))
      
      deconv_celltype_df$Method = factor(deconv_celltype_df$Method, 
                                         levels =  c("STRIDE", "RCTD", "Cell2location",  "NMFreg", "SPOTlight", "CCA"))
      deconv_celltype_df$Celltype = gsub("CD4Tconv", "CD4T", deconv_celltype_df$Celltype)
      deconv_celltype_df$Celltype = gsub("CD8Tex", "CD8T", deconv_celltype_df$Celltype)
      deconv_celltype_df$Celltype = gsub("Mono.Macro", "Mono/Macro", deconv_celltype_df$Celltype)
      
      celltypes = c("Malignant", "CD4T", "CD8T", "Tprolif",
                    "B", "Plasma", "Mono/Macro", "Mast",
                    "Endothelial", "Fibroblasts", "Myofibroblasts")
      deconv_celltype_df$Celltype = factor(deconv_celltype_df$Celltype, levels = celltypes)
      
      mycolor = c("grey", brewer.pal(12, "Paired")[c(3:4)],"#01665e", brewer.pal(12, "Paired")[c(1:2, 9:10, 5,7:8)])
      
      
      p = ggplot(deconv_celltype_df, aes(x = Method, y = Corr)) + 
        geom_jitter(aes(colour = Celltype), width = 0.2) + 
        stat_summary(fun = median, fun.min = median, fun.max = median,
                     geom = "crossbar", width = 0.5, size = 0.3) + 
        scale_y_continuous(limits = c(0,1)) + 
        labs(y = "Pearson's correlation", title = "Celltype-level correlation\nwith simulated fraction") + 
        theme_bw(base_size = 14) + 
        scale_colour_manual(values = mycolor) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
              axis.title.x = element_blank(),
              plot.title = element_text(hjust = 0.5), legend.title = element_blank())
      ggsave(paste(data_prefix, "celltype_corr_compare_methods.pdf", sep = "_"), p, width = 5.9, height = 4.8)
      
    }
  }
}   
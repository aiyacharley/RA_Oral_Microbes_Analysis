setwd("D:/我的工作.省中医/Projects/刘院士团队舌苔微生态/results_20260615/")
library(ggplot2)
library(igraph)
library(dplyr)
library(maaslin3)
library(ANCOMBC)
library(readxl)
library(microeco)
source('00.func.R')

#------------------------ humann3 -------------------------------------
#ps <- readRDS(file = "00_data/ps_tongue_humann3ko.rds");outdir <- "output/Figure7/humann3ko/"
ps <- readRDS(file = "00_data/ps_tongue_humann3pathway.rds");outdir <- "output/Figure7/humann3pathway/"

ps <- filter_phyloseq(ps, detection=1, prevalence=0.1)
ps
#------------------------ analysis -------------------------------------
ps_HCV1 <- phyloseq::subset_samples(ps, Time %in% c("HC","V1"))
sam_df <- data.frame(sample_data(ps_HCV1))
otu_df <- data.frame(otu_table(ps_HCV1))

var <- 'Time'
sam_df$Time <- factor(sam_df$Time)
set.seed(42)
# maaslin3
fit_out1 <- maaslin3::maaslin3(input_data = t(otu_df),
                               input_metadata = sam_df,
                               output = paste0(outdir,var),
                               formula = paste0('~',var,"+Age"), # 矫正Age因素
                               normalization = 'TSS',
                               transform = 'LOG',
                               correction = 'BH',
                               median_comparison_abundance = TRUE,
                               max_significance = 0.05,
                               max_pngs = 200,
                               plot_summary_plot = T,
                               coef_plot_vars = var,
                               heatmap_vars = var,
                               save_plots_rds = T,
                               #save_models = T,
                               plot_associations = T)
res_maaslin3 <- read.csv(file = paste0(outdir,var,"/all_results.tsv"), sep = '\t',header = T)
diff_maaslin3 <- subset(res_maaslin3, qval_individual<0.05 & metadata=="Time")
write.csv(file = paste0(outdir,"/sig_maaslin3.csv"), diff_maaslin3, row.names = F)
diff_maaslin3_HC <- subset(diff_maaslin3, coef<0)
diff_maaslin3_V1 <- subset(diff_maaslin3, coef>0)

# ANCOMBC2
set.seed(42)
res_ancombc2 <- ancombc2(
  data = ps_HCV1,
  #prv_cut = 0.1,
  fix_formula = "Time + Age",          # 固定效应：关注Time差异,调整Age变量
  group = "Time",                  # 指定分组变量
  p_adj_method = "BH",
  neg_lb = TRUE,                       # 推荐启用
  struc_zero = TRUE,               # 检测结构零
  n_cl = 6,                        # 使用多核并行加速
  verbose = TRUE
)
## 同时考虑统计显著性和结果稳健性
diff_ancombc2 <- subset(res_ancombc2$res, q_TimeV1 < 0.05 & diff_robust_TimeV1 == "TRUE")
write.csv(file = paste0(outdir,"/sig_ancombc2.csv"), diff_ancombc2, row.names = F)
diff_ancombc2_HC <- subset(diff_ancombc2, lfc_TimeV1<0)
diff_ancombc2_V1 <- subset(diff_ancombc2, lfc_TimeV1>0)

# lefse (microeco)
mt_HCV1 <- file2meco::phyloseq2meco(ps_HCV1)
res_lefse <- trans_diff$new(dataset = mt_HCV1, 
                            method = 'lefse', 
                            group = 'Time',
                            p_adjust_method = 'BH',
                            taxa_level = 'Species',
                            alpha = 0.05)
res_lefse$res_diff$Species <- gsub("s__", "", res_lefse$res_diff$Taxa)
res_lefse$res_diff$Taxa <- gsub("s__", "", res_lefse$res_diff$Taxa)
#res_lefse$res_diff$Species <- stringr::str_split(res_lefse$res_diff$Taxa, pattern = ";|\\|", simplify = TRUE)[,7]
### LDA>2
diff_lefse <- subset(res_lefse$res_diff,LDA>2)
res_lefse$plot_diff_bar(use_number = 1:nrow(diff_lefse))
write.csv(file = paste0(outdir,"/sig_lefse.csv"), diff_lefse, row.names = F)
diff_lefse_HC <- subset(diff_lefse, Group=="HC")
diff_lefse_V1 <- subset(diff_lefse, Group=="V1")


save(res_maaslin3,res_ancombc2,res_lefse, 
     diff_maaslin3,diff_maaslin3_HC,diff_maaslin3_V1,
     diff_ancombc2,diff_ancombc2_HC,diff_ancombc2_V1,
     diff_lefse,diff_lefse_HC,diff_lefse_V1,
     file = paste0(outdir,"/diff_bacteria_HCV1.rdata"))


## venn
library(ggvenn)
lists <- list(MaAsLin3=diff_maaslin3$feature, `ANCOMBC-2`=diff_ancombc2$taxon, LEfSe=diff_lefse$Species)
diff_set <- ggvenn::ggvenn(lists, fill_color = c('#E69F00','#2CA02C','#1F77B4'), stroke_size = 0.5, set_name_size = 4)
diff_set
ggsave(filename = paste0(outdir,"/venn_HCV1.pdf"), diff_set, width = 5.5, height = 5.5)
diff_set@data$Hit <- diff_set@data$MaAsLin3+diff_set@data$`ANCOMBC-2`+diff_set@data$LEfSe
write.csv(diff_set@data, file = paste0(outdir,"/Fig7A_venn_HCV1.csv"), row.names = F)
### 至少两种方法检测出差异的菌
#diff_set2 <- subset(diff_set@data, Hit>=2)
### 三种方法并集
diff_set2 <- diff_set@data
diff_bacteria <- diff_set2$`_key`


lists <- list(MaAsLin3=diff_maaslin3_HC$feature, `ANCOMBC-2`=diff_ancombc2_HC$taxon, LEfSe=diff_lefse_HC$Species)
diff_set <- ggvenn::ggvenn(lists, fill_color = c('#E69F00','#2CA02C','#1F77B4'), stroke_size = 0.5, set_name_size = 4)
diff_set
ggsave(filename = paste0(outdir,"/venn_HC.pdf"), diff_set, width = 5.5, height = 5.5)
diff_set@data$Hit <- diff_set@data$MaAsLin3+diff_set@data$`ANCOMBC-2`+diff_set@data$LEfSe
write.csv(diff_set@data, file = paste0(outdir,"/venn_HC.csv"), row.names = F)
### 至少两种方法检测出差异的菌
#diff_set2 <- subset(diff_set@data, Hit>=2)
### 三种方法并集
diff_set2 <- diff_set@data
diff_bacteria <- diff_set2$`_key`

lists <- list(MaAsLin3=diff_maaslin3_V1$feature, `ANCOMBC-2`=diff_ancombc2_V1$taxon, LEfSe=diff_lefse_V1$Species)
diff_set <- ggvenn::ggvenn(lists, fill_color = c('#E69F00','#2CA02C','#1F77B4'), stroke_size = 0.5, set_name_size = 4)
diff_set
ggsave(filename = paste0(outdir,"/venn_V1.pdf"), diff_set, width = 5.5, height = 5.5)
diff_set@data$Hit <- diff_set@data$MaAsLin3+diff_set@data$`ANCOMBC-2`+diff_set@data$LEfSe
write.csv(diff_set@data, file = paste0(outdir,"/venn_V1.csv"), row.names = F)
### 至少两种方法检测出差异的菌
#diff_set2 <- subset(diff_set@data, Hit>=2)
### 三种方法并集
diff_set2 <- diff_set@data
diff_bacteria <- diff_set2$`_key`


### 对于差异pathway，绘制热图
# 07.humann3_plotHeat.R


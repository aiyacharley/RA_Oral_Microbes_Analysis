setwd("D:/我的工作.省中医/Projects/刘院士团队舌苔微生态/results_20260615/")
library(ggplot2)
library(igraph)
library(dplyr)
library(maaslin3)
library(ANCOMBC)
library(readxl)
source('00.func.R')


#------------------------ sylph (input1) -------------------------------------
outdir <- "output/Figure2/"
ps <- readRDS(file = "00_data/ps_tongue_sylph.rds"); oName <- 'maaslin3_Readbase'
#------------------------ MAGs (inout2) -------------------------------------
# outdir <- "output/Figure2_MAGs/"
# ps <- readRDS(file = "00_data/ps_tongue_MAGs.rds"); oName <- 'maaslin3_MAGs'


#------------------------ analysis -------------------------------------
ps_HCV1 <- phyloseq::subset_samples(ps, Time %in% c("HC","V1"))
sam_df <- data.frame(sample_data(ps_HCV1))
otu_df <- data.frame(otu_table(ps_HCV1))

var <- 'Time'
sam_df$Time <- factor(sam_df$Time)
set.seed(1234)
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
                               #summary_plot_first_n = 25,
                               save_plots_rds = T,
                               #save_models = T,
                               plot_associations = T)
res_maaslin3 <- read.csv(file = paste0(outdir,var,"/all_results.tsv"), sep = '\t',header = T)
diff_maaslin3 <- subset(res_maaslin3, qval_individual<0.05 & metadata=="Time") 

# ANCOMBC2
set.seed(1234)
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

# lefse (microeco)
mt_HCV1 <- file2meco::phyloseq2meco(ps_HCV1)
res_lefse <- trans_diff$new(dataset = mt_HCV1, 
                            method = 'lefse', 
                            group = 'Time',
                            p_adjust_method = 'BH',
                            taxa_level = 'Species',
                            alpha = 0.05)
res_lefse$res_diff$Species <- stringr::str_split(res_lefse$res_diff$Taxa, pattern = ";|\\|", simplify = TRUE)[,7]
### LDA>2
diff_lefse <- subset(res_lefse$res_diff,LDA>2)
res_lefse$plot_diff_bar(use_number = 1:nrow(diff_lefse))


## venn
library(ggvenn)
lists <- list(MaAsLin3=diff_maaslin3$feature, `ANCOMBC-2`=diff_ancombc2$taxon, LEfSe=diff_lefse$Species)
diff_set <- ggvenn::ggvenn(lists, fill_color = c('#E69F00','#2CA02C','#1F77B4'), stroke_size = 0.5, set_name_size = 4)
diff_set
ggsave(filename = paste0(outdir,"/FigS2_venn.pdf"), diff_set, width = 5.5, height = 5.5)
diff_set@data$Hit <- diff_set@data$MaAsLin3+diff_set@data$`ANCOMBC-2`+diff_set@data$LEfSe
write.csv(diff_set@data, file = paste0(outdir,"/FigS2A_venn.csv"), row.names = F)

### 至少两种方法检测出差异的菌
#diff_set2 <- subset(diff_set@data, Hit>=2)
### 三种方法并集
diff_set2 <- diff_set@data
diff_bacteria <- diff_set2$`_key`
save(res_maaslin3,res_ancombc2,res_lefse,diff_set2, file = paste0(outdir,"/diff_bacteria_HCV1.rdata"))


ps1 <- ps@tax_table[diff_bacteria,]
tax_df <- data.frame(tax_table(ps1))
tax_df$Kingdom <- "k__Bacteria"
z1 <- apply(tax_df, 1, paste, collapse = "|") %>% as.data.frame()
diff_bacteria_taxon <- z1[,1]

env_cols = c('Age','BMI','DAS28_ESR','DAS28_CRP','cDAI','sDAI','SJC','ESR','TJC','VAS','RF','CRP','HAQ','PaGADA','PhGADA','Stiffness_time')
t1 <- trans_env$new(dataset = mt_HCV1, env_cols = env_cols)

# use other_taxa to select taxa you need
t1$cal_cor(method = "spearman", use_data = "other", p_adjust_method = "BH", other_taxa = diff_bacteria_taxon)
fig2b <- t1$plot_cor(color_vector = c("#1187A8", "white", "#ED8A10"),
                     cluster_ggplot = "both",
                     xtext_angle = 60)
#ggsave(filename = paste0(outdir,"/Fig2B_corrHeat.pdf"), fig2b, width = 9, height = 12)

## 相关性热图+三种检验信息
df_corr <- fig2b$plotlist[[1]]@data

res_maaslin3_all <- read.csv(file = paste0(outdir,var,"/all_results.tsv"), sep = '\t',header = T)
diff_maaslin3_target1 <- subset(res_maaslin3_all, metadata=="Time" & feature %in% diff_bacteria & model=="abundance")
diff_maaslin3_target2 <- subset(res_maaslin3_all, metadata=="Time" & feature %in% diff_bacteria & model=="prevalence")

diff_ancombc2_target <- subset(res_ancombc2$res, taxon %in% diff_bacteria)

diff_lefse_target1 <- subset(res_lefse$res_diff, Species %in% diff_bacteria & Group=="HC")
diff_lefse_target2 <- subset(res_lefse$res_diff, Species %in% diff_bacteria & Group=="V1")
diff_lefse_target1$LDA <- 0 - diff_lefse_target1$LDA
diff_lefse_target <- rbind(diff_lefse_target1, diff_lefse_target2)

h1 <- diff_maaslin3_target1[,c('feature','coef','qval_individual')]
h2 <- diff_maaslin3_target2[,c('feature','coef','qval_individual')]
h3 <- diff_ancombc2_target[,c('taxon','lfc_TimeV1','q_TimeV1')]
h4 <- diff_lefse_target[,c('Species','LDA','P.adj')]

h12 <- merge(h1,h2,by.x = 'feature',by.y = 'feature', all = T)
h123 <- merge(h12,h3,by.x = 'feature',by.y = 'taxon', all = T)
h1234 <- merge(h123,h4,by.x = 'feature',by.y = 'Species', all = T)
h1234coef <- h1234[,c('coef.x','coef.y','lfc_TimeV1','LDA')]
h1234qval <- h1234[,c('qval_individual.x','qval_individual.y','q_TimeV1','P.adj')]
colnames(h1234coef) <- c('MaAsLin3_abundance_coef','MaAsLin3_prevalence_coef','ANCOMBC2_lfc','LEfSe_LDA')
colnames(h1234qval) <- c('MaAsLin3_abundance_qval','MaAsLin3_prevalence_qval','ANCOMBC2_qval','LEfSe_Padj')
rownames(h1234coef) <- h1234$feature
rownames(h1234qval) <- h1234$feature

write.csv(df_corr, file = paste0(outdir,"/Fig2_corr.csv"), row.names = F)
write.csv(h1234coef, file = paste0(outdir,"/Fig2_coef_lfc_lda.csv"), row.names = T)
write.csv(h1234qval, file = paste0(outdir,"/Fig2_qval.csv"), row.names = T)

## 执行02.plotHeat.R绘制联合热图









setwd("D:/我的工作.省中医/Projects/刘院士团队舌苔微生态/results_20260615/")
library(ggplot2)
library(igraph)
library(dplyr)
library(readxl)
library(ComplexHeatmap)
library(circlize)
source('00.func.R')
#------------------------ sylph -------------------------------------
outdir <- "output/Figure3/"
ps <- readRDS(file = "00_data/ps_tongue_sylph.rds"); oName <- 'V1V4'
ps_RA <- phyloseq::subset_samples(ps, Group_tongue %in% c("V1N","V1R","V4N","V4R"))

sam_df <- data.frame(sample_data(ps_RA))
otu_df <- data.frame(otu_table(ps_RA))

var <- 'Time'
sam_df$Time <- factor(sam_df$Time)
set.seed(42)
# maaslin3
fit_out1 <- maaslin3::maaslin3(input_data = t(otu_df),
                               input_metadata = sam_df,
                               output = paste0(outdir,var),
                               formula = ~ Time + (1|PatientID),  # 配对
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
                               small_random_effects = TRUE,
                               plot_associations = T)
res_maaslin3 <- read.csv(file = paste0(outdir,var,"/all_results.tsv"), sep = '\t',header = T)
diff_maaslin3 <- subset(res_maaslin3, qval_individual<0.05 & metadata=="Time") 

# V1-V4: 治疗前后，配对样本，矫正PatientID因素
set.seed(42)
res_ancombc2 <- ancombc2(
  data = ps_RA,
  fix_formula = "Time",          # 固定效应：关注Time差异
  rand_formula = "(1|PatientID)", # 随机效应：控制个体配对效应
  group = "Time",                  # 指定分组变量
  p_adj_method = "BH",
  neg_lb = TRUE,                       # 推荐启用
  struc_zero = TRUE,               # 检测结构零
  n_cl = 8,                        # 使用多核并行加速
  verbose = TRUE
)
## 同时考虑统计显著性和结果稳健性
diff_ancombc2 <- subset(res_ancombc2$res, q_TimeV4 < 0.05 & diff_robust_TimeV4 == "TRUE")

# lefse (microeco)
mt_RA <- file2meco::phyloseq2meco(ps_RA)
res_lefse <- trans_diff$new(dataset = mt_RA, 
                            method = 'lefse', 
                            group = 'Time',
                            p_adjust_method = 'none', #'BH',
                            taxa_level = 'Species',
                            alpha = 0.05)
res_lefse$res_diff$Species <- stringr::str_split(res_lefse$res_diff$Taxa, pattern = ";|\\|", simplify = TRUE)[,7]
### LDA>2
diff_lefse <- subset(res_lefse$res_diff,LDA>2)
nrow(diff_lefse)
res_lefse$plot_diff_bar(use_number = 1:nrow(diff_lefse))
res_lefse$plot_diff_abund(use_number = 1:nrow(diff_lefse))

## venn
library(ggvenn)
lists <- list(MaAsLin3=diff_maaslin3$feature, `ANCOMBC-2`=diff_ancombc2$taxon, LEfSe=diff_lefse$Species)
diff_set <- ggvenn::ggvenn(lists, fill_color = c('#E69F00','#2CA02C','#1F77B4'), stroke_size = 0.5, set_name_size = 4)
diff_set
ggsave(filename = paste0(outdir,"/FigS3_venn.pdf"), diff_set, width = 5.5, height = 5.5)
diff_set@data$Hit <- diff_set@data$MaAsLin3+diff_set@data$`ANCOMBC-2`+diff_set@data$LEfSe
write.csv(diff_set@data, file = paste0(outdir,"/FigS2B_venn.csv"), row.names = F)

### 至少两种方法检测出差异的菌
#diff_set2 <- subset(diff_set@data, Hit>=2)
### 三种方法并集
diff_set2 <- diff_set@data
diff_bacteria <- diff_set2$`_key`
save(res_maaslin3,res_ancombc2,res_lefse,diff_set2, file = paste0(outdir,"/diff_bacteria_V1V4.rdata"))



diff_HCV1 <- read.csv(file = "output/Figure2/FigS2A_venn.csv", row.names = 1)
diff_V1V4 <- read.csv(file = "output/Figure3/FigS2B_venn.csv", row.names = 1)

target_bacteria <- intersect(rownames(diff_HCV1), rownames(diff_V1V4)) 
lists2 <- list(`HC vs V1`=rownames(diff_HCV1), `V1 vs V4`=rownames(diff_V1V4))
diff_set <- ggvenn::ggvenn(lists2, fill_color = c('#E69F00','#2CA02C'), stroke_size = 0.5, set_name_size = 4)
diff_set
ggsave(filename = paste0(outdir,"/FigS2C_venn2.pdf"), diff_set, width = 5.5, height = 5.5)

# 桑基图
target_bacteria_taxa <- apply(ps@tax_table[target_bacteria,c(2:7)], 1, paste, collapse = "|") %>% as.data.frame()
target_bacteria_taxa[,1]
plot_sankey(target_bacteria_taxa[,1], split = '|', saveWidget = paste0(outdir,"/FigS2D_sankeyNetwork.html"), height = 500, width = 800)

## 关注MAGs的丰度
mt <- phyloseq2meco(ps)

tmp <- trans_norm$new(dataset = mt)
norm_method <- "TSS"
mt2 <- tmp$norm(method = norm_method)
mt2$otu_table <- log(1e6*mt2$otu_table+1)


plot_list <- list()
for (f in target_bacteria){
  fig <- plot_box(dataset = mt2, 
                  y = f, 
                  group = "Time", 
                  testMethod = 'wilcox.test', 
                  signif_label = FALSE, 
                  ysqrt = FALSE,
                  paired_groups = c("V1", "V4"),  # 指定要配对的组
                  paired_id = "PatientID",         # 患者ID列名
                  show_paired_pval = TRUE          # 显示配对检验p值
  ) +
    labs(x = "", y = paste0("log1p(1e6*", norm_method, " Abundance)")) +
    scale_color_manual(values = c("HC" = "#44AF8F", "V1" = "#E07C30", "V4" = "#756FB3")) +
    scale_fill_manual(values = c("HC" = "#44AF8F", "V1" = "#E07C30", "V4" = "#756FB3")) +
    theme(legend.position = 'none')
  plot_list[[f]] <- fig
}
# 使用wrap_plots进行4×4布局
combined_plot <- wrap_plots(plot_list, ncol = 4, nrow = 4)
# 保存
ggsave(filename = paste0(outdir,"/Fig3_boxplot_target_HCV1V4.pdf"), 
       plot = combined_plot, width = 14, height = 14)








# extracting OTU data by manipulating otu_table
test <- clone(mt2)
test$otu_table <- test$otu_table[MAGs_target, ]
test$tidy_dataset()

res_lefse0 <- read.csv(file = "02_diff_V1V4/diff_lefse_tongue_species.csv", row.names = 1)
res_lefse <- res_lefse0[MAGs_target,c("LDA","Phylum","Family","Genus")]

res_ancombc2 <- read.csv(file = "02_diff_V1V4/diff_ancombc2paired_tongue_species.csv", row.names = 1)
res_ancombc <- res_ancombc2[MAGs_target,c("lfc_TimeV4","p_TimeV4","q_TimeV4")]

res_rf0 <- read.csv(file = "02_diff_V1V4/diff_rf_tongue_species.csv", row.names = 1)
res_rf <- res_rf0[MAGs_target,c("MeanDecreaseGini","P.unadj","P.adj")]

top_annotation = HeatmapAnnotation(Time=test$sample_table$Time, Age=test$sample_table$Age,BMI=test$sample_table$BMI,
                                   col = list(Time=c("HC" = "#44AF8F", "V1" = "#E07C30", "V4" = "#756FB3"),
                                              Age = colorRamp2(c(18, 75), c("white", "#608C1B")),
                                              BMI = colorRamp2(c(15, 40), c("white", "#2EBFBF"))))

left_annotation = rowAnnotation(Phylum=res_lefse$Phylum, 
                                Genus=res_lefse$Genus, 
                                LDA=anno_barplot(res_lefse$LDA,gp = gpar(fill = "#FF7F00"),
                                                 axis = TRUE, # 确保显示坐标轴
                                                 axis_param = list(
                                                   side = "top", # 将轴置于条形图底部
                                                   facing = "outside", # 标签朝外
                                                   direction = "reverse",
                                                   labels_rot = 0 # 标签水平放置
                                                 )),
                                col = list(Phylum = c("p__Pseudomonadota" = "#E41A1C", "p__Bacillota" = "#377EB8", "p__Actinomycetota" = "#4DAF4A"),
                                           Genus = c("g__Aggregatibacter" = "#984EA3", "g__Butyrivibrio" = "#FF7F00", "g__Pauljensenia" = "#FFFF33", "g__Eubacterium_B" = "#A65628")))
right_annotation = rowAnnotation(ANCOMBC2.lfc=res_ancombc$lfc_TimeV4, 
                                 ANCOMBC2.pval=(res_ancombc$p_TimeV4),
                                 MeanDecreaseGini=anno_barplot(res_rf$MeanDecreaseGini,gp = gpar(fill = "#1E78C0"),
                                                               axis_param = list(
                                                                 side = "top", # 将轴置于条形图底部
                                                                 facing = "outside", # 标签朝外
                                                                 labels_rot = 90 # 标签水平放置
                                                               )),
                                 col = list(ANCOMBC2.lfc = colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C")),
                                            ANCOMBC2.pval = colorRamp2(c(0.001, 0.05, 0.06), c("red","white","#377EB8"))))

ht <- ComplexHeatmap::Heatmap(test$otu_table, 
                              name = "Scale",
                        top_annotation = top_annotation,
                        left_annotation = left_annotation,
                        right_annotation = right_annotation,
                        cluster_rows = F, cluster_columns = F, 
                        show_column_names = F, 
                        col = col_fun, 
                        column_split = test$sample_table$Time,
                        row_split = c("up_down_up","down_up_down","down_up_down","down_up_down"),
                        row_gap = unit(2, "mm"), column_gap = unit(2, "mm"))

# 2. 绘制热图，并将图例置于底部
pdf(file = "02_diff_V1V4/heatmap_MAGs_target_HC_V1_V4.pdf", width = 12, height = 6)
ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
dev.off()


df1 <- read.csv("01_Maaslin3_Tongue/DAS28_ESR/all_results.tsv", sep = "\t", header = T)
#subset(df1, feature %in% "s__Oribacterium_parvum_MAG169" & model=="abundance" & metadata!="Age")
df1a <- subset(df1, feature %in% MAGs_target & model=="abundance" & metadata!="Age")
z1 <- df1a[,c('feature','metadata','coef','pval_individual','qval_individual')]

df2 <- read.csv("01_Maaslin3_Tongue/DAS28_CRP/all_results.tsv", sep = "\t", header = T)
df2a <- subset(df2, feature %in% MAGs_target & model=="abundance" & metadata!="Age")
z2 <- df2a[,c('feature','metadata','coef','pval_individual','qval_individual')]

df3 <- read.csv("01_Maaslin3_Tongue/sDAI/all_results.tsv", sep = "\t", header = T)
df3a <- subset(df3, feature %in% MAGs_target & model=="abundance" & metadata!="Age")
z3 <- df3a[,c('feature','metadata','coef','pval_individual','qval_individual')]

df4 <- read.csv("01_Maaslin3_Tongue/cDAI/all_results.tsv", sep = "\t", header = T)
df4a <- subset(df4, feature %in% MAGs_target & model=="abundance" & metadata!="Age")
z4 <- df4a[,c('feature','metadata','coef','pval_individual','qval_individual')]


z0 <- rbind(z1,z2,z3,z4)
# 1. 加载必要的R包
library(pheatmap)
library(tidyr) # 用于数据重塑
library(dplyr) # 用于数据处理

# 假设您的数据框名为 z0
# 2. 数据准备：将长格式数据重塑为热图矩阵
# 提取用于绘图的列：feature, metadata, coef
heatmap_data <- z0 %>%
  dplyr::select(feature, metadata, coef) %>%
  pivot_wider(names_from = metadata, values_from = coef) %>%
  as.data.frame()

# 将第一列（feature）设为行名，并转换为纯数值矩阵
row_names <- heatmap_data$feature
matrix_for_plot <- as.matrix(heatmap_data[, -1])
rownames(matrix_for_plot) <- row_names

# 3. 创建显著性标记矩阵
# 提取用于标记的列：feature, metadata, qval_individual
sig_data <- z0 %>%
  dplyr::select(feature, metadata, qval_individual) %>%
  pivot_wider(names_from = metadata, values_from = qval_individual) %>%
  as.data.frame()

# 构建逻辑矩阵：qval < 0.05 的位置标记为 TRUE (即需要标注*)
sig_matrix <- as.matrix(sig_data[, -1])
rownames(sig_matrix) <- sig_data$feature
# 确保显著性矩阵的行列顺序与绘图矩阵完全一致
sig_matrix <- sig_matrix[rownames(matrix_for_plot), colnames(matrix_for_plot)]
# 基于显著性阈值创建星号标记矩阵
display_matrix <- ifelse(sig_matrix < 0.001, "***",
                         ifelse(sig_matrix < 0.01, "**",
                                ifelse(sig_matrix < 0.05, "*", "")))
matrix_for_plot <- matrix_for_plot[c('s__Butyrivibrio_sp015258065_MAG137','s__Pauljensenia_sp900554605_MAG143','s__Eubacterium_B_sulci_MAG204','s__Aggregatibacter_sp000466335_MAG074'),]
display_matrix <- display_matrix[c('s__Butyrivibrio_sp015258065_MAG137','s__Pauljensenia_sp900554605_MAG143','s__Eubacterium_B_sulci_MAG204','s__Aggregatibacter_sp000466335_MAG074'),]

# 4. 绘制热图
pdf(file = "02_diff_V1V4/heatmap_Maaslin3_MAGs_target.pdf", width = 2, height = 6)
pheatmap(matrix_for_plot,
         show_rownames = F,
         # --- 核心参数：基于coef值着色 ---
         color = colorRampPalette(c("#94A4BD", "white", "#E58679"))(50), # 定义从蓝到红的渐变色
         scale = "none", # 不对数据进行行或列标准化，直接使用原始coef值
         # --- 显著性标注参数 ---
         display_numbers = display_matrix, # 传入自定义的标记矩阵，在指定位置显示"*"
         number_color = "black", # 设置星号为黑色
         fontsize_number = 12,   # 设置星号字体大小
         # --- 图形美化参数（可选）---
         cluster_rows = FALSE,   # 行不聚类，保持原始顺序
         cluster_cols = FALSE,   # 列不聚类，保持原始顺序
         border_color = NA,      # 去除单元格边框
         #main = "Coefficient Heatmap with Significance (*: qval < 0.05)", # 添加标题
         angle_col = 90,         # 倾斜列标签，避免重叠
         gaps_row = c(3), # 在行之间添加间隔，突出显示不同的feature
         # --- 图例与颜色条控制（可选）---
         # legend_breaks = c(-1, -0.5, 0, 0.5, 1), # 自定义图例断点
         # legend_labels = c("-1.0", "-0.5", "0", "0.5", "1.0") # 自定义图例标签
)
dev.off()

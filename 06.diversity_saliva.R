setwd("D:/我的工作.省中医/Projects/刘院士团队舌苔微生态/results_20260615/")
library(ggplot2)
library(patchwork)
library(ggplot2)
library(ggsignif)
library(ggsci)
library(dplyr)
library(microeco)
source('00.func.R')

#------------------------ sylph -------------------------------------
outdir <- "output/Figure6"
ps1 <- readRDS(file = "00_data/ps_tongue_sylph.rds")
ps1@sam_data$Time <- paste0(ps1@sam_data$Time,"_TC")
ps2 <- readRDS(file = "00_data/ps_saliva_sylph.rds")
ps2@sam_data$Time <- paste0(ps2@sam_data$Time,"_Sa")
ps <- merge_phyloseq(ps1, ps2)

alpha_diversity <- cal_alpha(ps,group = 'Time',signif_label=T)
write.csv(alpha_diversity$table_index, file = paste0(outdir,"/Fig6A_alpha_diversity.csv"), row.names = T)

group_colors <- c("HC_TC" = "#44AF8F", "V1_TC" = "#E07C30", "V4_TC" = "#756FB3", 
                  "HC_Sa" = "#8FCFBC", "V1_Sa" = "#FFCC99", "V4_Sa" = "#CCCCFF")

fig1a <- alpha_diversity$figure_diversity_shannon+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(x="")+
  scale_color_manual(values = group_colors)+
  scale_fill_manual(values = group_colors)

beta_diversity <- cal_ordination(ps, method_dist = 'bray', method = 'PCoA', group = 'Time')
sink(paste0(outdir,"/Figure6B_anosim.txt"), split = TRUE)
print(beta_diversity$rds_anosim)
sink()

write.csv(beta_diversity$table_PERMANOVA, file = paste0(outdir,"/Figure6C_PERMANOVA.csv"), row.names = F)
figure_anosim <- beta_diversity$figure_anosim

fig1b <- ggplot(figure_anosim@data, aes(x=class, y=rank, fill = class))+
  geom_boxplot(outliers = TRUE, alpha = 0.7,width=1,
               notch = TRUE, notchwidth = 0.5, staplewidth = 0, varwidth = TRUE, 
               na.rm = FALSE, 
               orientation = NA)+
  labs(title = figure_anosim@labels$title,
       subtitle = figure_anosim@labels$caption,
       x= "", y = "") +
  theme_classic() + 
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("Between" = "gray",group_colors))

p2 <- beta_diversity$figure_ord +
  scale_color_manual(values = group_colors)+
  scale_fill_manual(values = group_colors)

pairwise_results <- beta_diversity$table_PERMANOVA %>% filter(pairs != "Total") %>% select(pairs, R2, p.value)
# 构建要显示的文字
text_permanova <- paste0(
  "HC_TC vs HC_Sa: R² = ", pairwise_results$R2[pairwise_results$pairs == "HC_TC vs HC_Sa"], ", p = ", pairwise_results$p.value[pairwise_results$pairs == "HC_TC vs HC_Sa"], "\n",
  "V1_TC vs V1_Sa: R² = ", pairwise_results$R2[pairwise_results$pairs == "V1_TC vs V1_Sa"], ", p = ", pairwise_results$p.value[pairwise_results$pairs == "V1_TC vs V1_Sa"], "\n",
  "V4_TC vs V4_Sa: R² = ", pairwise_results$R2[pairwise_results$pairs == "V4_TC vs V4_Sa"], ", p = ", pairwise_results$p.value[pairwise_results$pairs == "V4_TC vs V4_Sa"], "\n"
)
# 添加到图中（放在右上角空白区域）
fig1c <- p2 + annotate(
  "text",
  x =  -0.6,  # 根据实际情况调整位置
  y =  0.6,
  label = text_permanova,
  size = 4,
  hjust = 0,  # 左对齐
  vjust = 1,  # 顶部对齐
  fontface = "plain"
)
fig1abc <- fig1a+fig1b+fig1c+plot_layout(width = c(1, 1.5, 2))
ggsave(paste0(outdir,"/Fig6ABC_alpha_beta_diversity.pdf"), fig1abc, width = 11, height = 4.5)


for (tax in c('Phylum','Family','Genus','Species')){
  print(tax)
  p1a <- plot_comp(ps, taxrank = tax, strata = 'Time') + theme(legend.position = "none")
  p1b <- plot_comp(ps, taxrank = tax, groupmean = 'Time', use_alluvium = T, sort_bacteria = F) + labs(y="")
  p1 <- p1a+p1b+plot_layout(width = c(15, 1))
  ggsave(filename = paste0(outdir,"/FigS6_comp_sylph_",tax,".pdf"), p1, width = 25, height = 6.5)
}
p1_phylum <- plot_comp(ps, taxrank = 'Phylum', groupmean = 'Time', use_alluvium = T, sort_bacteria = F)
p1_family <- plot_comp(ps, taxrank = 'Family', groupmean = 'Time', use_alluvium = T, sort_bacteria = F)
p1_genus <- plot_comp(ps, taxrank = 'Genus', groupmean = 'Time', use_alluvium = T, sort_bacteria = F)
p1_species <- plot_comp(ps, taxrank = 'Species', groupmean = 'Time', use_alluvium = T, sort_bacteria = F)
fig1def <- p1_phylum + p1_genus + p1_species
fig1def <- p1_family + p1_genus + p1_species
ggsave(filename = paste0(outdir,"/Fig6DEF_comp_sylph.pdf"), fig1def, width = 12, height = 5)

write.csv(p1_phylum@data, file = paste0(outdir,"/FigureS6_Phylum.csv"), row.names = F)
write.csv(p1_family@data, file = paste0(outdir,"/FigureS6_Family.csv"), row.names = F)
write.csv(p1_genus@data, file = paste0(outdir,"/FigureS6_Genus.csv"), row.names = F)
write.csv(p1_species@data, file = paste0(outdir,"/FigureS6_Species.csv"), row.names = F)


ps.rel <- filter_phyloseq(ps, rel = "total")
fig6g <- plot_detection_prevalence(ps.rel, min_prevalence = 0.1, min_detections = 0.005, step_detection = 100)
ggsave(filename = paste0(outdir,"/Fig6G_core_sylph.pdf"), fig6g, width = 11, height = 6)
df_wide <- fig1g@data %>%
  pivot_wider(
    names_from = DetectionThreshold,  # 将 DetectionThreshold 的值作为新列名
    values_from = Prevalence,          # 将 Prevalence 的值填充到新列中
    names_prefix = "Detection_"
  )
write.csv(df_wide, file = paste0(outdir,"/Figure6G_core.csv"), row.names = F)










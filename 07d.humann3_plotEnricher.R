# ====================================================================
# 07.humann3_plotEnricher.R
# KO / Module 富集气泡图 — HC vs V1 特有通路
# 筛选 qvalue < 0.05 的条目，去掉两组共有的通路
# ====================================================================

library(ggplot2)
library(dplyr)

# ---- 绘图函数 ----
plot_enrich_bubble <- function(file_hc, file_v1, out_pdf, out_csv, width = 7.2, height = 8) {
  hc <- read.csv(file_hc)
  v1 <- read.csv(file_v1)

  # 筛选
  hc_sig <- filter(hc, qvalue < 0.05)
  v1_sig <- filter(v1, qvalue < 0.05)

  # 去共有
  common_desc <- intersect(hc_sig$Description, v1_sig$Description)
  hc_uniq <- hc_sig %>% filter(!Description %in% common_desc) %>% mutate(Group = "HC")
  v1_uniq <- v1_sig %>% filter(!Description %in% common_desc) %>% mutate(Group = "V1")

  cat(basename(file_hc), "显著条目数:", nrow(hc_sig), "\n")
  cat(basename(file_v1), "显著条目数:", nrow(v1_sig), "\n")
  cat("共有通路数:", length(common_desc), "\n")
  cat("HC 特有通路数:", nrow(hc_uniq), "\n")
  cat("V1 特有通路数:", nrow(v1_uniq), "\n\n")

  # 合并
  dat <- bind_rows(hc_uniq, v1_uniq)

  # GeneRatio 转数值
  dat$GeneRatio <- sapply(strsplit(dat$GeneRatio, "/"), function(x) {
    as.numeric(x[1]) / as.numeric(x[2])
  })

  # 排序 (每组内 qvalue 大到小)
  dat <- dat %>%
    arrange(Group, desc(qvalue)) %>%
    mutate(Description = factor(Description, levels = unique(Description)))

  # 气泡图
  p <- ggplot(dat, aes(x = FoldEnrichment, y = Description,
                       size = RichFactor, colour = qvalue)) +
    geom_point() +
    scale_size_continuous(range = c(3, 10)) +
    scale_colour_gradient(low = "#DF6664", high = "#3A7EB9") +
    facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "FoldEnrichment", y = NULL,
         size = "RichFactor", colour = "q-value") +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "white", colour = "black"),
      strip.text = element_text(size = 11, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  # 输出绘图数据
  dat_out <- dat %>% select(-Description)  # 移除 factor 列，用原始 Description
  dat_out$Description <- as.character(dat$Description)
  write.csv(dat_out, out_csv, row.names = FALSE)

  ggsave(out_pdf, p, width = width, height = height)
  cat("→ 输出:", out_pdf, "\n")
  cat("→ 输出:", out_csv, "\n\n")
}

# ---- KO 气泡图 -> Fig7C ----
plot_enrich_bubble(
  file_hc  = "output/Figure7/humann3ko/enrichKO_HC.csv",
  file_v1  = "output/Figure7/humann3ko/enrichKO_V1.csv",
  out_pdf  = "output/Figure7/humann3ko/Fig7E_enrichKO_bubble.pdf",
  out_csv  = "output/Figure7/humann3ko/Fig7E_enrichKO_data.csv",
  height   = 7
)

# ---- Module 气泡图 -> Fig7D ----
plot_enrich_bubble(
  file_hc  = "output/Figure7/humann3ko/enrichModule_HC.csv",
  file_v1  = "output/Figure7/humann3ko/enrichModule_V1.csv",
  out_pdf  = "output/Figure7/humann3ko/Fig7F_enrichModule_bubble.pdf",
  out_csv  = "output/Figure7/humann3ko/Fig7F_enrichModule_data.csv",
  width    = 9.36,
  height   = 4.48
)


## V1组富集的KO可视化（可以选择特定的showCategory）
df_enrichKO <- read.csv(file = paste0(outdir,"Fig7E_enrichKO_data.csv"),sep = ',')
df_enrichKO_V1 <- subset(df_enrichKO, Group=="V1")
fig.cnet.ego <- enrichplot::cnetplot(resKO$obj, size_item = 0.5, hilight_alpha=0.3, node_label="share", color_edge = "category", showCategory = df_enrichKO_V1$Description)
ego.pairwise_termsim <- enrichplot::pairwise_termsim(resKO$obj)
fig.emap.ego <- enrichplot::emapplot(ego.pairwise_termsim, node_label_size=3, showCategory = df_enrichKO_V1$Description)
#fig.tree.ego <- enrichplot::treeplot(ego.pairwise_termsim, showCategory = df_enrichKO_V1$Description)

ggsave(filename = paste0(outdir,"/Fig7G_enrichKO_V1_emapplot.pdf"), fig.emap.ego, width = 6, height = 5)
ggsave(filename = paste0(outdir,"/enrichKO_V1_cnetplot.pdf"), fig.cnet.ego, width = 15, height = 8)


df_enrichModule <- read.csv(file = paste0(outdir,"Fig7F_enrichModule_data.csv"),sep = ',')
df_enrichModule_V1 <- subset(df_enrichModule, Group=="V1")
fig.cnet.ego <- enrichplot::cnetplot(resModule$obj,size_item = 0.5, hilight_alpha=0.3, node_label="share", color_edge = "category", showCategory = df_enrichModule_V1$Description)
ego.pairwise_termsim <- enrichplot::pairwise_termsim(resModule$obj)
fig.emap.ego <- enrichplot::emapplot(ego.pairwise_termsim, node_label_size=3) #, showCategory = df_enrichModule_V1$Description)
#fig.tree.ego <- enrichplot::treeplot(ego.pairwise_termsim, showCategory = df_enrichModule_V1$Description)

#ggsave(filename = paste0(outdir,"/Fig7E_enrichModule_V1_emapplot.pdf"), fig.emap.ego, width = 6, height = 5)
#ggsave(filename = paste0(outdir,"/Fig7F_enrichModule_V1_cnetplot.pdf"), fig.cnet.ego, width = 8, height = 6)

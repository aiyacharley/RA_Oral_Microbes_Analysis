# 02.plot.R
# Figure 2: Integrated clinical correlation heatmap
#   - Main panel: Taxa x Env correlation heatmap with significance stars
#   - Right panels: MaAsLin3 (abundance + prevalence), ANCOMBC2 LFC (heatmap)
#                   and LEfSe LDA (barplot) with q-value significance
# ============================================================================

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

# ---- Paths (outdir二选一，以sylph为主，MAGs补充)----
outdir <- "output/Figure2/"      # sylph
#outdir <- "output/Figure2_MAGs/"  # MAGs

corr_file <- paste0(outdir,"Fig2_corr.csv")
coef_file <- paste0(outdir,"Fig2_coef_lfc_lda.csv")
qval_file <- paste0(outdir,"Fig2_qval.csv")
out_pdf   <- paste0(outdir,"Fig2_corr.pdf")

# ---- Functions ----
sig_star <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01,  "**",
    ifelse(p < 0.05,   "*", ""))))
}

clean_taxa <- function(x) {
  x <- gsub("^s__", "", x)
  ifelse(nchar(x) > 45, paste0(substr(x, 1, 42), "..."), x)
}

# ---- Read data ----
corr_df <- read.csv(corr_file, stringsAsFactors = FALSE)
coef_df <- read.csv(coef_file, row.names = 1, stringsAsFactors = FALSE)
qval_df <- read.csv(qval_file, row.names = 1, stringsAsFactors = FALSE)

# ---- Reshape correlation data to wide format ----
corr_mat <- corr_df %>%
  select(Taxa, Env, Correlation) %>%
  pivot_wider(names_from = Env, values_from = Correlation) %>%
  column_to_rownames("Taxa") %>%
  as.matrix()

sig_mat <- corr_df %>%
  select(Taxa, Env, Significance) %>%
  pivot_wider(names_from = Env, values_from = Significance) %>%
  column_to_rownames("Taxa") %>%
  as.matrix()

# ---- Match species between datasets ----
common <- intersect(rownames(corr_mat), rownames(coef_df))
corr_mat <- corr_mat[common, , drop = FALSE]
sig_mat  <- sig_mat[common, , drop = FALSE]
coef_df  <- coef_df[common, , drop = FALSE]
qval_df  <- qval_df[common, , drop = FALSE]

# ---- Column (Env) order and grouping ----
env_order <- c("TJC", "SJC", "Stiffness_time", "VAS",
               "CRP", "ESR",
               "RF", "HAQ", "PaGADA", "PhGADA",
               "DAS28_ESR", "DAS28_CRP", "cDAI", "sDAI",
               "Age", "BMI")
env_order <- intersect(env_order, colnames(corr_mat))

# Apply ordering
corr_mat <- corr_mat[, env_order, drop = FALSE]
sig_mat  <- sig_mat[, env_order, drop = FALSE]

# Env group annotation
env_groups <- c("TJC"="Joint", "SJC"="Joint", "Stiffness_time"="Joint", "VAS"="Joint",
                "CRP"="Inflammation", "ESR"="Inflammation",
                "RF"="Immune", "HAQ"="Immune", "PaGADA"="Immune", "PhGADA"="Immune",
                "DAS28_ESR"="Activity", "DAS28_CRP"="Activity", "cDAI"="Activity", "sDAI"="Activity",
                "Age"="Demographics", "BMI"="Demographics")
col_group <- factor(env_groups[env_order], levels = c("Joint", "Inflammation", "Immune", "Activity", "Demographics"))

# ---- Clean row names ----
taxa_clean <- clean_taxa(rownames(corr_mat))

# ---- Color functions ----
cor_col_fun  <- colorRamp2(c(-0.4, 0, 0.4), c("#1187A8", "white", "#ED8A10"))
coef_col_fun <- colorRamp2(c(-3, 0, 3), c("#4575B4", "#F7F7F7", "#D73027"))

# ---- Prepare right-side annotations ----
abund_sig <- sig_star(qval_df$MaAsLin3_abundance_qval)
prev_sig  <- sig_star(qval_df$MaAsLin3_prevalence_qval)
ancom_sig <- sig_star(qval_df$ANCOMBC2_qval)
lefse_sig <- sig_star(qval_df$LEfSe_Padj)

# ---- Row clustering ----
set.seed(42)
row_hclust <- hclust(dist(corr_mat), method = "ward.D2")

# ---- Column top annotation (env groups) ----
top_ha <- HeatmapAnnotation(
  Group = col_group,
  col = list(Group = c(
    "Joint" = "#E41A1C",
    "Inflammation" = "#377EB8",
    "Immune" = "#4DAF4A",
    "Activity" = "#984EA3",
    "Demographics" = "#FF7F00"
  )),
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(3, "mm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 7, lineheight = 0.9),
    labels_gp = gpar(fontsize = 6),
    grid_height = unit(3, "mm"),
    grid_width = unit(3, "mm")
  )
)

# ---- Main heatmap ----
ht_main <- Heatmap(
  corr_mat,
  name = "Correlation",
  col = cor_col_fun,
  cluster_rows = row_hclust,
  cluster_columns = FALSE,
  column_order = env_order,
  top_annotation = top_ha,
  row_labels = taxa_clean,
  row_names_gp = gpar(fontsize = 7),
  row_names_side = "left",
  row_dend_side = "left",
  row_dend_width = unit(20, "mm"),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  column_names_centered = FALSE,
  heatmap_legend_param = list(
    title = "Spearman",
    title_gp = gpar(fontsize = 7, lineheight = 0.9),
    labels_gp = gpar(fontsize = 6),
    grid_height = unit(3, "mm"),
    grid_width = unit(3, "mm"),
    at = c(-0.4, -0.2, 0, 0.2, 0.4)
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    s <- sig_mat[i, j]
    if (s != "" && !is.na(s)) {
      grid.text(s, x, y, gp = gpar(fontsize = 7, col = "grey0"))
    }
  },
  na_col = "grey95"
)

# ---- Helper: single-column right-side heatmap ----
make_right_ht <- function(col_data, col_sig, ht_name, legend_title) {
  m <- as.matrix(col_data)
  #title_text <- gsub("\n", " ", ht_name)
  title_text <- ht_name
  colnames(m) <- title_text
  Heatmap(m,
    col = coef_col_fun,
    show_row_names = FALSE,
    show_column_names = TRUE,
    show_row_dend = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    width = unit(0.84, "cm"),
    column_names_rot = 45,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 6),
    heatmap_legend_param = list(
      title = legend_title,
      title_gp = gpar(fontsize = 7, lineheight = 0.9),
      labels_gp = gpar(fontsize = 6),
      grid_height = unit(3, "mm"),
      grid_width = unit(3, "mm")
    ),
    cell_fun = function(j, i, x, y, w, h, fill) {
      s <- col_sig[i]
      if (s != "" && !is.na(s)) {
        grid.text(s, x, y, gp = gpar(fontsize = 7, col = "grey0"))
      }
    },
    na_col = "grey95")
}

ht_abund <- make_right_ht(
  coef_df[, "MaAsLin3_abundance_coef", drop = FALSE],
  abund_sig, "MaAsLin3\nAbundance.Coef", "MaAsLin3\nAbundance.Coef")
ht_prev  <- make_right_ht(
  coef_df[, "MaAsLin3_prevalence_coef", drop = FALSE],
  prev_sig, "MaAsLin3\nPrevalence.Coef", "MaAsLin3\nPrevalence.Coef")
ht_ancom <- make_right_ht(
  coef_df[, "ANCOMBC2_lfc", drop = FALSE],
  ancom_sig, "ANCOMBC2\nLogFC", "ANCOMBC2\nLogFC")

# ---- Right side: LEfSe LDA (bar plot) ----
lefse_vals <- coef_df$LEfSe_LDA
bar_colors <- ifelse(lefse_vals > 0, "#D73027", "#4575B4")
bar_colors[is.na(lefse_vals)] <- "grey90"

# Prepare significance labels for LEfSe (displayed as text next to bars)
# Set the sign of LDA values so positive = enriched, negative = depleted

ha_lefse <- rowAnnotation(
  `LEfSe LDA` = anno_barplot(
    lefse_vals,
    gp = gpar(fill = bar_colors),
    bar_width = 0.65,
    axis = TRUE,
    axis_param = list(gp = gpar(fontsize = 5)),
    width = unit(1.8, "cm")
  ),
  ` ` = anno_text(lefse_sig, gp = gpar(fontsize = 7), just = "left"),
  show_annotation_name = c(TRUE, FALSE),
  annotation_name_gp = gpar(fontsize = 7),
  gap = unit(1, "mm")
)

# ---- Combine ----
ht_list <- ht_main + ht_abund + ht_prev + ht_ancom + ha_lefse

# ---- Save ----
dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)

pdf(out_pdf, width = 8, height = nrow(corr_mat)/9)
draw(ht_list,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE,
     row_title = NULL,
     column_title = NULL)
dev.off()


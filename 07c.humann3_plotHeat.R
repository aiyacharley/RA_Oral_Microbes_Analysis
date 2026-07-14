# 07.humann3_plotHeat.R
# Figure 7: HUMAnN3 pathway differential abundance heatmap
#   - Main panel: (intersection of MaAsLin3∩ANCOMBC2∩LEfSe) pathway
#                 abundance z-score heatmap. Columns = samples, ordered by
#                 Time (HC / V1 / V4). Top annotation = Time, Age, BMI.
#   - Right panels: MaAsLin3 abundance coefficient, ANCOMBC2 LFC,
#                   LEfSe LDA with significance stars
# ============================================================================

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(scales))

# ---- Paths ----
outdir <- "output/Figure7/humann3pathway/"
ps_file       <- "00_data/ps_tongue_humann3pathway.rds"
maaslin3_file <- paste0(outdir, "sig_maaslin3.csv")
ancombc2_file <- paste0(outdir, "sig_ancombc2.csv")
lefse_file    <- paste0(outdir, "sig_lefse.csv")
out_pdf       <- paste0(outdir, "Fig7D_humann3pathway_heat.pdf")

# ---- Functions ----
sig_star <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
    ifelse(p < 0.01,  "**",
    ifelse(p < 0.05,   "*", ""))))
}

clean_pathway <- function(x) {
  y <- sub("^[^:]+: ", "", x)
  y <- gsub("&beta;", "beta", y, fixed = TRUE)
  y <- gsub("&alpha;", "alpha", y, fixed = TRUE)
  y <- gsub("&gamma;", "gamma", y, fixed = TRUE)
  ifelse(nchar(y) > 55, paste0(substr(y, 1, 52), "..."), y)
}

# ======================================================================
# 1. Read data
# ======================================================================
cat("Reading phyloseq object ...\n")
ps <- readRDS(ps_file)
otu_mat <- as.matrix(otu_table(ps))
sam_df  <- as(sample_data(ps), "data.frame")
cat("  OTU table:", nrow(otu_mat), "pathways x", ncol(otu_mat), "samples\n")
cat("  Time levels:", paste(unique(sam_df$Time), collapse = ", "), "\n")

cat("Reading differential analysis results ...\n")
maaslin3 <- read.csv(maaslin3_file, stringsAsFactors = FALSE)
ancombc2 <- read.csv(ancombc2_file, stringsAsFactors = FALSE)
lefse    <- read.csv(lefse_file, stringsAsFactors = FALSE)

# ======================================================================
# 2. Extract significant features from each method (for Time)
# ======================================================================
# MaAsLin3: already filtered for Time — all entries significant
m3_feat <- unique(maaslin3$feature)
cat("  MaAsLin3 significant features:", length(m3_feat), "\n")

# ANCOMBC2: TimeV1 passed differential test
ancom_feat <- unique(ancombc2$taxon[ancombc2$diff_TimeV1 == TRUE])
cat("  ANCOMBC2 significant (TimeV1):", length(ancom_feat), "\n")

# LEfSe: already filtered — all entries significant
lefse_feat <- unique(lefse$Taxa)
cat("  LEfSe significant features:", length(lefse_feat), "\n")

# ---- Intersection of all three methods ----
sig_pathways <- Reduce(intersect, list(m3_feat, ancom_feat, lefse_feat))
cat("  Intersection (MaAsLin3 & ANCOMBC2 & LEfSe):", length(sig_pathways), "\n")

if (length(sig_pathways) == 0) {
  stop("No pathways are significant in all three methods — cannot draw heatmap.")
}

# ======================================================================
# 3. Build abundance matrix (rows = intersection pathways)
# ======================================================================
# Match pathway names to OTU table
avail <- intersect(sig_pathways, rownames(otu_mat))
cat("  Available in OTU table:", length(avail), "\n")
if (length(avail) == 0) stop("No intersection pathways found in OTU table.")

abund_raw <- otu_mat[avail, , drop = FALSE]

# Z-score each row
abund_z <- t(scale(t(abund_raw)))
abund_z[is.nan(abund_z)] <- 0

# ======================================================================
# 4. Prepare sample annotation (top) and column ordering
# ======================================================================
# Time column: make sure it's a factor with desired order
time_levels <- intersect(c("HC", "V1", "V4"), unique(sam_df$Time))
sam_df$Time <- factor(sam_df$Time, levels = time_levels)

# Sort samples by Time (then optionally by Age / BMI within group)
sam_ord <- sam_df[order(sam_df$Time), , drop = FALSE]
sample_order <- rownames(sam_ord)

# Align abundance matrix columns
abund_z <- abund_z[, sample_order, drop = FALSE]

# ---- Top annotation ----
# Age/BMI use colour gradient (columns coloured by value)
age_fun  <- colorRamp2(c(18, 75), c("white", "#608C1B"))
bmi_fun  <- colorRamp2(c(15, 40), c("white", "#2EBFBF"))

top_ha <- HeatmapAnnotation(
  Time = sam_ord$Time,
  Age  = sam_ord$Age,
  BMI  = sam_ord$BMI,
  col = list(
    Time = c("HC" = "#44AF8F", "V1" = "#E07C30", "V4" = "#756FB3"),
    Age  = age_fun,
    BMI  = bmi_fun
  ),
  annotation_name_gp = gpar(fontsize = 11),
  gap = unit(1, "mm"),
  show_legend = FALSE,
  annotation_legend_param = list(
    Time = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 9),
                grid_height = unit(4, "mm"), grid_width = unit(4, "mm"))
  )
)

# Column split by Time
col_split <- sam_ord$Time

# ======================================================================
# 5. Build coef_df & qval_df (rows = intersection pathways)
# ======================================================================
# Separate MaAsLin3 abundance / prevalence
maaslin3_abund <- maaslin3 %>% filter(model == "abundance") %>%
  distinct(feature, .keep_all = TRUE)
maaslin3_prev  <- maaslin3 %>% filter(model == "prevalence") %>%
  distinct(feature, .keep_all = TRUE)

# LEfSe signed LDA
lefse$LDA_signed <- ifelse(lefse$Group == "V1", lefse$LDA, -lefse$LDA)

fill_vec <- function(features, values, pool) {
  v <- setNames(rep(NA_real_, length(pool)), pool)
  m <- match(features, pool)
  v[m[!is.na(m)]] <- values[!is.na(m)]
  v
}

coef_df <- data.frame(
  MaAsLin3_abundance_coef  = fill_vec(maaslin3_abund$feature, maaslin3_abund$coef, avail),
  MaAsLin3_prevalence_coef = fill_vec(maaslin3_prev$feature,  maaslin3_prev$coef,  avail),
  ANCOMBC2_lfc             = fill_vec(ancombc2$taxon,         ancombc2$lfc_TimeV1, avail),
  LEfSe_LDA                = fill_vec(lefse$Taxa,             lefse$LDA_signed,    avail),
  row.names = avail,
  stringsAsFactors = FALSE
)

qval_df <- data.frame(
  MaAsLin3_abundance_qval  = fill_vec(maaslin3_abund$feature, maaslin3_abund$qval_individual, avail),
  MaAsLin3_prevalence_qval = fill_vec(maaslin3_prev$feature,  maaslin3_prev$qval_individual,  avail),
  ANCOMBC2_qval            = fill_vec(ancombc2$taxon,         ancombc2$q_TimeV1,             avail),
  LEfSe_Padj               = fill_vec(lefse$Taxa,             lefse$P.adj,                   avail),
  row.names = avail,
  stringsAsFactors = FALSE
)

# ======================================================================
# 6. Colour functions
# ======================================================================
# For z-score abundance heatmap
z_range <- max(abs(abund_z), na.rm = TRUE)
abund_col_fun <- colorRamp2(c(-z_range, 0, z_range),
                             c("#1187A8", "white", "#ED8A10"))

# For coefficient heatmaps (shared scale)
coef_all <- na.omit(c(coef_df$MaAsLin3_abundance_coef,
                      coef_df$MaAsLin3_prevalence_coef,
                      coef_df$ANCOMBC2_lfc))
coef_lim <- max(abs(coef_all), na.rm = TRUE)
coef_col_fun <- colorRamp2(c(-coef_lim, 0, coef_lim),
                            c("#4575B4", "#F7F7F7", "#D73027"))

# ======================================================================
# 7. Right-side significance stars
# ======================================================================
abund_sig <- sig_star(qval_df$MaAsLin3_abundance_qval)
prev_sig  <- sig_star(qval_df$MaAsLin3_prevalence_qval)
ancom_sig <- sig_star(qval_df$ANCOMBC2_qval)
lefse_sig <- sig_star(qval_df$LEfSe_Padj)

# ======================================================================
# 8. Row clustering
# ======================================================================
set.seed(42)
row_hclust <- hclust(dist(abund_z), method = "ward.D2")
# column clustering within each Time group uses Euclidean + ward.D2

# ======================================================================
# 9. Clean pathway display names
# ======================================================================
pathway_display <- clean_pathway(rownames(abund_z))

# ======================================================================
# 10. Main heatmap (abundance z-score)
# ======================================================================
ht_main <- Heatmap(
  abund_z,
  name = "Z-score",
  col = abund_col_fun,
  cluster_rows = row_hclust,
  cluster_columns = TRUE,
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2",
  column_split = col_split,
  column_gap = unit(1.5, "mm"),
  top_annotation = top_ha,
  show_column_names = FALSE,
  row_labels = pathway_display,
  row_names_gp = gpar(fontsize = 10),
  show_row_dend = FALSE,
  row_names_side = "left",
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 10, lineheight = 0.9),
    labels_gp = gpar(fontsize = 9),
    grid_height = unit(4, "mm"),
    grid_width = unit(4, "mm")
  ),
  na_col = "grey95"
)

# ======================================================================
# 11. Helper: single-column right-side heatmap
# ======================================================================
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
    width = unit(1.2, "cm"),
    column_names_rot = 45,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 9),
    heatmap_legend_param = list(
      title = legend_title,
      title_gp = gpar(fontsize = 10, lineheight = 0.9),
      labels_gp = gpar(fontsize = 9),
      grid_height = unit(4, "mm"),
      grid_width = unit(4, "mm")
    ),
    cell_fun = function(j, i, x, y, w, h, fill) {
      s <- col_sig[i]
      if (s != "" && !is.na(s)) {
        grid.text(s, x, y, gp = gpar(fontsize = 9, col = "grey0"))
      }
    },
    na_col = "grey95")
}

# ---- Right panels ----
ht_abund <- make_right_ht(
  coef_df[, "MaAsLin3_abundance_coef", drop = FALSE],
  abund_sig, "MaAsLin3\nAbundance.Coef", "MaAsLin3\nAbundance.Coef")
#ht_prev  <- make_right_ht(   # removed
#  coef_df[, "MaAsLin3_prevalence_coef", drop = FALSE],
#  prev_sig, "MaAsLin3\nPrevalence.Coef", "MaAsLin3\nPrevalence.Coef")
ht_ancom <- make_right_ht(
  coef_df[, "ANCOMBC2_lfc", drop = FALSE],
  ancom_sig, "ANCOMBC2\nLogFC", "ANCOMBC2\nLogFC")

# ======================================================================
# 12. LEfSe LDA bar plot (right side)
# ======================================================================
lefse_vals    <- coef_df$LEfSe_LDA
bar_colors    <- ifelse(lefse_vals > 0, "#D73027", "#4575B4")
bar_colors[is.na(lefse_vals)] <- "grey90"

ha_lefse <- rowAnnotation(
  `LEfSe LDA` = anno_barplot(
    lefse_vals,
    gp = gpar(fill = bar_colors),
    bar_width = 0.65,
    axis = TRUE,
    axis_param = list(gp = gpar(fontsize = 7)),
    width = unit(2.7, "cm")
  ),
  ` ` = anno_text(lefse_sig, gp = gpar(fontsize = 9), just = "left"),
  show_annotation_name = c(TRUE, FALSE),
  annotation_name_gp = gpar(fontsize = 9),
  gap = unit(1, "mm")
)

# ======================================================================
# 13. Combine & save
# ======================================================================
ht_list <- ht_main + ht_abund + ht_ancom + ha_lefse

dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)

n_rows <- nrow(abund_z)
pdf_height <- max(3.5, n_rows / 7 + 0.8)   # extra headroom for top annotation
pdf_width  <- 15

pdf(out_pdf, width = pdf_width, height = pdf_height)
draw(ht_list,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE,
     padding = unit(c(1, 6, 1, 1), "cm"),
     row_title = NULL,
     column_title = NULL)
dev.off()

cat("\nDone! PDF saved to:", out_pdf, "\n")
cat("  Pathways (intersection):", n_rows, "\n")
cat("  Samples:", ncol(abund_z), "\n")
cat("  Time groups:", paste(levels(sam_ord$Time), collapse = ", "), "\n")

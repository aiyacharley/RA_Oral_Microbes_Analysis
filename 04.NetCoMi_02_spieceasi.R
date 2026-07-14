# C:\Program Files\R\R-4.5.2
setwd("D:/我的工作.省中医/Projects/刘院士团队舌苔微生态/results_20260615")
library(NetCoMi)
library(igraph)
library(readxl)
library(ggvenn)
library(patchwork)
library(dplyr)
library(tidyr)
library(RColorBrewer)
source('00.func.R')

#outdir <- "output/Figure4/netConstruct_spring/"
outdir <- "output/Figure4/netConstruct_spieceasi/"
#outdir <- "output/Figure4/netConstruct_spearman/"
#outdir <- "output/Figure4/netConstruct_pearson/"

load(file = paste0(outdir,"/netConstruct.Rdata"))

# diff bacteria
diff_HCV1 <- read.csv(file = "output/Figure2/Fig3A_venn.csv", row.names = 1)
diff_V1V4 <- read.csv(file = "output/Figure3/Fig3B_venn.csv", row.names = 1)

target_bacteria2 <- intersect(rownames(diff_HCV1), rownames(diff_V1V4)) 
target_bacteria <- unique(c(rownames(diff_HCV1), rownames(diff_V1V4))) 
nodeMAGs <- intersect(rownames(ps.filter@otu_table), target_bacteria)
nodeMAGs2 <- intersect(rownames(ps.filter@otu_table), target_bacteria2)
print(nodeMAGs)

#-------------------------------- netAnalyze ----------------
run_netAnalyze <- function(net){
  netAnalyze(
    net,
    centrLCC = TRUE,
    clustMethod = "cluster_fast_greedy",
    hubPar = c("eigenvector","betweenness","closeness","strength"),
    hubQuant = 0.95,
    weightDeg = TRUE,
    normDeg = FALSE,
    normBetw = TRUE,
    normClose = TRUE,
    normEigen = TRUE,
    avDissIgnoreInf = TRUE
  )
}
ana.HC   <- run_netAnalyze(net.HC)
ana.V1   <- run_netAnalyze(net.V1)
ana.V4   <- run_netAnalyze(net.V4)
ana.HCV1 <- run_netAnalyze(net.HCV1)
ana.V1V4 <- run_netAnalyze(net.V1V4)
ana.HCV4 <- run_netAnalyze(net.HCV4)

#----------------------- netCompare --------------------
#netCompare.HCV1 <- netCompare(ana.HCV1, permTest = TRUE, nPerm = 1000L,adjust = "adaptBH", cores = 8, seed = 42);saveRDS(netCompare.HCV1, file = paste0(outdir,"/netCompare.HCV1.rds"))
#netCompare.V1V4 <- netCompare(ana.V1V4, permTest = TRUE, nPerm = 1000L,adjust = "adaptBH", cores = 8, seed = 42);saveRDS(netCompare.V1V4, file = paste0(outdir,"/netCompare.V1V4.rds"))
#netCompare.HCV4 <- netCompare(ana.HCV4, permTest = TRUE, nPerm = 1000L,adjust = "adaptBH", cores = 8, seed = 42);saveRDS(netCompare.HCV4, file = paste0(outdir,"/netCompare.HCV4.rds"))
## adjust = "adaptBH" 校正后p值均>0.05
netCompare.HCV1 <- readRDS(file = paste0(outdir,"/netCompare.HCV1.rds"))
netCompare.V1V4 <- readRDS(file = paste0(outdir,"/netCompare.V1V4.rds"))
netCompare.HCV4 <- readRDS(file = paste0(outdir,"/netCompare.HCV4.rds"))
comp_net_summary.HCV1 <- summary(netCompare.HCV1, groupNames = c("HC", "V1"), showCentr = "all", numbNodes = 300)
comp_net_summary.V1V4 <- summary(netCompare.V1V4, groupNames = c("V1", "V4"), showCentr = "all", numbNodes = 300)
comp_net_summary.HCV4 <- summary(netCompare.HCV4, groupNames = c("HC", "V4"), showCentr = "all", numbNodes = 300)


tax_table <- data.frame(ps.filter@tax_table@.Data)
all_phyla <- table(tax_table$Phylum) %>% sort(decreasing = T) %>% names()
phylum_colors <- setNames(
  RColorBrewer::brewer.pal(length(all_phyla), "Dark2"),
  all_phyla
)
c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
tax_table$color <- phylum_colors[tax_table$Phylum]
tmp_col <- tax_table[,8]
names(tmp_col) <- tax_table$Species

pdf(file = paste0(outdir,"/network_legend.pdf"), width = 3, height = 5)
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
# 添加图例
legend(x = 0, y = 0.8,
       legend = c("Negative", "Positive"),  # 图例标签
       col = c("#F18C8D", "#B2B2B2"),       # 线条颜色
       lty = c(1, 1),                       # 线条类型（1=实线，2=虚线）
       lwd = 2,                             # 线条宽度
       title = "Interaction",               # 图例标题
       title.adj = 0.2,                    # 标题居中
       title.font = 2,                     # 标题加粗
       cex = 0.9,                          # 文字大小
       bty = "n",                          # 显示边框
       bg = "white")                       # 白色背景
# 添加专业风格的图例
legend(x = 0, y = 0.6,
       legend = gsub("p__", "", names(phylum_colors)),  # 去掉了"p__"前缀
       col = phylum_colors,
       pch = 19,           # 更大更饱满的点
       pt.cex = 2,
       cex = 0.9,
       title = "Phylum",
       title.adj = 0.1,
       title.font = 2,     # 标题加粗
       bty = "n",          # 无边框，更清爽
       bg = "white")
dev.off()

pdf(file = paste0(outdir,"/network_net2HCV1.pdf"), width = 10, height = 5)
plot(ana.HCV1,groupNames = c('HC','V1'), 
     #nodeFilter = "names", nodeFilterPar = nodeMAGs,
     layoutGroup = "union", sameLayout = TRUE, layout = "spring",
     edgeFilter = "threshold", edgeFilterPar = 0,
     nodeSize = "mclr", 
     nodeSizeSpread = 2, 
     nodeTransp = 0,
     nodeColor = "colorVec",colorVec = tmp_col,
     cexTitle = 1, cexLabels = 0, labelScale=F,
     posCol = "#B2B2B2", negCol = "red",
     edgeTranspLow = 60, edgeTranspHigh = 60,
     mar = c(1,2,3,2),
     rmSingles = "none")
# plot(ana.HCV1,groupNames = c('HC','V1'),
#      nodeFilter = "names", nodeFilterPar = nodeMAGs,
#      layoutGroup = "union", sameLayout = TRUE, layout = "circle",
#      edgeFilter = "threshold", edgeFilterPar = 0,
#      nodeSize = "mclr",
#      nodeSizeSpread = 2,
#      nodeTransp = 0,
#      nodeColor = "colorVec",colorVec = tmp_col,
#      cexTitle = 1, cexLabels = 0.3, labelScale=F,
#      posCol = "#B2B2B2", negCol = "red",
#      edgeTranspLow = 60, edgeTranspHigh = 60,
#      mar = c(1,2,3,2),
#      rmSingles = "none")
dev.off()
pdf(file = paste0(outdir,"/network_net2V1V4.pdf"), width = 10, height = 5)
plot(ana.V1V4,groupNames = c('V1 (paired)','V4'), 
     #nodeFilter = "names", nodeFilterPar = nodeMAGs,
     layoutGroup = "union", sameLayout = TRUE, layout = "spring",
     edgeFilter = "threshold", edgeFilterPar = 0,
     nodeSize = "mclr", 
     nodeSizeSpread = 2, 
     nodeTransp = 0,
     nodeColor = "colorVec",colorVec = tmp_col,
     cexTitle = 1, cexLabels = 0, labelScale=F,
     posCol = "#B2B2B2", negCol = "red",
     edgeTranspLow = 60, edgeTranspHigh = 60,
     mar = c(1,2,3,2),
     rmSingles = "none")
# plot(ana.V1V4,groupNames = c('V1 (paired)','V4'), 
#      nodeFilter = "names", nodeFilterPar = nodeMAGs,
#      layoutGroup = "union", sameLayout = TRUE, layout = "circle",
#      edgeFilter = "threshold", edgeFilterPar = 0,
#      nodeSize = "mclr", 
#      nodeSizeSpread = 2, 
#      nodeTransp = 0,
#      nodeColor = "colorVec",colorVec = tmp_col,
#      cexTitle = 1, cexLabels = 0.3, labelScale=F,
#      posCol = "#B2B2B2", negCol = "red",
#      edgeTranspLow = 60, edgeTranspHigh = 60,
#      mar = c(1,2,3,2),
#      rmSingles = "none")
dev.off()


## 比较两个GCM
gobj1 <- calcGCM(net.HCV1$adjaMat1)
gobj2 <- calcGCM(net.HCV1$adjaMat2)
gcmtest <- testGCM(gobj1, gobj2, verbose = FALSE)
pdf(file = paste0(outdir,"/network_net2HCV1_GCMs.pdf"), width = 6, height = 6)
plotHeat(mat = gcmtest$diff, pmat = gcmtest$pAdjustDiff,
         type = "mixed", title = "Differences between GCs (HC-V1)",
         mar = c(0, 0, 2, 0))
dev.off()
## 比较两个GCM
gobj1 <- calcGCM(net.V1V4$adjaMat1)
gobj2 <- calcGCM(net.V1V4$adjaMat2)
gcmtest <- testGCM(gobj1, gobj2, verbose = FALSE)
pdf(file = paste0(outdir,"/network_net2V1V4_GCMs.pdf"), width = 6, height = 6)
plotHeat(mat = gcmtest$diff, pmat = gcmtest$pAdjustDiff,
         type = "mixed", title = "Differences between GCs (V1-V4, paired)",
         mar = c(0, 0, 2, 0))
dev.off()


library(igraph)
g0 <- graph_from_adjacency_matrix(net.HCV1$adjaMat1, 
                                 mode = "undirected", 
                                 weighted = TRUE,
                                 diag = FALSE)
g1 <- graph_from_adjacency_matrix(net.HCV1$adjaMat2, 
                                  mode = "undirected", 
                                  weighted = TRUE,
                                  diag = FALSE)
g1p <- graph_from_adjacency_matrix(net.V1V4$adjaMat1, 
                                  mode = "undirected", 
                                  weighted = TRUE,
                                  diag = FALSE)
g4 <- graph_from_adjacency_matrix(net.V1V4$adjaMat2, 
                                  mode = "undirected", 
                                  weighted = TRUE,
                                  diag = FALSE)
hist(net.HC$edgelist1$asso)
hist(net.V1$edgelist1$asso)
hist(net.V4$edgelist1$asso)
g0_coreness <- igraph::coreness(g0)
g1_coreness <- igraph::coreness(g1)
g1p_coreness <- igraph::coreness(g1p)
g4_coreness <- igraph::coreness(g4)
hist(g0_coreness)
hist(g1_coreness)
hist(g4_coreness)

kcore.HC <- 100*table(g0_coreness)/sum(table(g0_coreness))
kcore.V1 <- 100*table(g1_coreness)/sum(table(g1_coreness))
kcore.V4 <- 100*table(g4_coreness)/sum(table(g4_coreness))

# 将每个 kcore 数据转换为数据框，并添加 group 列
df_hc <- data.frame(kcore = names(kcore.HC), percentage = as.numeric(kcore.HC), group = "HC")
df_v1 <- data.frame(kcore = names(kcore.V1), percentage = as.numeric(kcore.V1), group = "V1")
df_v4 <- data.frame(kcore = names(kcore.V4), percentage = as.numeric(kcore.V4), group = "V4")

# 使用 bind_rows() 纵向合并，自动处理列名
merged_df <- bind_rows(df_hc, df_v1, df_v4, .id = "source")
merged_df$kcore <- as.numeric(merged_df$kcore)

merged_df$kcore <- factor(merged_df$kcore, levels = c(2,3,4,5))

fig_kcore <- ggplot(merged_df, aes(x = kcore, y = percentage, fill = group)) +
  geom_col(width = 0.9) +  # geom_col() 更直接
  # 新增：在每个分面的最上方显示该 group 的总和
  stat_summary(
    aes(label = after_stat(round(y, 1)), group = 1),
    fun = sum,
    geom = "text",
    vjust = -0.5,      # 调整垂直位置，使其在所有柱子之上
    size = 3,
    color = "black",
    fontface = "bold"
  ) +
  facet_wrap(~ group, ncol = 3) +  # 一列三行，竖向排列
  labs(
    x = "K-core Value",
    y = "Percentage (%)"
  ) +
  theme_test() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.text = element_text(face = "bold", size = 12),  # 分面标题样式
    plot.title = element_text(hjust = 0.5)
  )+
  scale_color_manual(values = c("HC" = "#44AF8F", "V1" = "#E07C30", "V4" = "#756FB3")) +
  scale_fill_manual(values = c("HC" = "#44AF8F", "V1" = "#E07C30", "V4" = "#756FB3"))

write.csv(merged_df, file = paste0(outdir,"/kcore_distribution.csv"),row.names = T)
ggsave(filename = paste0(outdir,'/kcore_distribution.pdf'),width = 10, height = 3)


# Zi-Pi plot
## 标签太多，设置show_text = F隐藏节点文本
abs_r <- 0
res.HC <- plot_zipi(ps = ps.HC, net = net.HC, abs_r = abs_r, title="HC", show_text = F)
res.V1 <- plot_zipi(ps = ps.V1, net = net.V1, abs_r = abs_r, title="V1", show_text = F)
res.V1p <- plot_zipi(ps = ps.V1p, net = net.V1p, abs_r = abs_r, title="V1 (paired)", show_text = F)
res.V4 <- plot_zipi(ps = ps.V4, net = net.V4, abs_r = abs_r, title="V4", show_text = F)
fig_ZiPi <- res.HC$figure + res.V1$figure + res.V4$figure

ggsave(filename = paste0(outdir,"/ZiPi.pdf"), plot = fig_ZiPi, width = 20, height = 5)
write.csv(res.HC$table, file = paste0(outdir,"/ZiPi_HC.csv"),row.names = F)
write.csv(res.V1$table, file = paste0(outdir,"/ZiPi_V1.csv"),row.names = F)
write.csv(res.V4$table, file = paste0(outdir,"/ZiPi_V4.csv"),row.names = F)

# 网络属性热图
network_attr <- as.data.frame(list(HC=res.HC$network_attr,V1=res.V1$network_attr,V4=res.V4$network_attr))
write.csv(network_attr, file = paste0(outdir,"/network_attr_heatmap.csv"),row.names = T)
pdf(file = paste0(outdir,"/network_attr_heatmap.pdf"), width = 4, height = 3)
pheatmap::pheatmap(network_attr[-c(1),], 
                   scale = 'row', cluster_rows = T, cluster_cols = T,
                   display_numbers = round(network_attr[-c(1),],2))
dev.off()


#-------------------------------- centralities -----------------------------
# 同样假设列表名为 my_list
null2zero <- function(my_list){
  for (i in 1:length(my_list)) {      # 遍历列表的每个索引
    if (is.null(my_list[[i]])) {      # 判断当前元素是否为NULL
      my_list[[i]] <- rep(NA, 283)     # 如果是，则替换为长度为10的NA向量
    }
  }
  return(my_list)
}
df_centralities_HC <- as.data.frame(null2zero(ana.HC$centralities))
df_centralities_V1 <- as.data.frame(null2zero(ana.V1$centralities))
df_centralities_V4 <- as.data.frame(null2zero(ana.V4$centralities))

sheet1_centralities <- list('HC' = df_centralities_HC[nodeMAGs,c("degree1","between1","close1","eigenv1")],
                           'V1' = df_centralities_V1[nodeMAGs,c("degree1","between1","close1","eigenv1")],
                           'V4' = df_centralities_V4[nodeMAGs,c("degree1","between1","close1","eigenv1")])
openxlsx::write.xlsx(sheet1_centralities, file = paste0(outdir,'/centralities_targetMAGs.xlsx'),row.names=T)

df_centralities_HCV1 <- as.data.frame(ana.HCV1$centralities)
df_centralities_V1V4 <- as.data.frame(ana.V1V4$centralities)

sheet_centralities <- list('HC' = df_centralities_HC,
                           'V1' = df_centralities_V1,
                           'V4' = df_centralities_V4,
                           'HCV1' = df_centralities_HCV1,
                           'V1V4' = df_centralities_V1V4)
openxlsx::write.xlsx(sheet_centralities, file = paste0(outdir,'/centralities.xlsx'),row.names=T)


# network property barplot
plot_network_index <- function(centralities, network_index=c("degree","between","close","eigenv")){
  #network_index <- "degree" # 选择网络属性，例如 "degree1", "between1", "close1", "eigenv1"
  network_index1 <- paste0(network_index,'1')
  
  xx <- cbind(centralities$HC[,network_index1,drop=F],
              centralities$V1[,network_index1,drop=F],
              centralities$V4[,network_index1,drop=F])
  colnames(xx) <- c("HC","V1","V4")
  # 转换为长格式（物种、分组、丰度）
  df_long <- xx %>%
    tibble::rownames_to_column("Species") %>% 
    pivot_longer(cols = c("HC","V1","V4"),
                 names_to = "Group",
                 values_to = "Value")
  
  df_long$Species <- factor(df_long$Species, levels = rev(nodeMAGs))
  
  fig <- ggplot(df_long, aes(x=Species,y=Value, fill=Group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    facet_grid(~Group, scales = "free_y") +
    theme_classic() +
    coord_flip()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "", x = "", y = network_index) +
    scale_fill_manual(values = c("HC" = "#44AF8F", "V1" = "#E07C30", "V4" = "#756FB3")) + 
    theme(legend.position = "none")
  return(fig)
}

fig_degree <- plot_network_index(sheet1_centralities, network_index = "degree")
fig_between <- plot_network_index(sheet1_centralities, network_index = "between")
fig_close <- plot_network_index(sheet1_centralities, network_index = "close")
fig_eigenv <- plot_network_index(sheet1_centralities, network_index = "eigenv")

fig_network_index1 <- fig_degree+
  fig_between+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+
  fig_close+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+
  fig_eigenv + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+
  plot_layout(ncol = 4)

ggsave(filename = paste0(outdir,"/network_index_targetMAGs.pdf"), plot = fig_network_index1, width = 15, height = 10)


library(ggtern)
tmpHC <- net.HC$edgelist1
tmpV1 <- net.V1$edgelist1
tmpV4 <- net.V4$edgelist1
tmpHC$edges <- paste0(tmpHC$v1,";",tmpHC$v2)
tmpV1$edges <- paste0(tmpV1$v1,";",tmpV1$v2)
tmpV4$edges <- paste0(tmpV4$v1,";",tmpV4$v2)
tmpHC <- tmpHC[,c("edges","asso")]
tmpV1 <- tmpV1[,c("edges","asso")]
tmpV4 <- tmpV4[,c("edges","asso")]

df1 <- merge(tmpHC, tmpV1, by = 'edges', all = T)
data <- merge(df1, tmpV4, by = 'edges', all = T)
colnames(data) <- c("edges","HC.assoc","V1.assoc","V4.assoc")
head(data)
data <- data %>%
  mutate(
    labels = paste0(
      case_when(
        is.na(HC.assoc) ~ "X",
        HC.assoc > 0 ~ "P",
        HC.assoc < 0 ~ "N",
        TRUE ~ "X"  # 处理可能的0值（根据sign函数，0返回0，但这里按需求归为X或可自定义）
      ),
      case_when(
        is.na(V1.assoc) ~ "X",
        V1.assoc > 0 ~ "P",
        V1.assoc < 0 ~ "N",
        TRUE ~ "X"
      ),
      case_when(
        is.na(V4.assoc) ~ "X",
        V4.assoc > 0 ~ "P",
        V4.assoc < 0 ~ "N",
        TRUE ~ "X"
      )
    )
  )
write.csv(data, file = paste0(outdir,"/Assoc3group.csv"),row.names = F)
data[is.na(data)] <- 0
data$HC.assoc <- abs(data$HC.assoc)
data$V1.assoc <- abs(data$V1.assoc)
data$V4.assoc <- abs(data$V4.assoc)

zz <- table(data$labels) %>% sort(decreasing = T)
zz2 <- paste(names(zz),"(",zz,")", sep = "")
data$labels <- factor(data$labels, levels = names(zz), labels = zz2)
col <- expand_colors(n=30)
# 基础三元图
p <- ggtern(data = data, aes(x = HC.assoc, y = V1.assoc, z = V4.assoc)) +
  geom_mask() +  # 添加三元图边界
  geom_point(aes(color = labels), alpha = 0.6) +  # 固定气泡大小
  scale_size(range = c(3, 3), guide = "none") +  # 固定大小为3mm，隐藏图例
  theme_tropical()+
  scale_color_manual(values = col)
p
ggsave(filename = paste0(outdir,"Assoc3group.pdf"),plot = p,width = 10, height = 5)






#----------------------------- no run ---------------------------------------
# 假设 nodeMAGsP 和 nodeMAGsN 分别是阳性/阴性节点集
netAnalyze.tmp <- netAnalyze.HCV1
nodeColors <- c(
  setNames(rep("#E5000C", length(nodeMAGsP)), nodeMAGsP),
  setNames(rep("#4477AA", length(nodeMAGsN)), nodeMAGsN)
)
xx <- setdiff(rownames(netAnalyze.tmp$input$assoMat1), nodeMAGs)
desired_order <- c(nodeMAGsP, nodeMAGsN, xx)

generate_two_circle_layout <- function(n_left = 11, n_right = 17, radius_left = 1, radius_right=1, center_dist = 2.5, nodeNames=NULL){
  
  # 1. 计算左侧圆圈上节点的坐标
  # 原理：将圆周等分，计算每个等分点上的正弦和余弦值作为坐标[1](@ref)
  angles_left <- seq(0, 2 * pi, length.out = n_left + 1)[-(n_left + 1)] # 生成n_left个等分角度
  x_left <- sin(angles_left) * radius_left - center_dist / 2 # X坐标：正弦值乘以半径，并向左偏移
  y_left <- cos(angles_left) * radius_left # Y坐标：余弦值乘以半径
  
  # 2. 计算右侧圆圈上节点的坐标
  angles_right <- seq(0, 2 * pi, length.out = n_right + 1)[-(n_right + 1)] # 生成n_right个等分角度
  x_right <- sin(angles_right) * radius_right + center_dist / 2 # X坐标：正弦值乘以半径，并向右偏移
  y_right <- cos(angles_right) * radius_right # Y坐标：余弦值乘以半径
  
  # 3. 合并坐标，并确保左侧节点在前，右侧节点在后
  layout_matrix <- matrix(c(x_left, x_right, y_left, y_right), ncol = 2, byrow = FALSE)
  
  # 4. 为行命名（可选，便于识别）
  rownames(layout_matrix) <- nodeNames
  colnames(layout_matrix) <- c("x", "y")
  
  return(layout_matrix)
}

# 生成你需要的布局 (11个节点在左，17个节点在右)
my_layout <- generate_two_circle_layout(n_left = 11, n_right = 17, radius_left = 0.5, radius_right = 1, center_dist = 2.5,nodeNames = nodeMAGs)

# 空跑一次提取图里的实际节点名
base_plot <- plot(netAnalyze.tmp, 
                  layout = "circle", # 随便给个默认布局先跑着
                  groupNames = c('HC','V1'), 
                  sameLayout = TRUE, 
                  layoutGroup = "union",
                  nodeFilter = "names", 
                  nodeFilterPar = nodeMAGs, 
                  doPlot = FALSE)

actual_nodes <- base_plot$labels[[1]] 

# 按照底层实际顺序，提取我们刚刚算好的XY坐标变成矩阵
final_layout_matrix <- as.matrix(my_layout[actual_nodes, c("x", "y")])

pdf(file = paste0(outdir,"/network_targetMAGs.pdf"), width = 11, height = 6)
fig_netHC <- plot(netAnalyze.HC,
                  title1 = "HC", showTitle = TRUE, 
                  layout = final_layout_matrix,
                  borderCol = "gray70", 
                  nodeSize = "mclr", 
                  nodeSizeSpread = 3, 
                  cexTitle = 2, cexLabels = 7,labelScale=TRUE,
                  mar = c(1,2,3,2), 
                  rmSingles = "none",
                  nodeFilter = "names", nodeFilterPar = nodeMAGs,
                  nodeTransp = 10,hubBorderCol = "gray",
                  nodeColor = "colorVec",colorVec = nodeColors,
                  edgeFilter = "threshold", edgeFilterPar = 0.6)
fig_netV1 <- plot(netAnalyze.V1,
                  title1 = "V1", showTitle = TRUE, 
                  layout = final_layout_matrix,
                  borderCol = "gray70", 
                  nodeSize = "mclr", 
                  nodeSizeSpread = 3, 
                  cexTitle = 2, cexLabels = 7,labelScale=TRUE,
                  mar = c(1,2,3,2), 
                  rmSingles = "none",
                  nodeFilter = "names", nodeFilterPar = nodeMAGs,
                  nodeTransp = 10,hubBorderCol = "gray",
                  nodeColor = "colorVec",colorVec = nodeColors,
                  edgeFilter = "threshold", edgeFilterPar = 0.6)
fig_netV4 <- plot(netAnalyze.V4,
                  title1 = "V4", showTitle = TRUE, 
                  layout = final_layout_matrix,
                  borderCol = "gray70", 
                  nodeSize = "mclr", 
                  nodeSizeSpread = 3, 
                  cexTitle = 2, cexLabels = 7,
                  mar = c(1,2,3,2), 
                  rmSingles = "none",
                  nodeFilter = "names", nodeFilterPar = nodeMAGs,
                  nodeTransp = 10,hubBorderCol = "gray",
                  nodeColor = "colorVec",colorVec = nodeColors,
                  edgeFilter = "threshold", edgeFilterPar = 0.6)
dev.off()

pdf(file = paste0(outdir,"/network_allMAGs.pdf"), width = 8, height = 8)
plot(netAnalyze.HC,
     title1 = "HC", showTitle = TRUE, 
     borderCol = "gray70", 
     nodeSize = "mclr", 
     nodeSizeSpread = 3, 
     cexTitle = 2, cexLabels = 4,labelScale=TRUE,
     mar = c(1,2,3,2), 
     rmSingles = T, #"none",
     nodeTransp = 10,hubBorderCol = "gray",
     colorVec = col,nodeColor = "cluster",
     edgeFilter = "threshold", edgeFilterPar = 0.6)
plot(netAnalyze.V1,
     title1 = "V1", showTitle = TRUE, 
     borderCol = "gray70", 
     nodeSize = "mclr", 
     nodeSizeSpread = 3, 
     cexTitle = 2, cexLabels = 4,labelScale=TRUE,
     mar = c(1,2,3,2), 
     rmSingles = T, #"none",
     nodeTransp = 10,hubBorderCol = "gray",
     colorVec = col,
     edgeFilter = "threshold", edgeFilterPar = 0.6)
plot(netAnalyze.V4,
     title1 = "V4", showTitle = TRUE, 
     borderCol = "gray70", 
     nodeSize = "mclr", 
     nodeSizeSpread = 3, 
     cexTitle = 2, cexLabels = 4,labelScale=TRUE,
     mar = c(1,2,3,2), 
     rmSingles = T, #"none",
     nodeTransp = 10,hubBorderCol = "gray",
     colorVec = col,
     edgeFilter = "threshold", edgeFilterPar = 0.6)
dev.off()

#----------------------- netCompare summary--------------------
comp_net_summary.HCV1 <- summary(netCompare.HCV1, groupNames = c("HC", "V1"), showCentr = "all", numbNodes = 300)
comp_net_summary.V1V4 <- summary(netCompare.V1V4, groupNames = c("V1", "V4"), showCentr = "all", numbNodes = 300)

capture.output(print(comp_net_summary.HCV1), file = paste0(outdir,"/netCompare.HCV1.txt"))
capture.output(print(comp_net_summary.V1V4), file = paste0(outdir,"/netCompare.V1V4.txt"))



# Differential network construction
diff_net.V4N_V4R <- diffnet(net.V4N_V4R, diffMethod = "fisherTest", adjust = "lfdr")
plot(diff_net.V4N_V4R,cexNodes = 1,cexLegend = 0.8,cexTitle = 2,mar = c(2,2,8,5),legendGroupnames = c("V4N","V4R"))
diffmat_sums <- rowSums(diff_net.V4N_V4R$diffAdjustMat)
diff_names <- names(diffmat_sums[diffmat_sums > 0])



plot_netcomi <- function(net, groupNames = "", outPDF = "./prefix", targetMAGs = NULL, alpha = 0.05){
  col <- expand_colors(n=20)
  pdf(file = paste0(outPDF,'.alpha',alpha,'.pdf'), width = 12, height = 8, onefile = TRUE)
  net.Analyze <- netAnalyze(net, clustMethod = "cluster_fast_greedy", hubPar = c("degree", "between", "closeness"),hubQuant = 0.9)
  
  # Differential network construction
  #diff_net <- diffnet(net, diffMethod = "fisherTest", alpha = 0.05, adjust = 'lfdr', lfdrThresh = 0.2, seed = 1234)
  diff_net <- diffnet(net, diffMethod = "fisherTest", alpha = alpha, adjust = 'none', lfdrThresh = 0.2, seed = 1234)
  plot(diff_net,cexNodes = 1,cexLegend = 0.8,cexTitle = 2,mar = c(2,2,5,2),legendGroupnames = groupNames)
  diffmat_sums <- rowSums(diff_net$diffAdjustMat)
  diff_names <- names(diffmat_sums[diffmat_sums > 0])

  fig_diff_net <- plot(net.Analyze,
                       groupNames = groupNames, 
                       sameLayout = TRUE, 
                       layoutGroup = "union", 
                       colorVec = col,
                       borderCol = "gray70", 
                       nodeSize = "mclr", 
                       nodeSizeSpread = 3, 
                       showTitle = TRUE, 
                       cexTitle = 2, , cexLabels = 2,
                       mar = c(1,2,3,1), 
                       repulsion = 1, 
                       labels = TRUE, 
                       rmSingles = "inboth",
                       nodeFilter = "names", nodeFilterPar = diff_names,
                       nodeTransp = 30, 
                       hubTransp = 10,
                       hubBorderCol = "black",
                       curve = 0,
                       curveAll = TRUE)
  if (!is.null(targetMAGs)){
    fig_target_net <- plot(net.Analyze,
                           groupNames = groupNames, 
                           sameLayout = TRUE, 
                           layoutGroup = "union", 
                           colorVec = col,
                           borderCol = "gray70", 
                           nodeSize = "mclr", 
                           nodeSizeSpread = 3, 
                           showTitle = TRUE, 
                           cexTitle = 2, , cexLabels = 2,
                           mar = c(1,2,3,1), 
                           repulsion = 1, 
                           labels = TRUE, 
                           rmSingles = "inboth",
                           nodeFilter = "names", nodeFilterPar = targetMAGs,
                           nodeTransp = 30, 
                           hubTransp = 10,
                           hubBorderCol = "black",
                           curve = 0,
                           curveAll = TRUE)
    dfvenn <- list(Maaslin3=targetMAGs, diffnet=diff_names)
    print(ggvenn(dfvenn))
    MAGs_intersect <- intersect(targetMAGs, diff_names)
    if (length(MAGs_intersect)>1){
      fig_target_net <- plot(net.Analyze,
                             groupNames = groupNames, 
                             sameLayout = TRUE, 
                             layoutGroup = "union", 
                             colorVec = col,
                             borderCol = "gray70", 
                             nodeSize = "mclr", 
                             nodeSizeSpread = 3, 
                             showTitle = TRUE, 
                             cexTitle = 2, , cexLabels = 2,
                             mar = c(1,2,3,1), 
                             repulsion = 1, 
                             labels = TRUE, 
                             rmSingles = "inboth",
                             nodeFilter = "names", nodeFilterPar = MAGs_intersect,
                             nodeTransp = 30, 
                             hubTransp = 10,
                             hubBorderCol = "black",
                             curve = 0,
                             curveAll = TRUE)
    }
  }
  dev.off()
}
plot_netcomi(net = net.HCV1, groupNames = c("HC","V1"), targetMAGs = nodeMAGs, outPDF = paste0(outdir,'/netcomi_',paste0(c("HC","V1"),collapse = '_')), alpha = 0.01)
plot_netcomi(net = net.V1V4, groupNames = c("V1","V4"), targetMAGs = nodeMAGs, outPDF = paste0(outdir,'/netcomi_',paste0(c("V1","V4"),collapse = '_')), alpha = 0.01)

diff_netN <- diffnet(net.V1N_V4N, diffMethod = "fisherTest", alpha = 0.01, adjust = 'none', lfdrThresh = 0.2, seed = 1234)
diff_netR <- diffnet(net.V1R_V4R, diffMethod = "fisherTest", alpha = 0.01, adjust = 'none', lfdrThresh = 0.2, seed = 1234)
diff_net1 <- diffnet(net.V1N_V1R, diffMethod = "fisherTest", alpha = 0.01, adjust = 'none', lfdrThresh = 0.2, seed = 1234)
diff_net4 <- diffnet(net.V4N_V4R, diffMethod = "fisherTest", alpha = 0.01, adjust = 'none', lfdrThresh = 0.2, seed = 1234)
diffmat_sumsN <- rowSums(diff_netN$diffAdjustMat)
diffmat_sumsR <- rowSums(diff_netR$diffAdjustMat)
diffmat_sums1 <- rowSums(diff_net1$diffAdjustMat)
diffmat_sums4 <- rowSums(diff_net4$diffAdjustMat)
diff_namesN <- names(diffmat_sumsN[diffmat_sumsN > 0])
diff_namesR <- names(diffmat_sumsR[diffmat_sumsR > 0])
diff_names1 <- names(diffmat_sums1[diffmat_sums1 > 0])
diff_names4 <- names(diffmat_sums4[diffmat_sums4 > 0])
dfvenn <- list(V1N_V4N=diff_namesN, V1R_V4R=diff_namesR, V1N_V1R = diff_names1, V4N_V4R = diff_names4)
fig_venn <- ggvenn(dfvenn)
pdf(file = paste0(outdir,'/diffnet_venn.pdf'), width = 8, height = 6)
print(fig_venn)
dev.off()




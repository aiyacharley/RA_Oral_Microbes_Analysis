library(phyloseq)
library(microeco)
library(microbiome)
library(ggplot2)

##### create phyloseq-class object; specify group (g1, g2) variables and subset specified levels (g1.level, g2.level); filter taxonomy by detection and prevalence cutoff;
create_phyloseq <- function(otu_table, tax_table=NULL, sam_data=NULL, file_tree=NULL, file_refseq=NULL, tidy_taxonomy=FALSE){
  # otu_table: OTU table
  # tax_table: taxonomy table
  # sam_data: samples metadata
  # file_tree: phylogenetic tree
  # file_refseq: reference sequences
  # tidy_taxonomy: Clean up the taxonomic table to make taxonomic assignments consistent
  
  if (is.null(tax_table)){
    # If no taxonomy file, tax_table is generated with the bacteria name, if the bacteria name is in mpa format, the tax_table is separated, otherwise only a list of Species is generated
    tax_table <- stringr::str_split(rownames(otu_table), pattern = '\\|', simplify = TRUE) %>% as.data.frame()
    # if (mpa_format) {
    #   tax_table <- stringr::str_split(rownames(otu_table), pattern = '\\|', simplify = TRUE) %>% as.data.frame()
    # }else{
    #   tax_table <- as.data.frame(row.names = rownames(otu_table),list(Species=rownames(otu_table)))
    # }
    
    if (ncol(tax_table)==1){
      colnames(tax_table) <- c('Species')
      f <- 'Species'
    }else{
      colnames(tax_table) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')[1:ncol(tax_table)]
      f <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')[ncol(tax_table)]
    }
    #print(tax_table[[f]])
    rownames(tax_table) <- tax_table[[f]]
    rownames(otu_table) <- tax_table[[f]]
    #rownames(tax_table) <- rownames(otu_table)
    #rownames(otu_table) <- rownames(otu_table)
  }else{
    tax_table <- tax_table[rownames(otu_table),]
  }
  # Clean up the taxonomic table to make taxonomic assignments consistent.
  if (tidy_taxonomy) tax_table <- microeco::tidy_taxonomy(tax_table)
  
  if (is.null(sam_data)){
    sam_data <- as.data.frame(row.names = colnames(otu_table),list(SampleID=colnames(otu_table)))
  }
  # Remove the sample whose sum is 0 and command sampleIDs in otu_table and sam_data
  TargetSamples <- intersect(colnames(otu_table)[colSums(otu_table)>0], rownames(sam_data))
  otu_table <- otu_table[,TargetSamples,drop=F]
  sam_data <- sam_data[TargetSamples,,drop=F]
  # create phyloseq
  OTU <- phyloseq::otu_table(as.matrix(otu_table), taxa_are_rows = TRUE)
  TAX <- phyloseq::tax_table(as.matrix(tax_table))
  META <- phyloseq::sample_data(sam_data)
  
  if (!is.null(file_tree)) {
    TREE <- phyloseq::phy_tree(file_tree)
  }else{
    TREE <- NULL
  }
  if (!is.null(file_refseq)) {
    REP <- phyloseq::refseq(file_refseq)
  }else{
    REP <- NULL
  }
  ps <- phyloseq::phyloseq(OTU, TAX, META, TREE, REP)
  return(ps)
}
##### filter phyloseq
filter_phyloseq <- function(ps, tax_level=NULL, rel='identity',
                            select.sample=NULL,select.taxa=NULL,
                            g1=NULL, g1.level=NULL, detection=0, prevalence=0,
                            delZeroVar=FALSE, delcorr=FALSE, collinearity=FALSE){
  # ps: phyloseq-class object
  # g1: group1 in sample metadata
  # g1.level: subgroups level included in group1
  # detection: detection cutoff
  # prevalence: prevalence cutoff
  # rel: Transformation to apply. The options include: 'CPM' (pre-sample normalization of the sum of the values to 1e+06), 'compositional' (ie relative abundance), 'Z', 'log10', 'log10p', 'hellinger', 'identity', 'clr', 'alr'
  
  # Agglomerate taxa of the same type
  if (!is.null(tax_level)){
    ps <- phyloseq::tax_glom(ps, taxrank = tax_level, NArm = TRUE, bad_empty = c('s__','g__','f__','o__','c__','p__',NA,'\t','',' ','Unassigned','Unclassified'))
    rownames(ps@otu_table) <- ps@tax_table@.Data[,tax_level]
    rownames(ps@tax_table) <- ps@tax_table@.Data[,tax_level]
    ps@phy_tree <- NULL
    ps@refseq <- NULL
  }
  # OTU table transformation to apply
  if (!is.null(rel)){
    if (rel=='CPM'){
      ps <- microbiome::transform(ps, 'compositional')
      ps@otu_table <- round(1000000*ps@otu_table,0)
    }else if (rel=='quantile'){
      otu <- preprocessCore::normalize.quantiles(ps@otu_table, keep.names = T)
      ps@otu_table <- otu_table(otu, taxa_are_rows = T)
    }else{
      # “total”, “max”, “frequency”, “normalize”, “range”, “rank”, “rrank”,
      # “standardize”, “pa”, “chi.square”, “hellinger”, “log”, “clr”, “rclr”, “alr”
      ps <- microbiome::transform(ps, rel)
    }
  }
  # select taxonomy
  if (!is.null(select.taxa)){
    ps <- phyloseq::prune_taxa(select.taxa, ps)
  }
  # select samples
  if (!is.null(select.sample)){
    ps <- phyloseq::prune_samples(select.sample, ps)
  }
  if (!is.null(g1) & !is.null(g1.level)){
    IDs <- rownames(subset(ps@sam_data, eval(parse(text = g1)) %in% g1.level))
    ps <- phyloseq::prune_samples(IDs, ps)
    ps@sam_data[[g1]] <- factor(ps@sam_data[[g1]], levels = g1.level)
  }
  # Filter the phyloseq object by detection, prevalence cutoff
  if (detection > 0 | prevalence > 0){
    ps <- microbiome::core(ps, detection = detection, prevalence = prevalence)
  }
  
  dataX <- as.data.frame(as.matrix(t(ps@otu_table)))
  #删除方差为0的变量
  if (delZeroVar){
    zerovar <- caret::nearZeroVar(dataX)
    if(length(zerovar)>0) dataX <- dataX[,-zerovar]
  }
  #删除强相关的变量，默认是0.9
  if(delcorr){
    highCorr = caret::findCorrelation(cor(dataX, method = 'pearson'), cutoff = 0.9)
    if(length(highCorr)>0) dataX <-  dataX[, -highCorr]
  }
  #解决多重共线性
  if(collinearity){
    comboInfo <- caret::findLinearCombos(dataX)
    if(length(comboInfo)>0) dataX <- dataX[, -comboInfo$remove]
  }
  Features <- colnames(dataX)
  ps <- phyloseq::prune_taxa(Features, ps)
  
  return(ps)
}
trans_abund_object <- function(dataset, taxrank='Phylum', groupmean = NULL, ntaxa = 15, sort_bacteria = FALSE){
  groupmean <- intersect(groupmean, colnames(dataset$sample_table))
  if (length(groupmean)==0) groupmean <- NULL
  if (length(groupmean)>1){
    # 使用apply函数和paste函数合并目标列
    tmplist <- list()
    for (f in groupmean){
      ls <- levels(factor(dataset$sample_table[[f]]))
      tmplist[[f]] <- ls
    }
    new_ls <- apply(expand.grid(tmplist), 1, function(x) paste(x, collapse = "_"))
    
    dataset$sample_table[['combinedxyz']] <- apply(dataset$sample_table[,groupmean], 1, function(x) paste(x, collapse = "_")) %>% factor(levels = new_ls)
    groupmean <- 'combinedxyz'
  }
  
  t1 <- trans_abund$new(dataset = dataset,
                        groupmean = groupmean,
                        taxrank = taxrank,
                        delete_taxonomy_prefix = TRUE,
                        ntaxa = ntaxa)
  # sort the samples by top bacteria
  if (sort_bacteria){
    tmp <- t1$data_abund[,c('Taxonomy','Sample','Abundance')]
    tmp2 <- tmp[tmp$Taxonomy==t1$data_taxanames[1],]
    tmp2_sorted <- tmp2[order(tmp2$Abundance, decreasing = TRUE), ]
    t1$data_abund$Sample <- factor(t1$data_abund$Sample, levels = tmp2_sorted$Sample)
  }
  return(t1)
}
plot_comp <- function(ps, taxrank='Phylum', strata = NULL, groupmean = NULL, ntaxa = 15, sort_bacteria = TRUE, 
                      clustering = FALSE, clustering_plot = FALSE, use_alluvium = FALSE){
  # ps: phyloseq object store and manage all the basic files.
  # taxrank: select specified taxonomy level
  # strata: Specifies a faceted variable
  # groupmean: calculate mean abundance for each group
  # ntaxa: how many taxa are selected to show
  # sort_bacteria: Whether sort the top1 bacteria to arrange the sample
  # clustering: whether order samples by the clustering
  # clustering_plot: whether add clustering plot
  # use_alluvium: whether add alluvium plot
  require(file2meco)
  require(microeco)
  dataset <- file2meco::phyloseq2meco(ps)
  
  if (length(strata)==0) strata <- NULL
  t1 <- trans_abund_object(dataset, taxrank=taxrank, groupmean = groupmean, ntaxa = ntaxa, sort_bacteria = sort_bacteria)
  
  if (!is.null(groupmean)){
    strata <- NULL
  }
  if (clustering_plot){
    Pic1 <- t1$plot_bar(others_color = 'grey80',
                        facet = strata,
                        clustering = clustering,
                        clustering_plot = clustering_plot, # clustering_plot=T时不能添加theme设置
                        xtext_keep = T,
                        xtext_angle = 90,
                        use_alluvium = use_alluvium)
  }else{
    Pic1 <- t1$plot_bar(others_color = 'grey80',
                        facet = strata,
                        clustering = clustering,
                        clustering_plot = clustering_plot, # clustering_plot=T时不能添加theme设置
                        xtext_keep = T,
                        xtext_angle = 90,
                        use_alluvium = use_alluvium)+
      #guides(fill = guide_legend(title = taxrank,nrow = 21))+
      theme(text = element_text(size = 12), axis.text.x= element_text(angle = 90, hjust = 1, vjust = .5))
  }
  Pic1
}
##### calculate alpha diversity
cal_alpha <- function(ps, group=NULL, strata=NULL, method_test='wilcox.test', signif_label=FALSE, adj.vars=NULL){
  # ps: phyloseq object
  # group: A string indicating the variable for group identifiers
  # strata: A string indicating the variable for strata identifiers
  # method_test: the name of the statistical test (e.g. t.test, wilcox.test etc.)
  # signif_label: whether show the p-value as labels: c('***'=0.001, '**'=0.01, '*'=0.05)
  # adj.vars: Specify adjust variables of sample metadata
  require(ggsci)
  res_alpha <- microbiome::alpha(ps, index = 'all')
  res_alpha <- dplyr::select(res_alpha, where(~ !all(is.na(.)))) # 去除全是NA的列
  alpha_index <- colnames(res_alpha)
  
  if (!is.null(adj.vars)){
    res_alpha <- adj_covariate_lm(data=res_alpha, covariates = as.matrix(ps@sam_data[,adj.vars]), vars=NULL)
  }
  res <- list(table_index=res_alpha)
  # plot
  if(!is.null(group)){
    res_alpha$group <- ps@sam_data[[group]]
    if(!is.null(strata)) res_alpha$strata <- ps@sam_data[[strata]]
    
    comps <- combn(unique(as.character(res_alpha$group)),m=2,simplify=F)
    res[['table_index']] <- res_alpha
  }
  ggparam1 <- list(
    theme_classic(),
    theme(text = element_text(size = 15)),
    geom_violin(alpha=.6, width=0.5),
    geom_boxplot(width = 0.05, alpha=.9),
    #geom_jitter(height = 0,width = 0.2, shape=1),
    scale_fill_d3(),
    scale_color_d3(),
    if(!is.null(group)) geom_signif(comparisons = comps,
                                    test = method_test,
                                    map_signif_level = signif_label,
                                    textsize = 4,
                                    step_increase = 0.1),
    if(!is.null(strata)) facet_wrap(strata~., nrow = 1)
  )
  for (p in alpha_index){
    if(!is.null(group)){
      fig <- ggplot(res_alpha,aes_string(x=rlang::sym('group'),y=rlang::sym(p),fill=rlang::sym('group'))) + ggparam1
    }else{
      res_alpha$group <- 'All'
      fig <- ggplot(res_alpha,aes_string(x=rlang::sym('group'),y=rlang::sym(p),fill=rlang::sym('group'))) + ggparam1
    }
    res[[paste0("figure_",p)]] <- fig
  }
  return(res)
}
##### Calculate dissimilarity matrix, anosim and PERMANOVA analysis
cal_dist <- function(ps, method='bray', group=NULL, adj.vars=NULL){
  # ps: phyloseq object
  # method: Dissimilarity index
  # group: A string indicating the variable for group identifiers for anosim and PERMANOVA analysis
  # adj.vars: adjust variables of sample metadata. Note: This feature is yet to be verified
  
  method <- match.arg(method, choices = c('manhattan', 'euclidean', 'canberra', 'clark', 'bray',
                                          'kulczynski', 'jaccard', 'gower', 'altGower', 'morisita',
                                          'horn', 'mountford', 'raup', 'binomial', 'chao', 'cao',
                                          'mahalanobis', 'chisq', 'chord', 'hellinger', 'aitchison','robust.aitchison'))
  dist_veg <- vegan::vegdist(t(ps@otu_table), method = method)
  # Dissimilarity adjusted by adj.vars
  if (!is.null(adj.vars)) dist_veg <- adj_distance(dist_veg, ps@sam_data[,adj.vars])
  anosim <- NULL
  table_adonis <- NULL
  if (!is.null(group)) {
    set.seed(1234)
    # Analysis of similarities (anosim)
    anosim <- vegan::anosim(dist_veg, grouping = ps@sam_data[[group]], permutations = 999)
    # PERMANOVA analysis
    table_adonis <- cal_PERMANOVA(dist_veg, group=ps@sam_data[[group]])
  }
  return(list(dist=dist_veg,anosim=anosim,PERMANOVA=table_adonis))
}
##### Permutational Multivariate Analysis of Variance Using Distance Matrices (PERMANOVA): Analysis of variance using distance matrices
cal_PERMANOVA <- function(dist, group=NULL){
  # dist: dissimilarity matrix
  # group: A string indicating the variable for group identifiers
  
  set.seed(1234)
  # Total
  adonis2res <- vegan::adonis2(dist ~ group, permutations = 999, by='margin')
  R2 <- round(adonis2res[1,3],6)
  pvalue <- adonis2res[1,5]
  # pairwise
  dune.pairwise.adonis <- pairwiseAdonis::pairwise.adonis(x=dist,
                                                          factors=group,
                                                          p.adjust.m = 'BH',
                                                          reduce = NULL,
                                                          perm = 999)
  dune.pairwise.adonis$R2 <- round(dune.pairwise.adonis$R2, 6)
  table_adonis <- rbind(dune.pairwise.adonis[,c('pairs','R2','p.value','p.adjusted')], c('Total',R2,pvalue,''))
  return(table_adonis)
}
##### Perform an ordination on phyloseq data
cal_ordination <- function(ps, method_dist='bray', method = 'PCoA', type='samples',
                           group=NULL, shape = NULL, label = NULL, ellipse_level=0.95, add_cent=TRUE){
  # ps: phyloseq object
  # method: several commonly-used ordination methods. Currently supported method options are: c('DCA', 'CCA', 'RDA', 'NMDS', 'MDS', 'PCoA')
  # method_dist: Dissimilarity index
  # type:  The plot type. Default is 'samples'. The currently supported options are c('samples', 'taxa', 'biplot', 'split', 'scree')
  # group: The name of the variable to map to colors in the plot
  # shape: The name of the variable to map to different shapes on the plot
  # label: The name of the variable to map to text labels on the plot
  # ellipse_level: The level at which to draw an ellipse
  # add_cent: plot central point
  if (!is.null(group)){
    # 使用apply函数和paste函数合并目标列
    if (length(group)>1){
      ps@sam_data[['Groupxyz']] <- apply(ps@sam_data[,group], 1, function(x) paste(x, collapse = "_"))
    }else{
      ps@sam_data[['Groupxyz']] <- ps@sam_data[,group] %>% unlist()
    }
    
  }else{
    ps@sam_data[['Groupxyz']] <- 'All'
  }
  group <- 'Groupxyz'
  
  res_dist <- cal_dist(ps, method = method_dist, group = group, adj.vars = NULL)
  
  method <- match.arg(method, choices = c('DCA', 'CCA', 'RDA', 'DPCoA', 'NMDS', 'MDS', 'PCoA')) # CAP, MDS==PCoA
  type <- match.arg(type, choices = c('samples', 'taxa', 'biplot', 'split', 'scree'))
  ord <- phyloseq::ordinate(ps, method = method, distance = res_dist$dist)
  #mtitle <- paste("ANOSIM: R = ",round(as.numeric(res_dist$anosim$statistic),3),", p = ", res_dist$anosim$signif, "\nPERMANOVA: R2 = ",round(as.numeric(res_dist$PERMANOVA$R2[nrow(res_dist$PERMANOVA)]),3),", p = ", res_dist$PERMANOVA$p.value[nrow(res_dist$PERMANOVA)])
  mtitle <- paste("PERMANOVA: R2 = ",round(as.numeric(res_dist$PERMANOVA$R2[nrow(res_dist$PERMANOVA)]),3),", p = ", res_dist$PERMANOVA$p.value[nrow(res_dist$PERMANOVA)])
  
  fig_ord1 <- phyloseq::plot_ordination(ps, ord, type = type, color = group, shape = shape, label = label, axes = c(1,2), title = mtitle) +
    geom_point(size=1)+
    theme_classic()+
    labs(color="Group")+
    #guides(color = guide_legend(title = 'Group'))+
    scale_color_d3()+scale_fill_d3()+
    list(
      if (ellipse_level!=0) stat_ellipse(linetype=5,level=ellipse_level)
    )
  
  plotdata <- fig_ord1$data
  if (type!='taxa' & type!='scree'){
    # fig_ord1 <- fig_ord1 + stat_ellipse(data = plotdata, aes_string(color=rlang::sym(group)),linetype=5,level=ellipse_level)
    # 求均值
    cent <- aggregate(cbind(plotdata[,1],plotdata[,2]) ~ plotdata[[group]], FUN = mean)
    names(cent) <- c(group,"x","y")
    # 合并到样本坐标数据中
    segs <- merge(plotdata, setNames(cent, c(group,'o1','o2')), by = group, sort = FALSE)
    fig_ord2 = fig_ord1 +
      geom_segment(data = segs, mapping = aes_string(xend = rlang::sym('o1'), yend = rlang::sym('o2'),color = rlang::sym(group)),show.legend=F) +
      geom_point(mapping = aes_string(x = rlang::sym('x'), y = rlang::sym('y'),color = rlang::sym(group)),data = cent, size = 5,pch = 24,fill = "white",show.legend=F)
    fig_ord2$data <- segs
    if (add_cent){
      pic <- fig_ord2
    }else{
      pic <- fig_ord1
    }
  }else{
    pic <- fig_ord1
  }
  require(ggvegan)
  fig_anosim <- autoplot(res_dist$anosim, notch = TRUE)+
    theme_classic()+
    theme(panel.grid = element_blank(),
          legend.position="top",
          text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  return(list(rds_dist=res_dist$dist,rds_anosim=res_dist$anosim,rds_ordination=ord,table_PERMANOVA=res_dist$PERMANOVA,figure_ord=pic, figure_anosim=fig_anosim))
}
##### Visualize OTU Core heatmap
plot_detection_prevalence <- function(ps, min_prevalence = 0.2, min_detections = 0, step_detection=100){
  # ps: phyloseq object
  # min_prevalence: minimum prevalence cutoff
  # step_detection: Set the number of steps for detection
  
  maxNum <- max(ps@otu_table)
  if (maxNum>1){
    digit <- 0
  }else{
    digit <- 3
  }
  detections <- round(seq(min_detections, maxNum, maxNum/step_detection),digit)
  prevalences <- seq(0.05, 1, 0.05)
  
  Pic <- microbiome::plot_core(ps, plot.type = 'heatmap', # heatmap, lineplot
                               prevalences = prevalences,
                               detections = detections,
                               colours = rev(brewer.pal(5, 'Spectral')),
                               min.prevalence = min_prevalence, horizontal = TRUE) +
    theme(axis.text.x= element_text(angle = 90, hjust = 1, vjust = .5))
  return(Pic)
}
##### sankey
plot_sankey <- function(taxa, split = "|", nodes_name = NULL, 
                        sankeyNetwork_iterations = 30, sankeyNetwork_fontSize = 10, 
                        saveWidget=NULL,height = NULL, width = NULL,
                        sinksRight = FALSE){
  require(networkD3)
  split_list <- strsplit(taxa, split = split, fixed = TRUE)
  # 创建一个空列表来存储每个元素生成的配对数据框
  pair_list <- list()
  for(i in seq_along(split_list)) {
    current_vec <- split_list[[i]]
    # 如果拆分后的元素数量大于1，才能生成配对
    if(length(current_vec) > 1) {
      # 生成 source 和 target 列
      source_col <- current_vec[1:(length(current_vec)-1)]
      target_col <- current_vec[2:length(current_vec)]
      # 创建临时数据框并存入列表
      pair_list[[i]] <- data.frame(source = source_col, target = target_col, stringsAsFactors = FALSE)
    }
  }
  # 将所有数据框行合并
  final_df <- do.call(rbind, pair_list)
  
  links <- final_df %>%
    dplyr::group_by(source, target) %>%
    dplyr::summarise(value = n()) %>% as.data.frame()
  
  # 创建节点，并调整顺序
  s0 <- unique(c(final_df$target,final_df$source))
  if (!is.null(nodes_name)){
    s1 <- intersect(nodes_name,s0)
    s2 <- setdiff(s0,nodes_name)
    s0 <- c(s1,s2)
  }
  nodes <- data.frame(
    list(name = s0)
  )
  # 使用 networkD3，连接必须使用 id 进行提供
  links$IDsource <- match(links$source, nodes$name)-1
  links$IDtarget <- match(links$target, nodes$name)-1
  # 绘制桑基图
  fig <- networkD3::sankeyNetwork(Links = links, Nodes = nodes,iterations = sankeyNetwork_iterations,
                                  Source = "IDsource", Target = "IDtarget",
                                  Value = "value", NodeID = "name",
                                  nodeWidth = 15, nodePadding = 10,
                                  height = height, width = width,
                                  fontSize = sankeyNetwork_fontSize, sinksRight = sinksRight)
  if (!is.null(saveWidget)){
    htmlwidgets::saveWidget(widget = fig, file = saveWidget)
  }
  return(fig)
}
# boxplot
plot_box <- function(dataset, y='', group = NULL, strata = NULL, 
                     testMethod = c('wilcox.test','t.test')[1],
                     signif_label = TRUE, ysqrt = FALSE,
                     paired_groups = NULL, # 新增：指定配对的组，如 c("V1","V4")
                     paired_id = NULL,     # 新增：患者ID列名
                     show_paired_pval = FALSE # 新增：是否显示配对检验p值
){
  require(ggsignif)
  require(dplyr)
  
  df <- as.data.frame(list(SID = rownames(dataset$sample_table)))
  df$y <- dataset$otu_table[y,] %>% unlist() %>% as.numeric()
  
  group <- intersect(group, colnames(dataset$sample_table))
  strata <- intersect(strata, colnames(dataset$sample_table))
  
  if (length(strata) > 0){
    df <- cbind(df, dataset$sample_table[,strata, drop=F])
  }
  if (length(group) > 0){
    V_tmp <- unique(as.character(dataset$sample_table[[group]]))
    if(length(V_tmp) >= 2) comps <- combn(V_tmp, m=2, simplify=F)
    
    df[['group']] <- factor((dataset$sample_table[,group]))
  } else {
    df$group <- 'all'
  }
  
  # 添加配对ID（如果有）
  if (!is.null(paired_id) && paired_id %in% colnames(dataset$sample_table)) {
    df[['paired_id']] <- dataset$sample_table[[paired_id]]
  }
  
  # 计算每组的中位数和平均值
  summary_stats <- df %>%
    group_by(group) %>%
    summarise(
      median_y = median(y, na.rm = TRUE),
      mean_y = mean(y, na.rm = TRUE),
      .groups = 'drop'
    )
  
  Pic <- ggplot(df, aes(x=group, y=y)) +
    geom_violin(aes(fill=group), width=0.5, alpha=0.6) +
    geom_boxplot(aes(color=group), width=0.05, outliers = TRUE) +
    # 添加中位数连接线（实线）
    geom_line(data = summary_stats, 
              aes(x = group, y = median_y, group = 1),
              color = "black", linetype = "solid", size = 1) +
    # 添加平均值连接线（虚线）
    geom_line(data = summary_stats, 
              aes(x = group, y = mean_y, group = 1),
              color = "red", linetype = "dashed", size = 1) +
    # 添加中位数点
    geom_point(data = summary_stats, 
               aes(x = group, y = median_y, group = 1),
               color = "black", size = 3, shape = 16) +
    # 添加平均值点
    geom_point(data = summary_stats, 
               aes(x = group, y = mean_y, group = 1),
               color = "red", size = 3, shape = 17) +
    theme_classic() +
    labs(y = "Abundance") +
    ggtitle(y) +
    theme(plot.title = element_text(size = 11)) +
    ggsci::scale_color_d3(palette = 'category20') +
    ggsci::scale_fill_d3(palette = 'category20') +
    list(
      # 如果指定了配对的组，添加配对样本连线
      if (!is.null(paired_groups) && length(paired_groups) == 2 && 
          !is.null(paired_id) && 'paired_id' %in% names(df)) {
        # 筛选出属于配对组的数据
        df_paired <- df[df$group %in% paired_groups, ]
        # 确保每个ID在两个组中都有数据
        paired_ids <- intersect(
          df_paired$paired_id[df_paired$group == paired_groups[1]],
          df_paired$paired_id[df_paired$group == paired_groups[2]]
        )
        df_paired <- df_paired[df_paired$paired_id %in% paired_ids, ]
        
        # 添加配对连线
        list(
          geom_line(data = df_paired, 
                    aes(x = group, y = y, group = paired_id),
                    color = "grey50", alpha = 0.5, size = 0.5,
                    position = position_dodge(0.2)),
          # 如果要求显示配对检验p值
          if (show_paired_pval && length(paired_ids) > 0) {
            # 准备配对数据
            df_wide <- df_paired %>%
              select(paired_id, group, y) %>%
              tidyr::pivot_wider(names_from = group, values_from = y)
            
            # 执行配对检验
            pval <- tryCatch({
              test_result <- wilcox.test(
                df_wide[[paired_groups[1]]], 
                df_wide[[paired_groups[2]]], 
                paired = TRUE
              )
              test_result$p.value
            }, error = function(e) NA)
            
            # 添加p值标注
            annotate("text", 
                     x = mean(as.numeric(factor(paired_groups, levels = levels(df$group)))),
                     y = max(df_paired$y, na.rm = TRUE) * 1.1,
                     label = ifelse(!is.na(pval) && pval < 0.001, "p < 0.001",
                                    ifelse(!is.na(pval), paste0("p = ", round(pval, 3)), "NA")),
                     size = 4, color = "darkred")
          }
        )
      },
      if (length(strata) > 0) facet_grid(as.formula(paste0("~", paste0(strata, collapse = "+"))), scale = "free"),
      if (length(group) > 0) geom_signif(comparisons = comps,
                                         test = testMethod, 
                                         map_signif_level = signif_label,
                                         textsize = 4,
                                         step_increase = 0.1),
      if (ysqrt) scale_y_sqrt()
    )
  
  Pic
}
# 网络节点角色
plot_zipi <- function(ps,net,abs_r=0,cut_z=2.5,cut_p=0.62, title="", show_text=TRUE){
  mt <- file2meco::phyloseq2meco(ps)
  adjacency_unweight <- net$assoMat1
  #随便使用 microeco 包构建一个网络
  t1 <- trans_network$new(dataset = mt, cor_method = 'spearman', filter_thres = 0.001)
  t1$cal_network(COR_p_thres = 0.01, COR_cut = abs_r)
  
  #然后用现成的网络将上述网络替换掉，这样就把我们的网络数据也添加至 microeco 对象里面了
  adjacency_unweight[abs(adjacency_unweight)<=abs_r] <- 0
  adjacency_unweight[abs(adjacency_unweight)>abs_r] <- 1
  g <- igraph::graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)  #邻接矩阵 -> igraph 的邻接列表，获得非含权的无向网络
  t1$res_network <- g
  
  #计算网络属性，并划分网络模块
  t1$cal_network_attr()
  t1$cal_module()
  
  #提取各节点的拓扑属性，包括 zi 和 pi 值等
  t1$get_node_table(node_roles = TRUE)
  result <- t1$res_node_table
  result$Abundance <- rowSums(mt$otu_table)
  
  # 根据 zi 和 pi 的值划分节点的角色
  # cut_z <- 2.5    #2.5
  # cut_p <- 0.62  # 0.62
  result$taxa_roles <- dplyr::case_when(
    result$z >= cut_z & result$p < cut_p ~ 'Module hubs',
    result$z >= cut_z & result$p >= cut_p ~ 'Network hubs',
    result$z < cut_z & result$p >= cut_p ~ 'Connectors',
    TRUE ~ 'Peripheral nodes'
  )
  z1 <- subset(result, taxa_roles %in% c('Module hubs','Network hubs','Connectors'))
  
  #作图展示
  fig <- ggplot(result,aes(x=p,y=z, color=taxa_roles)) +
    geom_point(aes(size=log(Abundance+1))) +
    scale_color_manual(values = c('Network hubs'='red','Module hubs'='#7570B3','Connectors'='#D95F02','Peripheral nodes'='#1B9E77')) +
    list(
      if (show_text){
        # 标记差异基因
        ggrepel::geom_text_repel(
          data = z1,
          aes(x = p, y = z, label=name),        
          min.segment.length = 0.1,
          max.overlaps = 10000,                    # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之。
          size=3,                                  # 字体大小
          box.padding=unit(0.5,'lines'),           # 标记的边距
          point.padding=unit(0.1, 'lines'), 
          segment.color='black',                   # 标记线条的颜色
          show.legend=F)
      }
      
    ) +
    # 横线和竖线分别为 zi=2.5 和 pi=0.62 的阈值线
    geom_hline(yintercept = cut_z, linetype='dashed')+
    geom_vline(xintercept = cut_p, linetype='dashed')+
    xlim(0,0.85)+
    labs(x='Among-module connectivity (Pi)', y='Within-module connectivity (Zi)', color='Taxa roles') +
    ggtitle(title)+
    theme_bw() +
    theme(legend.position = 'right',axis.title = element_text(size = 16))
  return(list(table=result, figure=fig, network_attr=t1$res_network_attr))
  # library(ggplot2)
  # t1$plot_taxa_roles(use_type = 1,add_label = T,add_label_text = "name")
  # t1$plot_taxa_roles(use_type = 2)
}
expand_colors <- function(color_values = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#326530",
                                           "#E6AB02","#A6761D","#ff523f","#2baeb5","#6E4821","#8e3af4",
                                           "#FF0000","#0000FF","#00FF00","#FFFF00","#FF00FF","#00FFFF"),
                          n = 99){
  # color_values = RColorBrewer::brewer.pal(8, "Dark2")
  options(warn = -1)
  
  
  if(n <= length(color_values)){
    total_colors <- color_values[1:n]
  }else{
    # message("Input colors are not enough to use. Add more colors automatically via color interpolation ...")
    ceiling_cycle_times <- ceiling(n/length(color_values))
    total_cycle_times <- lapply(seq_along(color_values), function(x){
      if((ceiling_cycle_times - 1) * length(color_values) + x <= n){
        ceiling_cycle_times
      }else{
        ceiling_cycle_times - 1
      }
    })
    total_cycle_times <- unlist(total_cycle_times)
    total_color_list <- lapply(seq_along(color_values), function(x){
      colorRampPalette(c(color_values[x], "white"))(total_cycle_times[x] + 1)
    })
    total_colors <- lapply(seq_len(ceiling_cycle_times), function(x){
      unlist(lapply(total_color_list, function(y){
        if((x + 1) <= length(y)) y[x]
      }))
    })
    total_colors <- unlist(total_colors)
  }
  total_colors
}
# MicrobiomeProfiler
enrichMicrobiome <- function(gene, type = c('KO','Module','COG','HMDB','MBKEGG','SMPDB','MDA')[1], 
                             pAdjustMethod = c("BH","fdr","BY","holm","hochberg","hommel","bonferroni","none")[1],
                             universe = NULL,
                             showCategory = 20,
                             orderBy = c('qvalue','pvalue')[1],
                             COG_dtype = c('category','pathway')[1],
                             title=""){
  require(MicrobiomeProfiler)
  title <- paste0("MicrobiomeProfiler enrich",type," ",title)
  para <- list(pAdjustMethod = pAdjustMethod, pvalueCutoff = 1, qvalueCutoff = 1, universe = universe)
  if (type == 'KO'){
    enrich_fun <- enrichKO
    para$gene <- gene
  }else if (type == 'Module'){
    enrich_fun <- enrichModule
    para$gene <- gene
  }else if (type == 'COG'){
    enrich_fun <- enrichCOG
    para$gene <- gene
    para$dtype <- COG_dtype
  }else if (type == 'HMDB'){
    enrich_fun <- enrichHMDB
    para$metabo_list <- gene
  }else if (type == 'MBKEGG'){
    enrich_fun <- enrichMBKEGG
    para$metabo_list <- gene
  }else if (type == 'SMPDB'){
    enrich_fun <- enrichSMPDB
    para$metabo_list <- gene
  }else if (type == 'MDA'){
    enrich_fun <- enrichMDA
    para$microbe_list <- gene
  }
  kk <- do.call(enrich_fun, para)
  fig.bar.ego <- barplot(kk, showCategory=showCategory,color=orderBy, orderBy=orderBy,label_format=150,decreasing=FALSE,font.size=10) + 
    ggtitle(title)+
    theme(legend.position = 'none')
  fig.ego <- enrichplot::dotplot(kk, showCategory=showCategory,color=orderBy, orderBy=orderBy,label_format=150,decreasing=FALSE,font.size=10) + 
    ggtitle("")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  fig <- fig.bar.ego + fig.ego
  return(list(obj=kk, figure=fig))
}

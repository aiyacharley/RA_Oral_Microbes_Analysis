# C:\Program Files\R\R-4.5.2
setwd("D:/我的工作.省中医/Projects/刘院士团队舌苔微生态/results_20260615")
library(NetCoMi)
library(readxl)
library(permute)
library(LaplacesDemon)
library(dplyr)
source('00.func.R')
outdir <- "output/Figure4/"

ps <- readRDS(file = "00_data/ps_tongue_sylph.rds")
# 转换成TSS，总和为1
ps <- microbiome::transform(ps, "total")
#colSums(ps@otu_table)

ps.HC <- phyloseq::subset_samples(ps, Time %in% c('HC'))
ps.V1 <- phyloseq::subset_samples(ps, Time %in% c('V1'))
ps.V4 <- phyloseq::subset_samples(ps, Time %in% c('V4'))

# 每组至少10%样本中丰度大于0.01%的OTU
ps.HC <- microbiome::core(ps.HC, detection = 1e-4, prevalence = 0.2)
ps.V1 <- microbiome::core(ps.V1, detection = 1e-4, prevalence = 0.2)
ps.V4 <- microbiome::core(ps.V4, detection = 1e-4, prevalence = 0.2)
taxa <- intersect(row.names(ps.HC@tax_table), intersect(row.names(ps.V1@tax_table), row.names(ps.V4@tax_table)))

ps.filter <- phyloseq::prune_taxa(taxa, ps)
ps.HC <- phyloseq::subset_samples(ps.filter, Time %in% c('HC'))
ps.V1 <- phyloseq::subset_samples(ps.filter, Time %in% c('V1'))
ps.V4 <- phyloseq::subset_samples(ps.filter, Time %in% c('V4'))
ps.V1p <- phyloseq::subset_samples(ps.filter, Group_tongue %in% c('V1N','V1R'))

ntaxa(ps.filter)
nsamples(ps.HC)
nsamples(ps.V1)
nsamples(ps.V4)
mean(otu_table(ps.filter) == 0)
summary(sample_sums(ps.filter))


#-------------------------------- netConstruct ----------------
### SPRING
# run_netConstruct1 <- function(ps1, ps2 = NULL, matchDesign = NULL){
#   netConstruct(data = ps1, data2 = ps2, matchDesign = matchDesign,
#                measure = "spring", measurePar = list(Rmethod = "approx",nlambda = 30),  
#                normMethod = "none", zeroMethod = "none", sparsMethod = "none",
#                weighted = TRUE, seed = 42)
# }
# net.HC <- run_netConstruct1(ps.HC);saveRDS(net.HC, file = paste0(outdir,"/netConstruct_spring/net.HC.rds"))
# net.V1 <- run_netConstruct1(ps.V1);saveRDS(net.V1, file = paste0(outdir,"/netConstruct_spring/net.V1.rds"))
# net.V1p <- run_netConstruct1(ps.V1p);saveRDS(net.V1p, file = paste0(outdir,"/netConstruct_spring/net.V1p.rds"))
# net.V4 <- run_netConstruct1(ps.V4);saveRDS(net.V4, file = paste0(outdir,"/netConstruct_spring/net.V4.rds"))
# net.HCV1 <- run_netConstruct1(ps.HC, ps.V1);saveRDS(net.HCV1, file = paste0(outdir,"/netConstruct_spring/net.HCV1.rds"))
# net.V1V4 <- run_netConstruct1(ps.V1p, ps.V4, matchDesign = c(1,1));saveRDS(net.V1V4, file = paste0(outdir,"/netConstruct_spring/net.V1V4.rds"))
# net.HCV4 <- run_netConstruct1(ps.HC, ps.V4);saveRDS(net.HCV4, file = paste0(outdir,"/netConstruct_spring/net.HCV4.rds"))
# save(ps.filter, ps.HC, ps.V1, ps.V4, ps.V1p, 
#      net.HC, net.V1, net.V1p, net.V4, net.HCV1, net.V1V4, net.HCV4, 
#      file = paste0(outdir,"/netConstruct_spring/netConstruct.Rdata"))

### spieceasi
run_netConstruct2 <- function(ps1, ps2 = NULL, matchDesign = NULL){
  netConstruct(data = ps1, data2 = ps2, matchDesign = matchDesign,
               measure = "spieceasi", measurePar = list(method = "mb",nlambda = 30,lambda.min.ratio = 1e-2),
               normMethod = "none", zeroMethod = "none", sparsMethod = "none",
               weighted = TRUE, seed = 42)
}
net.HC <- run_netConstruct2(ps.HC);saveRDS(net.HC, file = paste0(outdir,"/netConstruct_spieceasi/net.HC.rds"))
net.V1 <- run_netConstruct2(ps.V1);saveRDS(net.V1, file = paste0(outdir,"/netConstruct_spieceasi/net.V1.rds"))
net.V1p <- run_netConstruct2(ps.V1p);saveRDS(net.V1p, file = paste0(outdir,"/netConstruct_spieceasi/net.V1p.rds"))
net.V4 <- run_netConstruct2(ps.V4);saveRDS(net.V4, file = paste0(outdir,"/netConstruct_spieceasi/net.V4.rds"))
net.HCV1 <- run_netConstruct2(ps.HC, ps.V1);saveRDS(net.HCV1, file = paste0(outdir,"/netConstruct_spieceasi/net.HCV1.rds"))
net.V1V4 <- run_netConstruct2(ps.V1p, ps.V4, matchDesign = c(1,1));saveRDS(net.V1V4, file = paste0(outdir,"/netConstruct_spieceasi/net.V1V4.rds"))
net.HCV4 <- run_netConstruct2(ps.HC, ps.V4);saveRDS(net.HCV4, file = paste0(outdir,"/netConstruct_spieceasi/net.HCV4.rds"))
save(ps.filter, ps.HC, ps.V1, ps.V4, ps.V1p, 
     net.HC, net.V1, net.V1p, net.V4, net.HCV1, net.V1V4, net.HCV4, 
     file = paste0(outdir,"/netConstruct_spieceasi/netConstruct.Rdata"))
### spearman
# run_netConstruct3 <- function(ps1, ps2 = NULL, matchDesign = NULL){
#   netConstruct(data = ps1, data2 = ps2, matchDesign = matchDesign,
#                measure = "spearman",
#                sparsMethod = "t-test", alpha = 0.05, adjust = "adaptBH",
#                normMethod = "clr", zeroMethod = "multRepl",
#                weighted = TRUE, seed = 42)
# }
# net.HC <- run_netConstruct3(ps.HC);saveRDS(net.HC, file = paste0(outdir,"/netConstruct_spearman/net.HC.rds"))
# net.V1 <- run_netConstruct3(ps.V1);saveRDS(net.V1, file = paste0(outdir,"/netConstruct_spearman/net.V1.rds"))
# net.V1p <- run_netConstruct3(ps.V1p);saveRDS(net.V1p, file = paste0(outdir,"/netConstruct_spearman/net.V1p.rds"))
# net.V4 <- run_netConstruct3(ps.V4);saveRDS(net.V4, file = paste0(outdir,"/netConstruct_spearman/net.V4.rds"))
# net.HCV1 <- run_netConstruct3(ps.HC, ps.V1);saveRDS(net.HCV1, file = paste0(outdir,"/netConstruct_spearman/net.HCV1.rds"))
# net.V1V4 <- run_netConstruct3(ps.V1p, ps.V4, matchDesign = c(1,1));saveRDS(net.V1V4, file = paste0(outdir,"/netConstruct_spearman/net.V1V4.rds"))
# net.HCV4 <- run_netConstruct3(ps.HC, ps.V4);saveRDS(net.HCV4, file = paste0(outdir,"/netConstruct_spearman/net.HCV4.rds"))
# save(ps.filter, ps.HC, ps.V1, ps.V4, ps.V1p, 
#      net.HC, net.V1, net.V1p, net.V4, net.HCV1, net.V1V4, net.HCV4, 
#      file = paste0(outdir,"/netConstruct_spearman/netConstruct.Rdata"))
### pearson
run_netConstruct4 <- function(ps1, ps2 = NULL, matchDesign = NULL){
  netConstruct(data = ps1, data2 = ps2, matchDesign = matchDesign,
               measure = "pearson",
               sparsMethod = "t-test", alpha = 0.05, adjust = "adaptBH",
               normMethod = "clr", zeroMethod = "multRepl",
               weighted = TRUE, seed = 42)
}
net.HC <- run_netConstruct4(ps.HC);saveRDS(net.HC, file = paste0(outdir,"/netConstruct_pearson/net.HC.rds"))
net.V1 <- run_netConstruct4(ps.V1);saveRDS(net.V1, file = paste0(outdir,"/netConstruct_pearson/net.V1.rds"))
net.V1p <- run_netConstruct4(ps.V1p);saveRDS(net.V1p, file = paste0(outdir,"/netConstruct_pearson/net.V1p.rds"))
net.V4 <- run_netConstruct4(ps.V4);saveRDS(net.V4, file = paste0(outdir,"/netConstruct_pearson/net.V4.rds"))
net.HCV1 <- run_netConstruct4(ps.HC, ps.V1);saveRDS(net.HCV1, file = paste0(outdir,"/netConstruct_pearson/net.HCV1.rds"))
net.V1V4 <- run_netConstruct4(ps.V1p, ps.V4, matchDesign = c(1,1));saveRDS(net.V1V4, file = paste0(outdir,"/netConstruct_pearson/net.V1V4.rds"))
net.HCV4 <- run_netConstruct4(ps.HC, ps.V4);saveRDS(net.HCV4, file = paste0(outdir,"/netConstruct_pearson/net.HCV4.rds"))
save(ps.filter, ps.HC, ps.V1, ps.V4, ps.V1p, 
     net.HC, net.V1, net.V1p, net.V4, net.HCV1, net.V1V4, net.HCV4, 
     file = paste0(outdir,"/netConstruct_pearson/netConstruct.Rdata"))




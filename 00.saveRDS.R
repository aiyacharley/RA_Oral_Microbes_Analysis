setwd("D:/我的工作.省中医/Projects/刘院士团队舌苔微生态/results_20260615/")
library(file2meco)
library(phyloseq)
library(microeco)
library(readxl)
library(dplyr)
source('00.func.R')

#----------------------------- tongue --------------------------------------
meta0 <- read_xlsx("./00_data/00_metadata.xlsx",sheet = "meta_female") |> as.data.frame()
rownames(meta0) <- meta0$SampleID_Coating
otu <- read.delim("./00_data/tongue_abundance.tsv",sep = '\t', header = T,row.names = 1,check.names = F)

ps <- create_phyloseq(otu_table = otu, sam_data = meta0)
saveRDS(ps, file = "00_data/ps_tongue.rds")

meta0 <- read.csv("./00_data/tongue_meta.txt", sep = '\t', row.names = 1)
otu <- read.delim("./00_data/RA_tongue_relative_abundance_clean.tsv",sep = '\t', header = T,row.names = 1,check.names = F)
ps1 <- create_phyloseq(otu_table = otu, sam_data = meta0)
ps1@otu_table
saveRDS(ps1, file = "00_data/ps_tongue_sylph.rds")

#----------------------------- saliva --------------------------------------
meta0 <- read_xlsx("./00_data/00_metadata.xlsx",sheet = "meta_female") |> as.data.frame()
rownames(meta0) <- meta0$SampleID_Saliva
otu <- read.delim("./00_data/saliva_abundance.tsv",sep = '\t', header = T,row.names = 1,check.names = F)

ps <- create_phyloseq(otu_table = otu, sam_data = meta0)
saveRDS(ps, file = "00_data/ps_saliva_MAGs.rds")

meta0 <- read.csv("./00_data/saliva_meta.txt", sep = '\t', row.names = 1)
otu <- read.delim("./00_data/RA_saliva_relative_abundance_clean.tsv",sep = '\t', header = T,row.names = 1,check.names = F)
ps1 <- create_phyloseq(otu_table = otu, sam_data = meta0)
saveRDS(ps1, file = "00_data/ps_saliva_sylph.rds")

#----------------------------- tongue humann3 ko --------------------------------------
meta0 <- read.csv("./00_data/tongue_meta.txt", sep = '\t', row.names = 1)
otu <- read.delim("./00_data/humann3/tongue_genefamilies_merged_cpm_ko_unstratified.tsv",sep = '\t', header = T,row.names = 1,check.names = F)

ps1 <- create_phyloseq(otu_table = otu, sam_data = meta0)
saveRDS(ps1, file = "00_data/ps_tongue_humann3ko.rds")

#----------------------------- saliva humann3 ko --------------------------------------
meta0 <- read.csv("./00_data/saliva_meta.txt", sep = '\t', row.names = 1)
otu <- read.delim("./00_data/humann3/saliva_genefamilies_merged_cpm_ko_unstratified.tsv",sep = '\t', header = T,row.names = 1,check.names = F)

ps1 <- create_phyloseq(otu_table = otu, sam_data = meta0)
saveRDS(ps1, file = "00_data/ps_saliva_humann3ko.rds")

#----------------------------- tongue humann3 pathabundance --------------------------------------
meta0 <- read.csv("./00_data/tongue_meta.txt", sep = '\t', row.names = 1)
otu <- read.delim("./00_data/humann3/tongue_pathabundance_merged_cpm_unstratified.tsv",sep = '\t', header = T,row.names = 1,check.names = F)

ps1 <- create_phyloseq(otu_table = otu, sam_data = meta0)
saveRDS(ps1, file = "00_data/ps_tongue_humann3pathway.rds")

#----------------------------- saliva humann3 pathabundance --------------------------------------
meta0 <- read.csv("./00_data/saliva_meta.txt", sep = '\t', row.names = 1)
otu <- read.delim("./00_data/humann3/saliva_pathabundance_merged_cpm_unstratified.tsv",sep = '\t', header = T,row.names = 1,check.names = F)

ps1 <- create_phyloseq(otu_table = otu, sam_data = meta0)
saveRDS(ps1, file = "00_data/ps_saliva_humann3pathway.rds")


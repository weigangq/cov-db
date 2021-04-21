library(tidyverse)
library(pheatmap)
setwd("c:/Users/lai/Dropbox/cov2")

#####################
# Trace freqs
#######################
# global
x1 <- read_tsv("db-out-Areas/trace_N_America.tsv")
x2 <- read_tsv("db-out-Areas/trace_S_America.tsv")
x3 <- read_tsv("db-out-Areas/trace_Europe.tsv")
x4 <- read_tsv("db-out-Areas/trace_Asia.tsv")
x5 <- read_tsv("db-out-Areas/trace_Africa.tsv")
x6 <- read_tsv("db-out-Areas/trace_Oceania.tsv")

# US States
x7 <- read_tsv("db-out-States//trace_USA.tsv")
x8 <- read_tsv("db-out-States//trace_New_York.tsv")
x9 <- read_tsv("db-out-States//trace_California.tsv")
x10 <- read_tsv("db-out-States//trace_Washington.tsv")
x11 <- read_tsv("db-out-States//trace_Texas.tsv")
x12 <- read_tsv("db-out-States//trace_Michigan.tsv")


# pick only missense & high monthly freq
filter.trace <- function(x, cutoff=20) {
  out <- list("vector")
  count <- 1
  for(i in 1:nrow(x)){
    if(!is.na(x[i,'aaID'])) {
      max.row <- max(x[i,7:19])
      if ( max.row >= cutoff) { # only missense
        out[[count]] <- x[i,c(1,2,4,5,7:19)]
        count <- count + 1
      }
    }
  }
  out.df <- bind_rows(out)
  return(out.df)
}

##########################
# regional heatmap
#######################

x7a <- x7 %>% filter.trace()
x8a <- x8 %>% filter.trace()
x9a <- x9 %>% filter.trace()
x10a <- x10 %>% filter.trace()
x11a <- x11 %>% filter.trace()
x12a <- x12 %>% filter.trace()

x.us <- x7a %>% full_join(x10a, c("varID", "aaID", "locus")) %>% mutate(geo = "US") %>% select(-c("geo.x", "geo.y")) # add Washington
x.us <- x.us %>% full_join(x9a, c("varID", "aaID", "locus")) %>% mutate(geo = "US") %>% select(-c("geo.x", "geo.y")) # add california
x.us <- x.us %>% full_join(x8a, c("varID", "aaID", "locus")) %>% mutate(geo = "US") %>% select(-c("geo.x", "geo.y")) # add NY
x.us <- x.us %>% full_join(x11a, c("varID", "aaID", "locus")) %>% mutate(geo = "US") %>% select(-c("geo.x", "geo.y")) # add Texas
x.us <- x.us %>% full_join(x12a, c("varID", "aaID", "locus")) %>% mutate(geo = "US") %>% select(-c("geo.x", "geo.y")) # add Michigan

trace.states.by.heatmap <- function(x, n=6) {
  tag <- as.character(x[1,'geo'])
  x.mat <- as.matrix(x[,4:(13*n+3)])
  x.mat[is.na(x.mat)] <- 0
  rownames(x.mat) <- paste(x$varID, x$aaID, x$locus, sep="_")
  colnames(x.mat) <- colnames(x)[4:(13*n+3)]
  coul <- colorRampPalette(brewer.pal(9, "Blues"))(100)
  anno_column <- data.frame(Geo = factor(c(rep("USA",13), rep("Washington",13), rep("California", 13), rep("New_York", 13), rep("Texas", 13), rep("Michigan",13))), row.names = colnames(x.mat))
  gene.names.sorted <- names(sort(table(x$locus), decreasing = T))
  gene.colors <- c(1:7, rep(8, length(gene.names.sorted)-7)) # top 7
  names(gene.colors) <- gene.names.sorted
  anno_row <- data.frame(Locus = x$locus, row.names = rownames(x.mat))
  anno_colors <- list(Geo = c(USA = "#ff0000", Washington = "#ff9933", California = "#0099ff", New_York = "#33cc33", Texas = "#ff66ff", Michigan = "#cc0099"), Locus =  gene.colors)
  pheatmap(log10(x.mat+1), cluster_cols = F, scale = "none", clustering_distance_rows = "correlation", angle_col = 45, color = coul, cellwidth = 10, cellheight = 12, fontsize = 8, annotation_col = anno_column, annotation_colors = anno_colors, annotation_row = anno_row, filename = paste(tag, ".pdf", sep=""), main = tag)
}

trace.states.by.heatmap(x.us)


#############################
# global heatmap
################################
x1a <- x1 %>% filter.trace()
x2a <- x2 %>% filter.trace()
x3a <- x3 %>% filter.trace()
x4a <- x4 %>% filter.trace()
x5a <- x5 %>% filter.trace()
x6a <- x6 %>% filter.trace()

x.global <- x4a %>% full_join(x3a, c("varID", "aaID", "locus")) %>% mutate(geo = "global") %>% select(-c("geo.x", "geo.y")) # asia + europe

x.global <- x.global %>% full_join(x7a, c("varID", "aaID", "locus")) %>% mutate(geo = "global") %>% select(-c("geo.x", "geo.y")) # add N.America

x.global <- x.global %>% full_join(x2a, c("varID", "aaID", "locus")) %>% mutate(geo = "global") %>% select(-c("geo.x", "geo.y")) # add S. America

x.global <- x.global %>% full_join(x5a, c("varID", "aaID", "locus")) %>% mutate(geo = "global") %>% select(-c("geo.x", "geo.y")) # add Afica

x.global <- x.global %>% full_join(x6a, c("varID", "aaID", "locus")) %>% mutate(geo = "global") %>% select(-c("geo.x", "geo.y")) # add Oceania

trace.global.by.heatmap <- function(x, n=6) {
  tag <- as.character(x[1,'geo'])
  x.mat <- as.matrix(x[,4:(13*n+3)])
  x.mat[is.na(x.mat)] <- 0
  rownames(x.mat) <- paste(x$varID, x$aaID, x$locus, sep="_")
  colnames(x.mat) <- colnames(x)[4:(13*n+3)]
  coul <- colorRampPalette(brewer.pal(9, "Blues"))(100)
  anno_column <- data.frame(Geo = factor(c(rep("Asia",13), rep("Europe",13), rep("N_America", 13), rep("S_America", 13), rep("Africa", 13), rep("Oceania",13))), row.names = colnames(x.mat))
  gene.names.sorted <- names(sort(table(x$locus), decreasing = T))
  gene.colors <- c(1:7, rep(8, length(gene.names.sorted)-7)) # top 7
  names(gene.colors) <- gene.names.sorted
  anno_row <- data.frame(Locus = x$locus, row.names = rownames(x.mat))
  anno_colors <- list(Geo = c(N_America = "#ff0000", S_America = "#ff9933", Asia = "#0099ff", Europe = "#33cc33", Africa = "#ff66ff", Oceania = "#cc0099"), Locus =  gene.colors)
  pheatmap(log10(x.mat+1), cluster_cols = F, scale = "none", clustering_distance_rows = "correlation", angle_col = 45, color = coul, cellwidth = 10, cellheight = 12, fontsize = 8, annotation_col = anno_column, annotation_colors = anno_colors, annotation_row = anno_row, filename = paste(tag, ".pdf", sep=""), main = tag)
}

trace.global.by.heatmap(x.global)

#####################
# Do Not Use:
# Regional heatmap, NAm/USA/NY
##########################
x7a <- x7 %>% filter.trace()
x8a <- x8 %>% filter.trace()

x9 <- x7a %>% full_join(x8a, c("varID", "aaID", "locus")) %>% mutate(geo = "USA_NY") %>% select(-c("geo.x", "geo.y"))
x10 <- x9 %>% full_join(x1a, c("varID", "aaID", "locus")) %>% mutate(geo = "N_America_USA_NY") %>% select(-c("geo.x", "geo.y"))

x.region <- x10
# heatmaps
# plot three pops (regional)
trace.multi.geo.by.heatmap <- function(x, n) {
  tag <- as.character(x[1,'geo'])
  x.mat <- as.matrix(x[,4:(13*n+3)])
  x.mat[is.na(x.mat)] <- 0
  rownames(x.mat) <- paste(x$varID, x$aaID, x$locus, sep="_")
  colnames(x.mat) <- colnames(x)[4:(13*n+3)]
  coul <- colorRampPalette(brewer.pal(9, "Blues"))(100)
  anno_column <- data.frame(Geo = factor(c(rep("N_America",13), rep("USA",13), rep("New_York", 13))), row.names = colnames(x.mat))
  gene.names.sorted <- names(sort(table(x$locus), decreasing = T))
  gene.colors <- c(1:7, rep(8, length(gene.names.sorted)-7)) # top 7
  names(gene.colors) <- gene.names.sorted
  anno_row <- data.frame(Locus = x$locus, row.names = rownames(x.mat))
  anno_colors <- list(Geo = c(N_America = "red", USA = "green", New_York = "blue"), Locus =  gene.colors)
  pheatmap(log10(x.mat+1), cluster_cols = F, scale = "none", clustering_distance_rows = "correlation", angle_col = 45, color = coul, cellwidth = 15, cellheight = 12, fontsize = 8, annotation_col = anno_column, annotation_colors = anno_colors, annotation_row = anno_row, filename = paste(tag, ".pdf", sep=""), main = tag)
}

trace.multi.geo.by.heatmap(x10,3)

#############################################
# heatmaps
# plot single population
############################################

trace.by.heatmap <- function(x) {
  tag <- as.character(x[1,1])
  out <- list("vector")
  count <- 1
  for(i in 1:nrow(x)){
    if(!is.na(x[i,'conseq'])) {
      max.row <- max(x[i,7:19])
      conseq <- as.character(x[i,'conseq']) 
      if ( max.row >= 10 & conseq == 'missense') { # only missense
        out[[count]] <- x[i,]
        count <- count + 1
      }
    }
  }
  out.df <- bind_rows(out)
  x.mat <- as.matrix(out.df[,7:19])
  rownames(x.mat) <- paste(out.df$varID, out.df$aaID, out.df$locus, sep="_")
  colnames(x.mat) <- colnames(x)[7:19]
  sideCol <- if_else(out.df$conseq  == 'missense' , "red", if_else(out.df$conseq == 'synonymous', "green", "gray"))
  coul <- colorRampPalette(brewer.pal(9, "Blues"))(100)
  anno.snp <- data.frame(row.names = rownames(x.mat), conseq = out.df$conseq)
  ann_colors = list(conseq = c(synonymous = 'green', missense = 'red', noncoding = "gray"))
  #pheatmap(log10(x.mat+1), cluster_cols = F, scale = "none", clustering_distance_rows = "correlation", angle_col = 45, color = coul, annotation_row = anno.snp, annotation_colors = ann_colors, cellwidth = 15, cellheight = 12, fontsize = 8, filename = paste(tag, ".pdf", sep=""), main = tag)
  pheatmap(log10(x.mat+1), cluster_cols = F, scale = "none", clustering_distance_rows = "correlation", angle_col = 45, color = coul, cellwidth = 15, cellheight = 12, fontsize = 8, filename = paste(tag, ".pdf", sep=""), main = tag, na_col = 'gray')
  
  trace.by.heatmap(x1)
  trace.by.heatmap(x2)
  trace.by.heatmap(x3)
  trace.by.heatmap(x4)
  trace.by.heatmap(x5)
  trace.by.heatmap(x6)
  trace.by.heatmap(x7)
  trace.by.heatmap(x8)
  

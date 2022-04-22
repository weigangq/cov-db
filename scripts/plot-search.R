setwd("C:/Users/weiga/Dropbox/cov-db/scripts/")
library(tidyverse)
#library(stringdist)
library(interp)
library(ggrepel)
library(plotly)
library(threejs)
#########################
# plot landscapes

# recipe: https://rviews.rstudio.com/2020/12/14/plotting-surfaces-with-r/
plotly.land <- function(proj_file) {
  fn <- read_csv(proj_file)
  origdata <- as.data.frame(list(x = fn$Dim_1, y = fn$Dim_2, z = fn$fit))
  grid <- with(origdata, interp::interp(x, y, z))
  plot_ly(x = grid$x, y = grid$y, z = grid$z, name = paste("landscape", proj_file)) %>% 
    add_surface()
}

plotly.land("monotonic-pca.csv")
plotly.land("normal-pca.csv")
fn <- read_csv("additive-pca.csv")

#origdata <- as.data.frame(list(x = fn$Dim_1, y = fn$Dim_2, z = fn$fit))
#grid <- with(origdata, interp::interp(x, y, z))
#plot_ly(x = grid$x, y = grid$y, z = grid$z, name = ) %>% add_surface()

# recipe: https://ggplot2.tidyverse.org/reference/geom_contour.html
plot.land <- function(proj_file, df_end) {
  fn <- read_csv(proj_file)
  origdata <- as.data.frame(list(x = fn$Dim_1, y = fn$Dim_2, z = fn$fit))
  grid <- with(origdata, interp::interp(x, y, z))
  griddf <- subset(data.frame(x = rep(grid$x, nrow(grid$z)),
                              y = rep(grid$y, each = ncol(grid$z)),
                              z = as.numeric(grid$z)), !is.na(z))
  ggplot(griddf, aes(x, y, z = z)) +
    geom_contour_filled(bins = 10, alpha = 0.8) +
    geom_point(data = df_end, size = 3, color = "red", shape = 1) +
    geom_text(data = df_end, aes(label = id), vjust = 0, nudge_y = 0.1, color = "red", size = 6) +
    labs(title = paste("Fitness landscape: ", proj_file)) +
    theme_bw()  
}

#plot.land("normal-tsne.csv")
#plot.land("monotonic-tsne.csv")
#plot.land("normal-pca.csv")
#plot.land("monotonic-pca.csv")

############################
## Objective search outputs
x <- read_tsv("objective-search-2-landscapes.tsv", col_names = F)
colnames(x) <- c("tag", "gen", "close_id", "close_fit", "elite_id", "diff_close", "diff_fit", "landscape", "search_algorithm")

x.max.fit <- x %>% group_by(landscape) %>% summarise(max = max(close_fit))
x.end <- x %>% filter(close_id == 'H000')
x.last <- x %>% group_by(tag) %>% 
  summarise(gen = max(gen)) %>% 
  left_join(x, c("tag", "gen")) %>%
  filter(close_id != 'H000')

x %>% ggplot(aes(x = gen, y = close_fit, group = tag )) + 
  geom_line() + 
  geom_point(data = x.end, aes(gen, close_fit), shape =1, color = 2) +
  geom_point(data = x.last, aes(gen, close_fit), shape = 1, color = 1) +
  facet_wrap(~landscape) + 
  theme_bw() +
  geom_hline(data = x.max.fit, aes(yintercept = max), linetype = 2, color = 2) + 
  geom_hline(yintercept = 0, linetype = 2) + 
  labs(title = "Objective search: random fitness lanscape is deceptive", subtitle = "pop size:100; genome size: 50 bits; fitness: Normal (0,1); evolve for max 100 generations; 20 runs on each landscape")

x %>% group_by(tag, landscape) %>% 
  summarise(gen = max(gen)) %>% 
  ggplot(aes(x = landscape, y = gen)) +
#  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape=1, width = 0.25) + 
  theme_bw()

x %>% group_by(tag, landscape) %>% 
  summarise(gen = max(gen)) %>% 
  t.test(data = ., gen ~ landscape)

pca.norm <- read_csv("normal-pca.csv")
x.last.norm <- x %>% filter(landscape == 'normal') %>% 
  group_by(tag) %>% 
  summarise(gen = max(gen)) %>% 
  left_join(x, c("tag", "gen")) %>% 
  dplyr::select(tag, close_id, close_fit) %>% 
  left_join(pca.norm, c("close_id" = "label")) %>% 
  dplyr::select(2, 4:6) %>% 
  dplyr::rename(x = Dim_1, y = Dim_2, z = fit, id = close_id)

pca.mono <- read_csv("monotonic-pca.csv")
x.last.mono <- x %>% filter(landscape == 'monotonic') %>% 
  group_by(tag) %>% 
  summarise(gen = max(gen)) %>% 
  left_join(x, c("tag", "gen")) %>% 
  dplyr::select(tag, close_id, close_fit) %>% 
  left_join(pca.norm, c("close_id" = "label")) %>% 
  dplyr::select(2,4:6) %>% 
  dplyr::rename(x = Dim_1, y = Dim_2, z = fit, id = close_id)

plot.land("normal-pca.csv", x.last.norm)
plot.land("monotonic-pca.csv", x.last.mono)

plotly.land("normal-pca.csv")

setwd("~/GitHub/cov-db/sim-data/binary-sim")
library(tidyverse)
library(stringdist)
library(ggrepel)
library(akima)
library(lattice)s

## Landscape
land <- read_tsv("landscapeAdd.tsv", col_types = "cic")
# Calculate differences between strings as distance matrix
mat <- stringdistmatrix(land$String)
# Put the distance matrix into a 2-coordinate point for each string
cmd0 <- cmdscale(mat)
cmd0.df <- as.data.frame(cmd0)
fit <- land$Fitness
cmd.df <- cmd0.df %>% mutate(fit = fit )
colnames(cmd.df) <- c("x", "y", "fit")
# String data points with fitness values for landscape
fld <- with(cmd.df, interp(x = x, y = y, z = fit))

## Elite
elite <- read_tsv("eliteAddObj.tsv", col_types = "icicc")
head(elite)
# Calculate differences between strings as distance matrix
mat1 <- stringdistmatrix(elite$String)
# Put the distance matrix into a 2-coordinate point for each string
cmd1 <- cmdscale(mat1)
cmd1.df <- as.data.frame(cmd1)
cmd1.df <- cbind(cmd1.df, elite$Generation)
colnames(cmd1.df) <- c("x", "y", "Generation")

# Combined strings from landscape and elite
land_string <- as.data.frame(land$String)
colnames(land_string) <- c("String")
elite_string <- as.data.frame(elite$String)
colnames(elite_string) <- c("String")
strings <- rbind(land_string, elite_string)
mat2 <- stringdistmatrix(strings$String)
# Turn distance matrix into a data frame
mat2.dict <- as.matrix(mat2)
mat2.df <- as.data.frame(mat2.dict)

## Find closest string from landscape to elite strings
# select only elite as rows and landscape as columns
mat2.df <- mat2.df %>% slice(201:300)
mat2.df <- mat2.df %>% select(1:200)
# Find string in landscape that has the closest distance for each string in elite
closest <- as.matrix(apply( mat2.df, 1, which.min))
closest <- as.data.frame(closest)
closest <- closest %>% mutate(elite_string_index= 200 + 1:n())
colnames(closest) <- c("land_string_index", "elite_string_index")
# Minimum distance of elite to landscape point
#min_dist <- apply(mat2.df, MARGIN =  1, FUN = min, na.rm = T)
#min_dist <- as.data.frame(min_dist)
# Generate string_index for cmd0.df
cmd0.df <- cmd0.df %>% mutate(land_string_index= 1:n())
# Assign coordinates of closest landscape string as coordinates of elite string
elite_coor <- merge(closest,cmd0.df)
elite_coor <- elite_coor %>% select(V1, V2)
colnames(elite_coor) <- c("x", "y")

## Fitness landscape with elite points
filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               plot.title = title(main = "Additive-Objective",
                                  xlab = "x",
                                  ylab = "y"),
               key.title = title(main="Fitness ", cex.main = 1),
               plot.axes={points(elite_coor$x,elite_coor$y,
                                 pch = 19,
                                 col = "blue");
                 axis(1);axis(2); contour(fld, add = TRUE, lwd = 0.5)})

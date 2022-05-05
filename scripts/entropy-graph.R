library(tidyverse)
setwd("c:/users/winston/cov-db/scripts")

neut <- "../sim-data/entropy-neut1.tsv"
adpt <- "../sim-data/entropy-adpt.tsv"
bkg <- "../sim-data/entropy-bkg.tsv"
mix <- "../sim-data/entropy-mix.tsv"

entropy <- read_tsv(neut)
glimpse(entropy)

ggplot(entropy, aes(x=Generation, y=Entropy, color=Type)) + 
  geom_line() +
  theme_bw() + 
  xlab("Generation") + 
  ylab("Total Entropy (bits)") +
  labs(title = "Total Entropy Over Evolutionary Time")


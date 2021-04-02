library(tidyverse)
library(broom)
library(ggrepel)
library(ape)
setwd("../Dropbox/cov2/wright-fisher-sim/")

##################
# temporal analysis
#####################

x <- read_tsv("neut-temporal.tsv", col_names = F)
colnames(x) <- c("ref", "gen", "diff")
x %>% filter(gen <= 50) %>% ggplot(aes(x = gen, y = diff)) + geom_jitter(shape=1, color="gray") + facet_wrap(~ref) + theme_bw() + geom_smooth(method = "lm")

x %>% filter(gen <= 50) %>% ggplot(aes(x = diff)) + geom_density(fill="gray") + facet_wrap(~ref) + theme_bw() 

x1 <- x %>% filter(gen <=50)
x1.lm <- x1 %>% group_by(ref) %>% do(tidy(lm(diff ~ gen, data = .)))
x1.lm %>% filter(term == 'gen') %>% ggplot(aes(x=ref, y=estimate, label = ref)) + geom_point() + theme_bw() + geom_text_repel() + geom_smooth()


##################
# muller plot
#######################
x <- read_tsv("adpt-hi-fit-sites-cts.tsv", col_names = F)
colnames(x) <- c("snp", "gen", "ct")

x.wide <- x %>% pivot_wider(id_cols = c("snp", "gen"), names_from = "gen", values_from = "ct")

snps <- names(table(x$snp))
out <- list("vector", length(snps))
for(i in seq_along(snps)) {
  first_gen = 0
  xsub <- x %>% filter(snp == snps[i])
  for(g in 1:500){
    if(xsub[g, 'ct'] > 0) {
      first_gen = g
      break()
    }
  }
  out[[i]] <- tibble(snp = snps[i], first = first_gen)
}
first.gen <- bind_rows(out) %>% arrange(first)
x$snp <- factor(x$snp, levels = first.gen$snp)

x %>% ggplot(aes(x = gen, y = ct)) + geom_line() + theme_bw() +  facet_wrap(~snp, nrow=4) 

#x %>% ggplot(aes(x = gen, y = ct)) + geom_line() + theme_bw() +  geom_vline(data = first.gen, aes(xintercept = first), color = "red") + facet_wrap(~snp) 

##################
# site analysis
##################
x1 <- read_tsv("neut-pi.tsv", col_names = F)
#x1 <- x1 %>% mutate(model = "neut")
x2 <- read_tsv("bkg-pi.tsv", col_names = F)
#x2 <- x2 %>% mutate(model = "bkg")
x3 <- read_tsv("adpt-pi.tsv", col_names = F)
#x3 <- x3 %>% mutate(model = "adpt")

x <- bind_rows(x1, x2, x3)
colnames(x) <- c("model", "gen", "all", "syn", "mis")
x.long <- x %>% pivot_longer(3:5, names_to = "conseq", values_to = "pi")
x.long %>% ggplot(aes(x = gen, y = pi, color = conseq)) + geom_point(shape=1, alpha=0.5) + theme_bw() + facet_wrap(~model) + geom_smooth() + geom_hline(yintercept = 40, linetype = "dashed")

# sliding-window avgs:
# Don't use. Eough to use geo_smooth()
out <- vector("list", length(breaks))
win.size = 20
for(i in 1:480) {
  x.br <- x %>% filter(gen < i + win.size & gen >= i)
  x.mean <- x.br %>% group_by(model) %>% summarise(mean.pi = mean(pi))
  out[[i]] <- x.mean %>% mutate(gen.end = win.size + i)
}
x.win <- bind_rows(out)
x.win$model <- factor(x.win$model, levels = c("adpt", "neut", "bkg"))
x.win %>% ggplot(aes(x = gen.end, y = mean.pi, color = model)) + theme_bw() + geom_line(size = 2) + geom_hline(yintercept = 40, linetype= "dashed")


###############
# test runs
################
x.sample <- read_tsv("test-samples.tsv", col_names = F)
colnames(x.sample) <- c("tag", "gen", "id", "fit", "igs", "syn", "mis")    
x.sample <- x.sample %>% mutate(mut = igs + syn + mis)
x.long <- x.sample %>% pivot_longer(5:8, names_to = "mut_type", values_to = "num_mut")
x.long %>% filter(mut_type != 'mut' & mut_type != 'igs') %>% ggplot(aes(x = gen, y = num_mut, color = mut_type)) + geom_jitter(shape = 1, alpha = 0.2) + geom_smooth(method = "lm") + theme_bw() + geom_abline(intercept = 0, slope = 0.1, size=1.5) + scale_color_manual(values=c("magenta", "cyan"))
p = x.long %>% filter(mut_type == 'mut') %>% ggplot(aes(x = gen, y = num_mut)) + geom_jitter(size = 2, alpha=0.1) +  facet_wrap(~tag)+ theme_bw() + geom_smooth(method = "lm") 

d <- seq(0, 200) # generations
x.exp <- d * 0.1 # expected per site per day
x.sd <- sqrt(x.exp) # sd
x.err <- data.frame(day = d, exp.site = x.exp, x.hi = x.exp + 2*x.sd, x.lo = x.exp - 2*x.sd )
p + geom_line(data = x.err, aes(x = day, y = x.hi), linetype = "dashed", color="black", size = 1.5) + geom_line(data = x.err, aes(x = day, y = x.lo), linetype = "dashed", color="black", size=1.5) + geom_line(data = x.err, aes(x = day, y = exp.site), color="black", size=1.5)

##############
# sample analysis
##################
x.neut <- read_tsv("neut-samples.tsv", col_names = F)
x.adap <- read_tsv("adpt-samples.tsv", col_names = F)
x.bkg <- read_tsv("bkg-samples.tsv", col_names = F)
#x.mix <- read_tsv("mixed-samples.tsv", col_names = F)


x.sample <- bind_rows(x.neut, x.adap, x.bkg)
colnames(x.sample) <- c("tag", "gen", "id", "fit", "igs", "syn", "mis")    

x.sample <- x.sample %>% mutate(mut = igs + syn + mis)

x.fit <- x.sample %>% filter(fit > 0) %>% group_by(gen, tag) %>% summarise(mean.fit = mean(fit))

x.fit$tag <- factor(x.fit$tag, levels = c('adpt', 'neut', 'bkg'))

x.long <- x.sample %>% pivot_longer(5:8, names_to = "mut_type", values_to = "num_mut")

#x.long <- x.long %>% filter(tag != 'mixed')

x.long$tag <- factor(x.long$tag, levels = c("neut", "bkg", "adpt"))

x.long <- x.long %>% mutate(fit2 = if_else(fit > 1, 'pos', if_else(fit < 1, 'neg', 'neut')))

x.long$fit2 <- factor(x.long$fit2, levels = c('pos', 'neut', 'neg'))

# fit Poisson
d <- seq(0, 500) # generations
x.exp <- d * 0.1 # expected per site per day
x.sd <- sqrt(x.exp) # sd

x.err <- data.frame(day = d, exp.site = x.exp, x.hi = x.exp + 2*x.sd, x.lo = x.exp - 2*x.sd )


# facet by model
x.long %>% filter(mut_type != 'mut' & mut_type != 'igs') %>% ggplot(aes(x = gen, y = num_mut, color = mut_type)) + geom_jitter(shape = 1, alpha = 0.2) + geom_smooth(method = "lm") + theme_bw() + geom_abline(intercept = 0, slope = 0.1, size=1.5) + facet_wrap(~tag) +  scale_color_manual(values=c("magenta", "cyan"))

# total
p = x.long %>% filter(mut_type == 'mut') %>% ggplot(aes(x = gen, y = num_mut)) + geom_jitter(size = 2, aes(color=fit2), alpha=0.1) +  facet_wrap(~tag)+ theme_bw() + geom_smooth(method = "lm") 
 
p + geom_line(data = x.err, aes(x = day, y = x.hi), linetype = "dashed", color="black", size = 1.5) + geom_line(data = x.err, aes(x = day, y = x.lo), linetype = "dashed", color="black", size=1.5) + geom_line(data = x.err, aes(x = day, y = exp.site), color="black", size=1.5)

# plot fitness

x.fit %>% ggplot(aes(x = gen, y = mean.fit, color = tag)) + geom_line(size=2) + theme_bw()

# facet by mut_type
# x.long %>% filter(mut_type != 'igs') %>% ggplot(aes(x = gen, y = num_mut, color = model)) + geom_jitter(shape = 1, alpha = 0.2) + geom_smooth(method = "lm") + theme_bw() + geom_abline(intercept = 0, slope = 0.1) + facet_wrap(~mut_type)


#################################
# gene analysis
######################

gene.name <- read_tsv("gene-len-name.txt")
obs.gene <- read_csv("obs-genes.csv")

x.neut <- read_tsv("neut-genes.tsv", col_names = F)
x.bkg <- read_tsv("bkg-genes.tsv", col_names = F)
x.adapt <- read_tsv("adpt-genes.tsv", col_names = F)

x.gene = bind_rows(x.neut, x.bkg, x.adapt)
colnames(x.gene) <- c("model", "gid", "len", "prod", "syn", "mis")                                  

x.gene <- x.gene %>% mutate(syn2 = syn/len, mis2 = mis/len)

x.gene$model <- factor(x.gene$model, levels = c('adpt', 'neut', 'bkg'))

x.gene <- x.gene %>% left_join(gene.name, c("gid" = "id"))

obs.gene <- obs.gene %>% left_join(gene.name, c("Gene" = "name"))

obs.gene <- obs.gene %>% mutate(syn2 = Synonymous/len, mis2 = Missense/len)

g1 <- x.gene %>% select(model, syn2, mis2, name)
g2 <- obs.gene %>% select(syn2, mis2, Gene)
g2 <- g2 %>% mutate(model = rep("obs", 25)) 
g2 <- tibble(model= g2$model, syn2= g2$syn2, mis2 = g2$mis2, name = g2$Gene)

x.gene2 <- bind_rows(g1, g2)

p <- ggplot(data = x.gene, aes(x = syn2, y = mis2, color = model)) + geom_point() + geom_smooth(method = "lm") + theme_bw() + geom_abline(intercept = 0, slope = 1:3, linetype = "dashed") + geom_text_repel(aes(label = name)) 
  
#ggplot(data = x.gene2, aes(x = syn2, y = mis2, color = model)) + geom_point() + geom_smooth(method = "lm") + theme_bw() + geom_abline(intercept = 0, slope = 1:3, linetype = "dashed") + geom_text_repel(aes(label = name)) 

x.gene2$model <- factor(x.gene2$model, levels = c('adpt', 'neut', 'bkg', 'obs'))
x.obs <- x.gene2 %>% filter(model == 'obs')
x.gene2 %>% filter(model != 'obs') %>% ggplot(aes(x = syn2, y = mis2, color = model)) + geom_point(alpha=0.3) + geom_smooth(method = "lm") + theme_bw() + geom_abline(intercept = 0, slope = 1:3, linetype = "dashed") + geom_point(data = x.obs, aes(x=syn2, y=mis2), size=5, alpha = 0.5, color = "orange") + geom_text_repel(data = x.obs, aes(x=syn2, y=mis2, label=name), color = "orange")


x.gene %>% group_by(model) %>% do(tidy(lm(data = ., mis2 ~ syn2)))


###########################
# variant analysis
#########################

x.neut <- read_tsv("neut-vars.tsv", col_names = F)
x.bkg <- read_tsv("bkg-vars.tsv", col_names = F)
x.adapt <- read_tsv("adpt-vars.tsv", col_names = F)

x.var <- bind_rows(x.neut, x.bkg, x.adapt)
colnames(x.var) <- c("model", "pos", "snp.id", "gene.id", "ref.nt", "alt.nt", "ref.aa", "al.aa", "ref.codon", "alt.codon", "conseq", "fit", "count", "mult")
x.var$model <- factor(x.var$model, levels = c('adpt', 'neut', 'bkg'))

x.var2 <- x.var %>% mutate(fit2 = if_else(fit > 1, "adaptive", if_else(fit == 1, "neutral", "delet")))

x.var2$fit2 <- factor(x.var2$fit2, levels = c('adaptive', 'neutral', 'delet'))

x.var2 %>% filter(conseq == 'missense') %>% ggplot(aes(x = count, fill = fit2)) + geom_density(alpha=0.5) + theme_bw() + facet_wrap(~model) + scale_x_log10()


#################
# recurrence
##################

x.var2 %>% ggplot(aes(x = mult, fill = fit2)) + geom_bar(alpha=0.5) + theme_bw() + facet_wrap(~model)

x.var2 %>% ggplot(aes(x = mult, fill = conseq)) + geom_bar(alpha=0.5) + theme_bw() + facet_wrap(~model)

##################
# Coalescence
#################

t.neut <- read.tree("neut.dnd")
t.bkg <- read.tree("bkg.dnd")
t.adapt <- read.tree("adpt.dnd")

par(mfrow = c(4,1))
plot(t.neut, font = 1, edge.width = 2, no.margin = T, edge.color = "green", show.tip.label = F)
abline(v=c(0, 100, 200, 300, 400, 500), col = "gray")
plot(t.bkg, font = 1, edge.width = 2, no.margin = T, edge.color = "blue", show.tip.label = F)
abline(v=c(0, 100, 200, 300, 400, 500), col= "gray")
plot(t.adapt, font = 1, edge.width = 2, no.margin = T, edge.color = "red", show.tip.label = F)
axisPhylo(cex=2)
abline(v=c(0, 100, 200, 300, 400, 500), col= "gray")
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

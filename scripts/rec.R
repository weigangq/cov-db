################
# recombiantion analysis
#################

library(lubridate)
library(ggpubr)
library(grid)
library(gridExtra)
library(tidyverse)
setwd("C:/Users/lai/Dropbox/cov-browser/")
#library(tidyverse)

###############
# recurrence by periods (not cumulative)
################

x <- read_tsv("freq-sam.tsv", col_names = c("gen", "site", "freq", "recur", "gs"))
x %>% ggplot(aes(x=recur)) + geom_bar() + facet_wrap(~gen)

###################
# prob of recurrent mutation over time
##################

t1 <- read_tsv("freq-sam10.tsv", col_names = c("site", "freq", "recur", "gen"))
t1 <- t1 %>% mutate(sample.size=10)

t2 <- read_tsv("freq-sam20.tsv", col_names = c("site", "freq", "recur", "gen"))
t2 <- t2 %>% mutate(sample.size=20)

t3 <- read_tsv("freq-sam50.tsv", col_names = c("site", "freq", "recur", "gen"))
t3 <- t3 %>% mutate(sample.size=50)

t4 <- read_tsv("freq-sam100.tsv", col_names = c("site", "freq", "recur", "gen"))
t4 <- t4 %>% mutate(sample.size=100)

t5 <- read_tsv("freq-sam200.tsv", col_names = c("site", "freq", "recur", "gen"))
t5 <- t5 %>% mutate(sample.size=200)

t6 <- read_tsv("freq-sam400.tsv", col_names = c("site", "freq", "recur", "gen"))
t6 <- t6 %>% mutate(sample.size=400)

t.gen <- t # varying number of generations
t.sam <- bind_rows(t1, t2, t3, t4, t5, t6)

t.sam %>% ggplot(aes(x=recur, color = as.character(sample.size))) + geom_freqpoly(bins=30) + theme_bw()

t.sam %>% ggplot(aes(x=recur,  after_stat(density), color = as.character(sample.size))) + geom_freqpoly(bins=30) + theme_bw()

t %>% ggplot(aes(x=total.gen)) + geom_bar()
t.gr <- t %>% group_by(total.gen) %>% count()
t %>% group_by(total.gen) %>% summarise(single = sum())

t.sam %>% ggplot(aes(x=recur)) + geom_bar() + theme_bw() + facet_wrap(~sample.size)

##############
# fit poisson:
t.gen %>% ggplot(aes(x=recur)) + geom_bar() + theme_bw() + facet_wrap(~total.gen)
cts <- 1:30
t.gen.total <- t.gen %>% group_by(total.gen) %>% count()
t.gen.mean <- t.gen %>% group_by(total.gen) %>% summarise(m=mean(recur))
out <- vector("list", length = 6)
for(i in 1:length(t.gen.total$total.gen)) {
  x <- dpois(cts, lambda = as.numeric(t.gen.mean[i,2])) * as.numeric(t.gen.total[i,2])
  out[[i]] <- data.frame(total.gen=rep(as.numeric(t.gen.mean[i,1]),length(cts)), recur=cts, expected=x)
}
out.gen <- bind_rows(out)
t.gen.cts <- t.gen %>% group_by(total.gen, recur) %>% count()
gen <- t.gen.cts %>% left_join(out.gen, c("total.gen", "recur"))

gen %>% ggplot(aes(x=recur, y=n)) + geom_bar(stat = "identity", fill="gray") + theme_bw() + geom_line(aes(x=recur, y=expected, group=total.gen), color="blue") + geom_point(aes(x=recur, y=expected), color="blue", fill="white", shape=1) + facet_wrap(~total.gen) + ylab("num. sites") + xlab("num. mutation") + ggtitle("Number of recurent mutations increases with generation")

t.gen %>% ggplot(aes(y=freq, x=recur)) + geom_jitter(shape=1) + facet_wrap(~total.gen) + theme_bw() + ggtitle("High-freq SNPs emerge with high avg recurrence")

t.sam %>% ggplot(aes(y=freq, x=recur)) + geom_jitter(shape=1) + facet_wrap(~sample.size) + theme_bw() + ggtitle("Recurrent SNPs are low in frequency")
###########################

#################
# fit poisson for sample sizes
t.sam.total <- t.sam %>% group_by(sample.size) %>% count()
t.sam.mean <- t.sam %>% group_by(sample.size) %>% summarise(m=mean(recur))
out <- vector("list", length = 6)
for(i in 1:length(t.sam.total$sample.size)) {
  x <- dpois(cts, lambda = as.numeric(t.sam.mean[i,2])) * as.numeric(t.sam.total[i,2])
  out[[i]] <- data.frame(sample.size=rep(as.numeric(t.sam.mean[i,1]),length(cts)), recur=cts, expected=x)
}
out.sam <- bind_rows(out)
t.sam.cts <- t.sam %>% group_by(sample.size, recur) %>% count()
sam <- t.sam.cts %>% left_join(out.sam, c("sample.size", "recur"))

sam %>% ggplot(aes(x=recur, y=n)) + geom_bar(stat = "identity", fill="gray") + theme_bw() + geom_line(aes(x=recur, y=expected, group=sample.size), color="blue") + geom_point(aes(x=recur, y=expected), color="blue", fill="white", shape=1) + facet_wrap(~sample.size) + ylab("num. sites") + xlab("num. mutation") + ggtitle("Number of recurent mutations do not increase with sample size")



t.gr %>% ggplot(aes(x=total.gen, y=n)) + geom_point(shape=1, size=5) + geom_line() + theme_bw() + ylab("num. sites")


###################
# simulation tests: validation for equilibrium
# u=0.1, N=100, g=1000
########################

# check slope & variance

x <- read_tsv("sam.tsv", col_names = c("gen", "mut"))
x %>% ggplot(aes(x=gen, y=mut)) + geom_jitter(color="gray") + geom_abline(slope = 0.1, intercept = 0, color="red") + theme_bw() 

# check mutation-drif balance

theda <- read_tsv("theta.tsv", col_names = c("gen", "pi"))
theda %>% ggplot(aes(x=gen, y=pi)) + geom_point(shape=1) + geom_line(color="gray") + theme_bw() + geom_hline(yintercept = 20, color="red") + annotate(geom = "text", x=600, y=35, label="genome length = 29903\npop size = 100\nmut rate = 0.1 per genome per gen\nsample size = 10") + ggtitle("Wright-Fisher simulator of CoV genome evolution") + xlab("generation") + ylab("pairwise seq diff")

# Check seg sites and haplotypes

hap <- read_tsv("hap.tsv", col_names = c("gen", "segs", "haps", "cum.segs", "cum.haps"))
hap.long <- hap %>% gather(2:5, key = "type", value = "counts")
hap.long %>% ggplot(aes(x=gen, y=counts, color = type)) + geom_point(shape=1)  + theme_bw() + annotate(geom = "text", x=250, y=3000, label="genome length = 29903\npop size = 100\nmut rate = 0.1 per genome per gen\nsample size = 10") + ggtitle("Wright-Fisher simulator of CoV genome evolution") + xlab("generation") + scale_y_log10()

   
#######################
# recurrent mutation 
#######################

x <- read_tsv("freq.tsv", col_names = c("site", "freq", "recur", "gen"))
x %>% ggplot(aes(y=freq, x=recur)) + geom_jitter(shape=1, size=2) + theme_bw()

x %>% ggplot(aes(x=recur, after_stat(density))) + geom_histogram(binwidth = 1, fill="gray", color="blue") + theme_bw() + annotate(geom = "text", x=4, y=0.4, label="genome length = 29903\npop size = 100\nmut rate = 0.1 per genome per generation\nsample size = 10\ngen = 200\nmean = 1.63 mut/site") + ggtitle("Wright-Fisher simulator of CoV genome evolution") + xlab("num recurrent mutations")

table(x$recur)

intervals <- str_split(x$gen, "-") %>% lapply(function(a) {y <- as.integer(a); diff(y)}) %>% unlist()
stem(intervals)

##################
# trace lineages
###################

x <- read_tsv("samples.tsv", col_names = c("gen", "samp", "level", "anc"))

x.sam <- x %>% filter(gen==10)

p <- x.sam %>% ggplot(aes(x=anc, y=level)) + theme_bw()

for (i in 0:8) {
  df1 <- x.sam %>% filter(level == i)
  df2 <- x.sam %>% filter(level == i + 1)
  df <- data.frame(samp = df1$samp, lev1 = df1$level, lev2=df2$level, anc1=df1$anc, anc2=df2$anc)
  p <- p + geom_segment(data=df, aes(x=anc1, y=lev1, xend=anc2, yend=lev2), color="gray") + geom_point(data=df1, aes(x=anc, y=level), shape=1) + geom_point(data=df2, aes(x=anc, y=level), shape=1)
}

p
###################
# MLE of genome dates
#######################

x <- vector("list", 10)
for(k in 1:10) {
  t <- 1:200
  p <- dpois(k, lambda = t * 0.1, log=T)
  x[[k]]<- data.frame(num.mut =rep(k,100), num.day=t, prob=p)
}

x.df <- bind_rows(x)

ggplot(x.df, aes(x=num.day, y=prob, color=as.character(num.mut), group=as.character(num.mut))) + geom_line() + theme_bw()


####################
#  sim-wright-fisher
################

sam <- read_tsv("../../cov-browser/samples2.tsv", col_names = c("gen", "num.mut"))

p <- ggplot(sam, aes(x=gen, y=num.mut)) + geom_jitter(shape=1, color="gray") + theme_bw() + annotate(geom = "text", x=20, y=17, label="pop size = 1000\nmut rate = 0.1 per genome per generation\nsample size = 100") + ggtitle("Wright-Fisher simulator of CoV genome evolution")

d <- seq(0,100) # days
x.exp <- d * 0.1 # expected per site per day
x.sd <- sqrt(x.exp) # sd

x.err <- data.frame(day = d,exp.site = x.exp, x.hi = x.exp + 2*x.sd, x.lo = x.exp - 2*x.sd )

p + geom_line(data = x.err, aes(x = day, y = x.hi), linetype = "dashed", color="red") + geom_line(data = x.err, aes(x = day, y = x.lo), linetype = "dashed", color="red") + geom_line(data = x.err, aes(x = day, y = exp.site), color="red") 

#######################
frq <- read_tsv("../../cov-browser/sam", col_names = c("site", "freq"))

frq %>% ggplot(aes(x=freq)) + geom_histogram(binwidth = 5e-3) + theme_bw() + xlab("SNP frequency")


#####################
# simulating Poisson mutation
# prob of no recurrent mutation, modeled after the birthday problem
#########################

# first, plot birthday problem
n <- 1:100
y <- vector("double", length = max(n))
for(i in 1:50) {
  y[i] <- prod(seq(from=365-i+1, to=365))/365**i
}

plot(n, 1-y, type = "l", las=1 )
abline(h=0.5, lty=2)

# if lambda=1 mutation for each sample
n <-1:1000 # sample size
len.genome <- 29903
y <- vector("double", length = max(n))
for(i in 1:length(n)) {
  y[i] <- sum(log(seq(from=len.genome-i+1, to=len.genome))) - i*log(len.genome)
}

plot(n, 1-exp(y), type = "l", las=1, xlab = "sample size", ylab = "Prob(recurrent mutations)")
abline(h=0.5, lty=2)



###########################
mu <- 3e-6
days <- 100
num.sites <- function(num.days) {
  p <- rpois(n = len.genome, lambda = mu * num.days)
  sites <- which(p>0)
}

x <- read_tsv("../Dropbox/cov-browser/muts.txt", col_names = c("day", "mut"))

epi.lm <- lm(mut ~ day, data=x)

p <- ggplot(x, aes(x=day, y=mut)) + geom_jitter(shape=1, color="gray") + geom_smooth(method = "lm") + theme_bw()

d <- seq(10,150) # days
x.exp <- d * epi.lm$coefficients[2] + epi.lm$coefficients[1] # expected per site per day
x.sd <- sqrt(x.exp) # sd

x.err <- data.frame(day = d,exp.site = x.exp, x.hi = x.exp + 2*x.sd, x.lo = x.exp - 2*x.sd )

p + geom_line(data = x.err, aes(x = day, y = x.hi), linetype = "dashed") + geom_line(data = x.err, aes(x = day, y = x.lo), linetype = "dashed") + geom_line(data = x.err, aes(x = day, y = exp.site)) 

###########################
# 6/8/2020. jittering
#######################

y <- read_csv("../ycoord.csv")
y %>% ggplot(aes(x=st, y=nsnps)) + geom_point(shape=1) + geom_boxplot()

y <- read_csv("../yjitter.csv")
y %>% ggplot(aes(x=before, y=after)) + geom_point(shape=1, size=5) + geom_abline(intercept = 0, slope=1, color="red") + theme_bw()

####################
# 6/7/2020, data-5-12-2020 analysis
###########################
setwd("C:/Users/lai/Dropbox/cov-browser/data-05-12-2020/genome-analysis")
num.snp <- read_tsv("totalSNPs.txt2", col_names = c("acc", "no.snp"))
avg.id <- read_tsv("avgIdentity.txt2", col_names = c("acc", "avg.ref", "avg.acc"))
geo.date <- read_tsv("../acc-by-month/acc-date-geo.tsv", col_names = c("acc", "date", "geo.id", "area.id"))

acc.date <- geo.date %>% left_join(num.snp, "acc") %>% left_join(avg.id, "acc")

ref.length <- 29903;
# EPI_ISL_406030 (1/10 colletion date; 
# not EPI_ISL_403932, 1/14)
ref.date <- as.Date('2020-01-10')

acc.date <- acc.date %>% mutate(date2 = as.Date(date), diff.snp = (100-avg.ref)/100, diff.genome = no.snp/ref.length, time = date2 - ref.date)

acc.date %>% filter(diff.genome < 7.5e-4) %>% ggplot(aes(x=diff.genome)) + geom_histogram(bins = 100) # (remve top 30 rows as outliers)

acc.date <- acc.date %>% filter(diff.genome < 7.5e-4)

p1 = ggplot(data = acc.date, aes(x = date2, y = diff.genome)) + geom_jitter(aes(color = as.character(area.id)), shape = 1, size = 1) + theme_bw() + xlab("collection date") + ylab("seq.diff (per site)")+  theme(axis.text.x=element_text(angle=60, hjust=1), legend.position = "top")  + scale_x_date(date_breaks = "1 week")

# add poisson intervals
epi.lm <- lm(data = acc.date, diff.genome ~ time)
summary(epi.lm)

# add Poisson intervals
x <- seq(-30,120) # days
x.exp <- x * epi.lm$coefficients[2] + epi.lm$coefficients[1] # expected per site per day
x.exp2 <- x.exp * ref.length # expected per day
x.sd2 <- sqrt(x.exp2) # sd
x.sd <- x.sd2/ref.length # sd per site

x.err <- data.frame(day = x, date = as.Date('2020-01-10') + x , exp.site = x.exp, exp.day = x.exp2, sd.day = x.sd2, sd.site = x.sd, x.hi = x.exp + 2*x.sd, x.lo = x.exp - 2*x.sd )

scatterPlot <- p1 + geom_line(data = x.err, aes(x = date, y = x.hi), linetype = "dashed") + geom_line(data = x.err, aes(x = date, y = x.lo), linetype = "dashed") + geom_line(data = x.err, aes(x = date, y = exp.site)) + geom_vline(xintercept = as.Date(c('2020-01-23', '2020-02-02', '2020-03-13', '2020-03-11')))

# Marginal density plot of x (top panel)
xdensity <- ggplot(acc.date, aes(date2)) + geom_density(alpha=.5)
# +  scale_fill_manual(values = c('#999999','#E69F00')) + theme(legend.position = "none")

xdensity
# Marginal density plot of y (right panel)
ydensity <- ggplot(acc.date, aes(diff.genome)) + geom_density(alpha=.5) 

#+ scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
ydensity

# create blank slate
blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )

library("gridExtra")
grid.arrange(xdensity, blankPlot, scatterPlot, ydensity, ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))

###############
# binary string representation
##################
x <- read_tsv("st-index.tsv", col_names = c("id", "bin1", "bin2", "ind1", "ind2"))

# plot by decimal coords (not meaningful)
x %>% ggplot(aes(x=ind1, y=ind2)) + geom_jitter(size=3, shape=1) + theme_bw()

x.bin1 <- data.frame(bin1=sort(unique(x$bin1)), x1=1:72)
x.bin2 <- data.frame(bin2=sort(unique(x$bin2)), x2=1:57)
num.zero1 <- vector("integer", nrow(x.bin1)) 
for (i in seq_along(1:nrow(x.bin1))) {
  num.zero1[i] <- sum(str_detect(str_split(x.bin1$bin1,"")[[i]], "0"))
}
x.bin1 <- x.bin1 %>% mutate(num.zero1 = num.zero1)
num.zero2 <- vector("integer", nrow(x.bin2)) 
for (i in seq_along(1:nrow(x.bin2))) {
  num.zero2[i] <- sum(str_detect(str_split(x.bin2$bin2,"")[[i]], "0"))
}
x.bin2 <- x.bin2 %>% mutate(num.zero2 = num.zero2)


x1 <- x %>% left_join(x.bin1, "bin1")
x2 <- x1 %>% left_join(x.bin2, "bin2")

# Improvemnt 1. plot by num of zeros, equal-distant no matter position of the zeroes (better)
x2 %>% ggplot(aes(x=x1, y=x2))  + theme_minimal() + geom_vline(xintercept = 1:72, color="gray") + geom_hline(yintercept = 1:57, color = "gray")+ geom_jitter(size=3, color="red")

x.zero1 <- data.frame(num.zero1 = sort(unique(x2$num.zero1)), x3=1:8)
x.zero2 <- data.frame(num.zero2 = sort(unique(x2$num.zero2)), x4=1:8)

x2a <- x2 %>% arrange(num.zero1) %>% select(x1)
x2a <- data.frame(x1=unique(x2a$x1), x3=1:72)
x2b <- x2 %>% arrange(num.zero2) %>% select(x2)
x2b <- data.frame(x2=unique(x2b$x2), x4=1:57)

x3 <- x2 %>% left_join(x2a, "x1") %>% left_join(x2b, "x2")

# Improvement 2: sort by number of zeros as gradient
x3 %>% ggplot(aes(x=x3, y=x4, color = num.zero1+ num.zero2))  + theme_bw() + geom_point(size=5) + scale_color_gradient2(high="red", low="blue", mid="cyan", midpoint =6) + theme(legend.position = "bottom")

x3 %>% ggplot(aes(x=x3, y=x4, color = num.zero1+ num.zero2))  + theme_minimal() + geom_point(size=5) + scale_color_gradient(high="red", low="blue")

# Improvement 3: centerize ancestor
i <- 0
n <- nrow(x2a)
x2a.out = vector("list", nrow(x2a))
ct <- 1
while (i < n) {
  index <- ifelse(i %% 2 == 0, n-i/2, (i+1)/2)
  x2a.out[[ct]] <- data.frame(num=i, index=index)
  i <- i + 1
  ct <- ct + 1
}
x2a.df <- bind_rows(x2a.out)
x2a <- x2a %>% mutate(x5 = rev(x2a.df$index))

i <- 0
n <- nrow(x2b)
x2b.out = vector("list", nrow(x2b))
ct <- 1
while (i < n) {
  index <- ifelse(i %% 2 == 0, n-i/2, (i+1)/2)
  x2b.out[[ct]] <- data.frame(num=i, index=index)
  i <- i + 1
  ct <- ct + 1
}
x2b.df <- bind_rows(x2b.out)
x2b <- x2b %>% mutate(x6 = rev(x2b.df$index))

x4 <- x3 %>% left_join(x2a, "x1") %>% left_join(x2b, "x2")

p <- ggplot(data=x4, aes(x=x5, y=x6, color = num.zero1+ num.zero2))  + theme_bw() + geom_point(size=2) + scale_color_gradient2(high="red", low="blue", mid="cyan", midpoint =6) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Add rectangle grids
squares = data.frame(x=c(1,5,9,13,17,21,25,29,34),y=c(1,5,9,13,17,21,25,29,34))
x0 <- 36
y0 <- 29
p + geom_rect(data=squares, aes(NULL, NULL, xmin=x0-x, xmax=x0+x, ymin=y0-y, ymax=y0+y), color="gray", fill="gray", alpha = 0.1) 

# Add bands
df.out <- vector("list", 12)
for(i in 1:12) {
  df <- x4 %>% filter(num.zero1 + num.zero2==i) %>% select(x5,x6)
  df.out[[i]] <- data.frame(num.zero=i, xmin=range(df[[1]])[1], xmax=range(df[[1]])[2], ymin=range(df[[2]])[1], ymax=range(df[[2]])[2])
}

x.sqs <- bind_rows(df.out)
x0 <- 36
y0 <- 29
p + geom_rect(data=x.sqs, aes(NULL, NULL, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="gray", fill="gray", alpha = 0.1) 

# Iso-genetic lines
# Add path/polygon
# sort counter-clockwise using atan2
df.out <- vector("list", 11)
for(i in 1:11) {
  df <- x4 %>% filter(num.zero1 + num.zero2==i) %>% select(x5,x6)
  if (nrow(df) <3) { # need at least 3 points to make a polygon
    next
  }
  
  path <- vector("list", nrow(df)-1)
  for(j in 1:nrow(df)) { # sort clockwise
    df <- df %>% mutate(angle=atan2(x5-x0, x6-y0)) %>% arrange(angle)
  }  
  
#  df <- bind_rows(df, df[1,]) # to close the loop
  
#  for(j in 1:(nrow(df)-1)) {
#    path[[j]] <- data.frame(num.zero=i, x1=df$x5[j], y1=df$x6[j], x2=df$x5[j+1], y2=df$x6[j+1])
#  }
#  df.out[[i]] <- bind_rows(path)
  df.out[[i]] <- data.frame(id=i, x=df$x5, y= df$x6)
}
x.path <- bind_rows(df.out)
x.path <- x.path %>% mutate(id = paste("mut", id, sep=""))
shades <- data.frame(id = factor(paste("mut", 1:10, sep="")),
                      col = c("black", NA, "black", NA, "black", NA, "black", NA, "black", NA))
datapoly <- x.path %>% left_join(shades, "id")
root <- data.frame(x=x0, y=y0)
#p + geom_path(data=x.path, aes(x=x1, y=y1, xend=x2, yend=y2, color=num.zero)) + geom_point(data=root, aes(x=x, y=y), shape = 10, color="red", size=5)

x4 <- x4 %>% mutate(id = paste("mut", num.zero1 + num.zero2, sep=""))
x5 <- x4 %>% filter(id != 'mut0' & id != 'mut11' & id != 'mut12')

datapoly$id = factor(datapoly$id, levels = c(paste("mut", 1:10, sep = "")))
x5$id = factor(x5$id, levels = c(paste("mut", 1:10, sep = "")))

#ggplot(data=datapoly, aes(x=x, y=y)) + geom_polygon( aes(group=id, color=id), fill=NA) + theme_bw() + geom_point(data=root, aes(x=x, y=y), shape = 10, color="red", size=5) + theme_bw() + geom_point(data=x5, aes(x=x5, y=x6), size=2, color="magenta", shape=1)  

+ facet_wrap(~id, nrow = 5) + xlab("front seq") + ylab("back seq")

ggplot(data=datapoly, aes(x=x, y=y)) + geom_polygon( aes(group=id, color=id), fill=NA) + theme_bw() + geom_point(data=root, aes(x=x, y=y), shape = 10, color="red", size=5) + theme_bw() + geom_point(data=x5, aes(x=x5, y=x6), size=2, color="blue") + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# count major lineages
x4 %>% ggplot(aes(x=x5)) + geom_bar() # same front end
x4 %>% ggplot(aes(x=x6)) + geom_bar() # same backend
out <- vector("list", 72 * 57)
k <- 1
for (i in 1:72) {
  for (j in 1:57) {
    start = c(i,j)
    ct.up <- x4 %>% filter(x3.x == i & x4.x >= j) %>% count() %>% c()
    ct.right <- x4 %>% filter(x4.x == j & x3.x >= i) %>% count() %>% c()
    out[[k]] <- data.frame(front = i, back=j, ct=ct.up$n + ct.right$n)
    k <- k + 1
  }
}

out.df <- bind_rows(out)
out.df %>% ggplot(aes(x=ct)) + geom_histogram(binwidth = 3)
mel <- out.df %>% filter(ct > 20)
p <- ggplot(data=x4, aes(x=x3.x, y=x4.x, color = num.zero1+num.zero2))  + theme_bw() + geom_point(aes(size = num.zero1+num.zero2), shape=1) + scale_color_gradient2(high="red", low="blue", mid="cyan", midpoint =6) + theme(legend.position = "bottom") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Add bands
df.out <- vector("list", 12)
for(i in 1:12) {
  df <- x4 %>% filter(num.zero1 + num.zero2==i) %>% select(x3.x,x4.x)
  df.out[[i]] <- data.frame(num.zero=i, xmin=1, xmax=range(df[[1]])[2], ymin=1, ymax=range(df[[2]])[2])
}

x.sqs <- bind_rows(df.out)
p <- ggplot(data=x4, aes(x=x3.x, y=x4.x))  + theme_bw()  + geom_text(aes(label = num.zero1+ num.zero2)) + scale_x_continuous(breaks=x.lab$br, labels = x.lab$lab) + scale_y_continuous(breaks=y.lab$br, labels = y.lab$lab) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=90)) 

#+ geom_rect(data=x.sqs, aes(NULL, NULL, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="gray", fill=NA) + xlab("") + ylab("")

count.zero <- function(binStr){
  sum(str_detect(str_split(binStr,"")[[1]], "0"))
}

x.bins <- vector("list", 72)
for (i in 1:72) {
  df<- x4 %>% filter(x3.x == i) %>% select(bin1) %>% unique()
  x.bins[[i]] <- data.frame(br = i, lab = df$bin1, num.zero=count.zero(df$bin1))
}
x.lab <- bind_rows(x.bins)

y.bins <- vector("list", 57)
for (i in 1:57) {
  df<- x4 %>% filter(x4.x == i) %>% select(bin2) %>% unique()
  y.bins[[i]] <- data.frame(br = i, lab = df$bin2, num.zero=count.zero(df$bin2))
}
y.lab <- bind_rows(y.bins)

x.grid <- vector("integer", 8)
for(i in 1:8) {
  df <- x.lab %>% filter(num.zero == i-1) %>% summarise(max.br=max(br))
  x.grid[i] <- df$max.br
}

y.grid <- vector("integer", 8)
for(i in 1:8) {
  df <- y.lab %>% filter(num.zero == i-1) %>% summarise(max.br=max(br))
  y.grid[i] <- df$max.br
}

p + geom_vline(xintercept = x.grid, color = "gray") + geom_hline(yintercept = y.grid, color = "gray") + xlab("") + ylab("")

#######################
x <- read_tsv("bin.tsv", col_names = c('id', "snp1", "snp2", "snp3", "snp4", "snp5", "snp6", "snp7", "snp8"))
x2 <- read_tsv("bin.tsv2", col_names = c("id", "seq"))
x2 <- x2 %>% arrange(seq)

x.long$id_f = factor(x.long$id, levels = x2$id)

x.long <- x %>% gather(2:9, key = "snp", value = "pos")
x.long %>% ggplot(aes(x=snp, y=pos, group=id)) + geom_point() + geom_line() + theme_bw() 

#############
# clustering by t-sne
#################
x <- read_tsv("maf.pair-diff", col_names = c("a", "b", "len", "diff", "len2", "perc", "diff", "diff2"))
x.dist <- as.dist(xtabs(x$diff2 ~ x$b + x$a))
x.tsne <- tsne(x.dist, perplexity = 5, max_iter = 1000)
plot(x.tsne)

x.cmd <-cmdscale(x.dist)
plot(x.cmd)

###################
# 5/30/2020
# March genomes
######################
setwd("C:/Users/lai/Dropbox/cov-browser/data-05-12-2020")
x <- read_tsv("march.join2", col_names = c("acc", "date", "country.id", "area.id", "diff", "diff2", "diff3"))
x <- x %>% mutate(date2 = as.Date(date))
x %>% filter(diff<50) %>% ggplot(aes(x=date2, y=diff, color=as.character(area.id))) + geom_jitter(shape=1)
x %>% filter(diff < 50) %>% ggplot(aes(x=diff)) + geom_histogram(bins = 100)



###################
# 5/25/2020
# new genomes, N=17K
##############
setwd("C:/Users/lai/Dropbox/cov-browser/data-05-12-2020")
x <- read_tsv("snp-comp.txt", col_names = F)
names(x) <- c("pos", "acc", "num.allele", "sum", "ct1", "ct2", "freq1", "freq2", "share.type", "ct3", "mut", "mut.type")
x %>% group_by(mut.type) %>% count()
x %>% ggplot(aes(x=share.type, y= freq2, color=share.type, shape=mut.type)) + geom_jitter(size=5) + theme_bw() + scale_y_log10() + geom_hline(yintercept = 1, linetype="dashed") + ylab("alt. allele freq") + xlab("overlap with N=2345 samples")
x %>% group_by(type) %>% count()
x.ct <- x %>% group_by(type) %>% count()



######################
# analysis of origin
# 5/13/2020
######################

ref.length <- 29903;
epi <- read_csv("../network-all-sites/epi-db.csv", col_names = c("acc", "country", "geo", "date")) %>% filter(acc != 'EPI_ISL_402131')

st <- read_tsv("../network-maf2/no-impute.log", skip = 1, col_names = F)
names(st) <- c("st", "acc")

pid.st1 <- read_tsv("../network-all-sites/EPI_ISL_402123.pid", col_names = c("ref", "acc", "ct", "valid.len", "len", "pid", "diff"))
pid.st1 <- pid.st1 %>% mutate(seq.diff = ct/ref.length)

pid.st3 <- read_tsv("../network-maf2/EPI_ISL_406030.pid", col_names = c("ref", "acc", "ct", "valid.len", "len", "pid", "diff"))
pid.st3 <- pid.st3 %>% mutate(seq.diff = ct/ref.length)

pid.st4 <- read_tsv("../network-maf2/EPI_ISL_404895.pid", col_names = c("ref", "acc", "ct", "valid.len", "len", "pid", "diff"))
pid.st4 <- pid.st4 %>% mutate(seq.diff = ct/ref.length)

pid.st18 <- read_tsv("../network-maf2/EPI_ISL_413456.pid", col_names = c("ref", "acc", "ct", "valid.len", "len", "pid", "diff"))
pid.st18 <- pid.st4 %>% mutate(seq.diff = ct/ref.length)

st %>% filter(st=='ST3') %>% left_join(epi, "acc") %>% arrange(date)

st %>% filter(st=='ST18') %>% left_join(epi, "acc") %>% arrange(date)

#acc.early %>% filter(st %in% c("ST3", "ST9", "ST4", "ST41", "ST18")) %>% ggplot(aes(x=date, fill=st)) + geom_bar() + facet_wrap(~country) + theme_bw()

# ST3 vs ST4
st.key <- c("ST3", "ST9", "ST4", "ST41", "ST18", "ST7")
acc.key <- st %>% filter(st %in% st.key) %>% left_join(epi, "acc")

acc.st3 <- acc.key %>% left_join(pid.st3, "acc")
acc.st3 %>% ggplot(aes(x=date, y=seq.diff)) + geom_jitter(alpha=0.5, size=2, aes(color=st)) + theme_bw() + geom_smooth(method = "lm")

acc.st4 <- acc.key %>% left_join(pid.st4, "acc")
acc.st4 %>% ggplot(aes(x=date, y=seq.diff)) + geom_jitter(alpha=0.5, size=2, aes( color=st)) + theme_bw()  + geom_smooth(method = "lm")

df3 <- acc.st3 %>% select(st, date, seq.diff) %>% mutate(root="st3") %>% filter(!is.na(seq.diff))
df4 <- acc.st4 %>% select(st, date, seq.diff) %>% mutate(root="st4") %>% filter(!is.na(seq.diff))

df.key <- bind_rows(df3, df4)
df.key %>% ggplot(aes(x=date, y=seq.diff)) + geom_jitter(shape=1, size=1, aes(color=st)) + facet_wrap(~root) +  theme_bw() + geom_smooth(method = "lm") + xlab("collection date")

df3 <- acc.st3 %>% select(st, date, seq.diff) %>% mutate(root="st3") 
df4 <- acc.st4 %>% select(st, date, seq.diff) %>% mutate(root="st4")
lm.st3 <- lm(data = df3, seq.diff ~ date)
lm.st4 <- lm(data = df4, seq.diff ~ date)
anova(lm.st3, lm.st4)

# ST1 vs ST18
st.early <- c("ST1", "ST11", "ST12", "ST7", "ST3", "ST9", "ST6", "ST19", "ST33", "ST4", "ST41", "ST18", "ST50")

acc.early <- st %>% filter(st %in% st.early) %>% left_join(epi, "acc")


acc.st1 <- acc.early %>% left_join(pid.st1, "acc")
acc.st1 %>% ggplot(aes(x=date, y=seq.diff)) + geom_jitter(alpha=0.5, size=2, aes(color=st)) + theme_bw()  + geom_smooth(method = "lm")

acc.st18 <- acc.early %>% left_join(pid.st18, "acc")
acc.st18 %>% ggplot(aes(x=date, y=seq.diff)) + geom_jitter(alpha=0.5, size=2, aes(color=st)) + theme_bw() + geom_smooth(method = "lm")

df1 <- acc.st1 %>% select(st, date, seq.diff) %>% mutate(root="st1") %>% filter(!is.na(seq.diff))
df18 <- acc.st18 %>% select(st, date, seq.diff) %>% mutate(root="st18") %>% filter(!is.na(seq.diff))

df.early <- bind_rows(df1, df18)
df.early %>% ggplot(aes(x=date, y=seq.diff)) + geom_jitter(shape=1, size=1, aes(color=st)) + facet_wrap(~root) +  theme_bw() + geom_smooth(method = "lm") + xlab("collection date")
lm.st1 <- lm(data = df1, seq.diff ~ date)
lm.st18 <- lm(data = df18, seq.diff ~ date)
anova(lm.st1, lm.st18)
# combined plot

df <- bind_rows(df1, df18, df3, df4)
df %>% ggplot(aes(x=date, y=seq.diff)) + geom_jitter(shape=1, size=1, aes(color=st)) + facet_wrap(~root) +  theme_bw() + geom_smooth(method = "lm") + xlab("collection date")

#######################
# re-estimate rate using ST3
###################
epi2 <- epi %>% left_join(pid.st3, "acc") %>% mutate(date.ref = as.Date("2020-01-10"))
epi2 <- epi2 %>% mutate(time = date - date.ref)
geo <- read_tsv("../network-all-sites/epi.geo", col_names = c("acc", "area"))
epi2 <- epi2 %>% left_join(geo, "acc") %>% filter(!is.na(seq.diff))

p1 <- ggplot(data = epi2, aes(x = date, y = seq.diff)) + geom_jitter(aes(color = area), shape = 1, size = 1) + theme_bw() + xlab("collection date") + ylab("seq.diff (per site)") + scale_color_manual(breaks = c("Asia", "Africa", "Europe", "South America", "North America", "Oceania", "Cruise Ship"), values = c("#0099ff", "#ff66ff", "#33cc33", "#ff9933", "#ff0000", "#cc0099", "#00ffff")) +  theme(axis.text.x=element_text(angle=60, hjust=1))  + scale_x_date(date_breaks = "1 week")

# add poisson intervals
epi.lm <- lm(data = epi2, seq.diff ~ time) 
summary(epi.lm)

# add Poisson intervals
x <- seq(-18,80) # days
x.exp <- x * epi.lm$coefficients[2] + epi.lm$coefficients[1] # expected per site per day
x.exp2 <- x.exp * ref.length # expected per day
x.sd2 <- sqrt(x.exp2) # sd
x.sd <- x.sd2/ref.length # sd per site

x.err <- data.frame(day = x, date = as.Date('2020-01-10') + x , exp.site = x.exp, exp.day = x.exp2, sd.day = x.sd2, sd.site = x.sd, x.hi = x.exp + 2*x.sd, x.lo = x.exp - 2*x.sd )

p1 + geom_line(data = x.err, aes(x = date, y = x.hi), linetype = "dashed") + geom_line(data = x.err, aes(x = date, y = x.lo), linetype = "dashed") + geom_line(data = x.err, aes(x = date, y = exp.site)) 
+ geom_vline(xintercept = as.Date(c('2020-01-23', '2020-02-02', '2020-03-13')))




###################
# 5/10/2020 LD plot
#####################

x <- read_tsv("hapview-maf/maf.tsv")
names(x) <- c("s1", "s2", "D.prime", "lod", "r2", "lo", "hi", "dist", "tint")

x.sig <- x %>% filter(lod>=2)
# linkage decay over distance
#x.rec <- x.sig %>% filter(D.prime<1)
p1 <- x.sig %>% ggplot(aes(x=dist, y=r2)) + geom_point(shape=1) + theme_bw()

# LD decay: using moving average
breaks <- seq(100,29500, by=100)
out <- vector("list", length(breaks))
for(i in seq_along(breaks)) {
  x.br <- x.dist %>% filter(dist < breaks[i])
  out[[i]] <- data.frame(dist = breaks[i],
                         mean.r2 = mean(x.br$r2, na.rm = T),
                        sd.r2 = sd(x.br$r2, na.rm = T))
}
ma.r2 <- bind_rows(out)
p1 + geom_line(data = ma.r2, aes(x = dist, y=mean.r2), color="red") + geom_line(data = ma.r2, aes(x = dist, y=mean.r2+sd.r2), color="red", linetype="dashed") + geom_line(data = ma.r2, aes(x = dist, y=mean.r2-sd.r2), color="red", linetype="dashed")

# plot pairwise R2: x.sig only
site <- unique(c(x.sig$s1, x.sig$s2))
df.site <- data.frame(pos = site, site = as.character(site))
x.sig <- x.sig %>% mutate(pos1 = as.character(s1), pos2 = as.character(s2))
x.sig$pos1 <- factor(x.sig$pos1, levels = df.site$site)
x.sig$pos2 <- factor(x.sig$pos2, levels = df.site$site)


x.sig %>% ggplot(aes(x = pos2, y=pos1, fill=r2)) + geom_tile(colour="gray") + scale_fill_gradient(low = "white", high = "black") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_blank())

###############
# 5/5/2020
# Recombination rates (Saymon's code)
###############

rec1 <- read_tsv("../rec-LDjump/rec_regression_of_genome_window_ordered_by_time.txt")
rec1 <- rec1 %>% mutate(ratio = roh/theta )

s1=ggscatter(
  rec1, x = "date", y = "theta",
  color = "Geography", palette = "jco",
  add = "reg.line",
  conf.int = TRUE,
  size = 2.5
) + facet_wrap(~Geography, nrow = 3) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), color=Geography),
    label.y = 0.00039, size = 5, show.legend = F
  ) + xlab("") + ylab("Genome diversity (θ)") +
  scale_color_manual( values=c('#1A06EA', '#59AC06', '#EA1006'))+
  scale_fill_manual(values = c('#6675F5', '#3ADA07', '#F95E5D')) +
  theme(axis.text.y=element_text(size = 10, face = "bold"),
        axis.text.x=element_text(size = 8, face = "bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=20, face = "bold"), legend.text=element_text(size=16, face = "bold"),
        strip.text.x = element_text(size = 14, face = "bold"))


s2=ggscatter(
  rec1, x = "date", y = "roh",
  color = "Geography", palette = "jco",
  add = "reg.line",
  conf.int = TRUE,
  size = 2.5
) + facet_wrap(~Geography, nrow = 3) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), color=Geography),
    label.y = 0.073, size = 5, show.legend = F
  ) + xlab("") + ylab("Recombination rate (ρ)") +
  scale_color_manual( values=c('#1A06EA', '#59AC06', '#EA1006'))+
  scale_fill_manual(values = c('#6675F5', '#3ADA07', '#F95E5D')) +
  theme(axis.text.y=element_text(size = 10, face = "bold"),
        axis.text.x=element_text(size = 8, face = "bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=20, face = "bold"), legend.text=element_text(size=16, face = "bold"),
        strip.text.x = element_text(size = 14, face = "bold"))

s3=ggscatter(
  rec1, x = "date", y = "ratio",
  color = "Geography", palette = "jco",
  add = "reg.line",
  conf.int = TRUE,
  size = 2.5
) + facet_wrap(~Geography, nrow = 3) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), color=Geography),
    label.y = 450, size = 5, show.legend = F
  ) + xlab("") + ylab("Ratio (ρ/θ)") +
  scale_color_manual( values=c('#1A06EA', '#59AC06', '#EA1006'))+
  scale_fill_manual(values = c('#6675F5', '#3ADA07', '#F95E5D')) +
  theme(axis.text.y=element_text(size = 10, face = "bold"),
        axis.text.x=element_text(size = 8, face = "bold"),
        axis.title=element_text(size=15,face="bold"),
        legend.title=element_text(size=20, face = "bold"), legend.text=element_text(size=16, face = "bold"),
        strip.text.x = element_text(size = 14, face = "bold"))



Figure1 = ggarrange(s1, s2, s3,
                    labels = c("A", "B", "C"),
                    ncol = 3, legend = "none")

annotate_figure(Figure1, bottom = text_grob("Collection Date", face = "bold", size = 17))



###############
# 5/4/2020
# Evolution rate
###############
setwd("C:/Users/lai/Dropbox/cov-recom/scripts")
library(tidyverse)
ref.length <- 29903;


epi <- read_csv("../network-all-sites/epi-db.csv", col_names = c("acc", "country", "geo", "date")) %>% filter(acc != 'EPI_ISL_402131')
pid <- read_tsv("../network-all-sites/EPI_ISL_402123.pid", col_names = c("ref", "acc", "ct", "valid.len", "len", "pid", "diff"))
pid <- pid %>% mutate(seq.diff = ct/ref.length)
geo <- read_tsv("../network-all-sites/epi.geo", col_names = c("acc", "area"))

epi2 <- epi %>% left_join(pid, "acc") %>% mutate(date.ref = as.Date("2019-12-24"))
epi2 <- epi2 %>% mutate(time = date - date.ref)
epi2 <- epi2 %>% left_join(geo, "acc") %>% filter(!is.na(seq.diff))

p1 <- ggplot(data = epi2, aes(x = date, y = seq.diff)) + geom_jitter(aes(color = area), shape=1, size = 1) + theme_bw() + xlab("collection date") + ylab("seq.diff (per site)") + scale_color_manual(breaks = c("Asia", "Africa", "Europe", "South America", "North America", "Oceania", "Cruise Ship"), values = c("#0099ff", "#ff66ff", "#33cc33", "#ff9933", "#ff0000", "#cc0099", "#00ffff")) +  theme(axis.text.x=element_text(angle=60, hjust=1))  + scale_x_date(date_breaks = "1 week")

# + expand_limits(x= as.Date('2019-12-15'))
#+ geom_smooth(method = "lm") # doesn't extend

# add poisson intervals
#epi.lm <- lm(data = epi2, seq.diff ~ time + 0) # force intercept = 0
summary(epi.lm)

epi.lm <- lm(data = epi2, seq.diff ~ time) 
summary(epi.lm)

#p2 <- p1 + geom_hline(yintercept = epi.lm2$coefficients[1], linetype="dashed")

# Build prediction data frame
# pred_x = seq(-5,max(epi2$time))
# pred_lines = data.frame(x=pred_x,
                        date = pred_x + as.Date("2019-12-24"),
                        y=predict(epi.lm2, data.frame(time=pred_x)),
                        obs_Or_Pred=if_else(pred_x >=2, "Obs","Pred"))

#p3 <- p2 + geom_line(data = pred_lines, aes(x = date, y = y))

# add Poisson intervals
x <- seq(0,100) # days
x.exp <- x * epi.lm$coefficients[2] + epi.lm$coefficients[2] # expected per site per day
x.exp2 <- x.exp * ref.length # expected per day
x.sd2 <- sqrt(x.exp2) # sd
x.sd <- x.sd2/ref.length # sd per site

x.err <- data.frame(day = x, date = as.Date('2019-12-24') + x , exp.site = x.exp, exp.day = x.exp2, sd.day = x.sd2, sd.site = x.sd, x.hi = x.exp + 2*x.sd, x.lo = x.exp - 2*x.sd )

p1 + geom_line(data = x.err, aes(x = date, y = x.hi), linetype = "dashed") + geom_line(data = x.err, aes(x = date, y = x.lo), linetype = "dashed") + geom_line(data = x.err, aes(x = date, y = exp.site)) 
+ geom_vline(xintercept = as.Date(c('2020-01-23', '2020-02-02', '2020-03-13')))



##############
# 5/7/2020
# simulate origin
################

days <- 0:100
diff <- days * 1
days.worng <- days - 9
diff.wrong <- diff - diff[10]
df <- data.frame(days.real = days, diff.real = diff, days.delay = days.wrong, diff.delay = diff.wrong)
ggplot(data =df , aes(x = days.real, y = diff.real)) + geom_point() + geom_point(data = df, aes(x = days.delay, y = diff.real), color = "red")  

## 14th of May 2021

## test of specificity and sensitivity of each test and dataset

##############
# 1. datasets#
##############

# Sage et al. 1986
# nematodes
Nwormy = 15; Nnormal = 78; Nhybrid = 46; Npure = 21+26; Ntot = Nnormal+Nwormy
NwormyHybrid = 14; Nwormypure = 1; 
NnormalHybrid = Nhybrid - NwormyHybrid; NnormalPure = Npure -Nwormypure
SAGdata <- data.frame(wormy = c(hybrid=NwormyHybrid, pure=Nwormypure),
                          normal=c(hybrid=NnormalHybrid, pure=NnormalPure)) # total 78
sum(SAGdata) # N=93

# Moulia et al. 1993
MOUdata <- read.csv("data4fieldStudies/moulia1993.csv") # N=120
mean(MOUdata$load[MOUdata$genotype %in% "Mmm"])
mean(MOUdata$load[MOUdata$genotype %in% "Hybrid"]) # wrong in their table 2
mean(MOUdata$load[MOUdata$genotype %in% "Mmd"])

K <- kruskal.test(load ~ genotype, data = MOUdata)
K # ok

# Baird et al. 2012
BAIdata <- read.csv("https://raw.githubusercontent.com/alicebalard/parasiteLoad/master/data/WATWMdata.csv", na.strings = c("", " ", NA))
# Keep individuals with hybrid index, sex recorded and pinworms counted
BAIdata <- BAIdata[!is.na(BAIdata$HI) & !is.na(BAIdata$Sex) & !is.na(BAIdata$Aspiculuris.Syphacia),]
table(BAIdata$Aspiculuris.Syphacia >0)
length(BAIdata$Aspiculuris.Syphacia) # N=663
sum(BAIdata$Aspiculuris.Syphacia >0) / length(BAIdata$Aspiculuris.Syphacia) # N=585
# 71%

# Balard et al. 2020
BALdata <- read.csv("https://raw.githubusercontent.com/alicebalard/Article_IntensityEimeriaHMHZ/master/data/cleanedData.csv", na.strings = c("", " ", NA))
# Keep individuals with hybrid index, sex recorded and pinworms counted
BALdata <- BALdata[!is.na(BALdata$HI) & !is.na(BALdata$Sex) & !is.na(BALdata$Aspiculuris_Syphacia),]
length(BALdata$Aspiculuris_Syphacia) # N=585
# prevalence:
sum(BALdata$Aspiculuris_Syphacia >0) / length(BALdata$Aspiculuris_Syphacia) # N=585
# 52.4%

#################
# simulate data #
#################

N=100
myhybeffectsize = 0

MeanLoad <- function(L1, L2, hybeff, hybridIndex){
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  mean <- (L1 + (L2 - L1) * hybridIndex) * (1 - hybeff * heterozygoty)
  mean <- sapply(mean, function(x) {
    return(max(x, 0.01))
  })
  return(mean)
}

hybind = runif(n = N, min = 0, max = 1)
simdata = data.frame(HI =hybind, Sex = factor(rep(c("female", "male"), N/2)),
                     OPG= round(MeanLoad(
                       L1 = 50, L2 = 50, hybeff = myhybeffectsize,
                       hybridIndex = hybind)))


rnbinom()
## Functions


SizeNegBin <- function(A1, A2, Z, hybridIndex){
  aggregation <- Aggregation(A1, A2, Z, hybridIndex)
  aggregation <- sapply(aggregation, function(x) {
    return(max(x, 0.01))
  })
  size <- 1/aggregation
  return(size)
}

Aggregation <- function(A1, A2, Z, hybridIndex){
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  aggregation <- (A1 + (A2 - A1) * hybridIndex) + Z * heterozygoty
  return(aggregation)
}


Dlow 



######################################
# 2. calculate specificity (1-alpha) #
######################################
# Specificity=detecting H1


################################
# 3. calculate power (1 - beta)#
################################





##### previous
##################
# Sage et al. 1986
##################
Nwormy = 15; Nnormal = 78; Nhybrid = 46; Npure = 21+26; Ntot = Nnormal+Nwormy
NwormyHybrid = 14; Nwormypure = 1; 
NnormalHybrid = Nhybrid - NwormyHybrid; NnormalPure = Npure -Nwormypure

nematodesDF <- data.frame(wormy = c(hybrid=NwormyHybrid, pure=Nwormypure),
                          normal=c(hybrid=NnormalHybrid, pure=NnormalPure)) # total 78
chinem <- chisq.test(nematodesDF, correct = F) # they found 9.6, df=1, p <0.005
chinem

# Phi effect size chi-square
phi = sqrt(chinem$statistic/Ntot)
pwr.chisq.test(w = phi, N = Ntot, df = chinem$parameter, sig.level = 0.05)

# Chi squared power calculation 
# 
# w = 0.3847961
# N = 93
# df = 1
# sig.level = 0.05
# power = 0.9600163
# 
# NOTE: N is the number of observations

## cestodes
Nwormy = 9; Nnormal = 84; Nhybrid = 46; Npure = 21+26; Ntot = Nnormal+Nwormy
NwormyHybrid = 8; Nwormypure = 1; 
NnormalHybrid = Nhybrid - NwormyHybrid; NnormalPure = Npure -Nwormypure

cestodesDF <- data.frame(wormy = c(hybrid=NwormyHybrid, pure=Nwormypure),
                         normal=c(hybrid=NnormalHybrid, pure=NnormalPure)) # total 78
chices <- chisq.test(cestodesDF, correct = F) # they found 4.1, df=1, p <0.05
chices

# Phi effect size chi-square
phi = sqrt(chices$statistic/Ntot)
pwr.chisq.test(w = phi, N = Ntot, df = chices$parameter, sig.level = 0.05)

# Chi squared power calculation
# 
# w = 0.2581221
# N = 93
# df = 1
# sig.level = 0.05
# power = 0.7016971
# 
# NOTE: N is the number of observations

####################
# Moulia et al. 1993
####################
setwd("/home/alice/Desktop/Git/Article3_review/data4fieldStudies")
mouliadf <- read.csv("moulia1993.csv")

mean(mouliadf$load[mouliadf$genotype %in% "Mmm"])
mean(mouliadf$load[mouliadf$genotype %in% "Hybrid"]) # wrong in their table 2
mean(mouliadf$load[mouliadf$genotype %in% "Mmd"])

K <- kruskal.test(load ~ genotype, data = mouliadf)
K # ok
pairwise.wilcox.test(mouliadf$load,
                     mouliadf$genotype)

dunn.test(x=mouliadf$load, g = mouliadf$genotype)

# Kruskal-Wallis rank sum test
# 
# data: x and group
# Kruskal-Wallis chi-squared = 13.8175, df = 2, p-value = 0
# 
#           Comparison of x by group                            
#                (No adjustment)                                
#   Col Mean-|
#   Row Mean |     Hybrid        Mmd
#   ---------+----------------------
#        Mmd |   2.698917 
#            |    0.0035*
#            |
#        Mmm |   3.609349   0.642956
#            |    0.0002*     0.2601
# 
# alpha = 0.05
# Reject Ho if p <= alpha/2

# kwpower approximates power for the Kruskal-Wallis test, 
#  using a chi-square approximation under the null, and a non-central chi-square
#  approximation under the alternative. The noncentrality parameter is calculated 
#  using alternative means and the null variance structure.

kwpower(rep(10,3),c(0,1,2),"normal")


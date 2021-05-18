# devtools::install_github("alicebalard/parasiteLoad", ref = "v2.1")
library(parasiteLoad)
library(optimx)
library(fitdistrplus)

# for reproducibility
set.seed(8)

######################################
### Load data from the 4 field studies
### For Sage et al. and Moulia et al.: extract negbin parameters for simulation of bigger study

# Sage et al. 1986
# no easily interpretable individual data, so I took farm level data
SAGdata <- read.csv("../data4fieldStudies/Sage1986.csv")
SAGdata$HI <- (SAGdata$Hioriginal +8) /16
SAGdata$load <- SAGdata$mean_abundance

descdist(SAGdata$load)
hist(SAGdata$load)
fitSAG <- fitdist(data = round(SAGdata$load), distr = "nbinom")
fitSAG

# Moulia et al. 1993
MOUdata <- read.csv("../data4fieldStudies/moulia1993.csv") # N=120
MOUdata$HI <- 0.5
MOUdata$HI[MOUdata$genotype == "Mmm"] <- 1
MOUdata$HI[MOUdata$genotype == "Mmd"] <- 0

hist(MOUdata$load)
descdist(MOUdata$load)
fitMOU <- fitdist(data = round(MOUdata$load), distr = "nbinom")
fitMOU

# Baird et al. 2  012
BAIdata <- read.csv("https://raw.githubusercontent.com/alicebalard/parasiteLoad/master/data/WATWMdata.csv", na.strings = c("", " ", NA))
# Keep individuals with hybrid index, sex recorded and pinworms counted
BAIdata <- BAIdata[!is.na(BAIdata$HI) & !is.na(BAIdata$Sex) & !is.na(BAIdata$Aspiculuris.Syphacia),]

# Balard et al. 2020
BALdata <- read.csv("https://raw.githubusercontent.com/alicebalard/Article_IntensityEimeriaHMHZ/master/data/cleanedData.csv", na.strings = c("", " ", NA))
# Keep individuals with hybrid index, sex recorded and pinworms counted
BALdata <- BALdata[!is.na(BALdata$HI) & !is.na(BALdata$Sex) & !is.na(BALdata$Aspiculuris_Syphacia),]

######################################
### Simulate small & big datasets for each study 
### Shuffle the the nematode load to make it random along HI

### -> to do from there: BS

make_data_shuffle <- function(){
  ## Sage et al: simulate N=100, simulate N=600
  SAG100_df <- data.frame(HI = sample(seq(0,1, 0.001), size = 100),
                          load = rnegbin(n = 100, mu = fitSAG$estimate[["mu"]], theta = fitSAG$estimate[["size"]]),
                          Sex = as.factor(sample(rep(c("F", "M"), 1000), size = 100)))
  SAG600_df <- data.frame(HI = sample(seq(0,1, 0.001), size = 600),
                          load = rnegbin(n = 600, mu = fitSAG$estimate[["mu"]], theta = fitSAG$estimate[["size"]]),
                          Sex = as.factor(sample(rep(c("F", "M"), 1000), size = 600)))
  ## Moulia et al: simulate N=100 & N=600
  MOU100_df <- data.frame(HI = sample(seq(0,1, 0.001), size = 100),
                          load = rnegbin(n = 100, mu = fitMOU$estimate[["mu"]], theta = fitMOU$estimate[["size"]]),
                          Sex = as.factor(sample(rep(c("F", "M"), 1000), size = 100)))
  MOU600_df <- data.frame(HI = sample(seq(0,1, 0.001), size = 600),
                          load = rnegbin(n = 600, mu = fitMOU$estimate[["mu"]], theta = fitMOU$estimate[["size"]]),
                          Sex = as.factor(sample(rep(c("F", "M"), 1000), size = 600)))
  ## Baird et al: subsample N=100 & N=600
  BAI100_df <- data.frame(HI = sample(BAIdata$HI, size = 100), 
                          load = sample(BAIdata$Aspiculuris.Syphacia, size = 100),
                          Sex = as.factor(sample(factor(BAIdata$Sex), size = 100)))
  BAI600_df <- data.frame(HI = sample(BAIdata$HI, size = 600), 
                          load = sample(BAIdata$Aspiculuris.Syphacia, size = 600),
                          Sex = as.factor(sample(factor(BAIdata$Sex), size = 600)))
  ## Balard et al: subsample N=100 & N=600 (585 for BAL)
  BAL100_df <- data.frame(HI = sample(BALdata$HI, size = 100), 
                          load = sample(BALdata$Aspiculuris_Syphacia, size = 100),
                          Sex = as.factor(sample(factor(BALdata$Sex), size = 100)))
  BAL585_df <- data.frame(HI = sample(BALdata$HI, size = 585), 
                          load = sample(BALdata$Aspiculuris_Syphacia, size = 585),
                          Sex = as.factor(sample(factor(BALdata$Sex), size = 585)))
  return(list(SAG100_df = SAG100_df, SAG600_df = SAG600_df, 
              MOU100_df = MOU100_df, MOU600_df = MOU600_df, 
              BAI100_df = BAI100_df, BAI600_df = BAI600_df, 
              BAL100_df = BAL100_df, BAL585_df = BAL585_df))
}

######################################
# Apply the three test on each df:
runTests <- function(df){
  ## Test 1: chi2 wormy VS non wormy, cut at 250, hybrids 12.5% < HI < 87.5%
  chi2df <- data.frame(wormy = c(hybrid=length(df[df$load>250 & df$HI > 0.125 & df$HI < 0.875,]),
                                 pure=length(df[df$load>250 & (df$HI < 0.125 | df$HI > 0.875),])),
                       normal=c(hybrid=length(df[df$load<=250 & df$HI > 0.125 & df$HI < 0.875,]), 
                                pure=length(df[df$load<=250 & (df$HI < 0.125 | df$HI > 0.875),]))) # total 78
  chi2 <- chisq.test(chi2df) # they found 9.6, df=1, p <0.005
  
  ## Test 2: kruskal-Wallis test on load, 20% < HI < 60%
  df$genotype <- "Hybrid"
  df$genotype[df$HI < 0.2] <- "Mmd"
  df$genotype[df$HI > 0.6] <- "Mmm"
  KW <- kruskal.test(load ~ genotype, data = df)
  
  ## Test 3: maximum likelihood estimation and LRT, H0 (no diff in sex and both sides)
  ML <- parasiteLoad::analyse(data = df, response = "load", model = "negbin",
                              group = "Sex", hybridIndex = "HI")
  
  results = list(chi2 = chi2$p.value, KW = KW$p.value, ML = ML$H0$Gtest$pvalue)
  return(results)
}

mybootstrap <- function(numBS){
  # shuffled_data <- lapply(1:numBS, function(i) make_data_shuffle())
  ## BS function
  lapply(1:numBS, function(i) {
    d <- make_data_shuffle()
    lapply(d, runTests)
  })
}

numBS=1000

system.time(
try <- mybootstrap(numBS)
)
# user  system elapsed 
# 260.088   0.178 260.212 

vecSAG100chi2 <- vector(); vecSAG100KW <- vector(); vecSAG100ML <- vector()
vecSAG600chi2 <- vector(); vecSAG600KW <- vector(); vecSAG600ML <- vector()
vecMOU100chi2 <- vector(); vecMOU100KW <- vector(); vecMOU100ML <- vector()
vecMOU600chi2 <- vector(); vecMOU600KW <- vector(); vecMOU600ML <- vector()
vecBAI100chi2 <- vector(); vecBAI100KW <- vector(); vecBAI100ML <- vector()
vecBAI600chi2 <- vector(); vecBAI600KW <- vector(); vecBAI600ML <- vector()
vecBAL100chi2 <- vector(); vecBAL100KW <- vector(); vecBAL100ML <- vector()
vecBAL585chi2 <- vector(); vecBAL585KW <- vector(); vecBAL585ML <- vector()

for (i in 1:numBS){
  vecSAG100chi2 <- sum(c(vecSAG100chi2, try[[i]]$SAG100_df$chi2 < 0.05))/length(c(vecSAG100chi2, try[[i]]$SAG100_df$chi2 < 0.05))
  vecSAG100KW <- sum(c(vecSAG100KW, try[[i]]$SAG100_df$KW < 0.05))/length(c(vecSAG100KW, try[[i]]$SAG100_df$KW < 0.05))
  vecSAG100ML <- sum(c(vecSAG100ML, try[[i]]$SAG100_df$ML < 0.05))/length(c(vecSAG100ML, try[[i]]$SAG100_df$ML < 0.05))
  #
  vecSAG600chi2 <- sum(c(vecSAG600chi2, try[[i]]$SAG600_df$chi2 < 0.05))/length(c(vecSAG600chi2, try[[i]]$SAG600_df$chi2 < 0.05))
  vecSAG600KW <- sum(c(vecSAG600KW, try[[i]]$SAG600_df$KW < 0.05))/length(c(vecSAG600KW, try[[i]]$SAG600_df$KW < 0.05))
  vecSAG600ML <- sum(c(vecSAG600ML, try[[i]]$SAG600_df$ML < 0.05))/length(c(vecSAG600ML, try[[i]]$SAG600_df$ML < 0.05))
  #####
  vecMOU100chi2 <- sum(c(vecMOU100chi2, try[[i]]$MOU100_df$chi2 < 0.05))/length(c(vecMOU100chi2, try[[i]]$MOU100_df$chi2 < 0.05))
  vecMOU100KW <- sum(c(vecMOU100KW, try[[i]]$MOU100_df$KW < 0.05))/length(c(vecMOU100KW, try[[i]]$MOU100_df$KW < 0.05))
  vecMOU100ML <- sum(c(vecMOU100ML, try[[i]]$MOU100_df$ML < 0.05))/length(c(vecMOU100ML, try[[i]]$MOU100_df$ML < 0.05))
  #
  vecMOU600chi2 <- sum(c(vecMOU600chi2, try[[i]]$MOU600_df$chi2 < 0.05))/length(c(vecMOU600chi2, try[[i]]$MOU600_df$chi2 < 0.05))
  vecMOU600KW <- sum(c(vecMOU600KW, try[[i]]$MOU600_df$KW < 0.05))/length(c(vecMOU600KW, try[[i]]$MOU600_df$KW < 0.05))
  vecMOU600ML <- sum(c(vecMOU600ML, try[[i]]$MOU600_df$ML < 0.05))/length(c(vecMOU600ML, try[[i]]$MOU600_df$ML < 0.05))
  #####
  vecBAI100chi2 <- sum(c(vecBAI100chi2, try[[i]]$BAI100_df$chi2 < 0.05))/length(c(vecBAI100chi2, try[[i]]$BAI100_df$chi2 < 0.05))
  vecBAI100KW <- sum(c(vecBAI100KW, try[[i]]$BAI100_df$KW < 0.05))/length(c(vecBAI100KW, try[[i]]$BAI100_df$KW < 0.05))
  vecBAI100ML <- sum(c(vecBAI100ML, try[[i]]$BAI100_df$ML < 0.05))/length(c(vecBAI100ML, try[[i]]$BAI100_df$ML < 0.05))
  #
  vecBAI600chi2 <- sum(c(vecBAI600chi2, try[[i]]$BAI600_df$chi2 < 0.05))/length(c(vecBAI600chi2, try[[i]]$BAI600_df$chi2 < 0.05))
  vecBAI600KW <- sum(c(vecBAI600KW, try[[i]]$BAI600_df$KW < 0.05))/length(c(vecBAI600KW, try[[i]]$BAI600_df$KW < 0.05))
  vecBAI600ML <- sum(c(vecBAI600ML, try[[i]]$BAI600_df$ML < 0.05))/length(c(vecBAI600ML, try[[i]]$BAI600_df$ML < 0.05))
  #####
  vecBAL100chi2 <- sum(c(vecBAL100chi2, try[[i]]$BAL100_df$chi2 < 0.05))/length(c(vecBAL100chi2, try[[i]]$BAL100_df$chi2 < 0.05))
  vecBAL100KW <- sum(c(vecBAL100KW, try[[i]]$BAL100_df$KW < 0.05))/length(c(vecBAL100KW, try[[i]]$BAL100_df$KW < 0.05))
  vecBAL100ML <- sum(c(vecBAL100ML, try[[i]]$BAL100_df$ML < 0.05))/length(c(vecBAL100ML, try[[i]]$BAL100_df$ML < 0.05))
  #
  vecBAL585chi2 <- sum(c(vecBAL585chi2, try[[i]]$BAL585_df$chi2 < 0.05))/length(c(vecBAL585chi2, try[[i]]$BAL585_df$chi2 < 0.05))
  vecBAL585KW <- sum(c(vecBAL585KW, try[[i]]$BAL585_df$KW < 0.05))/length(c(vecBAL585KW, try[[i]]$BAL585_df$KW < 0.05))
  vecBAL585ML <- sum(c(vecBAL585ML, try[[i]]$BAL585_df$ML < 0.05))/length(c(vecBAL585ML, try[[i]]$BAL585_df$ML < 0.05))
}

resultsDF <- data.frame(chi2 = c(100 - vecSAG100chi2, 100 - vecSAG600chi2, 100 - vecMOU100chi2, 100 - vecMOU600chi2,
                                 100 - vecBAI100chi2, 100 - vecBAI600chi2, 100 - vecBAL100chi2, 100 - vecBAL585chi2), 
                        KW = c(100 - vecSAG100KW, 100 - vecSAG600KW, 100 - vecMOU100KW, 100 - vecMOU600KW,
                               100 - vecBAI100KW, 100 - vecBAI600KW, 100 - vecBAL100KW, 100 - vecBAL585KW),
                        ML = c(100 - vecSAG100ML, 100 - vecSAG600ML, 100 - vecMOU100ML, 100 - vecMOU600ML,
                               100 - vecBAI100ML, 100 - vecBAI600ML, 100 - vecBAL100ML, 100 - vecBAL585ML))

row.names(resultsDF) <- c("SAG100", "SAG600","MOU100", "MOU600","BAI100", "BAI600","BAL100", "BAL585" )

resultsDF

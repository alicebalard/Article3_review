devtools::install_github("alicebalard/parasiteLoad") # version with full Gtest Chisquare
library(parasiteLoad)
#library(optimx)

# Balard et al. 2020
BALdata <- read.csv("https://raw.githubusercontent.com/alicebalard/Article_IntensityEimeriaHMHZ/master/data/cleanedData.csv", na.strings = c("", " ", NA))
# Keep individuals with hybrid index, sex recorded and pinworms counted
BALdata <- BALdata[!is.na(BALdata$HI) & !is.na(BALdata$Sex) & !is.na(BALdata$Aspiculuris_Syphacia),]
length(BALdata$Aspiculuris_Syphacia) # N=585
# prevalence:
sum(BALdata$Aspiculuris_Syphacia >0) / length(BALdata$Aspiculuris_Syphacia) # N=585
# 52.4%

# for repetitivity
set.seed(8)

# Shuffle data:
shuffled_data <- data.frame(HI = sample(BALdata$HI), 
                            load = BALdata$Aspiculuris_Syphacia,
                            Sex = factor(BALdata$Sex))

# Apply test:
  T <- parasiteLoad::analyse(data = shuffled_data, response = "load", model = "negbin",
                      group = "Sex", hybridIndex = "HI")


### EXTRACT p-Value
T$H0

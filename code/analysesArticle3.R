### Code for figure making of review Res/tol 
### Feb 2022
### Alice Balard

#### Load data and functions within ####
source("dataPreparation_168mice.R")
source("myFunctions.R")
library(cowplot)
library(ggplot2)
library(dplyr)
library(lmtest)
library(pscl)
library(scales)
# mycolors <- scales::seq_gradient_pal("blue", "red", "Lab")(seq(0,1,length.out=8))
# New colors, scales too difficult to distinguish
mycolors <- c("blue", "blue", "blue","purple","purple", "red", "red", "red")

# what's my dataset?
table(art3SummaryDF$infection_isolate, art3SummaryDF$Mouse_genotype)
#                                 SCHUNT STRA SCHUNT-STRA STRA-BUSNA SCHUNT-PWD PWD-BUSNA BUSNA PWD
# Brandenburg64 (E. ferrisi)         10   10           2          6          7         6    10   8
# Brandenburg88 (E. falciformis)      4    7           3          6          5         6     7   6

######################
## KEEP ONLY OUTBRED!!
art3outbredDF <- art3FullDF[grep("-", art3FullDF$Mouse_genotype),]
# art3outbredDF$dpi <- as.factor(art3outbredDF$dpi) # important for violin

## summarise data per animal:
art3OutbreSummaryData <- makeSummaryTable(art3outbredDF)

## how many mice did not shed oocysts?
art3OutbreSummaryData[art3OutbreSummaryData$sumoocysts.per.tube == 0,]
ggplot(art3outbredDF[art3outbredDF$EH_ID %in% "LM0250",],
       aes(x=dpi, y=oocysts.per.tube)) + geom_line() + geom_point()

## Let's remove it:
art3outbredDF <- art3outbredDF[!art3outbredDF$EH_ID %in% "LM0250",]
art3OutbreSummaryData <- makeSummaryTable(art3outbredDF)

table(art3OutbreSummaryData$Mouse_subspecies, art3OutbreSummaryData$infection_isolate)
#                      Brandenburg64 (E. ferrisi) Brandenburg88 (E. falciformis)
# M.m.dom                                 2                              3
# Hybrid_mus_dom                         13                             11
# M.m.mus                                 6                              5

## what is the overall peak day for each parasite isolate? SAME as article 2 (as expected)
aggregate(art3OutbreSummaryData$dpi_max.OPG,
          list(art3OutbreSummaryData$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
# 1     Brandenburg64 (E. ferrisi) 21 6 0.46
# 2 Brandenburg88 (E. falciformis) 19 8 0.37
aggregate(art3OutbreSummaryData$dpi_minWeight,
          list(art3OutbreSummaryData$infection_isolate), 
          function(x) {paste(length(x), median(x), round(sd(x),2))})
# 1     Brandenburg64 (E. ferrisi) 21 5 0.94
# 2 Brandenburg88 (E. falciformis) 19 9 1.46

## Plot course of infection:
###### Course of infection FIGURE 2 ######
forplot <- art3outbredDF %>%
  group_by(infection_isolate, dpi) %>%
  summarise(mean = mean(OPG*10e-6, na.rm = TRUE),
            sd = sd(OPG*10e-6, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)

F2.1 <- ggplot(forplot, aes(dpi, mean, group = infection_isolate, col = infection_isolate)) + 
  geom_point(size = 3) +
  geom_line(aes(linetype=infection_isolate)) +
  scale_linetype_manual(values = c(1,2,1)) +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = .2)+
  ylab("million oocysts per gram of feces") +
  scale_x_continuous(breaks = 0:11, name = "days post infection") +
  scale_color_manual(values = c("darkgreen", "orange"))+
  theme(legend.position = c(0.25, 0.8)) +
  labs(color = "Eimeria isolate") 

forplot2 <-  art3outbredDF %>%
  group_by(infection_isolate, dpi) %>%
  summarise(mean = mean(relativeWeight, na.rm = TRUE),
            sd = sd(relativeWeight, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)

F2.2 <- ggplot(forplot2, aes(dpi, mean, group = infection_isolate, 
                             col = infection_isolate)) + 
  geom_point(size = 3) +
  geom_line(aes(linetype=infection_isolate)) +
  scale_linetype_manual(values = c(1,2,1)) +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = .2)+
  ylab("relative weight compared to day 0 (%)") +
  scale_x_continuous(breaks = 0:11, name = "days post infection") +
  scale_color_manual(values = c("darkgreen", "orange"))+
  theme(legend.position = c(0.25, 0.2)) +
  labs(color = "Eimeria isolate") 

Fig_courseOfInfection <- cowplot::plot_grid(F2.1, F2.2,
                                            labels=c("A", "B"), label_size = 20)

Fig_courseOfInfection

############# STATS
## LRT significance for each factor
myLRTsignificanceFactors <- function(modFull, modPar, modMouse, modInt){
  # print("significance of parasite:")
  return(list(signifParasite = lrtest(modFull, modPar),
              signifMouse = lrtest(modFull, modMouse),
              signifInter = lrtest(modFull, modInt)))
}

#############################
## RESISTANCE: inverse of OPG
modFULL <- glm.nb(max.OPG ~ infection_isolate*Mouse_subspecies, data = art3OutbreSummaryData)
modPara <- glm.nb(max.OPG ~ Mouse_subspecies, data = art3OutbreSummaryData)
modMous <- glm.nb(max.OPG ~ infection_isolate, data = art3OutbreSummaryData)
modinter <- glm.nb(max.OPG ~ infection_isolate+Mouse_subspecies, data = art3OutbreSummaryData)
myLRTsignificanceFactors(modFULL, modPara, modMous, modinter)
## NOTHING is significant
## predicted values:
predRes <- data.frame(ggpredict(modFULL, terms = c("infection_isolate", "Mouse_subspecies")))

## Within parasites:
data1 <- art3OutbreSummaryData[art3OutbreSummaryData$infection_isolate %in% "Brandenburg64 (E. ferrisi)",]
modFULL <- glm.nb(max.OPG ~ Mouse_subspecies, data = data1)
mod0 <- glm.nb(max.OPG ~ 1, data = data1)
lrtest(modFULL, mod0) # NOPE

data2 <- art3OutbreSummaryData[art3OutbreSummaryData$infection_isolate %in% "Brandenburg88 (E. falciformis)",]
modFULL <- glm.nb(max.OPG ~ Mouse_subspecies, data = data2)
mod0 <- glm.nb(max.OPG ~ 1, data = data2)
lrtest(modFULL, mod0) # NOPE

#############################
## WEIGHT LOSS
modFULL <- lm(relWL~infection_isolate*Mouse_subspecies, data = art3OutbreSummaryData)
modPara <- lm(relWL~Mouse_subspecies, data = art3OutbreSummaryData)
modMous <- lm(relWL~infection_isolate, data = art3OutbreSummaryData)
modinter <- lm(relWL~infection_isolate+Mouse_subspecies, data = art3OutbreSummaryData)
myLRTsignificanceFactors(modFULL, modPara, modMous, modinter)
## NOTHING is significant
## predicted values:
predImp <- data.frame(ggpredict(modFULL, terms = c("infection_isolate", "Mouse_subspecies")))

## Within parasites:
modFULL <- lm(relWL ~ Mouse_subspecies, data = data1)
mod0 <- lm(relWL ~ 1, data = data1)
lrtest(modFULL, mod0) # NOPE

modFULL <- lm(relWL ~ Mouse_subspecies, data = data2)
mod0 <- lm(relWL ~ 1, data = data2)
lrtest(modFULL, mod0) # NOPE

#############################
## TOLERANCE
modFULL <- lm(relWL ~ 0 + max.OPG : (infection_isolate * Mouse_subspecies), data = art3OutbreSummaryData)
modPara <- lm(relWL ~ 0 + max.OPG : (Mouse_subspecies), data = art3OutbreSummaryData)
modMous <- lm(relWL ~ 0 + max.OPG : (infection_isolate), data = art3OutbreSummaryData)
modinter <- lm(relWL ~ 0 + max.OPG : (infection_isolate + Mouse_subspecies), data = art3OutbreSummaryData)
myLRTsignificanceFactors(modFULL, modPara, modMous, modinter)
## MOUSE SUBSPECIES is significant
# p=0.02749
modFULL <- lm(relWL ~ 0 + max.OPG : (infection_isolate * Mouse_subspecies), data = art3OutbreSummaryData)
predTol <- data.frame(ggpredict(modFULL, terms = c("Mouse_subspecies", "infection_isolate"),
                                condition = c(max.OPG = 1000000)))  ## For a million OPG

## Within parasites:
modFULL <- lm(relWL ~ 0 + max.OPG : Mouse_subspecies, data = data1)
mod0 <- lm(relWL ~ 0 + max.OPG, data = data1)
lrtest(modFULL, mod0) # NOPE

modFULL <- lm(relWL ~ 0 + max.OPG : Mouse_subspecies, data = data2)
mod0 <- lm(relWL ~ 0 + max.OPG, data = data2)
lrtest(modFULL, mod0) # YASS in E88

###################################################################################
## Plot predicted values
cols <- c("blue", "purple", "red")

plotRES <- ggplot(data.frame(predRes), aes(x=group, y=predicted, col=group)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width=.1) +
  geom_point(size = 5) +
  scale_color_manual(values = cols)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_grid(.~x) +
  scale_y_continuous("(predicted) maximum million OPG \n(oocysts per gram of feces)")

plotIMP <- ggplot(data.frame(predImp), aes(x=group, y=predicted, col=group)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width=.1) +
  geom_point(size = 5) +
  scale_color_manual(values = cols)+
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_grid(.~x) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     "(predicted) maximum weight loss \nrelative to day of infection")

predTol$predicted <- predTol$predicted*4
predTol$relWL_OPGnull <- 0
names(predTol)[names(predTol) %in% c("predicted")] <- "relWL_NMOPG"
predTol <- melt(predTol, measure.vars = c("relWL_OPGnull", "relWL_NMOPG"))
names(predTol)[names(predTol) %in% c("variable", "value")] <- c("max.OPG", "relWL")
predTol$max.OPG <- as.character(predTol$max.OPG)
predTol$max.OPG[predTol$max.OPG %in% "relWL_OPGnull"] <- "0"
predTol$max.OPG[predTol$max.OPG %in% "relWL_NMOPG"] <- 4e6
predTol$max.OPG <- as.numeric(predTol$max.OPG)
names(predTol)[names(predTol) %in% "x"] <- "Mouse_subspecies"

# A data frame with labels for each facet
df_labels <- data.frame(group = c("Brandenburg64 (E. ferrisi)", "Brandenburg88 (E. falciformis)"),
                        label = c(NA, "Significant differences between subspecies\nLRT: mouse subspecies: G = 7, df = 2, p = 0.03"))

plotTOL <- ggplot(predTol, aes(x = max.OPG, y = relWL)) +
  geom_line(aes(group = Mouse_subspecies, col = Mouse_subspecies)) +
  # geom_label(aes(label = substring(label, 1, 1)), na.rm = T)+
  scale_color_manual(values = cols) +
  scale_x_continuous("maximum million OPG \n(oocysts per gram of feces)",
                     breaks = seq(0, 4000000, 1000000),
                     labels = seq(0, 4000000, 1000000)/1000000) +
  scale_y_continuous(name = "maximum weight loss \nrelative to day of infection",
                     breaks = seq(0,0.3, 0.05),
                     labels = scales::percent_format(accuracy = 5L))+
  facet_grid(.~group) +
  geom_point(data = art3OutbreSummaryData, size = 3, alpha = .5, aes(pch=Mouse_subspecies, col = Mouse_subspecies))+
  ggtitle("Tolerance \n(slope of B (max weight loss) on A (max parasite load), per genotype)") +
  coord_cartesian(ylim=c(0, 0.2)) +
  geom_text(x = 2e6, y = 0.19, aes(label = label), data = df_labels)

### For coupled plots:
names(predRes)[names(predRes) %in% "x"] <- "infection_isolate"
names(predRes)[names(predRes) %in% "group"] <- "Mouse_subspecies"
names(predImp)[names(predImp) %in% "x"] <- "infection_isolate"
names(predImp)[names(predImp) %in% "group"] <- "Mouse_subspecies"
## Re predTol
modFULL <- lm(relWL ~ 0 + max.OPG : (infection_isolate * Mouse_subspecies), data = art3OutbreSummaryData)
predTol <- data.frame(ggpredict(modFULL, terms = c("Mouse_subspecies", "infection_isolate"),
                                condition = c(max.OPG = 1000000)))  ## For a million OPG
names(predTol)[names(predTol) %in% "group"] <- "infection_isolate"
names(predTol)[names(predTol) %in% "x"] <- "Mouse_subspecies"

names(predRes)[names(predRes) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
  paste(names(predRes)[names(predRes) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "OPG", sep = "_")
names(predImp)[names(predImp) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
  paste(names(predImp)[names(predImp) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "relWL", sep = "_")
names(predTol)[names(predTol) %in% c("predicted", "conf.low",  "std.error", "conf.high")] <-
  paste(names(predTol)[names(predTol) %in% c("predicted", "conf.low",  "std.error", "conf.high")], "Tol", sep = "_")

finalplotDF <- merge(merge(predRes, predImp), predTol, by = c("Mouse_subspecies", "infection_isolate"))

# test correlations:
cor.test(finalplotDF$predicted_OPG[finalplotDF$infection_isolate %in% "Brandenburg64 (E. ferrisi)"],
         finalplotDF$predicted_relWL[finalplotDF$infection_isolate %in% "Brandenburg64 (E. ferrisi)"], method="spearman")
cor.test(finalplotDF$predicted_OPG[finalplotDF$infection_isolate %in% "Brandenburg64 (E. ferrisi)"],
         finalplotDF$predicted_Tol[finalplotDF$infection_isolate %in% "Brandenburg64 (E. ferrisi)"], method="spearman")

cor.test(finalplotDF$predicted_OPG[finalplotDF$infection_isolate %in% "Brandenburg88 (E. falciformis)"],
         finalplotDF$predicted_relWL[finalplotDF$infection_isolate %in% "Brandenburg88 (E. falciformis)"], method="spearman")
cor.test(finalplotDF$predicted_OPG[finalplotDF$infection_isolate %in% "Brandenburg88 (E. falciformis)"],
         finalplotDF$predicted_Tol[finalplotDF$infection_isolate %in% "Brandenburg88 (E. falciformis)"], method="spearman")

# Plot
## Plot1
PC1 <- ggplot(finalplotDF, aes(x=predicted_OPG, y=predicted_relWL)) +
  # geom_smooth(method = "lm", se = F, col = "black")+
  geom_errorbarh(aes(xmin = conf.low_OPG, xmax = conf.high_OPG), color = "grey") +
  geom_errorbar(aes(ymin = conf.low_relWL, ymax = conf.high_relWL), color = "grey") +
  geom_point(aes(col = Mouse_subspecies), size = 7)+
  scale_x_continuous(name = "Maximum million oocysts per gram \n (inverse of) RESISTANCE ",
                     breaks = seq(0.5,5,0.5)*1e6,
                     labels = seq(0.5,5,0.5)) +
  scale_y_continuous(name = "Maximum relative weight loss",
                     labels = scales::percent_format(accuracy = 1))+
  scale_color_manual(values = cols) +
  theme(legend.position = "top")+
  theme(text = element_text(size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(.~infection_isolate)

PC2 <- ggplot(finalplotDF, aes(x=predicted_OPG, y=-predicted_Tol)) +
  geom_errorbarh(aes(xmin = conf.low_OPG, xmax = conf.high_OPG), color = "grey") +
  geom_errorbar(aes(ymin = -conf.low_Tol, ymax = -conf.high_Tol), color = "grey") +
  geom_point(aes(col = Mouse_subspecies), size = 7)+
  scale_x_continuous(name = "Maximum million oocysts per gram \n (inverse of) RESISTANCE ",
                     breaks = seq(0.5,5,0.5)*1e6,
                     labels = seq(0.5,5,0.5)) +
  scale_y_continuous(name = "TOLERANCE (inverse of) slope of\n maximum weight loss on maximum oocysts per gram",
                     labels = scales::percent_format(accuracy = 1))+
  scale_color_manual(values = cols) +
  theme(legend.position = "top")+
  theme(text = element_text(size=15), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_grid(.~infection_isolate)

### Read plots
plotRES
plotIMP
plotTOL
pdf(file = "../figures/Figure6.pdf", width = 15)
plot_grid(PC1,PC2)
dev.off()
"Significant differences in tolerance between the three subspecies (LRT: mouse subspecies: G = 7, df = 2, p = 0.03) in E. falciformis"

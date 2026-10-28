#Figures

library(dplyr)
library(tidyverse)
library(qqplotr)
library(patchwork)

df <- read.csv("C:/Users/willi/Downloads/proteins_all.csv") 

df_avg <- read.csv("C:/Users/willi/Downloads/proteins_all.csv") %>% 
  mutate(Sample.ID = substr(Sample.ID, 1, nchar(Sample.ID)-4)) %>% unique()


#Boxplots for Early and Late Data
early_box <- ggplot(df_early, aes(x = Res_Lev, y = PR_Conc)) + geom_boxplot() + lims(y = c(0,2)) + theme_classic() + labs(title = "A", x = "Resistance Level", y = "Protein Concentration (mg protein/g dry weight)")
late_box <- ggplot(df_late, aes(x=Res_Lev, y = PR_Conc)) + geom_boxplot() + lims(y = c(0,2))+ theme_classic() + labs(title = "B", x = "Resistance Level", y = "")
early_box + late_box

#Boxplot for All
ggplot(df_avg, aes(x = Res_Lev, y = PR_Conc)) + geom_boxplot() + theme_classic() + labs(title = "A", x = "Resistance Level", y = "Protein Concentration (mg protein/g dry weight)")

#Exclude x121
df_RandM <- df %>% filter(Res_Lev != "Not Resistant")

#Evaluate normality
hist(df_RandM$PR_Conc)
ggplot(df_RandM, aes(sample = PR_Conc)) + qqplotr::stat_qq_band() + stat_qq_line() + stat_qq_point()

#Resistant vs Non-Resistant, Early
df_RM_June <- df %>% filter(Res_Lev != "Not Resistant", Month == "June")
ggplot(df_RM_June, aes(x=Res_Lev, y = PR_Conc)) + geom_boxplot()

#Resistant vs Non-Resistant, Late
df_RM_Aug <- df %>% filter(Res_Lev != "Not Resistant", Month == "August")
ggplot(df_RM_Aug, aes(x=Res_Lev, y=PR_Conc)) + geom_boxplot()


#------------------------------------------------------------------------------#
# Stats

library(tidyverse)
library(ggplot2)
library(car)
library(multcomp)
library(forcats)
library(ggsignif)
library(rstatix)

# Import data 
joint_df <- df_avg 
joint_df$Genotype <- NULL
ext_df <- joint_df |> filter(Res_Lev != "Medium")

# Exploratory data analysis
aug_df <- joint_df |> filter(Month == "August")
jun_df <- joint_df |> filter(Month == "June")
tapply(aug_df$PR_Conc, aug_df$Res_Lev, summary)
tapply(jun_df$PR_Conc, jun_df$Res_Lev, summary)

## Evaluating assumptions of each group to see if ANOVA is realistic

#Normality?

all_aov <- aov(PR_Conc ~ Res_Lev*Month, data = joint_df)
par(mfrow = c(1,2))
hist(all_aov$residuals)
qqPlot(all_aov$residuals, id = FALSE)
shapiro.test(all_aov$residuals)


res_aug_aov<- aov(PR_Conc ~ Res_Lev, data = aug_df)
par(mfrow = c(1,2))
hist(res_aug_aov$residuals)
qqPlot(res_aug_aov$residuals, id = FALSE)

res_jun_aov <- aov(PR_Conc ~ Res_Lev, data = jun_df)
par(mfrow = c(1,2))
hist(res_jun_aov$residuals)
qqPlot(res_jun_aov$residuals, id = FALSE)

# Equal Variance?
leveneTest(PR_Conc ~ Res_Lev*Month, data = joint_df)
plot(all_aov, which = 3)

leveneTest(PR_Conc ~ Res_Lev, data = aug_df)
leveneTest(PR_Conc ~ Res_Lev, data = jun_df)

# ANOVA

#August
aug_lm<- lm(PR_Conc ~ Res_Lev, data = aug_df)
aug_av<-aov(aug_lm)
summary(aug_av)
tukey_aug <- TukeyHSD(aug_av)
tukey_aug

#June
jun_lm <- lm(PR_Conc ~ Res_Lev, data = jun_df)
jun_av <- aov(jun_lm)
summary(jun_av)
tukey_jun <- TukeyHSD(jun_av)
tukey_jun

#Both (two way with interaction term)
summary(all_aov)
tukey_all <- TukeyHSD(all_aov)
tukey_all

#At Justin's recommendation, kruskall wallis and dunn post hoc
## All together
kruskal.test(PR_Conc ~ Res_Lev, data = joint_df)
prc <- joint_df |> dunn_test(PR_Conc ~ Res_Lev, p.adjust.method = "bonferroni")
prc

## June
kruskal.test(PR_Conc ~ Res_Lev, data = jun_df)
prc_j <- jun_df |> dunn_test(PR_Conc ~ Res_Lev, p.adjust.method = "bonferroni")
prc_j

## August
kruskal.test(PR_Conc ~ Res_Lev, data = aug_df)
prc_a <- aug_df |> dunn_test(PR_Conc ~ Res_Lev, p.adjust.method = "bonferroni")
prc_a







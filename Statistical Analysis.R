#Spectral Data Analysis
#------------------------------------------------------------SETUP----------------------------------------------------
#Installing necessary packages
library(ggplot2) #plotting
library(tidyverse) #dataframe handling
library(stringr) #you'll need this to mess with strings later
library(fdANOVA) #FANOVA for whole spectral signature comparison
library(mvnormtest) #assumption testing for MANOVA on important wavelengths
library(dplyr) #just good to have


set.seed(123)
#------------------------------------------------------------IMPORT DATA-----------------------------------------------------------
#Raw Data
df <- read.csv("C:/Users/willi/Desktop/Master's Thesis Workbook/Analyses/TBA_PhloemWL.csv") |> na.omit()

#Raw averaged data
df_EX <- read_csv("C:/Users/willi/Desktop/Master's Thesis Workbook/Analyses/TotalWL/df_EX.csv")

#Transformed average data
df_EX_2derv <- read_csv("C:/Users/willi/Desktop/Master's Thesis Workbook/Analyses/TotalWL/df_EX_2derv.csv")
#Widen raw data
df_EX_wide <- df_EX |> pivot_wider(names_from = Wavelength, values_from = Mean)

#Widen transformed data
df_EX_2derv_wide_all <- df_EX_2derv |> pivot_wider(names_from = Wavelength, values_from = Mean)


#10 wavelengths RAW
imp_wl <- c("X2400.000057", "X2141.108036", "X2153.665849", "X2160.000158", "X2134.883888",
            "X2166.371837", "X1462.948523", "X2147.368583", "X2192.238952", "X2122.543524")
df_EX_ten <- df_EX |> 
  filter(Wavelength %in% imp_wl)

df_EX_wide_10 <- df_EX |> filter(Wavelength %in% imp_wl) |> pivot_wider(names_from = Wavelength, values_from = Mean)

#10 WL transformed
df_EX_2derv_10 <- df_EX_2derv |> filter(Wavelength %in% imp_wl)

#10 WL transformed wider
df_EX_2derv_wide_10 <- df_EX_2derv |>
  filter(Wavelength %in% imp_wl) |> 
  pivot_wider(names_from = Wavelength, values_from = Mean)



#--------------------------------------------------------- IMPORTANT WAVELENGTH ANALYSIS ------------------------------------------------
#Assigning independent/dependent variables for MANOVA
dep_vars <- cbind(unlist(df_EX_2derv_wide_10$X1462.948523), unlist(df_EX_2derv_wide_10$X2122.543524), unlist(df_EX_2derv_wide_10$X2134.883888), unlist(df_EX_2derv_wide_10$X2141.108036), unlist(df_EX_2derv_wide_10$X2147.368583), unlist(df_EX_2derv_wide_10$X2153.665849),unlist(df_EX_2derv_wide_10$X2160.000158),unlist(df_EX_2derv_wide_10$X2166.371837),unlist(df_EX_2derv_wide_10$X2192.238952),unlist(df_EX_2derv_wide_10$X2400.000057))
ind_var <- df_EX_2derv_wide_10$Genotype

#Checking assumptions
mvnormtest::mshapiro.test(t(dep_vars)) #Shapiro-Wilk (Multivariate Normality Test)
biotools::boxM(dep_vars, ind_var) #Box's M-test (Homogeneity of Covariance Test)
distances <- mahalanobis(dep_vars, colMeans(dep_vars), cov(dep_vars)) #Mahalanobis distance for outliers
cutoff <- qchisq(0.999, df=ncol(dep_vars))
outliers <- distances > cutoff
sum(outliers)
heplots::cqplot(df_EX_2derv_wide_10[, imp_wl]) #Quantile plot

#MANOVA
manova_2derv <- manova(dep_vars ~ ind_var, data = df_EX_2derv_wide_10)
summary(manova_2derv,tol=0)

#post-hoc t tests
t_imp_wl_testList <- list()
p_vals <- numeric()

#Putting all t-tests in a list for later access
for(i in 1:length(imp_wl)){
  df_func <- df_EX_2derv_10 |> 
    filter(Wavelength == imp_wl[i])
  t_imp_wl_testList[[length(t_imp_wl_testList)+1]] <- t.test(df_func$Mean ~ df_func$Genotype)
  print(c("Student's t-test for Important Wavelength", imp_wl[i]))
  print(t_imp_wl_testList[[length(t_imp_wl_testList)]])
  p_vals[i] <- t_imp_wl_testList[[i]]$p.value
}

#Bonferroni correction for multiple comparisons
p_adjusted<- p.adjust(p_vals, method = "bonferroni")
print(p_adjusted)
#--------------------------------------------------------- WHOLE SPECTRAL SIGNATURE ANALYSIS -----------------------------------------------

#Wrangling data to fit the specs of the fanova.tests() function
x97 <- df_EX_2derv |> filter(Genotype == "X97")
x97$Sample.Name <- match(x97$Sample.Name, unique(x97$Sample.Name)) #reassigning sample id's to simple index for readability
x97 <- x97 |> pivot_wider(names_from = c(Genotype, Sample.Name), values_from = Mean)

x121 <- df_EX_2derv |> filter(Genotype == "X121")  
x121$Sample.Name <- match(x121$Sample.Name, unique(x121$Sample.Name)) #reassigning sample id's to simple index for readability
x121 <- x121 |> pivot_wider(names_from = c(Genotype, Sample.Name), values_from = Mean)

dat_fanova <- left_join(x97,x121) #joining dataframes
group.label <- c(rep("x97", times=34), rep("x121", times=33)) #making vector of genotypes that is going to align with new dataframe

dep_vars <- dat_fanova[,2:68] #identify response variables

plotFANOVA(x = dep_vars, group.label = group.label, means=TRUE)
fanova_results <- fanova.tests(dep_vars, group.label = group.label, test="FP") #conduct FANOVA
summary(fanova_results)

#mean difference curve
mean_curves <- df_EX_2derv |> 
  group_by(Genotype, Wavelength) |> 
  summarise(Mean = mean(Mean)) |> 
  pivot_wider(names_from = Genotype, values_from = Mean)

mean_curves <- mean_curves |> mutate(difference = X97-X121, Wavelength = as.numeric(str_remove_all(Wavelength, "X")))
ggplot(mean_curves, aes(x=Wavelength, y=difference)) + geom_line() + ylab("Difference in Transformed Reflectance") + xlab("Wavelength (nm)")

ggsave("Figure3.png", path = getwd(), units="mm", width = 85, height = 57, dpi = 300, device = "png", bg = "white")

#install.packages("optimx", version = "2021-10.12") # this package is required for 
#the parasite load package to work
require(devtools)

## install the pacakage of Alice Balard
devtools::install_github("alicebalard/parasiteLoad@v2.0", force = TRUE)
#force = TRUE)

library(parasiteLoad)
library(tidyverse)
library(tidyr)
library(dplyr)
library(cowplot)
library(randomForest)
library(ggplot2)
library(VIM) # visualizing missing data
library(mice) # imputing missing data without predictors
library(ggpubr)
library(optimx)
library(rfUtilities) # Implements a permutation test cross-validation for 
library(fitdistrplus) #testing distributions
library(logspline)
library(caret)
library(dplyr)
library(tidyr)

# read the data
hm <- read.csv("Data/Data_output/imputed_clean_data.csv")

# filter for the field mice
Field <- hm %>%
  filter(origin == "Field") %>%
    drop_na(HI)

# Create vectors for selecting relevant columns
EqPCR.cols      <- c("delta_ct_cewe_MminusE", "MC.Eimeria", "Ct.Eimeria")
#,"Ct.Mus")

Genes_v   <- c("IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10", 
                  "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                  "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                  "TICAM1", "TNF") #, "IL.12", "IRG6")

# select the gene columns
gene <-  Field %>%
  dplyr::select(c(Mouse_ID, "IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10", 
                  "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                  "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                  "TICAM1", "TNF"))

# data frame with only the genes
genes <- gene %>%
  dplyr::select(-Mouse_ID)


# load predicting weight loss model
weight_loss_predict <- readRDS("R/Models/predict_WL.rds")

set.seed(540)


#The predict() function in R is used to predict the values based on the input data.
predicted_WL <- predict(weight_loss_predict, genes)


# assign test.data to a new object, so that we can make changes
result_field <- genes

#add the new variable of predictions to the result object
result_field <- cbind(result_field, predicted_WL)

# add it to the field data 
Field <- cbind(Field, predicted_WL)

rm(gene,genes)

########## Analyzing the distribution of our data in order to 
# go on with the anaylsis 
Field %>% ggplot(aes(x = predicted_WL)) +
  geom_histogram(binwidth = 1.5)


##  predicted WL vs HI
Field %>%
    ggplot(aes(x = HI , y = predicted_WL , color = Sex)) +
    geom_smooth() +
    geom_point()

## body length vs predicted WL
Field %>%
    ggplot(aes(x = Body_Length , y = predicted_WL , color = Sex)) +
    geom_smooth() +
    geom_point()


## Let'S further analyse the distribution of WL
x <- Field$predicted_WL

descdist(data = x, discrete = FALSE)
descdist(data = x, discrete = FALSE, #data is continuous
         boot = 1000)


### quite close to normal distribution

## fitting different distributions
set.seed(10)
n = 25
size = 27
prob = .4
data = rbinom(x, size = size, prob = prob)
fit = fitdist(data = data, dist="binom", 
                   fix.arg=list(size = size), 
                   start=list(prob = 0.1))

summary(fit)


plot(fit)


## ---------------------------------------------------------------------------------------------------
normal_ <- fitdist(x, "norm")
weibull_ <- fitdist(x, "weibull")
gamma_ <- fitdist(x, "gamma")


# Define function to be used to test, get the log lik and aic
tryDistrib <- function(x, distrib){
  # deals with fitdistr error:
  fit <- tryCatch(MASS::fitdistr(x, distrib), error=function(err) "fit failed")
  return(list(fit = fit,
              loglik = tryCatch(fit$loglik, error=function(err) "no loglik computed"), 
              AIC = tryCatch(fit$aic, error=function(err) "no aic computed")))
}



findGoodDist <- function(x, distribs, distribs2){
  l =lapply(distribs, function(i) tryDistrib(x, i))
  names(l) <- distribs
  print(l)
  listDistr <- lapply(distribs2, function(i){
    if (i %in% "t"){
      fitdistrplus::fitdist(x, i, start = list(df =2))
    } else {
      fitdistrplus::fitdist(x,i)
    }}
  ) 
  par(mfrow=c(2,2))
  denscomp(listDistr, legendtext=distribs2)
  cdfcomp(listDistr, legendtext=distribs2)
  qqcomp(listDistr, legendtext=distribs2)
  ppcomp(listDistr, legendtext=distribs2)
  par(mfrow=c(1,1))
}


## ---------------------------------------------------------------------------------------------------
tryDistrib(x, "normal")
tryDistrib(x, "binomial")
tryDistrib(x, "student")
tryDistrib(x, "weibull")
tryDistrib(x, "weibullshifted")



## ---------------------------------------------------------------------------------------------------
findGoodDist(x, "normal", "weibull")


## ----normal-----------------------------------------------------------------------------------------
plot(normal_)
summary(normal_)
plot(gamma_)
summary(gamma_)
plot(weibull_)
summary(weibull_)


## ---------------------------------------------------------------------------------------------------
Field$Sex <- as.factor(Field$Sex)


##All
fitWL_Sex <- parasiteLoad::analyse(data = Field,
                        response = "predicted_WL",
                        model = "normal",
                        group = "Sex")


plot_WL_Sex<- bananaPlot(mod = fitWL_Sex$H3,
             data = Field,
             response = "predicted_WL",
             group = "Sex") +
    scale_fill_manual(values = c("brown", "forestgreen")) +
  scale_color_manual(values = c("brown", "forestgreen")) +
  theme_bw() 

plot_WL_Sex

ggsave(plot = plot_WL_Sex, filename = "figures/hybrid_sex.jpeg", width = 10, 
       height = 8, dpi = 1000)

# Create HI bar
HIgradientBar <- ggplot(data.frame(hi = seq(0,1,0.0001)),
                        aes(x=hi, y=1, fill = hi)) +
  geom_tile() +
  theme_void() +
  scale_fill_gradient(low = "blue", high = "red")  + 
  scale_x_continuous(expand=c(.01,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position = 'none')

plot_WL_Sex <- 
    plot_grid(plot_WL_Sex, 
          HIgradientBar,
          nrow = 2,
          align = "v",
          axis = "tlr",
          rel_heights = c(13, 1))


ggsave(plot = plot_WL_Sex, filename = "figures/hybrid_sex.jpeg", width = 10, 
       height = 8, dpi = 1000)



###################### according to infection

mc_field <- Field %>%
    drop_na(MC.Eimeria) 

mc_field$MC.Eimeria <- as.factor(mc_field$MC.Eimeria)
 

##All
fitWL_mc <- parasiteLoad::analyse(data = mc_field,
                                   response = "predicted_WL",
                                   model = "normal",
                                   group = "MC.Eimeria")


plot_WL_MC_eimeria<- bananaPlot(mod = fitWL_mc$H3,
                         data = mc_field,
                         response = "predicted_WL",
                         group = "MC.Eimeria") +
    scale_fill_manual(values = c("blueviolet", "limegreen")) +
    scale_color_manual(values = c("blueviolet", "limegreen")) +
    theme_bw() 

plot_WL_MC_eimeria

plot_WL_MC_eimeria <- 
    plot_grid(plot_WL_MC_eimeria, 
              HIgradientBar,
              nrow = 2,
              align = "v",
              axis = "tlr",
              rel_heights = c(13, 1))

plot_WL_MC_eimeria

ggsave(plot = plot_WL_MC_eimeria, filename = "figures/hybrid_mc_eimeria.jpeg", 
       width = 10, 
       height = 8, dpi = 1000)

#### oocysts
speparam <- c(L1start = 10.098368660  ,
              L1LB = 4.363865741 ,
              L1UB = 19.383819138 ,
              L2start = 10.098368660  ,
              L2LB = 4.363865741 ,
              L2UB = 19.383819138  ,
              alphaStart = 0, alphaLB = -5, alphaUB = 5,
              myshapeStart = 1, myshapeLB = 0.000000001, myshapeUB = 10)



## oocysts
oo_Field <- Field %>%
    drop_na(OPG) %>%
    filter(!OPG == 0)

##All
fitWL_oo <- parasiteLoad::analyse(data = oo_Field,
                                  response = "predicted_WL",
                                  model = "normal",
                                  group = "Sex")


plot_WL_OOC<- bananaPlot(mod = fitWL_oo$H3,
                                data = oo_Field,
                                response = "predicted_WL",
                                group = "Sex") +
    scale_fill_manual(values = c("blueviolet", "limegreen")) +
    scale_color_manual(values = c("blueviolet", "limegreen")) +
    theme_bw() 

plot_WL_OOC



## ---------------------------------------------------------------------------------------------------


ggplot(data = Field, aes(x = delta_ct_cewe_MminusE, y = predicted_WL)) +
  geom_point() +
  stat_smooth(method= "lm") 

Field2 <- Field %>%
  drop_na(delta_ct_cewe_MminusE)

cor(Field2$predicted_WL, Field2$delta_ct_cewe_MminusE)


tolerance <- lm(predicted_WL ~  delta_ct_cewe_MminusE, data = Field)


summary(tolerance)

confint(tolerance)



## ---------------------------------------------------------------------------------------------------
ggplot(data = Field, aes(x = OPG, y = predicted_WL)) +
  geom_point() +
  stat_smooth(method= "lm") +
  scale_x_log10()

Field2 <- Field %>%
  drop_na(OPG)

cor(Field2$predicted_WL, Field2$OPG)


tolerance <- lm(predicted_WL ~  OPG, data = Field)


summary(tolerance)

confint(tolerance)



## ---------------------------------------------------------------------------------------------------

tolerance <- lm(predicted_WL ~  OPG * delta_ct_cewe_MminusE, data = Field)


summary(tolerance)

confint(tolerance)



## ---------------------------------------------------------------------------------------------------
Field <- Field %>%
  dplyr::mutate(BMI = Body_Weight / (Body_Length)) #^2) which is the correct
# way to calculatebmi?

ggplot(data = Field, aes(x = BMI, y = predicted_WL)) +
  geom_point() +
  stat_smooth(method= "lm") 

bmi <- lm(predicted_WL ~ BMI, data = Field)

cor(Field$BMI, Field$predicted_WL, use = "complete.obs")

summary(bmi)

confint(bmi)



## ---------------------------------------------------------------------------------------------------
# load predicting parasite model
predict_parasite <- readRDS("R/Models/predict_Eimeria.rds")

Field_parasite <- Field %>%
  dplyr::select(all_of(Genes_v), eimeriaSpecies) %>%
  dplyr::filter(!eimeriaSpecies == "NA") %>%
   dplyr::filter(!eimeriaSpecies == "E_falciformis")


# rename to match the model
Field_parasite <- Field_parasite %>%
  dplyr::rename(current_infection = eimeriaSpecies)

# current infection should be a factor
Field_parasite$current_infection <- as.factor(Field_parasite$current_infection)


#The predict() function in R is used to predict the values based on the input data.

predictions_parasite <- predict(predict_parasite, Field_parasite)

# assign test.data to a new object, so that we can make changes
result_parasite <- Field_parasite

#add the new variable of predictions to the result object
result_parasite <- cbind(result_parasite, predictions_parasite)




## ---------------------------------------------------------------------------------------------------

conf_matrix_parasite <- 
  confusionMatrix(
    result_parasite$predictions_parasite,
    reference = result_parasite$current_infection)

print(conf_matrix_parasite)

conf_matrix_parasite$table

plt <- as.data.frame(conf_matrix_parasite$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))


ggplot(plt, aes(x = Prediction, y = reorder(Reference, desc(Reference)))) +
    geom_tile(aes(fill = Freq), colour = "white") +
    geom_text(aes(label = sprintf("%d", Freq)), vjust = 1) +
    scale_fill_gradient(low = "white", high = "Steelblue")  +
    labs(x = 'Predicted', y = 'Actual', fill = 'Number of observations', 
         title = 'Confusion Matrix') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> confusion_plot

confusion_plot

ggsave(filename = "figures/confusion_matrix_predicting_ferrisi_field.jpeg", 
       plot = confusion_plot, width = 6, height = 4, dpi = 1000)

## ---------------------------------------------------------------------------------------------------
model_MC <- readRDS("R/Models/predict_MC_Eimeria.rds")

set.seed(597)


Field_mc <- Field %>%
  dplyr::select(all_of(Genes_v), MC.Eimeria) %>%
  dplyr::filter(!MC.Eimeria == "NA") 

Field_mc$MC.Eimeria <- as.factor(Field_mc$MC.Eimeria)

#The predict() function in R is used to predict the values based on the input 
# data.
predictions_MC <- predict(model_MC, Field_mc)


#add the new variable of predictions to the result object
result_MC <- cbind(Field_mc, predictions_MC)


## ---------------------------------------------------------------------------------------------------

conf_matrix_MC <- 
  confusionMatrix(result_MC$predictions_MC, reference = result_MC$MC.Eimeria)

print(conf_matrix_MC)

conf_matrix_MC$table

plt <- as.data.frame(conf_matrix_MC$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

ggplot(plt, aes(x = Prediction, y =  Reference, fill= Freq)) +
        geom_tile() + geom_text(aes(label=Freq)) +
        scale_fill_gradient(low="white", high="darkturquoise") +
        labs(x = "Predictions",y = "Reference") 



## ---------------------------------------------------------------------------------------------------
Field_tol <- Field %>%
    mutate(tolerance = predicted_WL / delta_ct_cewe_MminusE)


Field_tol <- Field_tol %>%
  filter(!is.na(tolerance), MC.Eimeria == TRUE)

summary(Field_tol$tolerance)

Field_tol <- Field_tol %>%
    filter(tolerance > -5, tolerance < 30)


summary(Field_tol$tolerance)

hist(Field_tol$tolerance)

Field_tol %>%
    ggplot(aes(tolerance)) +
    geom_histogram()

parasiteLoad::getParamBounds("normal", data = Field_tol, response = "tolerance")

x <- Field_tol$tolerance

tryDistrib(x, "normal")
tryDistrib(x, "binomial")
tryDistrib(x, "student")
tryDistrib(x, "weibull")
tryDistrib(x, "weibullshifted")



##All
fitWL_tol <- parasiteLoad::analyse(data = Field_tol,
                        response = "tolerance",
                        model = "normal",
                        group = "Sex")



plot_tolerance_Sex<- bananaPlot(mod = fitWL_tol$H3,
             data = Field_tol,
             response = "tolerance",
             group = "Sex") +
    scale_fill_manual(values = c("blueviolet", "limegreen")) +
  scale_color_manual(values = c("blueviolet", "limegreen")) +
  theme_bw() 


plot_tolerance_Sex

# Create HI bar
HIgradientBar <- ggplot(data.frame(hi = seq(0,1,0.0001)),
                        aes(x=hi, y=1, fill = hi)) +
  geom_tile() +
  theme_void() +
  scale_fill_gradient(low = "blue", high = "red")  + 
  scale_x_continuous(expand=c(.01,0)) + 
  scale_y_continuous(expand=c(0,0)) +
  theme(legend.position = 'none')

plot_grid(plot_WL_Sex, 
          HIgradientBar,
          nrow = 2,
          align = "v",
          axis = "tlr",
          rel_heights = c(13, 1))
plot_WL_Sex


################## hybrid effect
Field <- Field %>%
    mutate(HI_2 = 2*HI*(1-HI), #linearize HI
           tolerance = predicted_WL / delta_ct_cewe_MminusE) 
# tolerance = health impact / infection intensity

i <- Field %>%
    filter(Sex == "M") %>%
    drop_na(tolerance)

cor(i$HI, i$tolerance, method = "spearman")
cor(i$HI, i$predicted_WL, method = "spearman")

i <- Field %>%
    filter(Sex == "F") %>%
    drop_na(tolerance)

cor(i$HI, i$tolerance, method = "spearman")
cor(i$HI, i$predicted_WL, method = "spearman")


ggplot(Field, aes(x = HI, HI_2)) +
    geom_point() +
    geom_line()

ggplot(Field, aes(x = HI_2, predicted_WL, color = Sex)) +
    geom_smooth(method = lm, se = TRUE) 

ggplot(Field, aes(x = HI_2, tolerance, color = Sex)) +
    geom_jitter() +
    geom_smooth(method = lm, se = TRUE) 


ggplot(Field, aes(x = HI_2, OPG)) +
    geom_point() +
    geom_line()

lm(formula = predicted_WL ~ HI_2 * Sex, data = Field)

df <- Field %>%
    filter(MC.Eimeria == TRUE)

model_tolerance <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
                      data = df)

summary(model_tolerance)

# lm(tolerance ~ hi + hi_2 * Sex)

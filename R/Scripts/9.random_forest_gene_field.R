## ----setup, include=FALSE---------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----libraries, message = FALSE, warnings = FALSE---------------------------------------------------
#install.packages("optimx", version = "2021-10.12") # this package is required for 
#the parasite load package to work
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
# Random Forests models
library(mice) #imputations
library(fitdistrplus) #testing distributions
library(logspline)
library(caret)


## ---------------------------------------------------------------------------------------------------
hm <- read.csv("output_data/2.imputed_MICE_data_set.csv")



## ----summary_stats_field----------------------------------------------------------------------------
Field <- hm %>%
  filter(origin == "Field") %>%
    drop_na(HI)


## ----genes------------------------------------------------------------------------------------------
EqPCR.cols      <- c("delta_ct_cewe_MminusE", "MC.Eimeria", "Ct.Eimeria") #,"Ct.Mus""delta_ct_ilwe_MminusE", )

Genes_wild   <- c("IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10", 
                  "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                  "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                  "TICAM1", "TNF") #, "IL.12", "IRG6")




## ---------------------------------------------------------------------------------------------------
#select the imputed gene columns
gene <-  Field %>%
  dplyr::select(c(Mouse_ID, "IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10", 
                  "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                  "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                  "TICAM1", "TNF"))

genes <- gene %>%
  dplyr::select(-Mouse_ID)

#remove rows with only nas
genes <- genes[,colSums(is.na(genes))<nrow(genes)]

#remove colums with only nas 
genes <- genes[rowSums(is.na(genes)) != ncol(genes), ]

# select the same rows from the gene data
gene <- gene[row.names(genes),]

# select the same rows from the field data
Field <- Field[row.names(genes),]



## ----predicting_field-------------------------------------------------------------------------------


# load predicting weight loss model
weight_loss_predict <- readRDS("r_scripts/models/predict_WL.rds")

set.seed(540)


#The predict() function in R is used to predict the values based on the input data.
predictions_field <- predict(weight_loss_predict, genes)


# assign test.data to a new object, so that we can make changes
result_field <- genes

#add the new variable of predictions to the result object
result_field <- cbind(result_field, predictions_field)

# add it to the field data 
Field <- cbind(Field, predictions_field)

Field <- Field %>%
  dplyr::mutate(predictions_pos = (-1) * predictions_field)



## ---- warning=FALSE, echo=FALSE, message=FALSE------------------------------------------------------

require(devtools)

devtools::install_github("alicebalard/parasiteLoad@v2.0", force = TRUE)

#force = TRUE)

library(parasiteLoad)


## ---------------------------------------------------------------------------------------------------

Field %>% ggplot(aes(x = predictions_field)) +
  geom_histogram(binwidth = 1.5)



## ---------------------------------------------------------------------------------------------------
Field %>%
    ggplot(aes(x = HI , y = predictions_field , color = Sex)) +
    geom_smooth() +
    geom_point()


Field %>%
    ggplot(aes(x = Body_Length , y = predictions_field , color = Sex)) +
    geom_smooth() +
    geom_point()


## ---------------------------------------------------------------------------------------------------

x <- Field$predictions_pos

descdist(data = x, discrete = FALSE)
descdist(data = x, discrete = FALSE, #data is continuous
         boot = 1000)




## ---------------------------------------------------------------------------------------------------
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



parasiteLoad::getParamBounds("normal", data = Field, response = "predictions_pos")


speparam <- c(L1start = -10.510012081,
                     L1LB = -13.335976577,
                     L1UB = -4.233544127,
                     L2start = -10.510012081,
                     L2LB = -13.335976577,
                     L2UB = -4.233544127,
                     alphaStart = 0, alphaLB = -5, alphaUB = 5,
                     myshapeStart = 1, myshapeLB = 1e-9, myshapeUB = 5)

##All
fitWL_Sex <- parasiteLoad::analyse(data = Field,
                        response = "predictions_pos",
                        model = "normal",
                        group = "Sex")

Field$predictions_field

plot_WL_Sex<- bananaPlot(mod = fitWL_Sex$H3,
             data = Field,
             response = "predictions_pos",
             group = "Sex") +
    scale_fill_manual(values = c("blueviolet", "limegreen")) +
  scale_color_manual(values = c("blueviolet", "limegreen")) +
  theme_bw() 

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



## ---------------------------------------------------------------------------------------------------
Field <- Field %>%
  dplyr::rename(WL = predictions_pos)
ggplot(data = Field, aes(x = delta_ct_cewe_MminusE, y = WL)) +
  geom_point() +
  stat_smooth(method= "lm") 

Field2 <- Field %>%
  drop_na(delta_ct_cewe_MminusE)

cor(Field2$WL, Field2$delta_ct_cewe_MminusE)


tolerance <- lm(WL ~  delta_ct_cewe_MminusE, data = Field)


summary(tolerance)

confint(tolerance)



## ---------------------------------------------------------------------------------------------------
ggplot(data = Field, aes(x = OPG, y = WL)) +
  geom_point() +
  stat_smooth(method= "lm") +
  scale_x_log10()

Field2 <- Field %>%
  drop_na(OPG)

cor(Field2$WL, Field2$OPG)


tolerance <- lm(WL ~  OPG, data = Field)


summary(tolerance)

confint(tolerance)



## ---------------------------------------------------------------------------------------------------

tolerance <- lm(WL ~  OPG * delta_ct_cewe_MminusE, data = Field)


summary(tolerance)

confint(tolerance)



## ---------------------------------------------------------------------------------------------------
Field <- Field %>%
  dplyr::mutate(BMI = Body_Weight / (Body_Length)) #^2) which is the correct
# way to calculatebmi?

ggplot(data = Field, aes(x = BMI, y = predictions_field)) +
  geom_point() +
  stat_smooth(method= "lm") 

bmi <- lm(WL ~ BMI, data = Field)

cor(Field$BMI, Field$WL, use = "complete.obs")

summary(bmi)

confint(bmi)



## ---------------------------------------------------------------------------------------------------
# load predicting parasite model
predict_parasite <- readRDS("r_scripts/models/predict_infecting_parasite.rds")

Field_parasite <- Field %>%
  dplyr::select(all_of(Genes_wild), eimeriaSpecies) %>%
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

ggplot(plt, aes(x = Prediction, y =  Reference, fill= Freq)) +
        geom_tile() + geom_text(aes(label=Freq)) +
        scale_fill_gradient(low="white", high="darkturquoise") +
        labs(x = "Predictions",y = "Reference") 



## ---------------------------------------------------------------------------------------------------
model_MC <- readRDS("r_scripts/models/predict_MC.rds")

set.seed(597)


Field_mc <- Field %>%
  dplyr::select(all_of(Genes_wild), MC.Eimeria) %>%
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
Field <- Field %>%
  mutate(tolerance = WL / delta_ct_cewe_MminusE)

Field$tolerance

Field_tol <- Field %>%
  filter(!is.na(tolerance))

Field_tol <- Field_tol[-37, ]


hist(Field_tol$tolerance)

parasiteLoad::getParamBounds("normal", data = Field_tol, response = "tolerance")


speparam <- c(L1start = 10,
                     L1LB = 1e-9,
                     L1UB = 20,
                     L2start = 10,
                     L2LB = 1e-9,
                     L2UB = 20,
                     alphaStart = 0, alphaLB = -5, alphaUB = 5,
                     myshapeStart = 1, myshapeLB = 1e-9, myshapeUB = 5)

##All
fitWL_Sex <- parasiteLoad::analyse(data = Field_tol,
                        response = "tolerance",
                        model = "normal",
                        group = "Sex")



plot_WL_Sex<- bananaPlot(mod = fitWL_Sex$H3,
             data = Field_tol,
             response = "tolerance",
             group = "Sex") +
    scale_fill_manual(values = c("blueviolet", "limegreen")) +
  scale_color_manual(values = c("blueviolet", "limegreen")) +
  theme_bw() 

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

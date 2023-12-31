## ----setup, include=FALSE---------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---------------------------------------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(optimx)


## ---------------------------------------------------------------------------------------------------
hm <- read.csv("output_data/2.imputed_MICE_data_set.csv")



## ---------------------------------------------------------------------------------------------------
Gene_lab   <- c("IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF") #"IL.12", "IRG6")



Genes_wild   <- c("IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10", 
                  "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                  "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                  "TICAM1", "TNF") #, "IL.12", "IRG6")



## ---------------------------------------------------------------------------------------------------
# Selecting the field samples

field <- hm %>%
  dplyr::filter(origin == "Field") 

field <- unique(field)

#make a factor out of the melting curves (important for later visualization)
field <- field %>%
  dplyr::mutate(MC.Eimeria = as.factor(MC.Eimeria))

genes_mouse <- field %>% 
  dplyr::select(all_of(Genes_wild))

genes <- genes_mouse

#remove rows with only nas
genes <- genes[,colSums(is.na(genes))<nrow(genes)]

#remove colums with only nas 
genes <- genes[rowSums(is.na(genes)) != ncol(genes), ]


##select same rows in the first table
field <- field[row.names(genes), ]





## ---- warning=FALSE, echo=FALSE, message=FALSE------------------------------------------------------
require(devtools)

devtools::install_github("alicebalard/parasiteLoad@v2.0", force = TRUE)

#force = TRUE)

library(parasiteLoad)


## ---------------------------------------------------------------------------------------------------
x <- field$IDO1



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



## ----IFNy-------------------------------------------------------------------------------------------
# remove NA in HI
field <- field %>% 
  drop_na(HI)

field$Sex <- as.factor(field$Sex)

parasiteLoad::getParamBounds("weibull", data = field, response = "IFNy")



IFNy <- parasiteLoad::analyse(data = field,
                        response = "IFNy",
                        model = "weibull",
                        group = "Sex")


parasiteLoad::analyse(data = field,
                        response = "IFNy",
                        model = "weibull",
                        group = "Sex")

bananaPlot(mod = IFNy$H0,
             data = field,
             response = "IFNy",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----CXCR3 gene-------------------------------------------------------------------------------------
CXCR3 <- parasiteLoad::analyse(data = field,
                        response = "CXCR3",
                        model = "weibull",
                        group = "Sex")



bananaPlot(mod = CXCR3$H0,
             data = field,
             response = "CXCR3",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----IL.6-------------------------------------------------------------------------------------------

IL.6 <- parasiteLoad::analyse(data = field,
                        response = "IL.6",
                        model = "weibull",
                        group = "Sex")
##All
print(IL.6)


bananaPlot(mod = IL.6$H0,
             data = field,
             response = "IL.6",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----IL.10------------------------------------------------------------------------------------------

IL.10 <- parasiteLoad::analyse(data = field,
                        response = "IL.10",
                        model = "weibull",
                        group = "Sex")
##All
print(IL.10)


bananaPlot(mod = IL.10$H0,
             data = field,
             response = "IL.10",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----IL.13------------------------------------------------------------------------------------------


IL.13 <- parasiteLoad::analyse(data = field,
                        response = "IL.13",
                        model = "weibull",
                        group = "Sex")
##All
print(IL.13)


bananaPlot(mod = IL.13$H0,
             data = field,
             response = "IL.13",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----IL1RN------------------------------------------------------------------------------------------

IL1RN <- parasiteLoad::analyse(data = field,
                        response = "IL1RN",
                        model = "weibull",
                        group = "Sex")
##All
print(IL1RN)


bananaPlot(mod = IL1RN$H0,
             data = field,
             response = "IL1RN",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----CXCR3------------------------------------------------------------------------------------------

CXCR3 <- parasiteLoad::analyse(data = field,
                        response = "CXCR3",
                        model = "weibull",
                        group = "Sex")
##All
print(CXCR3)


bananaPlot(mod = CXCR3$H0,
             data = field,
             response = "CXCR3",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----CASP1------------------------------------------------------------------------------------------

CASP1 <- parasiteLoad::analyse(data = field,
                        response = "CASP1",
                        model = "weibull",
                        group = "Sex")
##All
print(CASP1)


bananaPlot(mod = CASP1$H0,
             data = field,
             response = "CASP1",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----CXCL9------------------------------------------------------------------------------------------
CXCL9 <- parasiteLoad::analyse(data = field,
                        response = "CXCL9",
                        model = "weibull",
                        group = "Sex")
##All
print(CXCL9)


bananaPlot(mod = CXCL9$H0,
             data = field,
             response = "CXCL9",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----IDO1-------------------------------------------------------------------------------------------
IDO1 <- parasiteLoad::analyse(data = field,
                        response = "IDO1",
                        model = "weibull",
                        group = "Sex")
##All
print(IDO1)


bananaPlot(mod = IDO1$H0,
             data = field,
             response = "IDO1",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----IRGM1------------------------------------------------------------------------------------------

IRGM1 <- parasiteLoad::analyse(data = field,
                        response = "IRGM1",
                        model = "weibull",
                        group = "Sex")
##All
print(IRGM1)


bananaPlot(mod = IRGM1$H0,
             data = field,
             response = "IRGM1",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----MPO--------------------------------------------------------------------------------------------

MPO <- parasiteLoad::analyse(data = field,
                        response = "MPO",
                        model = "weibull",
                        group = "Sex")
##All
print(MPO)


bananaPlot(mod = MPO$H0,
             data = field,
             response = "MPO",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----MUC2-------------------------------------------------------------------------------------------


MUC2 <- parasiteLoad::analyse(data = field,
                        response = "MUC2",
                        model = "weibull",
                        group = "Sex")
##All
print(MUC2)


bananaPlot(mod = MUC2$H0,
             data = field,
             response = "MUC2",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----MUC5AC-----------------------------------------------------------------------------------------


MUC5AC <- parasiteLoad::analyse(data = field,
                        response = "MUC5AC",
                        model = "weibull",
                        group = "Sex")
##All
print(MUC5AC)


bananaPlot(mod = MUC5AC$H0,
             data = field,
             response = "MUC5AC",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----MYD88------------------------------------------------------------------------------------------


MYD88 <- parasiteLoad::analyse(data = field,
                        response = "MYD88",
                        model = "weibull",
                        group = "Sex")
##All
print(MYD88)


bananaPlot(mod = MYD88$H0,
             data = field,
             response = "MYD88",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----NCR1-------------------------------------------------------------------------------------------


NCR1 <- parasiteLoad::analyse(data = field,
                        response = "NCR1",
                        model = "weibull",
                        group = "Sex")
##All
print(NCR1)


bananaPlot(mod = NCR1$H0,
             data = field,
             response = "NCR1",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----PRF1-------------------------------------------------------------------------------------------


PRF1 <- parasiteLoad::analyse(data = field,
                        response = "PRF1",
                        model = "weibull",
                        group = "Sex")
##All
print(PRF1)


bananaPlot(mod = PRF1$H0,
             data = field,
             response = "PRF1",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----RETNLB-----------------------------------------------------------------------------------------


RETNLB <- parasiteLoad::analyse(data = field,
                        response = "RETNLB",
                        model = "weibull",
                        group = "Sex")
##All
print(RETNLB)


bananaPlot(mod = RETNLB$H0,
             data = field,
             response = "RETNLB",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----SOCS1------------------------------------------------------------------------------------------

SOCS1 <- parasiteLoad::analyse(data = field,
                        response = "SOCS1",
                        model = "weibull",
                        group = "Sex")
##All
print(SOCS1)


bananaPlot(mod = SOCS1$H0,
             data = field,
             response = "SOCS1",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----TICAM1-----------------------------------------------------------------------------------------


TICAM1 <- parasiteLoad::analyse(data = field,
                        response = "TICAM1",
                        model = "weibull",
                        group = "Sex")
##All
print(TICAM1)


bananaPlot(mod = TICAM1$H0,
             data = field,
             response = "TICAM1",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


## ----TNF--------------------------------------------------------------------------------------------

TNF <- parasiteLoad::analyse(data = field,
                        response = "TNF",
                        model = "weibull",
                        group = "Sex")
##All
print(TNF)


bananaPlot(mod = TNF$H0,
             data = field,
             response = "TNF",
             group = "Sex") +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red")) +
  theme_bw()

  


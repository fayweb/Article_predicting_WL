
#install.packages("optimx", version = "2021-10.12") # this package is required for 
#the parasite load package to work
#require(devtools)

## install the pacakage of Alice Balard
#devtools::install_github("alicebalard/parasiteLoad@v2.0", force = TRUE)
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
library(ggeffects)
library(ggbeeswarm)
library(ggdist)
library(jtools)


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

## ----------------------------------------------------------------------------

# Can the predicted weight loss be predicted by infection intensities
Field2 <- Field %>%
  drop_na(delta_ct_cewe_MminusE) %>%
    filter(MC.Eimeria == TRUE)

ggplot(data = Field2, aes(x = delta_ct_cewe_MminusE, y = predicted_WL)) +
    geom_point() +
    stat_smooth(method= "lm") 


cor(Field2$predicted_WL, Field2$delta_ct_cewe_MminusE)


model_WL <- lm(predicted_WL ~  delta_ct_cewe_MminusE, data = Field2)

summary(model_WL)
confint(model_WL)


model_WL_infection <- 
    ggpredict(model_WL) %>% 
    plot(colors = "darkorange2") +   # Use a refined shade of blue
    labs(title = NULL) +  
    xlab("Infection intensity of Eimeria spp. in the caecum") +
    ylab("Predicted values of weight loss") +
    theme_bw()

model_WL_infection

ggsave(filename = "figures/linear_model_WL_infection_field.jpeg", 
       plot = model_WL_infection, 
       width = 6, height = 4, dpi = 1000)


## ---------------------------------------------------------------------------------------------------
ggplot(data = Field, aes(x = OPG, y = predicted_WL)) +
  geom_point() +
  stat_smooth(method= "lm") +
  scale_x_log10()

Field2 <- Field %>%
  drop_na(OPG)

ggplot(data = Field2, aes(x = OPG, y = predicted_WL)) +
    geom_point() +
    stat_smooth(method= "lm") +
    scale_x_log10()


cor(Field2$predicted_WL, Field2$OPG)
model <- lm(predicted_WL ~  OPG, data = Field)


summary(model)

confint(model)



## ---------------------------------------------------------------------------------------------------

model <- lm(predicted_WL ~  OPG * delta_ct_cewe_MminusE, data = Field)


summary(model)

confint(model)


## ---------------------------------------------------------------------------------------------------
Field <- Field %>%
  dplyr::mutate(BMI = Body_Weight / (Body_Length)) #^2) which is the correct
# way to calculatebmi?

ggplot(data = Field, aes(x = BMI, y = predicted_WL)) +
  geom_point() +
  stat_smooth(method= "lm") 

bmi <- lm(predicted_WL ~ BMI, data = Field)

bmi_plot <- 
    ggpredict(bmi) %>%
    plot(color = "purple") +
    labs(y = "Predicted health impact", x = "BMI") +
    theme(title = element_blank())


cor(Field$BMI, Field$predicted_WL, use = "complete.obs")

summary(bmi)

confint(bmi)

ggsave(filename = "figures/BMI_WL.jpeg", plot = bmi_plot, width = 5, height = 4, dpi = 1000)

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


# Define function to be used to test, get the log lik and aic
tryDistrib <- function(x, distrib){
    # deals with fitdistr error:
    fit <- 
        tryCatch(MASS::fitdistr(x, distrib), error=function(err) "fit failed")
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

tryDistrib(x, "normal") #-208.0489
tryDistrib(x, "binomial")
tryDistrib(x, "student")
tryDistrib(x, "weibull")
tryDistrib(x, "weibullshifted")

Field_tol$Sex <- as.factor(Field_tol$Sex)


##All
fitWL_tol <- parasiteLoad::analyse(data = Field_tol,
                        response = "tolerance",
                        model = "normal",
                        group = "Sex")



plot_tolerance_Sex <- 
    bananaPlot(mod = fitWL_tol$H3,
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

#plot_grid(plot_WL_Sex, 
#          HIgradientBar,
 #         nrow = 2,
  #        align = "v",
   #       axis = "tlr",
    #      rel_heights = c(13, 1))



################## hybrid effect
Field <- Field %>%
    mutate(HE = 2*HI*(1-HI), #linearize HI
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


ggplot(Field, aes(x = HI, HE)) +
    geom_point() +
    geom_line()

ggplot(Field, aes(x = HE, predicted_WL, color = Sex)) +
    geom_smooth(method = lm, se = TRUE) 

ggplot(Field, aes(x = HE, tolerance, color = Sex)) +
    geom_jitter() +
    geom_smooth(method = lm, se = TRUE) 


ggplot(Field, aes(x = HE, OPG)) +
    geom_point() +
    geom_line()

model <- lm(formula = predicted_WL ~ HE * Sex, data = Field)
summary(model)

df <- Field %>%
    filter(MC.Eimeria == TRUE)
###########################################################
#########################################################
model_tolerance <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
                      data = df)

summary(model_tolerance)

model <- lm(tolerance ~ HI + HE * Sex, df)
summary(model)


######################################################################
Field_asp <- Field %>%
    drop_na(Aspiculuris_sp) %>%
    filter(Aspiculuris_sp != 0) %>%
    mutate(tolerance = predicted_WL / Aspiculuris_sp) %>%
    drop_na(tolerance)

model <- lm(tolerance ~ Aspiculuris_sp, data = Field_asp)
summary(model)



model_tolerance <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
                      data = df)

summary(model_tolerance)

model <- lm(tolerance ~ HI + HE * Sex, df)
summary(model)
####################################################################
Field_asp <- Field %>%
    drop_na(Aspiculuris_sp) %>%
    filter(Aspiculuris_sp != 0) %>%
    mutate(tolerance = predicted_WL / Aspiculuris_sp) %>%
    drop_na(tolerance)

model <- lm(tolerance ~ Aspiculuris_sp, data = Field_asp)
summary(model)

#########################################################################
df <- Field %>%
    filter(eimeriaSpecies == "E_falciformis")


model_tolerance <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
                      data = df)

summary(model_tolerance)


#################################
df <- Field %>%
    filter(eimeriaSpecies == "E_ferrisi")


model_tolerance <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
                      data = df)

summary(model_tolerance)

#################
df <- Field %>%
    filter(eimeriaSpecies == "E_falciformis")


model_tolerance <- lm(predicted_WL ~ delta_ct_cewe_MminusE, 
                      data = df)

summary(model_tolerance)


#################################
df <- Field %>%
    drop_na(ILWE_Crypto_Ct) %>%
    filter(ILWE_Crypto_Ct != 0) %>%
    mutate(tolerance = predicted_WL / -ILWE_Crypto_Ct) %>%
    drop_na(tolerance)


model_tolerance <- lm(predicted_WL ~ ILWE_Crypto_Ct, 
                      data = df)

summary(model_tolerance)

###################################
df <- Field %>%
    drop_na(Syphacia_sp) %>%
    filter(Syphacia_sp != 0) %>%
    mutate(tolerance = predicted_WL / Syphacia_sp) %>%
    drop_na(tolerance)


model_tolerance <- lm(predicted_WL ~ Syphacia_sp, 
                      data = df)

summary(model_tolerance)

#####################################
###################################
df <- Field %>%
    drop_na(Heterakis_sp) %>%
    filter(Heterakis_sp != 0) %>%
    mutate(tolerance = predicted_WL / Heterakis_sp) %>%
    drop_na(tolerance)


model_tolerance <- lm(predicted_WL ~ Heterakis_sp, 
                      data = df)

summary(model_tolerance)


################################
df <- Field %>%
    drop_na(Trichuris_muris) %>%
    filter(Trichuris_muris != 0) %>%
    mutate(tolerance = predicted_WL / Trichuris_muris) %>%
    drop_na(tolerance)


model_tolerance <- lm(predicted_WL ~ Trichuris_muris, 
                      data = df)

summary(model_tolerance)

##################################################################
##
Field <- Field %>%
    mutate(infection_intensity_Eim = delta_ct_cewe_MminusE)
Field_par <- Field %>%
    mutate(
        infected_Aspiculuris = Aspiculuris_sp != 0,
        infected_syphasia = Syphacia_sp != 0,
        infected_crypto = ILWE_Crypto_Ct != 0
    )

Field_par <- Field_par %>%
    mutate(
        infected_Aspiculuris = as.factor(infected_Aspiculuris),
        infected_syphasia = as.factor(infected_syphasia),
        infected_crypto = as.factor(infected_crypto)
    )

model1 <- lm(predicted_WL ~ MC.Eimeria + infected_Aspiculuris + HI + HE
            + infected_syphasia + infected_crypto, data = Field_par)

summary(model1)

plot_summs(model1, plot.distributions = TRUE, robust = TRUE, scale = TRUE,
           colors = "mediumblue") -> plot1
plot1

ggsave(filename = "figures/coefficient_plot_model1.jpeg", plot = plot1, 
       width = 5, height = 4, dpi = 300)

#######################################
model2 <- lm(predicted_WL ~ MC.Eimeria * HE, Field)
summary(model2)
plot_summs(model2, plot.distributions = TRUE) 


## hybrid index + infectopm
Field$MC.Eimeria <- as.factor(Field$MC.Eimeria)
model3 <- lm(predicted_WL ~ MC.Eimeria *infection_intensity_Eim * HE + 
                HI + HE, Field)
summary(model3)

plot_summs(model3, plot.distributions = TRUE, robust = TRUE, scale = TRUE,
           colors = "mediumblue") -> plot2

plot2 

ggsave(filename = "figures/coefficient_plot_model2.jpeg", plot = plot2, 
       width = 5, height = 4, dpi = 300)

ggpredict(model3, terms = c("MC.Eimeria", "delta_ct_cewe_MminusE")) %>%
    plot()

## hybrid index + infectopm
model4 <- lm(predicted_WL ~  HI + HE, Field)
summary(model4)
plot_summs(model4,  plot.distributions = TRUE, robust = TRUE, scale = TRUE,
           colors = "mediumblue")  -> plot_4

ggsave(filename = "figures/coefficient_plot_HE.jpeg", plot = plot_4, 
       width = 5, height = 4, dpi = 300)

plot_summs(model3, model4,  robust = TRUE, 
           scale = TRUE) -> model1_2
model1_2

ggsave(filename = "figures/coefficient_plot_model1_2.jpeg", plot = model1_2, 
       width = 8, height = 6, dpi = 300)

model4 <- lm(predicted_WL ~  MC.Eimeria, Field)
summary(model4)









##################################
#raincloud plots
ggplot(Field, aes(x = MC.Eimeria, y = predicted_WL, color = MC.Eimeria)) +
    geom_violin(trim = FALSE, alpha = 0.5) + 
    geom_quasirandom(aes(color = MC.Eimeria), size = 1, alpha = 0.8) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
    scale_color_brewer(palette = "Dark2") +
    theme_minimal() +
    labs(title = "Raincloud Plot of predicted_WL by MC.Eimeria Group",
         x = "MC.Eimeria Group",
         y = "Predicted WL")

# Filter out NA values in MC.Eimeria
Field_filtered <- Field %>% filter(!is.na(MC.Eimeria))




# Define colors
colors <- c("TRUE" = "forestgreen", "FALSE" = "purple")

ggplot(Field_filtered, aes(y = MC.Eimeria, x = predicted_WL, fill = MC.Eimeria)) + 
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        alpha = 0.5,
        .width = 0, 
        justification = -.2, 
        point_colour = NA,
        orientation = "y"  # Set orientation to y
    ) + 
    scale_fill_manual(values = colors) +
    geom_boxplot(
        width = .15, 
        outlier.shape = NA,
        orientation = "y"  # Set orientation to y
    ) +
    stat_dots(
        # ploting on left side
        side = "left",
        # adjusting position
        justification = 1.1,
        # adjust grouping (binning) of observations
        binwidth = 0.25,
        alpha = 0.5) +
    geom_point(
        shape = 95,
        size = 15,
        alpha = .2,
        color = "gray50",
        position = position_dodge(width = 0.75)
    ) +
    coord_cartesian(ylim = c(1.2, 2.9), clip = "off") +
    theme_minimal() +
    labs(y = "Infection status with Eimerai spp.", 
         x = "Predicted detrimental immune signature") -> raincloud_plots__eimeria

raincloud_plots__eimeria

ggsave(plot = raincloud_plots__eimeria, filename = "figures/raincloud_eimeria.jpeg", 
       width = 6, 
       height = 4, dpi = 1000)

###############################################
########################
# raincloud plots parasites
model_par <- lm(predicted_WL ~ ILWE_Crypto_Ct + Aspiculuris_sp + Syphacia_sp, data = Field)
summary(model_par)

Field_par <- Field_par %>%
    pivot_longer(cols = c("infected_Aspiculuris", "infected_syphasia", 
                          "infected_crypto"), names_to = "parasite",
                 values_to = "infection_status")  %>% 
    drop_na(infection_status)


ggplot(Field_par, aes(y = infection_status, x = predicted_WL, 
                  fill = infection_status)) + 
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        alpha = 0.7,
        .width = 0, 
        justification = -.2, 
        point_colour = NA,
        orientation = "y"  # Set orientation to y
    ) + 
    scale_fill_manual(values = colors) +
    geom_boxplot(
        width = .15, 
        outlier.shape = NA,
        orientation = "y"  # Set orientation to y
    ) +
    geom_point(
        shape = 95,
        size = 15,
        alpha = .2,
        color = "gray50",
        position = position_dodge(width = 0.75)
    ) +
    stat_dots(
        # ploting on left side
        side = "left",
        # adjusting position
        justification = 1.1,
        # adjust grouping (binning) of observations
        binwidth = 0.25) +
    coord_cartesian(ylim = c(1.2, 2.9), clip = "off") +
    theme_minimal() +
    labs(y = "Infection status with Eimerai spp.", 
         x = "Predicted detrimental immune signature") +
    facet_wrap(~parasite) -> parasites

parasites

ggsave(plot = parasites, filename = "figures/raincloud_parasites.jpeg", 
       width = 10, 
       height = 4, dpi = 1000)

###############################################################
################################################

model <- lm(predicted_WL ~ eimeriaSpecies, data = Field)
summary(model)



amplicon <- read.csv("https://raw.githubusercontent.com/ferreira-scm/Eimeria_AmpSeq/master/data/Wild/Sample_selection_Metabarcoding_Complete.csv")


glimpse(amplicon)
glimpse(Field)


eimer_sp <- amplicon %>%
    dplyr::select(c(Mouse_ID, Species)) %>%
    rename(amplicon_species = Species) %>%
    drop_na(amplicon_species)

eimer_sp$Mouse_ID <- gsub(pattern = "_", replacement = "", x = eimer_sp$Mouse_ID)


Field <- Field %>%
    left_join(eimer_sp, by = "Mouse_ID")




# create new variable depending on amplicon sequencing
Field <- Field %>%
    mutate(species_Eimeria = case_when(
        is.na(eimeriaSpecies) ~ amplicon_species,
        !is.na(eimeriaSpecies) ~ eimeriaSpecies))

Field$species_Eimeria <- gsub(pattern = "Negative", replacement = "uninfected", 
                              x = Field$species_Eimeria)

Field <- Field %>%
    mutate(infection_status = case_when(
        is.na(MC.Eimeria) & species_Eimeria == "uninfected" ~ "FALSE",
        is.na(MC.Eimeria) & species_Eimeria == "E_falciformis" ~ "TRUE",
        is.na(MC.Eimeria) & species_Eimeria == "E_ferrisi" ~ "TRUE",
        !is.na(MC.Eimeria) ~ MC.Eimeria
    ))

Field$infection_status <- as.factor(Field$infection_status)
############

infected <- parasiteLoad::analyse(data = Field,
                                            response = "predicted_WL",
                                            model = "normal",
                                            group = "infection_status")



plot_infected <- 
    bananaPlot(mod = infected$H3,
               data = Field,
               response = "predicted_WL",
               group = "infection_status") +
    scale_fill_manual(values = c("blueviolet", "limegreen")) +
    scale_color_manual(values = c("blueviolet", "limegreen")) +
    theme_bw() 

plot_infected

##############
model <- lm(predicted_WL ~ species_Eimeria, data = Field)
summary(model)


ggpredict(model = model) %>%
    plot()



##All
only_inf <- Field %>%
    dplyr::filter(species_Eimeria == c("E_falciformis", "E_ferrisi"))

only_inf$species_Eimeria <- as.factor(only_inf$species_Eimeria)

fitWL_species_Diff <- parasiteLoad::analyse(data = only_inf,
                                   response = "predicted_WL",
                                   model = "normal",
                                   group = "species_Eimeria")



plot_species <- 
    bananaPlot(mod = fitWL_species_Diff$H3,
               data = only_inf,
               response = "predicted_WL",
               group = "species_Eimeria") +
    scale_fill_manual(values = c("blueviolet", "limegreen")) +
    scale_color_manual(values = c("blueviolet", "limegreen")) +
    theme_bw() 


plot_species

A <- Field %>%
    dplyr::select(c(Mouse_ID, infection_status, species_Eimeria))


write.csv(A, "tables/table_for_saskia.csv", row.names = FALSE)


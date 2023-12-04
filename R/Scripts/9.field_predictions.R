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
library(leaflet)
library(webshot)
library(htmlwidgets)
library(cowplot)
library(gridExtra)
library(magick)

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
               "TICAM1", "TNF") #"IL.12", "IRG6")

# select the gene columns
gene <-  Field %>%
  dplyr::select(c(Mouse_ID, all_of(Genes_v)))

# data frame with only the genes
genes <- gene %>%
  dplyr::select(-Mouse_ID)


# load predicting weight loss model
weight_loss_predict <- readRDS("R/Models/predict_WL.rds")

set.seed(540)


#The predict() function in R is used to predict the values 
#based on the input data.
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

## 
normal_ <- fitdist(x, "norm")
weibull_ <- fitdist(x, "weibull")
gamma_ <- fitdist(x, "gamma")


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


## Now fit the distributions to the predicted weight loss data
tryDistrib(x, "normal") # -782.4131
tryDistrib(x, "binomial") #-784.7632
tryDistrib(x, "student")
tryDistrib(x, "weibull")
tryDistrib(x, "weibullshifted")


## Compare again between normal and weibull
findGoodDist(x, "normal", "weibull")


## plot the distributions
plot(normal_)
summary(normal_)
plot(gamma_)
summary(gamma_)
plot(weibull_)
summary(weibull_)


# Testing differences between female and male hybrids of M.m. musculus and 
#m.m.domesticus in predicted weight loss
Field$Sex <- as.factor(Field$Sex)


##All
fitWL_Sex <- parasiteLoad::analyse(data = Field,
                        response = "predicted_WL",
                        model = "normal",
                        group = "Sex")


plot_WL_Sex<- bananaPlot(mod = fitWL_Sex$H3,
             data = Field,
             response = "predicted_WL",
             group = "Sex",
  cols = c("white", "white")) +
    scale_fill_manual(values = c("orange", "forestgreen")) +
  scale_color_manual(values = c("orange", "forestgreen")) +
  theme_bw() +
    theme(legend.position="none",
         axis.title.x=element_blank(),
          axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_Sex

ggsave(plot = plot_WL_Sex, filename = "figures/hybrid_sex.jpeg", width = 10, 
       height = 8, dpi = 1000)


# Adjust the gradient bar plot to include axis labels and remove space
HIgradientBar <- ggplot(data.frame(hi = seq(0,1,0.0001)), 
                        aes(x=hi, y=1, fill = hi)) +
    geom_tile() +
    scale_x_continuous(breaks=seq(0, 1, by=0.25), 
                       labels=c("0", "0.25", "0.5", "0.75", "1")) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_void() +
    theme(legend.position = 'none', 
          plot.margin = unit(c(-1, 0, 0, 0), "npc"), 
          # This removes space around the plot
          axis.text.x = element_text(color = "black", 
                                     angle = 0, vjust = 0.5, hjust=0.5)) 
        # Adjust text vjust for positioning

HIgradientBar

# Create the combined plot with the gradient bar as the "axis"
plot_WL_Sex_combined <- 
    plot_grid(plot_WL_Sex ,
              HIgradientBar,  
              nrow = 2, 
              rel_heights = c(1.3, 1/8),
              align = "hv",
              axis = "tb",
              vjust = c(-1,2),
              scale = 1.1) 
   
# Display the combined plot
plot_WL_Sex_combined


ggsave(plot = plot_WL_Sex_combined, 
       filename = "figures/hybrid_sex.jpeg", width = 10, 
       height = 8, dpi = 1000)


####################### Mapping ######################################
leaflet(data = Field) %>%
    addTiles() %>%
    addMarkers(
        lng = ~Longitude, lat = ~Latitude, popup = ~as.character(Mouse_ID))

colorPalette <- colorRampPalette(c("blue", "red"))


colors <- colorPalette(100)[as.numeric(cut(Field$HI, breaks = 100))]

leaflet_map <-
    leaflet(data = Field) %>%
    addTiles() %>%
    addCircleMarkers(lng = ~Longitude, lat = ~Latitude, color = ~colors, 
                     radius = 5, fillOpacity = 0.8, stroke = FALSE, 
                     popup = ~as.character(HI))


#Read the Leaflet map image
leaflet_image <- magick::image_read("figures/Hybrid_map.jpeg")

# Convert to a raster for grid plotting
leaflet_raster <- rasterGrob(leaflet_image, interpolate = TRUE)

# Combine the ggplot and raster image
combined_plot <- grid.arrange(
    leaflet_raster,
    plot_WL_Sex_combined,
    ncol = 2, # Set the number of columns to 2 for horizontal alignment
    widths = c(1, 1))

# Add the annotations to the grob
combined_grob <- arrangeGrob(
    grobs = list(combined_plot,
                 textGrob("A", x = unit(0.1, "npc"), y = unit(0.95, "npc"), 
                          gp = gpar(fontface = "bold", cex = 1.5)),
                 textGrob("B", x = unit(0.51, "npc"), y = unit(0.95, "npc"), 
                          gp = gpar(fontface = "bold", cex = 1.5))),
    ncol = 3
)


# Add annotations
grid.text("A", x = unit(0.1, "npc"), y = unit(0.95, "npc"), 
          gp = gpar(fontface = "bold", cex = 1.5))
grid.text("B", x = unit(0.51, "npc"), y = unit(0.95, "npc"), 
          gp = gpar(fontface = "bold", cex = 1.5))


ggsave(plot = combined_plot, 
       filename = "figure_panels/banana_map_immune_signature.jpeg", width = 16, 
       height = 8, dpi = 1000)


##############################################################################
### Testing diffences between infected and uninfected hybrid mice

# Melting Curve analysis
Field_mc <- Field %>%
    drop_na(MC.Eimeria)

Field_mc$MC.Eimeria <- as.factor(Field_mc$MC.Eimeria)

fitWL_mc <- parasiteLoad::analyse(data = Field_mc,
                                   response = "predicted_WL",
                                   model = "normal",
                                   group = "MC.Eimeria")


plot_WL_mc <- 
    bananaPlot(mod = fitWL_mc$H3,
                         data = Field_mc,
                         response = "predicted_WL",
                         group = "MC.Eimeria",
                         cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3"), 
                       name = "Melting Curve analysis") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Melting Curve analysis") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_mc

# Create the combined plot with the gradient bar as the "axis"
plot_WL_mc_combined <- 
    plot_grid(plot_WL_mc ,
              HIgradientBar,  
              nrow = 2, 
              rel_heights = c(1.3, 1/8),
              align = "hv",
              axis = "tb",
              vjust = c(-1,2),
              scale = 1) 

# Display the combined plot
plot_WL_mc_combined


ggsave(plot = plot_WL_mc_combined, 
       filename = "figures/hybrid_infected.jpeg", width = 10, 
       height = 8, dpi = 1000)



############ Testing according to delta ct (not recommended, just validating)
Field_ct <-  Field %>%
    drop_na(delta_ct_cewe_MminusE) %>%
    dplyr::mutate(delta_infection = 
                      case_when(delta_ct_cewe_MminusE < -5 ~ "uninfected"),delta_ct_cewe_MminusE > -5 ~ "infected")

Field_ct$delta_infection <- as.factor(Field_ct$delta_infection)

##All
fitWL_ct<- parasiteLoad::analyse(data = Field_ct,
                                  response = "predicted_WL",
                                  model = "normal",
                                  group = "delta_infection")


plot_WL_ct <- 
    bananaPlot(mod = fitWL_ct$H3,
               data = Field_ct,
               response = "predicted_WL",
               group = "delta_infection",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("indianred3","steelblue1"),
                      name = "qPCR values - Eimeria detection in Caecum") +
    scale_color_manual(values = c("indianred3", "steelblue1"),
                       name = "qPCR values - Eimeria detection in Caecum") +
    theme_bw() +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_ct


# Create the combined plot with the gradient bar as the "axis"
combined_plot <- 
    plot_grid(plot_WL_ct ,
              HIgradientBar,  
              nrow = 2, 
              rel_heights = c(1.3, 1/8),
              align = "hv",
              axis = "tb",
              vjust = c(-1,2),
              scale = 1) 

# Display the combined plot
combined_plot


ggsave(plot = plot_WL_mc_combined, 
       filename = "figures/hybrid_infected.jpeg", width = 10, 
       height = 8, dpi = 1000)


###################################
#################################


# cryüp
Field <- Field %>%
    mutate(crypto_infected = 
               case_when(
                   Oocyst_Predict_Crypto == 0 ~ "uninfected",
                   Oocyst_Predict_Crypto > 0 ~ "infected",
                   ILWE_Crypto_Ct == NA ~ NA
               ))

Field_crypto <- Field %>%
    drop_na(crypto_infected)

Field_crypto$crypto_infected <- as.factor(Field_crypto$crypto_infected)

##All
fitWL_crypto <- parasiteLoad::analyse(data = Field_crypto,
                                     response = "predicted_WL",
                                     model = "normal",
                                     group = "crypto_infected")


# plot it
plot_WL_crypto <- 
    bananaPlot(mod = fitWL_crypto$H3,
               data = Field_crypto,
               response = "predicted_WL",
               group = "crypto_infected",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3")) +
    scale_color_manual(values = c("steelblue1", "indianred3")) +
    theme_bw()  +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_mc

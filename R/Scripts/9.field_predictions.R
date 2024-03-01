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
library(leaflet)
library(webshot)
library(htmlwidgets)
library(cowplot)
library(gridExtra)
library(magick)
library(grid) 
library(patchwork)

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
    labs(y = "Predicted detrimental health impact, 
         Immune signature") 

plot_WL_Sex

### create the hybrid bar

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

#HIgradientBar

plot_WL_Sex <- plot_WL_Sex +
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
         plot.margin = unit(c(0, 0, 0, 0), "cm"))

#plot_WL_Sex
#
ggsave(plot = plot_WL_Sex, filename = "figures/hybrid_sex.jpeg", width = 6, 
       height = 6, dpi = 1000)

# Ensure the bottom plot (HIgradientBar) is ready (assuming HIgradientBar is already defined)
HIgradientBar <- HIgradientBar + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Use patchwork to combine the plots without any space between them
combined_plot <- plot_WL_Sex / HIgradientBar + 
    plot_layout(heights = c(1, 0.1)) # Adjust the relative heights as needed

# Print the combined plot
#combined_plot


ggsave(plot = combined_plot, 
       filename = "figures/hybrid_sex.jpeg", width = 4, height = 4, dpi = 1000)


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

leaflet_raster <- ggplot() +
    background_image(leaflet_image) + coord_fixed()
###################################



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
    scale_fill_manual(values = c("forestgreen", "purple"), 
                       name = "Melting Curve analysis") +
    scale_color_manual(values = c("forestgreen", "purple"),
                       name = "Melting Curve analysis") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

#plot_WL_mc

# Create the combined plot with the gradient bar as the "axis"
plot_WL_mc_combined <-  (plot_WL_mc / HIgradientBar) + 
    plot_layout(heights = c(1, 0.1)) # Adjust the relative heights as needed


# Display the combined plot
#plot_WL_mc_combined


ggsave(plot = plot_WL_mc_combined, 
       filename = "figures/hybrid_meltingcurve.jpeg",  
       width = 7, height = 7, dpi = 1000)
####################
#######################
######################
# Use patchwork to combine the plots without any space between them
map_plot <- (combined_plot | plot_WL_mc_combined) +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
map_plot <- map_plot + 
    plot_annotation(title = 'Fig. 8', 
                    theme = theme(plot.title = element_text(size = 13, hjust = 0),
                                  legend.position = "none"))

# Display the panel figure
#print(map_plot)


ggsave(plot = map_plot, 
       filename = "figure_panels/banana_map_immune_signature.jpeg", width = 10, 
       height = 5, dpi = 1000)

############ Testing according to delta ct (not recommended, just validating)
Field_ct <-  Field %>%
    drop_na(delta_ct_cewe_MminusE) %>%
    dplyr::mutate(delta_infection = 
                      case_when(delta_ct_cewe_MminusE < -5 ~ "uninfected",
                                delta_ct_cewe_MminusE > -5 ~ "infected"))

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


ggsave(plot = combined_plot, 
       filename = "figures/hybrid_infected_ct.jpeg",  width = 6, height = 6, dpi = 1000)


##################### Eimeria oocysts
Field_ooc_eim <-  Field %>%
    drop_na(OPG) %>%
    dplyr::mutate(oocysts_present = 
                      case_when(OPG > 0 ~ "TRUE",
                                OPG == 0 ~ "FALSE"))

Field_ooc_eim$oocysts_present <- as.factor(Field_ooc_eim$oocysts_present)

##All
fitWL_ooc<- parasiteLoad::analyse(data = Field_ooc_eim,
                                 response = "predicted_WL",
                                 model = "normal",
                                 group = "oocysts_present")


plot_WL_ooc <- 
    bananaPlot(mod = fitWL_ooc$H3,
               data = Field_ooc_eim,
               response = "predicted_WL",
               group = "oocysts_present",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1","indianred3"),
                      name = "Eimeria oocysts present") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Eimeria oocysts present") +
    theme_bw() +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_ooc


# Create the combined plot with the gradient bar as the "axis"
combined_plot_ooc <- 
    plot_grid(plot_WL_ooc ,
              HIgradientBar,  
              nrow = 2, 
              rel_heights = c(1.3, 1/8),
              align = "hv",
              axis = "tb",
              vjust = c(-1,2),
              scale = 1) 

# Display the combined plot
combined_plot_ooc


ggsave(plot = combined_plot_ooc, 
       filename = "figures/hybrid_eimeria_ooc.jpeg",  width = 6, height = 6, dpi = 1000)

###################################


####################################################
#####  crpyto ct
Field_crypto <- Field %>%
    drop_na(ILWE_Crypto_Ct) %>%
    mutate(crypto_infected_ct = 
               case_when(
                   ILWE_Crypto_Ct > 0 ~ "TRUE",
                   ILWE_Crypto_Ct == 0 ~ "FALSE",
                   ILWE_Crypto_Ct == NA ~ NA
               ))

Field_crypto <- Field_crypto %>%
    drop_na(crypto_infected)

Field_crypto$crypto_infected_ct <- as.factor(Field_crypto$crypto_infected_ct)



fitWL_crypto_ct <- parasiteLoad::analyse(data = Field_crypto,
                                      response = "predicted_WL",
                                      model = "normal",
                                      group = "crypto_infected_ct")


# plot it
plot_WL_crypto_ct <- 
    bananaPlot(mod = fitWL_crypto_ct$H3,
               data = Field_crypto,
               response = "predicted_WL",
               group = "crypto_infected_ct",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3"),
                      name = "Cryptosporidum detected in qPCR") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Cryptosporidum detected in qPCR") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_crypto_ct

# Create the combined plot with the gradient bar as the "axis"
plot_WL_crypto_Ct <- 
    plot_grid(plot_WL_crypto_ct ,
              HIgradientBar,  
              nrow = 2, 
              rel_heights = c(1.3, 1/8),
              align = "hv",
              axis = "tb",
              vjust = c(-1,2),
              scale = 1) 

# Display the combined plot
plot_WL_crypto_Ct


ggsave(plot = plot_WL_crypto_Ct, 
       filename = "figures/hybrid_crypto_qpcr.jpeg",  width = 6, height = 6, dpi = 1000)

rm(Field_crypto, Field_ct, Field_mc, gamma_, normal_, weibull_)

################################################
##############
######Worms?
##################First aspiculuris
Field %>%
    pivot_longer(cols = c("Trichuris_muris", "Mastophorus_muris", 
                          "Catenotaenia_pusilla", "Heligmosomoides_polygurus",
                          "Heterakis_sp", "Aspiculuris_sp", "Syphacia_sp",
                          "Taenia_sp", "Hymenolepis_sp"), names_to = "Worms", 
                 values_to = "Counts_worms") %>%
    ggplot(aes(x = Worms, y = Counts_worms)) +
    geom_violin()




#####  Aspiculuris_sp
wormy_field <- Field %>%
    mutate(wormy_asp = 
               case_when(
                   Aspiculuris_sp > 0 ~ "TRUE",
                   Aspiculuris_sp == 0 ~ "FALSE",
                   Aspiculuris_sp == NA ~ NA
               ))

wormy_field <- wormy_field %>%
    drop_na(wormy_asp)

wormy_field$wormy_asp <- as.factor(wormy_field$wormy_asp)

fitWL_worm <- parasiteLoad::analyse(data = wormy_field,
                                         response = "predicted_WL",
                                         model = "normal",
                                         group = "wormy_asp")


# plot it
plot_WL_worms <- 
    bananaPlot(mod = fitWL_worm$H3,
               data = wormy_field,
               response = "predicted_WL",
               group = "wormy_asp",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3"),
                      name = "Detected Aspiculuris") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Detected Aspiculuris") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_worms

# Create the combined plot with the gradient bar as the "axis"
plot_WL_worms <- 
    plot_grid(plot_WL_worms ,
              HIgradientBar,  
              nrow = 2, 
              rel_heights = c(1.3, 1/8),
              align = "hv",
              axis = "tb",
              vjust = c(-1,2),
              scale = 1) 

# Display the combined plot
plot_WL_worms


ggsave(plot = plot_WL_worms, 
       filename = "figures/hybrid_aspiculuris.jpeg",  width = 6, height = 6, 
       dpi = 1000)

##############################################################
#####  Syphacia
wormy_field <- Field %>%
    mutate(wormy_asp = 
               case_when(
                   Syphacia_sp > 0 ~ "TRUE",
                   Syphacia_sp == 0 ~ "FALSE",
                   Syphacia_sp == NA ~ NA
               ))

wormy_field <- wormy_field %>%
    drop_na(wormy_asp)

wormy_field$wormy_asp <- as.factor(wormy_field$wormy_asp)

fitWL_worm <- parasiteLoad::analyse(data = wormy_field,
                                    response = "predicted_WL",
                                    model = "normal",
                                    group = "wormy_asp")


# plot it
plot_WL_syphacia <- 
    bananaPlot(mod = fitWL_worm$H3,
               data = wormy_field,
               response = "predicted_WL",
               group = "wormy_asp",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3"),
                      name = "Detected Syphacia") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Detected Syphacia") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_syphacia

# Create the combined plot with the gradient bar as the "axis"
plot_WL_syphacia <- 
    plot_grid(plot_WL_syphacia ,
              HIgradientBar,  
              nrow = 2, 
              rel_heights = c(1.3, 1/8),
              align = "hv",
              axis = "tb",
              vjust = c(-1,2),
              scale = 1) 

# Display the combined plot
plot_WL_syphacia


ggsave(plot = plot_WL_syphacia, 
       filename = "figures/hybrid_syphacia.jpeg",  width = 6, height = 6, 
       dpi = 1000)

##############  Trichuris
wormy_field <- Field %>%
    mutate(wormy_asp = 
               case_when(
                   Trichuris_muris > 0 ~ "TRUE",
                   Trichuris_muris == 0 ~ "FALSE",
                   Trichuris_muris == NA ~ NA
               ))

wormy_field <- wormy_field %>%
    drop_na(wormy_asp)

wormy_field$wormy_asp <- as.factor(wormy_field$wormy_asp)

fitWL_worm <- parasiteLoad::analyse(data = wormy_field,
                                    response = "predicted_WL",
                                    model = "normal",
                                    group = "wormy_asp")


# plot it
plot_WL_trichuris <- 
    bananaPlot(mod = fitWL_worm$H3,
               data = wormy_field,
               response = "predicted_WL",
               group = "wormy_asp",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3"),
                      name = "Detected Trichuris") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Detected Trichuris") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_trichuris

# Create the combined plot with the gradient bar as the "axis"
plot_WL_trichuris <- 
    plot_grid(plot_WL_trichuris ,
              HIgradientBar,  
              nrow = 2, 
              rel_heights = c(1.3, 1/8),
              align = "hv",
              axis = "tb",
              vjust = c(-1,2),
              scale = 1) 

# Display the combined plot
plot_WL_trichuris


ggsave(plot = plot_WL_trichuris, 
       filename = "figures/hybrid_trichuris.jpeg", width = 6, height = 6, dpi = 1000)


##############  Heterakis
wormy_field <- Field %>%
    mutate(wormy_asp = 
               case_when(
                   Heterakis_sp > 0 ~ "TRUE",
                   Heterakis_sp == 0 ~ "FALSE",
                   Heterakis_sp == NA ~ NA
               ))

wormy_field <- wormy_field %>%
    drop_na(wormy_asp)

wormy_field$wormy_asp <- as.factor(wormy_field$wormy_asp)

fitWL_worm <- parasiteLoad::analyse(data = wormy_field,
                                    response = "predicted_WL",
                                    model = "normal",
                                    group = "wormy_asp")


# plot it
plot_WL_heterakis <- 
    bananaPlot(mod = fitWL_worm$H3,
               data = wormy_field,
               response = "predicted_WL",
               group = "wormy_asp",
               cols = c("white", "white")) +
    scale_fill_manual(values = c("steelblue1", "indianred3"),
                      name = "Detected Heterakis") +
    scale_color_manual(values = c("steelblue1", "indianred3"),
                       name = "Detected Heterakis") +
    theme_bw()  +
    theme(legend.position = c(0.5, 0.05),
          legend.direction = "horizontal",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Predicted detrimental health impact, 
         Immune signature")

plot_WL_heterakis

# Create the combined plot with the gradient bar as the "axis"
plot_WL_heterakis <- 
    plot_grid(plot_WL_heterakis ,
              HIgradientBar,  
              nrow = 2, 
              rel_heights = c(1.3, 1/8),
              align = "hv",
              axis = "tb",
              vjust = c(-1,2),
              scale = 1) 

# Display the combined plot
plot_WL_heterakis


ggsave(plot = plot_WL_heterakis, 
       filename = "figures/hybrid_heterakis.jpeg",  width = 6, height = 6, dpi = 1000)



infection_plots <- 
    plot_grid(plot_WL_mc_combined, combined_plot, combined_plot_ooc, 
          plot_WL_crypto_Ct, plot_WL_worms, plot_WL_syphacia, plot_WL_trichuris,
          plot_WL_heterakis,
        nrow = 4, align = c("hv"), labels = LETTERS[1:8])

ggsave(infection_plots, filename = "figure_panels/infection_panel.jpeg",
       width = 12, height = 24, dpi = 1000)

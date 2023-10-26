#if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("org.Mm.eg.db")
# BiocManager::install("clusterProfiler")

library(FactoMineR)
library(reshape2)
library(corrplot)
library(factoextra)
library(lmtest)
library(ggpubr)
library(janitor)
library(pheatmap)
library(visdat)
library(scatterplot3d)
library(clusterProfiler) # gene enrichment analysis
library(org.Mm.eg.db) # gene ids identifiers Mus musculus
library(viridis)
library(tidyr)
library(dplyr)



hm <- read.csv("Data/Data_output/imputed_clean_data.csv")


# WOrking with laboratory data only
# Select genes
lab <- hm %>%
    dplyr::filter(origin == "Lab")


Genes_v   <- c("IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10",
               "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
               "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
               "TICAM1", "TNF") #"IL.12", "IRG6")


genes <- lab[ ,colnames(lab) %in% Genes_v]


# PCA
## we can now run a normal pca on the complete data set
res.pca <- PCA(genes)

## How much does each dimension contribute to variance?

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70), barfill = "seagreen2") -> 
    variance_contrib

variance_contrib


ggsave(filename = "figures/contributions_all_dimentsions.jpeg", plot = variance_contrib, 
       width = 6, height = 4, dpi = 1000)

fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#DB6212", "#CC8733", "#5f25e6", "#073DA8"),
             repel = TRUE, title = "") -> pca_col


ggsave(filename = "figures/pca_fviz_package.jpeg", plot = pca_col, 
       width = 10, height = 5, dpi = 300)

fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#DB6212", "#CC8733", "#5f25e6", "#073DA8"), 
             repel = TRUE, title = "")

## Description of the dimensions
## We get a correlation between each variable and the first dimension
dimdesc(res.pca)


# Convert mouse_id to a data frame
mouse <- data.frame(Mouse_ID = lab[,1])
mouse_id <- data.frame(Mouse_ID = lab[,1])

# Add the new column pc1 to the mouse_id data frame
mouse$pc1 <- res.pca$ind$coord[, 1]

mouse$pc2 <- res.pca$ind$coord[, 2]  # indexing the second column

lab <- lab %>% 
    left_join(mouse, by = "Mouse_ID")


## We also need to extract the data for the variable contributions to each of 
## the pc axes.
pca.vars <- res.pca$var$coord %>% data.frame

pca.vars$vars <- rownames(pca.vars)

pca.vars.m <- melt(pca.vars, id.vars = "vars")

source("R/Functions/circle_fun.R")

circ <- circleFun(c(0,0),2,npoints = 500)


#It’s possible to use the function corrplot() [corrplot package] to highlight 
#the most contributing variables for each dimension:
var.contrib <- as.data.frame(res.pca$var$contrib)
var.contrib.matrix <- data.matrix(var.contrib)
corrplot(var.contrib.matrix, is.corr=FALSE) 


pca_var <- as.data.frame(pca.vars)




fviz_pca_biplot(res.pca, 
                col.ind = lab$current_infection, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Infection groups",
                title = "") 

############################### res - tol? 

infecto <- lab %>% 
    filter(!current_infection == "uninfected", MC.Eimeria == "TRUE") %>%
    dplyr::select(all_of(Genes_v), Mouse_ID, delta_ct_cewe_MminusE,
                  MC.Eimeria, OOC, current_infection, pc1, pc2) %>%
    drop_na(delta_ct_cewe_MminusE)


genes <- infecto %>%
    dplyr::select(all_of(Genes_v)) 


# Perform PCA on cleaned data
res.pca <- PCA(genes)

infecto$current_infection <- as.factor(infecto$current_infection)

# Biplot: Color by 'current_infection' and adjust size by 'delta_ct_cewe_MminusE'
biplot_ooc <- fviz_pca_biplot(res.pca, 
                col.ind = infecto$current_infection, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Infection groups",
                title = "",
                pointsize = infecto$OOC) # adjust point size
#palette = c("E.ferrisi" = "blue", "E.falciformis" = "red"))
ggsave(plot = biplot_ooc, filename = "figures/biplot_ooc.jpeg", width = 10, height = 8,
       dpi = 1000)

# Biplot: Color by 'current_infection' and adjust size by 'delta_ct_cewe_MminusE'
biplot_delta <- fviz_pca_biplot(res.pca, 
                              col.ind = infecto$current_infection, palette = "jco", 
                              addEllipses = TRUE, label = "var",
                              col.var = "black", repel = TRUE,
                              legend.title = "Infection groups",
                              title = "",
                              pointsize = infecto$delta_ct_cewe_MminusE) # adjust point size
#palette = c("E.ferrisi" = "blue", "E.falciformis" = "red"))
ggsave(plot = biplot_delta, filename = "figures/biplot_delta.jpeg", width = 10, height = 8,
       dpi = 1000)


# For PC1
infecto_pc1_summary <- infecto %>%
    group_by(current_infection, pc1_sign = ifelse(pc1 > 0, "PC1 Positive", "PC1 Negative")) %>%
    summarise(
        mean_OOC = mean(OOC, na.rm = TRUE),
        sd_OOC = sd(OOC, na.rm = TRUE),
        count = n(),
        .groups = "drop"
    )

# For PC2
infecto_pc2_summary <- infecto %>%
    group_by(current_infection, pc2_sign = ifelse(pc2 > 0, "PC2 Positive", "PC2 Negative")) %>%
    summarise(
        mean_OOC = mean(OOC, na.rm = TRUE),
        sd_OOC = sd(OOC, na.rm = TRUE),
        count = n(),
        .groups = "drop"
    )

print(infecto_pc1_summary)
print(infecto_pc2_summary)

# Extract the OOC values for each group and condition
OOC_E_falciformis_PC1_Positive <- infecto[infecto$current_infection == "E_falciformis" & infecto$pc1 > 0, "OOC"]
OOC_E_falciformis_PC1_Negative <- infecto[infecto$current_infection == "E_falciformis" & infecto$pc1 <= 0, "OOC"]

OOC_E_ferrisi_PC1_Positive <- infecto[infecto$current_infection == "E_ferrisi" & infecto$pc1 > 0, "OOC"]
OOC_E_ferrisi_PC1_Negative <- infecto[infecto$current_infection == "E_ferrisi" & infecto$pc1 <= 0, "OOC"]

# Perform the t-tests
ttest_E_falciformis_PC1 <- t.test(OOC_E_falciformis_PC1_Positive, OOC_E_falciformis_PC1_Negative)
ttest_E_ferrisi_PC1 <- t.test(OOC_E_ferrisi_PC1_Positive, OOC_E_ferrisi_PC1_Negative)

# Print the results
cat("\nT-test results for E_falciformis (PC1 Positive vs PC1 Negative):\n")
print(ttest_E_falciformis_PC1)

cat("\nT-test results for E_ferrisi (PC1 Positive vs PC1 Negative):\n")
print(ttest_E_ferrisi_PC1)


# Extract the OOC values for each group and condition
OOC_E_falciformis_PC2_Positive <- infecto[infecto$current_infection == "E_falciformis" & infecto$pc2 > 0, "OOC"]
OOC_E_falciformis_PC2_Negative <- infecto[infecto$current_infection == "E_falciformis" & infecto$pc2 <= 0, "OOC"]

OOC_E_ferrisi_PC2_Positive <- infecto[infecto$current_infection == "E_ferrisi" & infecto$pc2 > 0, "OOC"]
OOC_E_ferrisi_PC2_Negative <- infecto[infecto$current_infection == "E_ferrisi" & infecto$pc2 <= 0, "OOC"]

# Perform the t-tests
ttest_E_falciformis_PC2 <- t.test(OOC_E_falciformis_PC2_Positive, OOC_E_falciformis_PC2_Negative)
ttest_E_ferrisi_PC2 <- t.test(OOC_E_ferrisi_PC2_Positive, OOC_E_ferrisi_PC2_Negative)

# Print the results
cat("\nT-test results for E_falciformis (PC2 Positive vs PC2 Negative):\n")
print(ttest_E_falciformis_PC2)

cat("\nT-test results for E_ferrisi (PC2 Positive vs PC2 Negative):\n")
print(ttest_E_ferrisi_PC2)


#################
# Extract the delta_ct_cewe_MminusE values for each group and condition
delta_E_falciformis_PC1_Positive <- infecto[infecto$current_infection == "E_falciformis" & infecto$pc1 > 0, "delta_ct_cewe_MminusE"]
delta_E_falciformis_PC1_Negative <- infecto[infecto$current_infection == "E_falciformis" & infecto$pc1 <= 0, "delta_ct_cewe_MminusE"]

delta_E_ferrisi_PC1_Positive <- infecto[infecto$current_infection == "E_ferrisi" & infecto$pc1 > 0, "delta_ct_cewe_MminusE"]
delta_E_ferrisi_PC1_Negative <- infecto[infecto$current_infection == "E_ferrisi" & infecto$pc1 <= 0, "delta_ct_cewe_MminusE"]

# Perform the t-tests
ttest_E_falciformis_PC1 <- t.test(delta_E_falciformis_PC1_Positive, delta_E_falciformis_PC1_Negative)
ttest_E_ferrisi_PC1 <- t.test(delta_E_ferrisi_PC1_Positive, delta_E_ferrisi_PC1_Negative)

# Print the results
cat("\nT-test results for E_falciformis (PC1 Positive vs PC1 Negative) with delta_ct_cewe_MminusE:\n")
print(ttest_E_falciformis_PC1)

cat("\nT-test results for E_ferrisi (PC1 Positive vs PC1 Negative) with delta_ct_cewe_MminusE:\n")
print(ttest_E_ferrisi_PC1)


# Extract the delta_ct_cewe_MminusE values for each group and condition
delta_E_falciformis_PC2_Positive <- infecto[infecto$current_infection == "E_falciformis" & infecto$pc2 > 0, "delta_ct_cewe_MminusE"]
delta_E_falciformis_PC2_Negative <- infecto[infecto$current_infection == "E_falciformis" & infecto$pc2 <= 0, "delta_ct_cewe_MminusE"]

delta_E_ferrisi_PC2_Positive <- infecto[infecto$current_infection == "E_ferrisi" & infecto$pc2 > 0, "delta_ct_cewe_MminusE"]
delta_E_ferrisi_PC2_Negative <- infecto[infecto$current_infection == "E_ferrisi" & infecto$pc2 <= 0, "delta_ct_cewe_MminusE"]

# Perform the t-tests
ttest_E_falciformis_PC2 <- t.test(delta_E_falciformis_PC2_Positive, delta_E_falciformis_PC2_Negative)
ttest_E_ferrisi_PC2 <- t.test(delta_E_ferrisi_PC2_Positive, delta_E_ferrisi_PC2_Negative)

# Print the results
cat("\nT-test results for E_falciformis (PC2 Positive vs PC2 Negative) with delta_ct_cewe_MminusE:\n")
print(ttest_E_falciformis_PC2)

cat("\nT-test results for E_ferrisi (PC2 Positive vs PC2 Negative) with delta_ct_cewe_MminusE:\n")
print(ttest_E_ferrisi_PC2)


library(ggplot2)

# 
primary <- full %>%
    filter(infection == "primary", origin == "Lab")

challenge <- full %>%
    filter(infection == "challenge", origin == "Lab")

ggplot(subset(challenge, dpi != 0), aes(x = dpi, y = delta_ct_cewe_MminusE, color = current_infection)) + 
    geom_point() +
    geom_smooth() +
    labs(title = "Oocyst shedding during the course of infection", 
         x = "Days post infection", 
         y = "Oocysts per gram", 
         color = "Treatment") +
    theme_minimal() +
    theme(
        text = element_text(family = "serif"),
        title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "gray85"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)
    ) +
    
    # Labels and scales
    labs(
        title = "Oocyst shedding during the course of infection", 
        x = "Days post infection", 
        y = "Oocysts per gram", 
        color = "Treatment",
        shape = "Treatment"
    ) +
    scale_color_brewer(palette = "Set1") +
    scale_shape_manual(values = c(16, 17, 18)) 

ggplot(subset(primary, dpi != 0), aes(x = dpi, y = OOC, color = current_infection)) + 
    # Points and smooth line
    geom_point(aes(shape = current_infection), size = 5, alpha = 05) +
    geom_smooth(se = FALSE, method = "loess", lwd = 1.5) +
    
    # Theme adjustments for publication quality
    

    
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

ggplot(Field, aes(x = HI, y = predicted_WL)) +
    geom_point() +
    geom_jitter() +
    geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE)

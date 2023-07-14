## ----message = FALSE, warnings = FALSE--------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(stringr)
library(FactoMineR)
library(reshape2)
library(corrplot)
library(factoextra)
library(lmtest)
library(ggpubr)
library(janitor)
library(pheatmap)
library(visdat)


## ---------------------------------------------------------------------------------------------------
hm <- read.csv("output_data/2.imputed_MICE_data_set.csv")


## ---------------------------------------------------------------------------------------------------
Gene_lab   <- c("IFNy", "CXCR3", "IL.6", "IL.13",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF") # "IL.12", "IRG6")



## ----black_and_white_pca----------------------------------------------------------------------------
#select the genes and lab muce
lab <- hm %>%
  dplyr::filter(origin == "Lab", Position == "mLN") #selecting for mln to avoid
# duplicates

lab <- unique(lab)

gene <- lab %>%
  dplyr::select(c(Mouse_ID, all_of(Gene_lab)))

genes <- unique(gene)

genes <- genes[, -1]

#remove rows with only nas
genes <- genes[,colSums(is.na(genes))<nrow(genes)]

#remove colums with only nas 
genes <- genes[rowSums(is.na(genes)) != ncol(genes), ]

#select same rows in the first table
gene <- gene[row.names(genes), ]

# we can now run a normal pca on the complete data set
res.pca <- PCA(genes)


## ----percentage_of_explained_variance_per_dimension-------------------------------------------------
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70))


## ----PCA_gene_lab_fig_1_1, resolution = 1000--------------------------------------------------------
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#DB6212", "#CC8733", "#5f25e6", "#073DA8"),
             repel = TRUE, title = "")


## ----individual_plot_pca, resolution = 1000---------------------------------------------------------
fviz_pca_ind(res.pca, col.ind = "cos2", 
                  gradient.cols = c("#DB6212", "#CC8733", "#5f25e6", "#073DA8"), 
                  repel = TRUE, title = "")


## ----dimensions, include = FALSE, echo = FALSE, warnings = FALSE------------------------------------
#Description of the dimensions
# We get a correlation between each variable and the first dimension
dimdesc(res.pca)


## ---- echo = FALSE----------------------------------------------------------------------------------

mouse_id <- gene %>%
  dplyr::select(Mouse_ID)



mouse_id$pc1 <- res.pca$ind$coord[, 1] # indexing the first column

mouse_id$pc2 <- res.pca$ind$coord[, 2]  # indexing the second column




lab <- lab %>% 
  left_join(mouse_id, by = "Mouse_ID")



#We also need to extract the data for the variable contributions to each of the pc axes.
pca.vars <- res.pca$var$coord %>% data.frame


pca.vars$vars <- rownames(pca.vars)

pca.vars.m <- melt(pca.vars, id.vars = "vars")

source("r_scripts/functions/circle_fun.R")

circ <- circleFun(c(0,0),2,npoints = 500)



## ----correlations_genes_dimensions, resolution = 1000-----------------------------------------------
#It’s possible to use the function corrplot() [corrplot package] to highlight 
#the most contributing variables for each dimension:
var.contrib <- res.pca$var$contrib
corrplot(var.contrib, is.corr=FALSE) 


## ----contr_var_pc_genes_dimension1_figure_1_2, resolution = 1000------------------------------------
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 18, 
             title = "Contribution of immune genes to the first dimension of the PCA")

# res.pca$var$contrib


## ----contr_var_pc_genes_dimension2_figure_1_3, resolution = 1000------------------------------------
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 18, 
             title = "Contribution of immune genes to the second dimension of the PCA")


## ----contr_var_pc1_2_genes, echo = FALSE------------------------------------------------------------
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 18)


## ----contr_individuals_genes, echo = FALSE----------------------------------------------------------
# Total contribution on PC1 and PC2
fviz_contrib(res.pca, choice = "ind", axes = 1:2)


## ----pca_biplot_genes_infecting_parasite_figure_1_4, resolution = 1000------------------------------
#select same rows in the first table
lab <- lab[row.names(genes), ]

fviz_pca_biplot(res.pca, 
                col.ind = lab$Parasite_challenge, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Infection groups",
                title = "") 



## ----Linear_model_1, echo = FALSE-------------------------------------------------------------------
# predicting weight loss with the pc1 and pc2
model_1_pc1_pc2 <- lm(WL_max ~ pc1 + pc2, data = lab)
summary(model_1_pc1_pc2)
AIC(model_1_pc1_pc2)


## ----linear_model_2, echo = FALSE-------------------------------------------------------------------

model_2_pc1_pc2_challenge <- lm(WL_max ~ pc1 + pc2 + Parasite_challenge, data = lab)
summary(model_2_pc1_pc2_challenge)
AIC(model_2_pc1_pc2_challenge)


## ----linear_model_3, echo = FALSE-------------------------------------------------------------------

model_3_Parasite_challenge_hybrid_status <- lm(WL_max ~ pc1 + pc2 + hybrid_status, 
                 data = lab)
summary(model_3_Parasite_challenge_hybrid_status)
AIC(model_3_Parasite_challenge_hybrid_status)


## ----LLR_test---------------------------------------------------------------------------------------
# Compare
llr_test <- anova(model_1_pc1_pc2, model_2_pc1_pc2_challenge)
print(llr_test)


## ----Likelihood_ratio_test--------------------------------------------------------------------------
# model_2_pc1_pc2_challenge <- lm(WL_max ~ pc1 + pc2 + Parasite_challenge, data = lab)

weight_no_pc1 <- lm(WL_max ~ pc2 + Parasite_challenge, data = lab)
weight_no_pc2 <- lm(WL_max ~ pc1  + Parasite_challenge, data = lab)
weight_no_Parasite_challenge <- lm(WL_max ~ pc1 + pc2, data = lab)
lrtest(model_2_pc1_pc2_challenge, weight_no_pc1)
lrtest(model_2_pc1_pc2_challenge, weight_no_pc2)
lrtest(model_2_pc1_pc2_challenge, weight_no_Parasite_challenge)
lrtest(weight_no_pc1, weight_no_pc2)



## ----scatter_plot_pc1_WL_linear_model, resolution = 1000--------------------------------------------
# for the model with pc1 and pc2
ggplot(lab, aes(x = pc1, y = WL_max)) +
  geom_point(aes(color = Parasite_challenge)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "PC1", y = "Maximum Weight Loss") +
  theme_bw()




## ----scatter_plot_pc2_WL_linear_model, resolution = 1000--------------------------------------------
# for the model with pc2 and Parasite_challenge
ggplot(lab, aes(x = pc2, y = WL_max)) +
  geom_point(aes(color = Parasite_challenge)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "PC2", y = "Maximum Weight Loss") +
  theme_bw()



## ----residual_plot_pc1, resolution = 1000-----------------------------------------------------------
# calculate residuals for the model with pc1 and pc2
lab$residuals_pc1_pc2 <- resid(model_1_pc1_pc2)

ggplot(lab, aes(x = pc1, y = residuals_pc1_pc2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "PC1", y = "Residuals") +
  theme_bw()

## ----residual_plot_pc2, resolution = 1000-----------------------------------------------------------
ggplot(lab, aes(x = pc2, y = residuals_pc1_pc2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "PC2", y = "Residuals") +
  theme_bw()

## ----residual_plot_parasite_Challenge, resolution = 1000--------------------------------------------
ggplot(lab, aes(x = Parasite_challenge, y = residuals_pc1_pc2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Infection group", y = "Residuals") +
  theme_bw()


## ----scatter_p1_infections_colours, resolution = 1000-----------------------------------------------

# First, make sure Parasite_challenge is a factor
lab$Parasite_challenge <- as.factor(lab$Parasite_challenge)

# Then, define the color for each level of Parasite_challenge
color_mapping <- c("E_falciformis" = "salmon", 
                   "E_ferrisi" = "green", 
                   "uninfected" = "blue")

# Now create the scatter plot using this color mapping
ggplot(lab, aes(x = pc1, y = WL_max, color = Parasite_challenge)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(color = Parasite_challenge)) +
  scale_color_manual(values = color_mapping) +
  labs(x = "PC1", y = "Maximum Weight Loss") +
  theme_bw()



## ----scatter_p2_infections_colours, resolution = 1000-----------------------------------------------
# Now create the scatter plot using this color mapping
ggplot(lab, aes(x = pc2, y = WL_max, color = Parasite_challenge)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(color = Parasite_challenge)) +
  scale_color_manual(values = color_mapping) +
  labs(x = "PC2", y = "Maximum Weight Loss") +
  theme_bw()


## ----scatterplot3d_pc1_pc2, resolution = 1000-------------------------------------------------------
library(scatterplot3d)


# Create a new color column in your dataframe by mapping 'Parasite_challenge' to your colors
lab$color <- color_mapping[lab$Parasite_challenge]

# 3D scatter plot
scatterplot3d(lab$pc1, lab$pc2, lab$WL_max, pch = 16, color = lab$color,
              xlab = "PC1", ylab = "PC2", zlab = "Maximum Weight Loss")



## ---------------------------------------------------------------------------------------------------
 # turn the data frame into a matrix and transpose it. We want to have each cell 
 # type as a row name 
 gene <- t(as.matrix(gene))
 
 # turn the first row into column names
 gene %>%
     row_to_names(row_number = 1) -> heatmap_data
 
 heatmap_data <- as.data.frame(heatmap_data)
 
 table(rowSums(is.na(heatmap_data)) == nrow(heatmap_data))

 
# turn the columns to numeric other wise the heatmap function will not work
 heatmap_data[] <- lapply(heatmap_data, function(x) as.numeric(as.character(x)))

 # remove columns with only NAs 
 heatmap_data <- Filter(function(x)!all(is.na(x)), heatmap_data) 
 
 #remove rows with only Nas
 heatmap_data <-  heatmap_data[, colSums(is.na(heatmap_data)) != 
                                   nrow(heatmap_data)]
 
  
#Prepare the annotation data frame
annotation_df <- as_tibble(lab) %>%
    dplyr::select(c("Mouse_ID",  "WL_max", "Parasite_challenge")) 
  
annotation_df <- unique(annotation_df) 

annotation_df <- as.data.frame(annotation_df)




### Prepare the annotation columns for the heatmap
rownames(annotation_df) <- annotation_df$Mouse_ID


# Match the row names to the heatmap data frame
rownames(annotation_df) <- colnames(heatmap_data)

#remove the unecessary column
annotation_df <- annotation_df %>% dplyr::select(-Mouse_ID, )




## ---- Heatmap_gene_expression, resolution = 1000----------------------------------------------------
# Define colors for each parasite
parasite_colors <- c("E_falciformis" = "coral2",
                     "E_ferrisi" = "chartreuse4",
                     "uninfected" = "cornflowerblue")

# Generate the heat map
pheatmap(heatmap_data, annotation_col = annotation_df,  scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_colors = list(Parasite_challenge = parasite_colors)) # use annotation_colors



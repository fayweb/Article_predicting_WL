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


hm <- read.csv("Data/Data_output/imputed_clean_data.csv")


field <- hm %>%
    filter(origin == "Field")

genes_mouse <- field %>%
  dplyr::select(c(Mouse_ID, "IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10", removed IL.10 
                  # as we had too few measurements to impute
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF")) 


genes <- genes_mouse %>%
  dplyr::select(-Mouse_ID)

# we can now run a normal pca on the complete data set
res.pca <- PCA(genes)

fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70))

fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)

fviz_pca_ind(res.pca, col.ind = "cos2", 
                  gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                  repel = TRUE)

#Description of the dimensions
# We get a correlation between each variable and the first dimension
dimdesc(res.pca)

str(res.pca)

field$pc1 <- res.pca$ind$coord[, 1] # indexing the first column

field$pc2 <- res.pca$ind$coord[, 2]  # indexing the second column



#We also need to extract the data for the variable contributions to each of the pc axes.
pca.vars <- res.pca$var$coord %>% data.frame


pca.vars$vars <- rownames(pca.vars)

pca.vars.m <- melt(pca.vars, id.vars = "vars")

source("R/Functions/circle_fun.R")

circ <- circleFun(c(0,0),2,npoints = 500)

#It’s possible to use the function corrplot() [corrplot package] to highlight 
#the most contributing variables for each dimension:
var.contrib <- res.pca$var$contrib
corrplot(var.contrib, is.corr=FALSE) 

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 18)

# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 18)


fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 18)

#The most important (or, contributing) variables can be highlighted on the 
#correlation plot as follow:
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
             )


# Total contribution on PC1 and PC2
fviz_contrib(res.pca, choice = "ind", axes = 1:2)


fviz_pca_biplot(res.pca, 
                col.ind = field$MC.Eimeria, palette = "jco", 
                addEllipses = TRUE, fieldel = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Infected with Eimeria") 

 # turn the data frame into a matrix and transpose it. We want to have each cell 
 # type as a row name 
 genes_mouse <- t(as.matrix(genes_mouse))
 
 # turn the first row into column names
 genes_mouse %>%
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
annotation_df <- as_tibble(field) %>%
    dplyr::select(c("Mouse_ID", "HI",
                    "Sex")) 
  
annotation_df <- unique(annotation_df) 

annotation_df <- as.data.frame(annotation_df)

 ### Prepare the annotation columns for the heatmap

rownames(annotation_df) <- annotation_df$Mouse_ID

# Match the row names to the heatmap data frame
rownames(annotation_df) <- colnames(heatmap_data)

#remove the unecessary column
annotation_df <- annotation_df %>% dplyr::select(-Mouse_ID)




## heatmap
pheatmap(heatmap_data, annotation_col = annotation_df, scale = "row")


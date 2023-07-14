## ---------------------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------------------
hm <- read.csv("output_data/2.imputed_MICE_data_set.csv")


## ---------------------------------------------------------------------------------------------------
Gene_field   <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF") #"IL.12", "IRG6")

#add a suffix to represent changes in data file
Gene_field_imp <- paste(Gene_field, "imp", sep = "_")

Genes_wild   <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10", 
                  "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                  "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                  "TICAM1", "TNF") #, "IL.12", "IRG6")

Genes_wild_imp <- paste(Genes_wild, "imp", sep = "_")

Facs_field <- c("Position", "CD4", "Treg", "Div_Treg", "Treg17", "Th1", 
                    "Div_Th1", "Th17", "Div_Th17", "CD8", "Act_CD8", 
                    "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8","Treg_prop", 
                    "IL17A_CD4")  

Facs_field_imp <- paste(Facs_field, "imp", sep = "_")

Facs_wild <- c( "Treg", "CD4", "Treg17", "Th1", "Th17", "CD8",
                     "Act_CD8", "IFNy_CD4", "IL17A_CD4", "IFNy_CD8")

Facs_wild_imp <- paste(Facs_wild, "imp", sep = "_")


## ---------------------------------------------------------------------------------------------------
# somehow the field samples have the origin na,
# fix that
hm$origin[is.na(hm$origin)] <- "Field"

field <- hm %>%
  dplyr::filter(origin == "Field") 

field <- unique(field)

#make a factor out of the melting curves (important for later visualization)
field <- field %>%
  dplyr::mutate(MC.Eimeria = as.factor(MC.Eimeria))

genes_mouse <- field %>%
  dplyr::select(c(Mouse_ID, "IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10", removed IL.10 
                  # as we had too few measurements to impute
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF")) 


genes <- genes_mouse %>%
  dplyr::select(-Mouse_ID)


#remove rows with only nas
genes <- genes[,colSums(is.na(genes))<nrow(genes)]

#remove colums with only nas 
genes <- genes[rowSums(is.na(genes)) != ncol(genes), ]


##select same rows in the first table
field <- field[row.names(genes), ]


# we can now run a normal pca on the complete data set
res.pca <- PCA(genes)


## ---------------------------------------------------------------------------------------------------
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70))


## ---------------------------------------------------------------------------------------------------
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)

## ---------------------------------------------------------------------------------------------------
fviz_pca_ind(res.pca, col.ind = "cos2", 
                  gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                  repel = TRUE)


## ----dimensions, include = FALSE, echo = FALSE, warnings = FALSE------------------------------------
#Description of the dimensions
# We get a correlation between each variable and the first dimension
dimdesc(res.pca)


## ---- echo = FALSE, include = FALSE-----------------------------------------------------------------
str(res.pca)


## ---- echo = FALSE----------------------------------------------------------------------------------


field$pc1 <- res.pca$ind$coord[, 1] # indexing the first column

field$pc2 <- res.pca$ind$coord[, 2]  # indexing the second column




#We also need to extract the data for the variable contributions to each of the pc axes.
pca.vars <- res.pca$var$coord %>% data.frame


pca.vars$vars <- rownames(pca.vars)

pca.vars.m <- melt(pca.vars, id.vars = "vars")

source("r_scripts/functions/circle_fun.R")

circ <- circleFun(c(0,0),2,npoints = 500)



## ----correlations_genes_dimensions, echo = FALSE----------------------------------------------------
#It’s possible to use the function corrplot() [corrplot package] to highlight 
#the most contributing variables for each dimension:
var.contrib <- res.pca$var$contrib
corrplot(var.contrib, is.corr=FALSE) 


## ----contr_var_pc_genes, echo = FALSE---------------------------------------------------------------
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 18)



## ---------------------------------------------------------------------------------------------------
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 18)


## ----contr_var_pc1_2_genes, echo = FALSE------------------------------------------------------------
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 18)


## ----pca_contribution_genes, echo = FALSE-----------------------------------------------------------

#The most important (or, contributing) variables can be highlighted on the 
#correlation plot as follow:
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
             )


## ----contr_individuals_genes, echo = FALSE----------------------------------------------------------
# Total contribution on PC1 and PC2
fviz_contrib(res.pca, choice = "ind", axes = 1:2)


## ----pca_biplot_genes, echo = FALSE-----------------------------------------------------------------


fviz_pca_biplot(res.pca, 
                col.ind = field$MC.Eimeria, palette = "jco", 
                addEllipses = TRUE, fieldel = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Infected with Eimeria") 



## ---------------------------------------------------------------------------------------------------
genes_mouse <- genes_mouse %>%
  rename_with(~str_remove(., '_imp'))

 #select same rows in the first table
genes_mouse <- genes_mouse[row.names(genes), ]
 
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




## ---- echo = FALSE----------------------------------------------------------------------------------
pheatmap(heatmap_data, annotation_col = annotation_df, scale = "row")


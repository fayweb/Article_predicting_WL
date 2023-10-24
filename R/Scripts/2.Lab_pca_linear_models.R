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


### Contributions to the first dimension

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 18, 
             title = "Contribution of immune genes to the first dimension of the PCA", 
             fill =  "seagreen2") -> contributions_pc1

contributions_pc1

ggsave(filename = "figures/contributions_pc1.jpeg", plot = contributions_pc1, 
       width = 6, height = 4, dpi = 1000)
# res.pca$var$contrib


### Contributions to the second dimension

## Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 18, 
             title = "Contribution of immune genes to the second dimension of the PCA",
             fill =  "seagreen2") -> contributions_pc2

contributions_pc2

ggsave(filename = "figures/contributions_pc2.jpeg", plot = contributions_pc2, 
       width = 6, height = 4, dpi = 1000)


fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 18)

# Total contribution on PC1 and PC2
fviz_contrib(res.pca, choice = "ind", axes = 1:2)

#select same rows in the first table
lab <- lab[row.names(genes), ]


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
fviz_pca_biplot(res.pca, 
                col.ind = infecto$current_infection, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Infection groups",
                title = "",
                pointsize = infecto$OOC) # adjust point size
                #palette = c("E.ferrisi" = "blue", "E.falciformis" = "red"))

    
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












################## Linear models: Predicting weight loss with the PCA eigenvectors
# predicting weight loss with the pc1 and pc2
model_1_pc1_pc2 <- lm(WL_max ~ pc1 + pc2, data = lab)
summary(model_1_pc1_pc2)
AIC(model_1_pc1_pc2)

### use the ggefects package
# 

# predicting weight loss with pc1
model_1_pc1 <- lm(WL_max ~ pc1, data = lab)
summary(model_1_pc1)
AIC(model_1_pc1)


# Here the base of comparison is E_ferrisi. I should change it to compare to the
# uninfected mice
model_2_pc1_pc2_challenge <- lm(WL_max ~ pc1 + pc2 + infection, data = lab)
summary(model_2_pc1_pc2_challenge)
AIC(model_2_pc1_pc2_challenge)


# homozygous / heterozygous infection


# correlations between infection and the immune responses
plot(model_2_pc1_pc2_challenge)

# covariance matrix of the fixed effects of the model


model_3_infection_hybrid_status <- lm(WL_max ~ pc1 + pc2 + hybrid_status, 
                                      data = lab)
summary(model_3_infection_hybrid_status)
AIC(model_3_infection_hybrid_status)

### information of previous infections, some mice are actually reacting for the 
# first time, as the



# Compare
llr_test <- anova(model_1_pc1_pc2, model_2_pc1_pc2_challenge)
print(llr_test)

# model_2_pc1_pc2_challenge <- lm(WL_max ~ pc1 + pc2 + current_infection, data = lab)

weight_no_pc1 <- lm(WL_max ~ pc2 + current_infection, data = lab)
weight_no_pc2 <- lm(WL_max ~ pc1  + current_infection, data = lab)
weight_no_infection <- lm(WL_max ~ pc1 + pc2, data = lab)
lrtest(model_2_pc1_pc2_challenge, weight_no_pc1)
lrtest(model_2_pc1_pc2_challenge, weight_no_pc2)
lrtest(model_2_pc1_pc2_challenge, weight_no_infection)
lrtest(weight_no_pc1, weight_no_pc2)


## Visualizing the regression models
### scatter plot with regression lines 

# for the model with pc1 and pc2
ggplot(lab, aes(x = pc1, y = WL_max)) +
  geom_point(aes(color = current_infection)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "PC1", y = "Maximum Weight Loss") +
  theme_bw()


# for the model with pc2 and infection
ggplot(lab, aes(x = pc2, y = WL_max)) +
  geom_point(aes(color = current_infection)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "PC2", y = "Maximum Weight Loss") +
  theme_bw()



### Residual Plots

# calculate residuals for the model with pc1 and pc2
lab$residuals_pc1_pc2 <- resid(model_1_pc1_pc2)

ggplot(lab, aes(x = pc1, y = residuals_pc1_pc2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "PC1", y = "Residuals") +
  theme_bw()


ggplot(lab, aes(x = pc2, y = residuals_pc1_pc2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "PC2", y = "Residuals") +
  theme_bw()


ggplot(lab, aes(x = current_infection, y = residuals_pc1_pc2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Infection group", y = "Residuals") +
  theme_bw()


### 3D plots



# First, make sure infection is a factor
lab$current_infection <- as.factor(lab$infection)

# Then, define the color for each level of infection
color_mapping <- c("E_falciformis" = "salmon", 
                   "E_ferrisi" = "green", 
                   "uninfected" = "blue")

# Now create the scatter plot using this color mapping
ggplot(lab, aes(x = pc1, y = WL_max, color = current_infection)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(color = current_infection)) +
  scale_color_manual(values = color_mapping) +
  labs(x = "PC1", y = "Maximum Weight Loss") +
  theme_bw()



# Now create the scatter plot using this color mapping
ggplot(lab, aes(x = pc2, y = WL_max, color = current_infection)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(color = current_infection)) +
  scale_color_manual(values = color_mapping) +
  labs(x = "PC2", y = "Maximum Weight Loss") +
  theme_bw()




# Heatmap
### repeating the heatmap on the now imputed data

# turn the data frame into a matrix and transpose it. We want to have each cell 
# type as a row name 
gene <- t(as.matrix(data.frame(mouse_id,genes)))

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


# Define colors for each parasite
parasite_colors <- c("E_falciformis" = "coral2",
                     "E_ferrisi" = "chartreuse4",
                     "uninfected" = "cornflowerblue")
# Define your own colors
#my_colors <- colorRampPalette(c("red", "white", "blue"))(100)

# Generate the heat map
pheatmap(heatmap_data, annotation_col = annotation_df,# color = my_colors,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_colors = list(current_infection = parasite_colors)) # use annotation_colors


####### Gene enrichment analysis
#create a new vector n to match the genes with gene ids in the package
# Enrichr

# Get the available keytypes in the database
keytypes <- keytypes(org.Mm.eg.db)

# View the list of keytypes
print(keytypes)

gene_ens <-  c("ENSMUSG00000055170", #IFNG, 
               "ENSMUSG00000050232", #CXCR3"
               "ENSMUSG00000025746", #IL6
               "ENSMUSG00000020383", #IL13",
               "ENSMUSG00000026981", #IL1RN
               "ENSMUSG00000025888", #CASP1
               "ENSMUSG00000029417", #CXCL9
               "ENSMUSG00000031551", #IDO1
               "ENSMUSG00000046879", #IRGM1
               "ENSMUSG00000009350", #MPO"
               "ENSMUSG00000025515", #MUC2
               "ENSMUSG00000037974", #MUC5AC
               "ENSMUSG00000032508", #MYD88
               "ENSMUSG00000062524", #NCR1
               "ENSMUSG00000037202.7", #PRF1
               "ENSMUSG00000022650", #RETNLB
               "ENSMUSG00000038037", #SOCS1
               "ENSMUSG00000047123", #TICAM1
               "ENSMUSG00000024401") #TNF             


#All genes
# Perform gene ontology enrichment analysis
enrich_result <- enrichGO(gene = gene_ens,
                          OrgDb = org.Mm.eg.db,
                          ont = "BP",
                          keyType = "ENSEMBL")



# Extract the relevant columns from the enrichment table
enriched_terms <- enrich_result$Description
p_values <- enrich_result$p.adjust
gene_ratio <- enrich_result$GeneRatio

enriched_terms_df <- data.frame(enrich_result$Description, 
                                enrich_result$p.adjust)

enriched_terms_df <- enriched_terms_df %>%
  dplyr::rename("GO_Term" = "enrich_result.Description", "p_value" = 
                  "enrich_result.p.adjust" )

write.csv(enriched_terms_df, "Data/Data_output/enriched_sorted_terms.csv", row.names = FALSE)

# Sort the enriched terms based on p-values
sorted_terms <- as.data.frame(enriched_terms[order(p_values)])


sorted_terms[1:10,]


# Create a data frame for the bar plot
bar_data <- data.frame(GO_Term = sorted_terms[1:35,], p_value = -log10(p_values[1:35]))

# Sort the data frame in ascending order
bar_data <- bar_data[order(bar_data$p_value), ]

# Create the bar plot using ggplot2
ggplot(bar_data, aes(x = GO_Term, y = p_value, fill = p_value)) +
  geom_segment(aes(xend = GO_Term, yend = 0), color = "mediumvioletred", size = 1.5) +
  geom_point(size = 3, shape = 19, color = "mediumvioletred", fill = "white") +
  coord_flip() +
  labs(x = "Enriched GO Terms", y = "-log10(p-value)",
       title = "Gene Ontology Enrichment Analysis") +
  theme_minimal()

########################################
# FInding out which genes are associated with specific GO:Terms
# make a data frame out of ENSEMBLID and the descriptions of the genes
genes_sel <- data.frame(geneID = enrich_result$geneID, 
                        Description = enrich_result$Description)

# Filter the data frame based on conditions
filter_genes <- function(x)genes_sel %>%
  filter(str_detect(geneID, paste(gene_ens, collapse = "|")), 
         Description == x)

filter_genes("positive regulation of inflammatory response")

# "positive regulation of inflammatory response"
# ENSMUSG00000055170/ #IFNG,
#ENSMUSG00000025746/ IL6
#ENSMUSG00000025888/ CASP1
#ENSMUSG00000031551/ IDO1
#ENSMUSG00000024401 TNF

filter_genes("regulation of inflammatory response")
#ENSMUSG00000020383/ IL.13
#ENSMUSG00000032508/MYD88


pca_var <- pca_var %>%
  dplyr::mutate(Pro_infl = case_when(
    vars == "IFNy" ~ "positive regulation of inflammatory response",
    vars == "IL.6" ~ "positive regulation of inflammatory response",
    vars == "CASP1" ~ "positive regulation of inflammatory response",
    vars == "IDO1" ~ "positive regulation of inflammatory response",
    vars == "TNF" ~ "positive regulation of inflammatory response",
    vars == "IL.13" ~ "regulation of inflammatory response",
    vars == "MYD88" ~ "regulation of inflammatory response",
    TRUE ~ "none"))

filter_genes("positive regulation of T cell activation")
# ENSMUSG00000055170/ IFNG
# ENSMUSG00000025746/ IL6
# ENSMUSG00000038037 SOCS1

filter_genes("negative regulation of T cell activation")
# ENSMUSG00000031551/ IDO1
#ENSMUSG00000038037 SOCS1
pca_var <- pca_var %>%
  dplyr::mutate(T_activ = case_when(
    vars == "IFNy" ~ "positive regulation of T cell activation",
    vars == "IL.6" ~ "positive regulation of T cell activation",
    vars == "SOCS1" ~ "positive / negative regulation of T cell activation",
    vars == "IDO1" ~ "negative regulation of T cell activation",
    TRUE ~ "none"))

filter_genes("negative regulation of response to cytokine stimulus")
# ENSMUSG00000025746 IL6
# ENSMUSG00000026981 IL1RN

filter_genes("positive regulation of response to cytokine stimulus")
#ENSMUSG00000025888/ CASP1
# ENSMUSG00000046879 IRGM1
pca_var <- pca_var %>%
  dplyr::mutate(Cytokine_response = case_when(
    vars == "IL.6" ~ "negative regulation of response to cytokine stimulus",
    vars == "IL1RN" ~ "negative regulation of response to cytokine stimulus",
    vars == "CASP1" ~ "positive regulation of response to cytokine stimulus",
    vars == "IRGM1" ~ "positive regulation of response to cytokine stimulus",
    TRUE ~ "none"))

filter_genes("positive regulation of immune effector process")
# ENSMUSG00000055170/ IFNG
#ENSMUSG00000025746/ IL6
#ENSMUSG00000020383/ IL13
#ENSMUSG00000047123/ TICAM1
#ENSMUSG00000024401 TNF

filter_genes("negative regulation of immune system process")
#ENSMUSG00000055170/ IFNg
#ENSMUSG00000031551/ IDO1
#ENSMUSG00000038037/ SOC1
#ENSMUSG00000024401 TNF


filter_genes("negative regulation of cytokine production")
#ENSMUSG00000055170/ IFNg
#ENSMUSG00000025746/ IL6
#ENSMUSG00000020383/ IL13
#ENSMUSG00000031551/ IDO1
#ENSMUSG00000024401 TNF

filter_genes("positive regulation of cytokine production involved in inflammatory response")
# ENSMUSG00000025746/ IL6
#ENSMUSG00000032508/MYD88
#ENSMUSG00000047123/TICAM1
#ENSMUSG00000024401 TNF
pca_var <- pca_var %>%
  dplyr::mutate(Cytokine_response = case_when(
    vars == "IFNg" ~ "negative regulation of cytokine production",
    vars == "IL6" ~ "regulation of cytokine production",
    vars == "IL13" ~ "negative regulation of cytokine production",
    vars == "IDO1" ~ "negative regulation of cytokine production",
    vars == "TNF" ~ "regulation of cytokine production",
    vars == "MYD88" ~ "positive regulation of cytokine production",
    vars == "TICAM1" ~ "positive regulation of cytokine production",
    TRUE ~ "none"))


## Now go on to select the interest groupings seen on the pca
############################# IL.13
# "positive regulation of mast cell activation involved in immune response"
#  "leukocyte activation involved in inflammatory response" 
#"positive regulation of B cell proliferation"
gene_ens <-    "ENSMUSG00000020383"#IL13"

# Perform gene ontology enrichment analysis
enrich_result <- enrichGO(gene = gene_ens,
                          OrgDb = org.Mm.eg.db,
                          ont = "BP",
                          keyType = "ENSEMBL")

# View the enrichment result
enrich_result

# Extract the relevant columns from the enrichment table
enriched_terms <- enrich_result$Description
p_values <- enrich_result$p.adjust
gene_ratio <- enrich_result$GeneRatio

# Sort the enriched terms based on p-values
sorted_terms <- as.data.frame(enriched_terms[order(p_values)])
sorted_terms[1:50,]


# Create a data frame for the bar plot
bar_data <- data.frame(GO_Term = sorted_terms[1:35,], p_value = -log10(p_values[1:35]))

# Sort the data frame in ascending order
bar_data <- bar_data[order(bar_data$p_value), ]

# Create the bar plot using ggplot2
ggplot(bar_data, aes(x = GO_Term, y = p_value, fill = p_value)) +
  geom_segment(aes(xend = GO_Term, yend = 0), color = "mediumvioletred", size = 1.5) +
  geom_point(size = 3, shape = 19, color = "mediumvioletred", fill = "white") +
  coord_flip() +
  labs(x = "Enriched GO Terms", y = "-log10(p-value)",
       title = "Gene Ontology Enrichment Analysis, IL13") +
  theme_minimal()


####################### TICAM1
# "MyD88-independent toll-like receptor signaling pathway" 
# "macrophage activation involved in immune response" 
#"positive regulation of B cell proliferation"                                             
#[14] "positive regulation of interferon-beta production" 
#  "positive regulation of cytokine production involved in immune response"  
# "positive regulation of interleukin-6 production"  
gene_ens <-  "ENSMUSG00000047123" #TICAM1

# Perform gene ontology enrichment analysis
enrich_result <- enrichGO(gene = gene_ens,
                          OrgDb = org.Mm.eg.db,
                          ont = "BP",
                          keyType = "ENSEMBL")

# View the enrichment result
enrich_result

# Extract the relevant columns from the enrichment table
enriched_terms <- enrich_result$Description
p_values <- enrich_result$p.adjust
gene_ratio <- enrich_result$GeneRatio

# Sort the enriched terms based on p-values
sorted_terms <- as.data.frame(enriched_terms[order(p_values)])
sorted_terms[1:50,]


# Create a data frame for the bar plot
bar_data <- data.frame(GO_Term = sorted_terms[1:35,], p_value = -log10(p_values[1:35]))

# Sort the data frame in ascending order
bar_data <- bar_data[order(bar_data$p_value), ]

# Create the bar plot using ggplot2
ggplot(bar_data, aes(x = GO_Term, y = p_value, fill = p_value)) +
  geom_segment(aes(xend = GO_Term, yend = 0), color = "mediumvioletred", size = 1.5) +
  geom_point(size = 3, shape = 19, color = "mediumvioletred", fill = "white") +
  coord_flip() +
  labs(x = "Enriched GO Terms", y = "-log10(p-value)",
       title = "Gene Ontology Enrichment Analysis, TICAM1") +
  theme_minimal()


###################################################################################

################### IRGM1, SOCS1, MUC2
gene_ens <- c("ENSMUSG00000046879", #IRGM1
              "ENSMUSG00000038037", #SOCS1
              "ENSMUSG00000025515") #MUC2

# Perform gene ontology enrichment analysis
enrich_result <- enrichGO(gene = gene_ens,
                          OrgDb = org.Mm.eg.db,
                          ont = "BP",
                          keyType = "ENSEMBL")

# View the enrichment result
enrich_result

# Extract the relevant columns from the enrichment table
enriched_terms <- enrich_result$Description
p_values <- enrich_result$p.adjust
gene_ratio <- enrich_result$GeneRatio

# Sort the enriched terms based on p-values
sorted_terms <- as.data.frame(enriched_terms[order(p_values)])
sorted_terms[1:10,]


# Create a data frame for the bar plot
bar_data <- data.frame(GO_Term = sorted_terms[1:10,], p_value = -log10(p_values[1:10]))

# Sort the data frame in ascending order
bar_data <- bar_data[order(bar_data$p_value), ]

# Create the bar plot using ggplot2
ggplot(bar_data, aes(x = reorder(GO_Term, p_value), y = p_value, fill = p_value)) +
  geom_segment(aes(xend = GO_Term, yend = 0), color = "mediumvioletred", size = 1.5) +
  geom_point(size = 3, shape = 19, color = "mediumvioletred", fill = "white") +
  coord_flip() +
  labs(x = "Enriched GO Terms", y = "-log10(p-value)",
       title = "Gene Ontology Enrichment Analysis, IRGM1, SOCS1, MUC2") +
  theme_minimal() -> plotA

plotA

ggsave(filename = "figures/IRGM1_SOCS1_MUC2.jpeg", plot = plotA, 
       width = 10, height = 5, dpi = 300)

####################  MUC5AC, IL1Rn, MPO 
gene_ens <- c("ENSMUSG00000026981", #IL1RN
              "ENSMUSG00000009350", #MPO"
              "ENSMUSG00000037974") #MUC5AC

# Perform gene ontology enrichment analysis
enrich_result <- enrichGO(gene = gene_ens,
                          OrgDb = org.Mm.eg.db,
                          ont = "BP",
                          keyType = "ENSEMBL")

# View the enrichment result
enrich_result

# Extract the relevant columns from the enrichment table
enriched_terms <- enrich_result$Description
p_values <- enrich_result$p.adjust
gene_ratio <- enrich_result$GeneRatio

# Sort the enriched terms based on p-values
sorted_terms <- as.data.frame(enriched_terms[order(p_values)])

# Define the vector of values for which you want to filter
filter_terms <- c("interleukin-1-mediated signaling pathway",
                   "negative regulation of cytokine-mediated signaling pathway",
                   "inflammatory response to antigenic stimulus",
                   "negative regulation of response to cytokine stimulus",
                   "acute inflammatory response")

# Filter the dataframe based on the column values
sorted_terms <- sorted_terms %>%
    filter(`enriched_terms[order(p_values)]` %in% filter_terms)


# Create a data frame for the bar plot
bar_data <- data.frame(GO_Term = sorted_terms[1:5,], p_value = -log10(p_values[1:5]))

# Sort the data frame in ascending order
bar_data <- bar_data[order(bar_data$p_value), ]

# Create the bar plot using ggplot2
ggplot(bar_data, aes(x = reorder(GO_Term, p_value), y = p_value, fill = p_value)) +
  geom_segment(aes(xend = GO_Term, yend = 0), color = "mediumvioletred", size = 1.5) +
  geom_point(size = 3, shape = 19, color = "mediumvioletred", fill = "white") +
  coord_flip() +
  labs(x = "Enriched GO Terms", y = "-log10(p-value)",
       title = "Gene Ontology Enrichment Analysis, MUC5AC, IL1RN, MPO") +
  theme_minimal() -> plotB


plotB
ggsave(filename = "figures/MUC5AC_IL1Rn_MPO.jpeg", plot = plotB, 
       width = 10, height = 5, dpi = 300)


############################ IL13
gene_ens <-  "ENSMUSG00000020383" #IL13"

# Perform gene ontology enrichment analysis
enrich_result <- enrichGO(gene = gene_ens,
                          OrgDb = org.Mm.eg.db,
                          ont = "BP",
                          keyType = "ENSEMBL")

# View the enrichment result
enrich_result

# Extract the relevant columns from the enrichment table
enriched_terms <- enrich_result$Description
p_values <- enrich_result$p.adjust
gene_ratio <- enrich_result$GeneRatio

# Sort the enriched terms based on p-values
sorted_terms <- as.data.frame(enriched_terms[order(p_values)])


# Define the vector of values for which you want to filter
filter_terms <- c("negative regulation of cytokine production",
                  "cell activation involved in immune response",
                  "regulation of inflammatory response",
                  "positive regulation of leukocyte activation")

# Filter the dataframe based on the column values
sorted_terms <- sorted_terms %>%
    filter(`enriched_terms[order(p_values)]` %in% filter_terms)

# Create a data frame for the bar plot
bar_data <- data.frame(GO_Term = sorted_terms[1:4,], p_value = -log10(p_values[1:4]))

# Sort the data frame in ascending order
bar_data <- bar_data[order(bar_data$p_value), ]

# Create the bar plot using ggplot2
ggplot(bar_data, aes(x = GO_Term, y = p_value, fill = p_value)) +
  geom_segment(aes(xend = GO_Term, yend = 0), color = "mediumvioletred", size = 1.5) +
  geom_point(size = 3, shape = 19, color = "mediumvioletred", fill = "white") +
  coord_flip() +
  labs(x = "Enriched GO Terms", y = "-log10(p-value)",
       title = "Gene Ontology Enrichment Analysis, IL-13") +
  theme_minimal() -> plotC

ggsave(filename = "figures/IL13.jpeg", plot = plotC, 
       width = 10, height = 5, dpi = 300)


################## TICAM1, NCR1, PRF1, CXCR3, RETNLB, IL.6, CXCL9, CASP1, MYD88, TNF
gene_ens <-  c("ENSMUSG00000047123", #TICAM1
               "ENSMUSG00000062524", #NCR1
               "ENSMUSG00000037202.7", #PRF1
               "ENSMUSG00000050232", #CXCR3"
               "ENSMUSG00000022650", #RETNLB
               "ENSMUSG00000025746", #IL6
               "ENSMUSG00000029417", #CXCL9
               "ENSMUSG00000025888", #CASP1
               "ENSMUSG00000032508", #MYD88
               "ENSMUSG00000024401") #TNF 


# Perform gene ontology enrichment analysis
enrich_result <- enrichGO(gene = gene_ens,
                          OrgDb = org.Mm.eg.db,
                          ont = "BP",
                          keyType = "ENSEMBL")

# View the enrichment result
enrich_result

# Extract the relevant columns from the enrichment table
enriched_terms <- enrich_result$Description
p_values <- enrich_result$p.adjust
gene_ratio <- enrich_result$GeneRatio

# Sort the enriched terms based on p-values
sorted_terms <- as.data.frame(enriched_terms[order(p_values)])
sorted_terms[1:50,]


# Create a data frame for the bar plot
bar_data <- data.frame(GO_Term = sorted_terms[1:15,], p_value = -log10(p_values[1:15]))

# Sort the data frame in ascending order
bar_data <- bar_data[order(bar_data$p_value), ]

# Create the bar plot using ggplot2
ggplot(bar_data, aes(x = reorder(GO_Term, p_value), y = p_value, fill = p_value)) +
  geom_segment(aes(xend = GO_Term, yend = 0), color = "mediumvioletred", size = 1.5) +
  geom_point(size = 3, shape = 19, color = "mediumvioletred", fill = "white") +
  coord_flip() +
  labs(x = "Enriched GO Terms", y = "-log10(p-value)",
       title = "Gene Ontology Enrichment Analysis: 
       TICAM1, NCR1, PRF1, CXCR3, RETNLB, IL.6, CXCL9, CASP1, MYD88, TNF") +
  theme_minimal() -> plotD

ggsave(filename = "figures/most_immune_genes.jpeg", plot = plotD, 
       width = 10, height = 5, dpi = 300)

## save the variance contribution of each gene 
##save the normalized data 
write.csv(pca_var, "Data/Data_output/variance_contr_gene_lab.csv", row.names = TRUE)

# save the lab data frame for figures
write.csv(lab, "Data/Data_output/lab_pca.csv", row.names = FALSE)


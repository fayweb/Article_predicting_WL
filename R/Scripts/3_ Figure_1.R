library(ggrepel)
library(tidyverse)
library(tidyr)
library(dplyr)
library(scales)
library(cowplot)
library(ggthemes)
library(grid)
library(ggplot2)
library(ggpmisc)
library(broom)
library(knitr)
library(kableExtra)
library(webshot)
library(RColorBrewer)
library(ggeffects)
library(pheatmap)
library(modelsummary)
library(pdp)
library(gt)
library(ggdist)
library(broom)
library(reshape2)
library(knitr)
library(stargazer)
library(kableExtra)
library(sjmisc)
library(sjlabelled)
library(jtools)
library(sjPlot)
library(FactoMineR)
library(ggridges)
library(Polychrome)
library(patchwork)
library(gridExtra)
library(cowplot)
library(patchwork)
library(ggpubr)
library(factoextra)
library(jtools)


# read the data
hm <- read.csv("Data/Data_output/imputed_clean_data.csv")


# Select laboratory data 
# Select genes
lab <- hm %>%
    dplyr::filter(origin == "Lab")

# create a vector to select genes
Genes_v   <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",
               "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
               "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
               "TICAM1", "TNF") #"IL.12", "IRG6")

# create a data frame containing the continuous gene expression variables
genes <- lab[ ,colnames(lab) %in% Genes_v]


# increase maximum overlaps
options(ggrepel.max.overlaps = Inf)

# PCA
## we can now run a normal pca on the complete data set
res.pca <- PCA(genes)

# Convert mouse_id to a data frame (to facilitate data joining)
mouse <- data.frame(Mouse_ID = lab[,1])

# Add the new column pc1 to the mouse_id data frame
mouse$PC1 <- res.pca$ind$coord[, 1]

mouse$PC2 <- res.pca$ind$coord[, 2]  # indexing the second column
mouse$PC3 <-  res.pca$ind$coord[, 3]
mouse$PC4 <-  res.pca$ind$coord[, 4]
mouse$PC5 <-  res.pca$ind$coord[, 5]
# join the coordinates
lab <- lab %>% 
    left_join(mouse, by = "Mouse_ID")

## We also need to extract the data for the variable contributions to each of 
## the pc axes

# read the variance explained by each gene for the pca 
vpg <- read.csv("Data/Data_output/variance_contr_gene_lab.csv")

# Change the first column of the variance contribution of variables to the gene
# names
vpg <- vpg %>%
  dplyr::rename(Variable = vars, PC1 = Dim.1, PC2 = Dim.2)

# add cos2 to lab
lab <- lab %>% mutate(cos2 = lab$PC1^2 + lab$PC2^2)

# Then, define the color for each level of infection
color_mapping <- c("E_falciformis" = "salmon", 
                   "E_ferrisi" = "forestgreen", 
                   "uninfected" = "cornflowerblue")

# PCA graph of individuals
pca_individuals <-
    ggplot(lab, aes(x = PC1, y = PC2, color = current_infection)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_point(size = 5, alpha = 0.5, color = "black",  shape = 21, aes(fill = current_infection)) +
    labs(x = "PC1 (34.37%)", y = "PC2 (16.03%)",# title = "PCA graph of individuals",
         colour = "Current infection") +
    theme_minimal() +
    theme(#plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "right") +
    scale_color_manual(values = color_mapping)# +
    #guides(color = guide_legend(override.aes = list(size = 4))) 

#pca_individuals

ggsave(filename = "figures/pca_individuals.jpeg", plot = pca_individuals, 
       width = 6, height = 4, dpi = 300)


####### PCA graph of variables


# Add cos2 variable to the dataframe
vpg$cos2 <- with(vpg, PC1^2 + PC2^2)


# Define custom gradient colors
gradient_colors <- c("#B00B69", "#A55EA7", "#1D1CC9")

# Define the breaks and labels for the color legend
breaks <- c(0, 50, 100, 150)
labels <- c("0", "50", "100", "150")

# Plotting the factor map 
pca_variables <-
  ggplot(vpg, aes(x = PC1, y = PC2, color = cos2)) +
  geom_segment(aes(xend = 0, yend = 0), color = "gray50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(size = 3) +
  geom_label_repel(aes(label = Variable), size = 2, box.padding = 0.5, max.overlaps = Inf) +
  coord_equal() +
  xlab("PC1 (34.37%)") +
  ylab("PC2 (16.03%") +
  #ggtitle("PCA Plot of Variables") +
  theme_minimal()  +
  scale_color_gradient(low = "blue", high = "orange")+
    labs(color = "Squared distance to origin")

#pca_variables

ggsave(filename = "figures/pca_variables.jpeg", plot = pca_variables, 
       width = 5, height = 6, dpi = 300)
################################################################
# Create a custom color palette for 19 genes
# build-in color palette
#display.brewer.all(colorblindFriendly = TRUE)

color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(19)

#pca_variables <-
    ggplot(vpg, aes(x = PC1, y = PC2, color = Variable)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # Segment and points
    geom_segment(aes(xend = 0, yend = 0), color = "gray50") +
    geom_point(size = 3,  color = "black", shape = 21, aes(fill = Variable)) +
    # Labels
    geom_label_repel(
        aes(label = Variable, fill = Variable), size = 3.5, box.padding = 0.5, 
        max.overlaps = Inf,
        color = "white",  # Color for the text inside the label
        segment.color = "black") +  # Color for the connecting lines
    # Axes and theme
    coord_equal() +
    xlab("PC1 (34.37%)") +
    ylab("PC2 (16.03%)") +
    theme_minimal() + 
    theme(legend.position = "none") +
   # labs(title = "PCA graph of variabes") +
    # Coloring for the 19 genes
    scale_color_manual(values = color_palette)

#print(pca_variables)


#ggsave(filename = "figures/pca_variables.jpeg", plot = pca_variables, 
#       width = 12, height = 6, dpi = 600)


fviz_pca_biplot(res.pca, 
                col.ind = lab$current_infection, palette = c("E_falciformis" = "salmon", 
                                                             "E_ferrisi" = "forestgreen", 
                                                             "uninfected" = "cornflowerblue"),
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Infection groups",
                title = "") -> biplot

biplot

ggsave(filename = "figures/biplot.jpeg", plot = biplot, 
       width = 12, height = 6, dpi = 600)



# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 18, 
             title = "PC1", 
             fill =  "seagreen2") -> contributions_pc1

contributions_pc1

ggsave(filename = "figures/contributions_pc1.jpeg", plot = contributions_pc1, 
       width = 6, height = 4, dpi = 1000)
# res.pca$var$contrib


### Contributions to the second dimension

## Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 18, 
             title = "PC2",
             fill =  "seagreen2") -> contributions_pc2

contributions_pc2

ggsave(filename = "figures/contributions_pc2.jpeg", plot = contributions_pc2, 
       width = 6, height = 4, dpi = 1000)


##################################################################################################
#################################################
#################################################
#################################################

# Load the required packages
###PC1 PC2 linear regression
lab$current_infection <- factor(lab$current_infection, 
        levels = c("uninfected", "E_falciformis", "E_ferrisi"))

lab$mouse_strain <- as.factor(lab$mouse_strain)
lab$immunization <- factor(lab$immunization, levels = 
                                    c("naive", "uninfected",
                                      "heterologous", "homologous"))

model_1 <- lm(WL_max ~ PC1 + PC2 +
               current_infection + delta_ct_cewe_MminusE +
                mouse_strain + immunization + 
                weight_dpi0, data = lab )

summary(model_1)

tab_model(model_1)

stargazer(model_1,
          type = "text", out = "tables/predictors_weightloss.txt", 
          title = "Linear models - Predicting maximum weight loss")

#, file = "tables/predicting_weight_loss_linear.doc")

# Extract the residuals from the model
residuals <- resid(model_1)

# Create a data frame with the residuals
residuals_df <- data.frame(residuals = residuals)

# Create the QQ plot
residuals_1 <-
    ggplot(residuals_df, aes(sample = residuals)) +
    stat_qq(color = "blue") +
    ggtitle("QQ Plot of Residuals") +
    xlab("Theoretical Quantiles") +
    ylab("Sample Quantiles") +
    theme_minimal()

residuals_1

ggsave(filename = "figures/residuals_model_1.jpeg", 
       plot = residuals_1, 
       width = 12, height = 6, dpi = 600)

# Extract the fitted values from the model
fitted_values <- fitted(model_1)

# Create a data frame with the residuals and the fitted values
data_df <- data.frame(residuals = residuals, fitted_values = fitted_values)

# Create the scatter plot
residuals_vs_fitted <-
    ggplot(data_df, aes(x = fitted_values, y = residuals)) +
    geom_point(color = "blue") +
    ggtitle("Residuals vs Fitted Values") +
    xlab("Fitted Values") +
    ylab("Residuals") +
    theme_minimal()

residuals_vs_fitted

ggsave(filename = "figures/residuals_vs_fitted.jpeg", 
       plot = residuals_vs_fitted, 
       width = 12, height = 6, dpi = 600)

#########
# without parasite data
#model_2 <- lm(WL_max ~ PC1 + PC2 + mouse_strain + weight_dpi0, data = lab)
#summary(model_2)

# without host data
model_2 <- lm(WL_max ~ PC1 + PC2 + current_infection + delta_ct_cewe_MminusE, data = lab)

summary(model_2)

# only pc1 + pc2
model_3 <- lm(WL_max ~ PC1 + PC2 , data = lab)

summary(model_3)



plot_coefs(model_1, model_2, model_3)


model_5 <- lm(WL_max ~ PC1 + PC2 + current_infection + delta_ct_cewe_MminusE + 
                  weight_dpi0, data = lab)

summary(model_5)
plot_coefs(model_1, model_2, model_3, model_4, model_5)

## Please cite as:
##  Hlavac, Marek (2018). stargazer: Well-Formatted Regression and Summary Statistics Tables.
stargazer(model_1, model_2, model_3,
          type = "text",
          out = "tables/stargazer.docx", 
          title = "Linear models - Predicting maximum weight loss",
          align = TRUE)

export_summs(model_1, model_2, model_3,
             scale = TRUE, to.file = "docx", 
             file.name = "tables/lab_model1_3.docx")


models <- list(
    "Model 1" = model_1,
    "Model 2" = model_2,
    "Model 3" = model_3)

modelsummary(models, stars = TRUE, gof_omit = "IC|Adj|F|RMSE|Log", 
             output = "tables/lab_model1_3.docx")

modelsummary(models = as.list(model_1, model_2, model_3), 
             output = "tables/lab_model1.docx")



model_4 <- lm(WL_max ~ PC1 * current_infection + PC2 *current_infection, 
              data = lab)

summary(model_4)
plot_summs(model_4) -> coefs4

coefs4

ggsave(filename = "figures/plot_sums_mix_PCA.jpeg", plot = coefs4, width = 6, height = 4)
#see the ggefects
effects <- ggpredict(model_4)

pc1_current_infection <- 
    ggpredict(model_4, terms = c("PC1")) %>% 
    plot(colors = "darkorchid") +   # Use a refined shade of blue
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = NULL) +  # This removes the title
  #  ggtitle("Effect of PC1 on Predicted Weight Loss") +
    xlab("Principal Component 1 (PC1)") +
    ylab("Predicted values of weight loss") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    )


pc1_current_infection

ggsave(filename = "figures/pc1_current_infection.jpeg", 
       plot = pc1_current_infection, 
       width = 6, height = 4, dpi = 1000)

pc2_current_infection <- 
    ggpredict(model_4, terms = c("PC2")) %>% 
    plot(colors = "darkorchid") +   # Use a refined shade of blue
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = NULL) +  # This removes the title
   # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Principal Component 2 (PC2)") +
    ylab("Predicted values of weight loss") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    )

pc2_current_infection

ggsave(filename = "figures/pc2_current_infection.jpeg", 
       plot = pc2_current_infection, 
       width = 6, height = 4, dpi = 1000)


plot_summs(model_4)
modelsummary(model_4, stars = TRUE, gof_omit = "IC|Adj|F|RMSE|Log", 
             output = "tables/mixed_effects_pca.docx")

# Now create the scatter plot using this color mapping
ggpredict(model_4, terms = c("PC1", "current_infection")) %>% 
    plot() +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Principal Component 1 (PC1)") +
    ylab("Predicted values of weight loss") +
    theme_minimal() +
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = color_mapping) +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    ) -> pc1_WL_current_infection

pc1_WL_current_infection

ggsave("figures/pc1_WL_current_infection.jpeg", pc1_WL_current_infection, width = 8, height = 6, dpi = 1000)


# Now create the scatter plot using this color mapping
# Now create the scatter plot using this color mapping
ggpredict(model_4, terms = c("PC2", "current_infection")) %>% 
    plot() +  
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Principal Component 2 (PC2)") +
    ylab("Predicted values of weight loss") +
    theme_minimal() +
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = color_mapping) +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    ) -> pc2_WL_current_infection

pc2_WL_current_infection

ggsave("figures/pc2_WL_current_infection.jpeg", pc2_WL_current_infection, width = 8, height = 6, dpi = 1000)



# Combine the figures
# Combine the plots
panel_figure <- 
    (pca_variables| biplot)  +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
panel_figure <- panel_figure + 
    plot_annotation(title = 'Fig. 3', 
                    theme = theme(plot.title = element_text(size = 20, hjust = 0))) #+
   # plot_layout(heights = c(10, 1, 1,1), 
       #         widths = c(1,1,1,1)) # Adjust according to your layout needs

# Display the panel figure
#print(panel_figure)

# Save the panel figure
ggsave("figure_panels/pca_biplot_variable.jpeg", 
       panel_figure, width = 12, height = 5, dpi = 300)

# Combine the plots
panel_figure <- 
    (pca_variables| biplot)  +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
panel_figure <- panel_figure + 
    plot_annotation(title = 'Fig. 3', 
                    theme = theme(plot.title = element_text(size = 20, hjust = 0))) #+
   # plot_layout(heights = c(10, 1, 1,1), 
       #         widths = c(1,1,1,1)) # Adjust according to your layout needs

# Display the panel figure
#print(panel_figure)

# Save the panel figure
ggsave("figure_panels/pca_biplot_variable.jpeg", 
       panel_figure, width = 12, height = 5, dpi = 300)



################################
############################
# Figure 4
# Combine the plots
panel_figure2 <- 
    (contributions_pc1 | contributions_pc2 ) /
    (pc1_WL_current_infection | pc2_WL_current_infection) /
    free(coefs4) +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
panel_figure2 <- panel_figure2 + 
    plot_annotation(title = 'Fig. 5', 
                    theme = theme(plot.title = element_text(size = 20, hjust = 0))) +
    plot_layout(heights = c(1, 1,1), 
                        widths = c(1,1,1))

# Save the panel figure
ggsave("figure_panels/mixed_models_PCA.jpeg", 
       panel_figure2, width = 13, height = 12, dpi = 300)

figure_panel_1 <- ggarrange(pca_variables, biplot, 
                            pc1_current_infection, pc2_current_infection,
                            pc1_WL_current_infection, pc2_WL_current_infection,
                            labels = c("A", "B", "C", "D", "E", "F"),
                            ncol = 2, nrow = 3,
                            )

# Adding the title "Figure 1" to the entire arrangement
figure_panel_1 <- annotate_figure(figure_panel_1, 
                                  top = text_grob("Fig. 2", size = 14, 
                                                  face = "bold"))

print(figure_panel_1)


ggsave("figure_panels/pca_panel.jpeg", figure_panel_1, 
       width = 11, height = 9, dpi = 300)


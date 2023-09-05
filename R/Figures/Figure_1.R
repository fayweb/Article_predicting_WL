library(ggrepel)
library(tidyverse)
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
library(pdp)
library(broom)
library(reshape2)
library(knitr)
library(stargazer)
library(kableExtra)


# read the lab data with pca vectors
lab <- read.csv("Data/Data_output/lab_pca.csv")

# change the labels pc1 and pc2 to PC1 / PC2
lab <- lab %>%
  dplyr::rename(PC1 = "pc1", PC2 = "pc2")

# read the variance explained by each gene for the pca 
vpg <- read.csv("Data/Data_output/variance_contr_gene_lab.csv")

# Change the first column of the variance contribution of variables to the gene
# names
vpg <- vpg %>%
  dplyr::rename(Variable = vars, PC1 = Dim.1, PC2 = Dim.2)

# add cos2 to lab
lab <- lab %>% mutate(cos2 = lab$PC1^2 + lab$PC2^2)

# Define color palette
color_palette <- c("E_ferrisi" = "#66C2A5", "uninfected" = "#8DA0CB", 
                   "E_falciformis" = "#FC8D62")

# PCA graph of individuals
pca_individuals <-
  ggplot(lab, aes(x = PC1, y = PC2, color = infection, shape = infection)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = "PC1 (52.32%)", y = "PC2 (11.79%)", #title = "PCA graph of individuals",
       colour = "Current infection", shape ="Current infection") +
  theme_minimal() +
  theme(#plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "right") +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = c("E_ferrisi" = 17, "uninfected" = 16, "E_falciformis" = 18)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

pca_individuals

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
  geom_point(size = 3) +
  geom_label_repel(aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf) +
  coord_equal() +
  xlab("PC1 (52.32%)") +
  ylab("PC2 (11.79%)") +
  #ggtitle("PCA Plot of Variables") +
  theme_minimal() + 
  #theme(legend.position = "right",
   #plot.title = element_text(size = 12, face = "bold")) +
  guides(color = guide_colorbar(title = "Squared Distance from Origin")) +
  scale_color_gradientn(colors = gradient_colors, guide = "none")  


ggsave(filename = "figures/pca_variables.jpeg", plot = pca_variables, 
       width = 12, height = 6, dpi = 600)


#pca_variables

######################## Enriched Terms data frame
enriched_terms_df <- read.csv("Data/Data_output/enriched_sorted_terms.csv")


# First, transform the p-values to -log10 scale
enriched_terms_df$p_value_log <- -log10(enriched_terms_df$p_value)


# Create the lollipop plot using ggplot2
enrichment_terms_plot <- 
  ggplot(enriched_terms_df[1:15,], aes(x = reorder(GO_Term, p_value_log), y = p_value_log)) +
  geom_segment(aes(xend = GO_Term, yend = 0, color = p_value_log), linewidth = 1.5) +
  geom_point(aes(fill = p_value_log), size = 3, shape = 21, color = "mediumvioletred") +
  scale_fill_gradientn(colours = rev(gradient_colors)) +
  scale_color_gradientn(colours = rev(gradient_colors)) +
  coord_flip() +
  labs(x = "Enriched GO Terms", y = "-log10(p-value)",
       title = "Gene Ontology Enrichment Analysis") +
  theme_minimal() +
    theme( plot.title = element_text(size = 12, face = "bold")) 

enrichment_terms_plot

ggsave(filename = "figures/enrichment_terms_plot.jpeg", plot = enrichment_terms_plot, 
       width = 12, height = 6, dpi = 600)


pca_individuals <- pca_individuals + coord_fixed(ratio = 1)



##################################################################################################
#################################################
#################################################
#################################################

# Load the required packages
###PC1 PC2 linear regression
lab <- lab %>%
    dplyr::mutate(current_infection = case_when(
        current_infection == "infected_eimeria" ~ "E_falciformis",
        TRUE ~ current_infection))

model_1 <- lm(WL_max ~ PC1 + PC2 + current_infection + delta_ct_cewe_MminusE +
                mouse_strain + immunization + 
                weight_dpi0, data = lab )

summary(model_1)

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
    ylab("Sample Quantiles")

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
    ylab("Residuals")

residuals_vs_fitted

ggsave(filename = "figures/residuals_vs_fitted.jpeg", 
       plot = residuals_vs_fitted, 
       width = 12, height = 6, dpi = 600)

#########
# without parasite data
model_2 <- lm(WL_max ~ PC1 + PC2 + mouse_strain + weight_dpi0, data = lab)
summary(model_2)

# without host data
model_3 <- lm(WL_max ~ PC1 + PC2 + current_infection + delta_ct_cewe_MminusE +
                  immunization + weight_dpi0, data = lab)

summary(model_3)

# only pc1 + pc2
model_4 <- lm(WL_max ~ PC1 + PC2 , data = lab)

summary(model_4)

## Please cite as:
##  Hlavac, Marek (2018). stargazer: Well-Formatted Regression and Summary Statistics Tables.
stargazer(model_1, model_2, model_3, model_4, 
          type = "html", out = "figures/predictors_weightloss.html", 
          title = "Linear models - Predicting maximum weight loss")

#correcting for nas in delta ct
model_2 <- lm(WL_max ~ PC1 + PC2 + mouse_strain + weight_dpi0, data = lab %>% 
                  drop_na(delta_ct_cewe_MminusE))

model_4 <- lm(WL_max ~ PC1 + PC2 , data = lab %>% 
                  drop_na(delta_ct_cewe_MminusE))

# Anova of different models
anova_mod <- anova(model_1, model_2, model_3, model_4)

stargazer(anova_mod, type = "html", out = "figures/anova_model.html", title = 
              "Analysis of Variance Table")

#see the ggefects
effects <- ggpredict(model_4)

pc1_current_infection <- 
    ggpredict(model_4, terms = c("PC1")) %>% 
    plot(colors = "blue")

pc1_current_infection

ggsave(filename = "figures/pc1_current_infection.jpeg", 
       plot = pc1_current_infection, 
       width = 12, height = 6, dpi = 600)

pc2_current_infection <- 
    ggpredict(model_4, terms = c("PC2")) %>% 
    plot(colors = "blue")

pc2_current_infection

ggsave(filename = "figures/pc2_current_infection.jpeg", 
       plot = pc2_current_infection, 
       width = 12, height = 6, dpi = 600)


# produce the table without levels (immunization and mouse_strains)
#not possible


# Combine the figures

figure_panel_1 <- 
    plot_grid(pca_variables, pca_individuals, 
              pc1_current_infection, pc2_current_infection,
              cols  = 2, rel_heights = c(3, 2), 
              labels = c("A", "B", "C", "D"))

ggsave("figure_panels/figure_panel_1.jpeg", figure_panel_1, 
       width = 16, height = 10, dpi = 300)


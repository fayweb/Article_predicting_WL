library(ggrepel)
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
#pca_individuals <- 
pca_individuals <-
  ggplot(lab, aes(x = PC1, y = PC2, color = infection, shape = infection)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = "PC1", y = "PC2", title = "PCA graph of individuals",
       colour = "Current infection", shape ="Current infection") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "right") +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = c("E_ferrisi" = 17, "uninfected" = 16, "E_falciformis" = 18)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

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
#pca_variables <-
  ggplot(vpg, aes(x = PC1, y = PC2, color = cos2)) +
  geom_segment(aes(xend = 0, yend = 0), color = "gray50") +
  geom_point(size = 3) +
  geom_label_repel(aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf) +
  coord_equal() +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA Plot of Variables") +
  theme_minimal() +
  theme(legend.position = "right") +
  guides(color = guide_colorbar(title = "Squared Distance from Origin")) +
  scale_color_gradientn(colors = gradient_colors, guide = "none")  +
  theme(plot.title = element_text(size = 12))


ggsave(filename = "figures/pca_variables.jpeg", plot = pca_variables, 
       width = 12, height = 6, dpi = 600)

enriched_terms_df <- read.csv("Data/Data_output/enriched_sorted_terms.csv")


# First, transform the p-values to -log10 scale
enriched_terms_df$p_value_log <- -log10(enriched_terms_df$p_value)


# Create the lollipop plot using ggplot2
enrichment_terms_plot <- 
ggplot(enriched_terms_df[1:30,], aes(x = reorder(GO_Term, p_value_log), y = p_value_log)) +
  geom_segment(aes(xend = GO_Term, yend = 0, color = p_value_log), size = 1.5) +
  geom_point(aes(fill = p_value_log), size = 3, shape = 21, color = "mediumvioletred") +
  scale_fill_gradientn(colours = rev(gradient_colors)) +
  scale_color_gradientn(colours = rev(gradient_colors)) +
  coord_flip() +
  labs(x = "Enriched GO Terms", y = "-log10(p-value)",
       title = "Gene Ontology Enrichment Analysis") +
  theme_minimal()

ggsave(filename = "figures/enrichment_terms_plot.jpeg", plot = enrichment_terms_plot, 
       width = 12, height = 6, dpi = 600)



#################################################
# Load the required packages
###PC1 PC2 linear regression
model <- lm(WL_max ~ PC1 + PC2, data = lab)
summary(model)
ggeffect(model)


# Generate predicted effects
effects <- ggeffects(model)

# Generate predicted effects
predicted_effects <- ggeffects(model)

# Plot for PC1
plot1 <- ggplot(predicted_effects[["PC1"]], aes(x = PC1, y = predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  labs(title = "Predicted Effects of PC1 on WL_max",
       x = "PC1", y = "Predicted WL_max")

print(plot1)

# Plot for PC2
plot2 <- ggplot(predicted_effects[["PC2"]], aes(x = PC2, y = predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  labs(title = "Predicted Effects of PC2 on WL_max",
       x = "PC2", y = "Predicted WL_max")

print(plot2)



model_1 <- lm(WL_max ~ PC1 + PC2 + primary_infection * challenge_infection, data = lab)
summary(model)

model_2 <- lm(WL_max ~ PC1 + PC2 + challenge_infection, data = lab)
summary(model)

anova(model_1, model_2)

model_3 <- lm(WL_max ~ PC1 + PC2, data = lab)
summary(model_3)

anova(model_1, model_2, model_3)

# Load the required packages
###PC1 PC2 linear regression
model <- lm(WL_max ~ PC1 + PC2, data = lab)
summary(model)
ggeffect(model)

plot(ggeffect(model))

# Generate equation text
eq_text <- paste("WL_max =", round(coef(model)[1], 2),
                 "+", round(coef(model)[2], 2), "PC1",
                 "+", round(coef(model)[3], 2), "PC2")

# Calculate R-squared
predicted <- predict(model)
r_squared <- summary(model)$r.squared
r_squared_text <- paste("R-squared =", round(r_squared, 2))

# Plot the data with the equation and R-squared
ggplot(lab, aes(x = predicted, y = WL_max)) +
    geom_point(color = "#336699", size = 3, alpha = 0.6) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#990000", 
                size = 0.8) +
    labs(x = "Predicted", y = "Observed") +
    ggtitle("PC1 and PC2 predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "none") +
    annotate("text", x = max(predicted), y = min(lab$WL_max),
             label = eq_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    stat_poly_eq(
        formula = y ~ x,
        label.x = "right", label.y = "bottom",
        label = paste("R^2 =", round(r_squared, 2)),
        parse = TRUE,
        size = 4,
        family = "serif",
        fontface = "bold",
        aes(label = paste("R^2 =", round(summary(model)$r.squared, 2))),
        label.x.npc = 0.70, label.y.npc = 0.4
    )# -> linear_pc1_pc2_WL


# Tidy model summary
tidy_model <- tidy(model)

# Create a table of model coefficients
coef_table <- tidy_model %>%
    mutate(
        term = ifelse(term == "(Intercept)", "Intercept", term),
        Estimate = round(estimate, 2),
        `Std. Error` = round(std.error, 2),
        `t value` = round(statistic, 2),
        `Pr(>|t|)` = p.value
    ) %>%
    select(term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`) %>%
    rename(`Pr(>|t|)` = `Pr(>|t|)`)

# Print the coefficient table
write.csv(coef_table, "tables/coefficient_table_pc1_pc2_linear.csv", 
          row.names = FALSE)

ggsave(filename = "figures/linear_pc1_pc2_WL.jpeg", plot = linear_pc1_pc2_WL, 
       width = 6, height = 4, dpi = 300)

# Calculate residuals
residuals <- lab$WL_max - predicted

# Generate equation text for residuals
eq_text_residuals <- paste("Residuals =", round(coef(model)[1], 2),
                           "+", round(coef(model)[2], 2), "PC1",
                           "+", round(coef(model)[3], 2), "PC2")

# Plot the residuals
residuals_pc1_pc2_WL <-
    ggplot(lab, aes(x = predicted, y = residuals)) +
    geom_point(color = "#336699", size = 3, alpha = 0.6) +
    geom_hline(yintercept=0, color = "#990000", 
               size = 0.8) +
    labs(x = "Predicted", y = "Residuals") +
    ggtitle("Residual plot of PC1 and PC2 predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "none") +
    annotate("text", x = max(predicted), y = min(residuals),
             label = eq_text_residuals, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") 


ggsave(filename = "figures/residuals_pc1_pc2_WL.jpeg", plot = residuals_pc1_pc2_WL, 
       width = 12, height = 6, dpi = 600)


#########################################################################
#### PC1 + PC2 + heter/hom infections predicting WL
# create a new variable that shows if an infection is primary or homologous
# or heterologous
lab <- lab %>% 
    dplyr::mutate(
        immunization = 
                      case_when(
                          infection_history == "falciformis_ferrisi" ~ "heterologous",
                          infection_history == "ferrisi_falciformis" ~ "heterologous",
                          infection_history == "falciformis_uninfected" ~ "uninfected",
                          infection_history == "ferrisi_uninfected" ~ "uninfected",
                          infection_history == "ferrisi_ferrisi" ~ "homologous",
                          infection_history == "falciformis_falciformis" ~ "homologous",
                          infection_history == "uninfected_falciformis" ~ "naive",
                          infection_history == "uninfected_ferrisi" ~ "naive",
                          infection_history == "uninfected" ~ "uninfected",
                          TRUE ~ "NA"
                          ))

# Perform linear regression
model <- lm(WL_max ~ PC1 + PC2 + immunization, data = lab)

# Generate equation text
eq_text <- paste("WL_max =", round(coef(model)[1], 2),
                 "+", round(coef(model)[2], 2), "PC1",
                 "+", round(coef(model)[3], 2), "PC2")

# Calculate R-squared
predicted <- predict(model)
r_squared <- summary(model)$r.squared
r_squared_text <- paste("R-squared =", round(r_squared, 2))

# Plot the data with the equation and R-squared
linear_pc1_pc2_immunization <-
    ggplot(lab, aes(x = predicted, y = WL_max, color = immunization)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#990000",
                size = 0.8) +
    labs(x = "Predicted", y = "Observed") +
    ggtitle("PC1, PC2, and immunization predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(lab$WL_max),
             label = eq_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    annotate("text", x = max(predicted), y = min(lab$WL_max) - 2,
             label = r_squared_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_brewer(palette = "Set1")  # Use a color palette from RColorBrewer

ggsave(filename = "figures/linear_pc1_pc2_immunization.jpeg", 
       plot = linear_pc1_pc2_immunization, 
       width = 12, height = 6, dpi = 600)

# Calculate residuals
residuals <- lab$WL_max - predicted

# Generate equation text for residuals
eq_text_residuals <- paste("Residuals =", round(coef(model)[1], 2),
                           "+", round(coef(model)[2], 2), "PC1",
                           "+", round(coef(model)[3], 2), "PC2")

# Plot the residuals
residuals_pc1_pc2_immunization <-
    ggplot(lab, aes(x = predicted, y = residuals, color = immunization)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_hline(yintercept=0, color = "#990000", 
               size = 0.8) +
    labs(x = "Predicted", y = "Residuals") +
    ggtitle("Residual plot of PC1, PC2, and Infection Type predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(residuals),
             label = eq_text_residuals, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_brewer(palette = "Set1")  # Use a color palette from RColorBrewer



ggsave(filename = "figures/residuals_pc1_pc2_immunization.jpeg", 
       plot = residuals_pc1_pc2_immunization, 
       width = 12, height = 6, dpi = 600)

######################################## just immunization
model_immunization <- lm(WL_max ~ immunization, data = lab)

# Generate equation text
eq_text <- paste("WL_max =", round(coef(model_immunization)[1], 2), " + ", 
                 "immunization")

# Calculate R-squared
predicted <- predict(model_immunization)
r_squared <- summary(model_immunization)$r.squared
r_squared_text <- paste("R-squared =", round(r_squared, 2))

# Plot the data with the equation and R-squared
linear_immunization <-
    ggplot(lab, aes(x = predicted, y = WL_max, color = immunization)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#990000",
                size = 0.8) +
    labs(x = "Predicted", y = "Observed") +
    ggtitle("Infection Type predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(lab$WL_max),
             label = eq_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    annotate("text", x = max(predicted), y = min(lab$WL_max) - 2,
             label = r_squared_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_brewer(palette = "Set1")  # Use a color palette from RColorBrewer

ggsave(filename = "figures/linear_immunization.jpeg", 
       plot = residuals_pc1_pc2_immunization, 
       width = 12, height = 6, dpi = 600)

# Calculate residuals
residuals <- lab$WL_max - predicted

# Generate equation text for residuals
eq_text_residuals <- paste("Residuals =", round(coef(model_immunization)[1], 2), " + ", 
                           "immunization")

# Plot the residuals
residuals_immunization <-
    ggplot(lab, aes(x = predicted, y = residuals, color = immunization)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_hline(yintercept=0, color = "#990000", 
               size = 0.8) +
    labs(x = "Predicted", y = "Residuals") +
    ggtitle("Residual plot of Infection Type predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(residuals),
             label = eq_text_residuals, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_brewer(palette = "Set1")  # Use a color palette from RColorBrewer


ggsave(filename = "figures/residuals_immunization.jpeg", 
       plot = residuals_immunization, 
       width = 12, height = 6, dpi = 600)

##################################### interaction of infection type 
model_interaction <- lm(WL_max ~ PC1*immunization + PC2*immunization, data = lab)

# Calculate R-squared
predicted <- predict(model_interaction)
r_squared <- summary(model_interaction)$r.squared
r_squared_text <- paste("R-squared =", round(r_squared, 2))

# Plot the data with the R-squared
ggplot(lab, aes(x = predicted, y = WL_max, color = immunization)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#990000",
                size = 0.8) +
    labs(x = "Predicted", y = "Observed") +
    ggtitle("Interaction between PC1, PC2 and Infection Type predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(lab$WL_max),
             label = r_squared_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_brewer(palette = "Set1")  # Use a color palette from RColorBrewer

# Calculate residuals
residuals <- lab$WL_max - predicted

# Plot the residuals
ggplot(lab, aes(x = predicted, y = residuals, color = immunization)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_hline(yintercept=0, color = "#990000", 
               size = 0.8) +
    labs(x = "Predicted", y = "Residuals") +
    ggtitle("Residual plot of Interaction between PC1, PC2, and Infection Type predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    scale_color_brewer(palette = "Set1")  # Use a color palette from RColorBrewer


step_model <- step(lm(WL_max ~ PC1 + PC2 + immunization, data = lab), direction="both")
summary(step_model)


##############################################immunization + current infection
###################################### interaction of immunization and Parasite_challenge
model_current_immunization <- lm(WL_max ~ immunization * Parasite_challenge, data = lab)
summary(model_current_immunization)

# Generate equation text
eq_text <- paste("WL_max =", round(coef(model_current_immunization)[1], 2), " + ", 
                 "immunization*Parasite_challenge")

# Calculate R-squared
predicted <- predict(model_current_immunization)
r_squared <- summary(model_current_immunization)$r.squared
r_squared_text <- paste("R-squared =", round(r_squared, 2))

# Plot the data with the equation and R-squared
immunization_parasite_chal_lm <-
  ggplot(lab, aes(x = predicted, y = WL_max, color = immunization)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#990000",
                size = 0.8) +
    labs(x = "Predicted", y = "Observed") +
    ggtitle("Infection Type and Parasite Challenge Interaction predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(lab$WL_max),
             label = eq_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    annotate("text", x = max(predicted), y = min(lab$WL_max) - 2,
             label = r_squared_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_brewer(palette = "Set1")  # Use a color palette from RColorBrewer 


# Calculate residuals
residuals <- lab$WL_max - predicted

# Generate equation text for residuals
eq_text_residuals <- paste("Residuals =", round(coef(model_current_immunization)[1], 2), " + ", 
                           "immunization*Parasite_challenge")

# Plot the residuals
ggplot(lab, aes(x = predicted, y = residuals, color = immunization)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_hline(yintercept=0, color = "#990000", 
               size = 0.8) +
    labs(x = "Predicted", y = "Residuals") +
    ggtitle("Residual plot of Infection Type and Parasite Challenge Interaction predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(residuals),
             label = eq_text_residuals, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_brewer(palette = "Set1")  # Use a color palette from RColorBrewer

########################### hybrid status
model_hybrid_status <- lm(WL_max ~ PC1 + PC2 + hybrid_status, data = lab)

# Calculate R-squared
predicted <- predict(model_hybrid_status)
r_squared <- summary(model_hybrid_status)$r.squared
r_squared_text <- paste("R-squared =", round(r_squared, 2))

# Generate equation text
eq_text <- paste("WL_max =", round(coef(model_hybrid_status)[1], 2),
                 "+", round(coef(model_hybrid_status)[2], 2), "*PC1",
                 "+", round(coef(model_hybrid_status)[3], 2), "*PC2")

# Plot the data with the equation and R-squared
ggplot(lab, aes(x = predicted, y = WL_max, color = hybrid_status)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#990000",
                size = 0.8) +
    labs(x = "Predicted", y = "Observed") +
    ggtitle("PC1, PC2, and Hybrid Status predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(lab$WL_max),
             label = eq_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    annotate("text", x = max(predicted), y = min(lab$WL_max) - 2,
             label = r_squared_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_brewer(palette = "Set1")  # Use a color palette from RColorBrewer


#######################################################hybrids
# Define the model
model_hybrid_only <- lm(WL_max ~ hybrid_status, data = lab)

# Generate equation text
eq_text <- paste("WL_max =", round(coef(model_hybrid_only)[1], 2),
                 "+", round(coef(model_hybrid_only)[2], 2), "*hybrid_status")

# Calculate R-squared
predicted <- predict(model_hybrid_only)
r_squared <- summary(model_hybrid_only)$r.squared
r_squared_text <- paste("R-squared =", round(r_squared, 2))

# Plot the data with the equation and R-squared
ggplot(lab, aes(x = predicted, y = WL_max, color = hybrid_status)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#990000",
                size = 0.8) +
    labs(x = "Predicted", y = "Observed") +
    ggtitle("Hybrid Status predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(lab$WL_max),
             label = eq_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    annotate("text", x = max(predicted), y = min(lab$WL_max) - 2,
             label = r_squared_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_brewer(palette = "Set1")  # Use a color palette from RColorBrewer

# Calculate residuals
residuals <- lab$WL_max - predicted

# Generate equation text for residuals
eq_text_residuals <- paste("Residuals =", round(coef(model_hybrid_only)[1], 2), " + ", 
                           round(coef(model_hybrid_only)[2], 2), "*hybrid_status")

# Plot the residuals
ggplot(lab, aes(x = predicted, y = residuals, color = hybrid_status)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_hline(yintercept=0, color = "#990000", 
               size = 0.8) +
    labs(x = "Predicted", y = "Residuals") +
    ggtitle("Residual plot of Hybrid Status predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(residuals),
             label = eq_text_residuals, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_brewer(palette = "Set1")  # Use a color palette from RColorBrewer


############### Parasite model + pc1 / pc2
# Define the model
model_pc1_pc2_parasite_challenge <- lm(WL_max ~ PC1 + PC2 + Parasite_challenge, data = lab)

# Generate equation text
eq_text <- paste("WL_max =", round(coef(model_pc1_pc2_parasite_challenge)[1], 2),
                 "+", round(coef(model_pc1_pc2_parasite_challenge)[2], 2), "*PC1",
                 "+", round(coef(model_pc1_pc2_parasite_challenge)[3], 2), "*PC2",
                 "+", round(coef(model_pc1_pc2_parasite_challenge)[4], 2), "*Parasite_challenge")

# Calculate R-squared
predicted <- predict(model_pc1_pc2_parasite_challenge)
r_squared <- summary(model_pc1_pc2_parasite_challenge)$r.squared
r_squared_text <- paste("R-squared =", round(r_squared, 2))


# Plot the data with the equation and R-squared
ggplot(lab, aes(x = predicted, y = WL_max, color = Parasite_challenge)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#990000",
                size = 0.8) +
    labs(x = "Predicted", y = "Observed") +
    ggtitle("PC1, PC2, and Parasite Challenge predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(lab$WL_max),
             label = eq_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    annotate("text", x = max(predicted), y = min(lab$WL_max) - 2,
             label = r_squared_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_manual(values = color_palette)  # Use the defined color palette


####################### only parasite
# Define the model
model_parasite_challenge <- lm(WL_max ~ Parasite_challenge, data = lab)

# Generate equation text
eq_text <- paste("WL_max =", round(coef(model_parasite_challenge)[1], 2),
                 "+", round(coef(model_parasite_challenge)[2], 2), "*Parasite_challenge")

# Calculate R-squared
predicted <- predict(model_parasite_challenge)
r_squared <- summary(model_parasite_challenge)$r.squared
r_squared_text <- paste("R-squared =", round(r_squared, 2))


# Plot the data with the equation and R-squared
parasite_challenge_lm <-
  ggplot(lab, aes(x = predicted, y = WL_max, color = Parasite_challenge)) +
    geom_point(size = 3, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#990000",
                size = 0.8) +
    labs(x = "Predicted", y = "Observed") +
    ggtitle("Parasite Challenge predicting weight loss") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right") +
    annotate("text", x = max(predicted), y = min(lab$WL_max),
             label = eq_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    annotate("text", x = max(predicted), y = min(lab$WL_max) - 2,
             label = r_squared_text, hjust = 1, vjust = -0.2,
             color = "black", size = 4, fontface = "bold") +
    scale_color_manual(values = color_palette)  # Use the defined color palette





################################Create the figure panel
parasite_challenge_lm
#figure_panel_1 <- 
  plot_grid(
    
    pca_individuals, pca_variables, linear_pc1_pc2_WL, immunization_parasite_chal_lm,
    labels = c("A", "B", "C", "D"),
    ncol = 2, align = "h", axis = "lr",,
    rel_widths = c(1, 1.2), rel_heights = c(1.2, 1))

ggsave("figure_panel_1.png", figure_panel_1, width = 10, height = 8, dpi = 300)


library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggridges)
library(tidyr)
library(rlang)
library(patchwork)


# read the data
hm <- read.csv("Data/Data_output/imputed_clean_data.csv")


# Select laboratory data 
# Select genes
lab <- hm %>%
    dplyr::filter(origin == "Lab")

Field <- hm %>%
    dplyr::filter(origin == "Field")

# create a vector to select genes
Genes_v   <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",
               "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
               "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
               "TICAM1", "TNF") #"IL.12", "IRG6")

# check the distributions of weight loss for strains


# Define colors
colors <- c("TRUE" = "firebrick3", "FALSE" = "steelblue")

# transform mouse strain into factor
lab$mouse_strain <- as.factor(lab$mouse_strain)

lab$mouse_strain <- gsub(pattern = "_", " ", lab$mouse_strain)

# order factor levels
lab$mouse_strain <- factor(lab$mouse_strain, 
                           levels = names(
                               sort(tapply(lab$WL_max, lab$mouse_strain, median))))

#Numbers of each mouse strain

lab %>%
    group_by(mouse_strain) %>%
    # Summarize the data to get counts for each mouse strain
    summarize(count = n()) %>%
    # Reorder mouse_strain by count
    mutate(mouse_strain = reorder(mouse_strain, count)) %>%
    # Plotting
    ggplot(aes(x = mouse_strain, y = count, fill = mouse_strain)) +
    geom_bar(stat = "identity") + 
    # Specify stat = "identity" for pre-summarized data
    geom_text(aes(label = count), vjust = -0.3) + 
    # Add count labels on top of bars
    scale_fill_viridis_d() + 
    # Use a nice color scale, like Viridis
    theme_minimal() + # Apply a minimal theme for a cleaner look
    theme(axis.text.x = element_text(angle = 50)) +
    labs(#title = "Mouse strains in experimental infections",
             x = "Mouse Strain", 
         y = "Number of mice") +# Add label
    guides(fill = "none") -> m_s

m_s

ggsave(filename = "~/GitHub/Article_predicting_WL/figures/mice_strains_n.jpeg",
       plot = m_s, width = 8, height = 6)




# Creating a density plot for the Hybrid Index (HI)
ggplot(Field, aes(x = HI)) + 
    geom_density(fill = "steelblue", alpha = 0.7) +
    geom_vline(aes(xintercept = mean(HI, na.rm = TRUE)),
               color = "red", linetype = "dashed", size = 1) +
    labs(title = "Distribution of Hybrid Index (HI) Among Wild Mice",
         x = "Hybrid Index (HI)",
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) -> h_w

h_w

ggsave(filename = "~/GitHub/Article_predicting_WL/figures/densityplot_HI.jpeg",
       plot = m_s, width = 8, height = 6)
# The red dashed line represents the mean HI value, providing a reference for 
#the central tendency of hybridization.



ggplot(lab %>%
           filter(infection == "challenge"), 
       aes(x = WL_max, y = mouse_strain, fill = mouse_strain)) + 
    geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(height = 0), 
                        scale = 0.9, alpha = 0.6, point_shape = 21, point_size = 2, point_alpha = 1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, position = position_nudge(x = 0.2)) +
    coord_flip() +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 70, vjust = 0.5, hjust=0.3)) +
    xlab("Maximum relative weight loss") +
    ylab("Mouse Strain") -> strains_weight

strains_weight

ggsave(filename = "figures/strains_weight.jpeg", plot = strains_weight, 
       width = 8, height = 6, dpi = 1000)





ggplot(lab, aes(x = WL_max, y = mouse_strain, fill = mouse_strain)) + 
    geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(height = 0), 
                        scale = 0.9, alpha = 0.6, point_shape = 21, point_size = 2, point_alpha = 1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, position = position_nudge(x = 0.2)) +
    coord_flip() +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 70, vjust = 0.5, hjust=0.3)) +
    xlab("Maximum relative weight loss") +
    ylab("Mouse Strain") +
    facet_grid(~infection) -> strains_weight_challenge

strains_weight_challenge

ggsave(filename = "figures/strains_weight_chalenge.jpeg", plot = strains_weight_challenge, 
       width = 10, height = 6, dpi = 1000)

###################### eimeria strains
# Then, define the color for each level of infection
color_mapping <- c("E_falciformis" = "salmon", 
                   "E_ferrisi" = "forestgreen", 
                   "uninfected" = "skyblue")

lab$current_infection <- gsub(pattern = "_", replacement = ". ", lab$current_infection)

ggplot(lab %>%
           filter(infection == "challenge"), 
       aes(x = WL_max, y = Parasite_challenge, fill = Parasite_challenge)) + 
    geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(height = 0), 
                        scale = 0.9, alpha = 0.6, point_shape = 21, point_size = 2, point_alpha = 1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, position = position_nudge(x = 0.2)) +
   # coord_flip() +
    theme_minimal() +
    scale_fill_manual(values = color_mapping) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.3)) +
    xlab("Maximum relative weight loss") +
    ylab("Parasite strain - challenge infections") -> eimeria_weight_chal

eimeria_weight_chal

ggsave(filename = "figures/eimeria_strains_weight_c.jpeg", plot = eimeria_weight_chal, 
       width = 8, height = 6, dpi = 1000)

ggplot(lab %>%
           filter(infection == "primary"), 
       aes(x = WL_max, y = Parasite_primary, fill = Parasite_primary)) + 
    geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(height = 0), 
                        scale = 0.9, alpha = 0.6, point_shape = 21, point_size = 2, point_alpha = 1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, position = position_nudge(x = 0.2)) +
    # coord_flip() +
    theme_minimal() +
    scale_fill_manual(values = c("E_falciformis" = "salmon", 
                                  "E_ferrisi" = "forestgreen")) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.3)) +
    xlab("Maximum relative weight loss") +
    ylab("Parasite strain - primary infections") -> eimeria_weight_prim

eimeria_weight_prim

ggsave(filename = "figures/eimeria_strains_weight_p.jpeg", plot = eimeria_weight_prim, 
       width = 8, height = 6, dpi = 1000)

# prim vs challenge
lab$current_infection <- gsub(pattern = "_", replacement = ". ", lab$current_infection)

ggplot(lab, 
       aes(x = WL_max, y = current_infection, fill = current_infection)) + 
    geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(height = 0), 
                        scale = 0.9, alpha = 0.6, point_shape = 21, point_size = 2, point_alpha = 1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, position = position_nudge(x = 0.2)) +
    # coord_flip() +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.3)) +
    xlab("Maximum relative weight loss") +
    ylab("Parasite strain") +
    scale_color_manual(values = color_mapping) +
    facet_wrap(~ infection, nrow = 2)-> eimeria_weight_challenge

eimeria_weight_challenge

ggsave(filename = "figures/eimeria_strains_weight_challenge.jpeg", plot = eimeria_weight_challenge, 
       width = 8, height = 10, dpi = 1000)


###### is the difference between parasite strains significant?
parasite_WL <- lm(formula = WL_max ~ current_infection, data = lab)
summary(parasite_WL)


##### is the differnece between mouse strains significant?`
mouse_WL <- lm(WL_max ~ mouse_strain, data = lab)
summary(mouse_WL)               


###################################################
### read uncleaend laboratory data
hm_full <- read.csv("Data/Data_output/full_data_prior_imputation.csv")

Challenge <- hm_full %>% 
    filter(origin == "Lab")

#relative weight loss per day - challenge
Challenge %>%
    filter(infection == "challenge") %>%
    ggplot(aes(x = dpi, y = relative_weight, color = Parasite_challenge, 
               fill = Parasite_challenge)) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.4, 
               shape = 21, stroke = 0.5, size = 3) + # Adjusted for outlines
    geom_smooth(aes(fill = Parasite_challenge), method = "loess", se = TRUE, alpha = 0.2) + # Add smooth line with confidence intervals
    scale_color_manual(values = color_mapping) + # Apply custom color mapping
    scale_fill_manual(values = color_mapping) + # Ensure fills match colors for confidence intervals
    labs(#title = "Relative Weight by Days Post Infection",
         x = "Days Post Infection (dpi)",
         y = "Relative weight, challenge infections",
         color = "Infection group",
         fill = "Infection group") + # Added for consistency with the legend
    theme_minimal() + # Use a minimal theme for a cleaner look
    theme(legend.position = "right", # Adjust legend position
          plot.title = element_text(hjust = 0.5), # Center the plot title
          legend.title.align = 0.5) -> Rwc

Rwc

ggsave(filename = "figures/relative_WL_challenge.jpeg", plot = Rwc, 
       width = 8, height = 10, dpi = 1000)

#relative weight loss per day - primary
Challenge %>%
    filter(infection == "primary") %>%
    ggplot(aes(x = dpi, y = relative_weight, color = Parasite_primary, 
               fill = Parasite_primary)) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.4, 
                shape = 21, stroke = 0.5, size = 3) + # Adjusted for outlines
    geom_smooth(aes(fill = Parasite_primary), method = "loess", se = TRUE, alpha = 0.2) + # Add smooth line with confidence intervals
    scale_color_manual(values = color_mapping) + # Apply custom color mapping
    scale_fill_manual(values = color_mapping) + # Ensure fills match colors for confidence intervals
    labs(#title = "Relative Weight by Days Post Infection",
         x = "Days Post Infection (dpi)",
         y = "Relative weight, primary infections",
         color = "Infection group",
         fill = "Infection group") + # Added for consistency with the legend
    theme_minimal() + # Use a minimal theme for a cleaner look
    theme(legend.position = "right", # Adjust legend position
          plot.title = element_text(hjust = 0.5), # Center the plot title
          legend.title.align = 0.5) -> Rwp

Rwp

ggsave(filename = "figures/relative_WL_primary.jpeg", plot = Rwp, 
       width = 8, height = 10, dpi = 1000)



# Adjusted ggplot2 code for OoC variable
Challenge %>%
    filter(infection == "challenge") %>%
    ggplot(aes(x = dpi, y = OO4sq, color = Parasite_challenge,
               fill = Parasite_challenge)) + # Adjust y-axis to OoC
    geom_jitter(width = 0.2, height = 0, alpha = 0.6, 
                shape = 21, stroke = 0.5, size = 3) + # Adjusted for outlines
    geom_smooth(aes(fill = Parasite_challenge), 
                method = "loess", se = TRUE, alpha = 0.2) + # Add smooth line with confidence intervals
    scale_color_manual(values = color_mapping) + # Apply custom color mapping
    scale_fill_manual(values = color_mapping) + # Ensure fills match colors for confidence intervals
    labs(#title = "Oocysts per Gram by Days Post Infection",
         x = "Days Post Infection (dpi)",
         y = "Oocysts per Gram, challenge infections",
         color = "Infection group",
         fill = "Infection group") + # Adjusted labels for clarity
    theme_minimal() + # Use a minimal theme for a cleaner look
    theme(legend.position = "right", # Adjust legend position
          plot.title = element_text(hjust = 0.5), # Center the plot title
          legend.title.align = 0.5) -> ooc_challenge

ooc_challenge

ggsave(filename = "figures/oocysts_challenge.jpeg", plot = ooc_challenge, 
       width = 8, height = 10, dpi = 1000)


# Adjusted ggplot2 code for OoC variable
Challenge %>%
    filter(infection == "primary", !dpi == "0") %>%
    ggplot(aes(x = dpi, y = OOC, color = Parasite_primary,
               fill = Parasite_primary)) + # Adjust y-axis to OoC
    geom_jitter(width = 0.2, height = 0, alpha = 0.6, 
                shape = 21, stroke = 0.5, size = 3) + # Adjusted for outlines
    geom_smooth(aes(fill = Parasite_primary), 
                method = "loess", se = TRUE, alpha = 0.2) + # Add smooth line with confidence intervals
    scale_color_manual(values = color_mapping) + # Apply custom color mapping
    scale_fill_manual(values = color_mapping) + # Ensure fills match colors for confidence intervals
    labs(#title = "Oocysts per Gram by Days Post Infection",
         x = "Days Post Infection (dpi)",
         y = "Oocysts per Gram, primary infections",
         color = "Infection group",
         fill = "Infection group") + # Adjusted labels for clarity
    theme_minimal() + # Use a minimal theme for a cleaner look
    theme(legend.position = "right", # Adjust legend position
          plot.title = element_text(hjust = 0.5), # Center the plot title
          legend.title.align = 0.5) -> ooc_primary

ooc_primary

ggsave(filename = "figures/oocysts_primary.jpeg", plot = ooc_primary, 
       width = 8, height = 10, dpi = 1000)

ooc_primary + ooc_challenge + Rwp + Rwc +
    m_s + strains_weight +
eimeria_weight_prim + eimeria_weight_chal + plot_layout(ncol = 2)



# Combine the plots
panel_figure <- #(ooc_primary | ooc_challenge) / # oocysts
                (Rwp | Rwc) /
                (strains_weight) / #| m_s  ) /
                (eimeria_weight_prim | eimeria_weight_chal) +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
panel_figure <- panel_figure + 
    plot_annotation(title = 'Fig. 1', 
                    theme = theme(plot.title = element_text(size = 20, hjust = 0)))

# Control sizes of each plot within the panel
# This is a generic example. You'll need to adjust the widths, heights, and layout design based on your specific needs.
panel_figure <- panel_figure + 
    plot_layout(heights = c(1, 1, 1), 
                widths = c(1, 1, 1)) # Adjust according to your layout needs

# Display the panel figure
print(panel_figure)

# Save the panel figure
ggsave('figure_panels/experimental_design_simple.jpeg', 
       panel_figure, width = 10, height = 12, dpi = 300)

    
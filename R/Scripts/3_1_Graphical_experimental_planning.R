library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggridges)
library(tidyr)
library(rlang)
library(patchwork)
library(gt)
library(stringr)
library(maps)
library(ggdist)


# read the data
hm <- read.csv("Data/Data_output/imputed_clean_data.csv")

### read uncleaend laboratory data
hm_full <- read.csv("Data/Data_output/full_data_prior_imputation.csv")

Challenge <- hm_full %>% 
    filter(origin == "Lab")

# Select laboratory data 
# Select genes
lab <- hm %>%
    dplyr::filter(origin == "Lab")

Field <- hm %>%
    dplyr::filter(origin == "Field")

# create a vector to select genes
Genes_v   <- c("IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10",
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
    dplyr::group_by(mouse_strain) %>%
    # Summarize the data to get counts for each mouse strain
    dplyr::summarize(count = n()) %>%
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
    labs(#title = "Distribution of Hybrid Index (HI) Among Wild Mice",
         x = "Hybrid Index (HI)",
         y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) -> h_w

h_w

ggsave(filename = "~/GitHub/Article_predicting_WL/figures/densityplot_HI.jpeg",
       plot = m_s, width = 8, height = 6)
# The red dashed line represents the mean HI value, providing a reference for 
#the central tendency of hybridization.


# Base world map
world_map <- map_data("world")

lon_range <- range(Field$Longitude) + c(-0.1, 0.1)  # Expanding the range a bit for padding
lat_range <- range(Field$Latitude) + c(-0.1, 0.1)   # Expanding the range a bit for padding

# Plot with zoom
map_hybrids <-
    ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "gray90", color = "white") +
    geom_point(data = Field, aes(x = Longitude, y = Latitude, color = HI), size = 3) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(#title = "Mouse HI Values Across Locations", 
         x = "Longitude", y = "Latitude") +
    theme_minimal() +
    coord_fixed(1.3) +  # Keep map aspect ratio
    xlim(lon_range) +   # Set x-axis limits to zoom in
    ylim(lat_range)     # Set y-axis limits to zoom in

map_hybrids 

ggsave(filename = "~/GitHub/Article_predicting_WL/figures/map_HI.jpeg",
       plot = map_hybrids, width = 8, height = 6)

###################################################
####################################################
###################################################
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



ggplot(Challenge %>%
           filter(infection == "challenge") %>%
           group_by(Mouse_ID), 
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

ggplot(Challenge %>%
           filter(infection == "primary") %>%
           group_by(Mouse_ID), 
       aes(x = WL_max, y = Parasite_primary, fill = Parasite_primary)) + 
    geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(height = 0), 
                        scale = 0.9, alpha = 0.6, point_shape = 21, point_size = 2, point_alpha = 1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, position = position_nudge(x = 0.2)) +
    # coord_flip() +
    theme_minimal() +
    scale_fill_manual(values = color_mapping) +
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

############################
#Hybrids


# Combine the plots
panel_hybr <- 
    (h_w | map_hybrids) +
    plot_layout(guides = 'collect') + # Collect all legends into a single legend
    plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)

# Add a figure title
panel_hybr <- panel_hybr + 
    plot_annotation(title = 'Fig. 2', 
                    theme = theme(plot.title = element_text(size = 13, hjust = 0)))

# Display the panel figure
print(panel_hybr)

# Save the panel figure
ggsave('figure_panels/hyb_plot_design.jpeg', 
       panel_hybr, width = 6, height = 3, dpi = 300)
#######################
############### make tables to present mouse strains and infection history
#eimeria strains 
Eim_strains <- 
    lab %>%
    group_by(infection_history) %>%
    dplyr::count() %>%
    mutate(infection_history = str_replace(infection_history, "falciformis", "E. falciformis"),
           infection_history = str_replace(infection_history, "ferrisi", "E. ferrisi"),
           infection_history = str_replace(infection_history, "uninfected", "Uninfected"),
           infection_history = str_replace(infection_history, "_", " X "))

Eim_strains                

as.data.frame(Eim_strains) %>%
    dplyr::rename("Infection history" = infection_history, 
                  "Number of mice" = n) %>%
    gt() %>%
    tab_header(title = md("Infection history of the mice used in experimental laboratory infections")) 


# prim 
#prim <- 
#count(lab$Parasite_primary)
x <- lab %>% filter(infection == "challenge")
count(x$Parasite_challenge)
model <- lm(WL_max ~ mouse_strain, data = lab)
summary(model)

model2 <- lm(WL_max ~ Parasite_primary, 
             data = Challenge %>% filter(infection == "primary"))
summary(model2)


#########################################################################
############################# genes
lab %>%
    pivot_longer(cols = all_of(Genes_v), 
                 names_to = "Genes", values_to = "Expression") %>%
    ggplot(aes(x = Parasite_primary, y = Expression)) +
    geom_density() +
    facet_wrap(~Genes)

color_mapping <- c("E. falciformis" = "salmon", 
                   "E. ferrisi" = "forestgreen", 
                   "uninfected" = "skyblue")

lab %>%
    pivot_longer(cols = all_of(Genes_v), 
                 names_to = "Genes", values_to = "Expression") %>%
    ggplot(aes(x = Expression, fill = current_infection)) +
    ggdist::stat_halfeye(
        adjust = .5, 
        width = .6, 
        alpha = 0.5,
        .width = 0, 
        justification = -.2, 
        point_colour = NA,
        orientation = "y"  # Set orientation to y
    ) +
    geom_boxplot(position = "dodge2",
        width = .5, 
        outlier.shape = NA,
        orientation = "y"  # Set orientation to y
    ) +
    facet_wrap(~Genes,  scales = 'free', ncol = 4) +
    labs(x = "Expression Level", y = "Density") +
    theme_minimal() +
    scale_fill_manual(values = color_mapping)  +
    theme_minimal() +
    labs(y = "Density", 
         x = "Immune gene expression level") -> density_imm

density_imm

ggsave(filename = "~/GitHub/Article_predicting_WL/figures/density_immune_genes.jpeg",
        plot = density_imm, width = 10, height = 8, dpi = 300)


# mean dpi at peak weight loss falciformis
s <- 
    Challenge %>%
    group_by(Mouse_ID, infection) %>% 
    filter(infection == "primary", Parasite_primary == "E_falciformis",
           relative_weight == min(relative_weight)) 
mean(s$dpi)

# mean dpi at peak weight loss falciformis - challenge
s <- 
    Challenge %>%
    group_by(Mouse_ID, infection) %>% 
    filter(infection == "challenge", Parasite_challenge == "E_falciformis",
           relative_weight == min(relative_weight)) 
mean(s$dpi)

# mean dpi at peak weight loss primary ferrisi
s <- 
    Challenge %>%
    group_by(Mouse_ID, infection) %>% 
    filter(infection == "primary", Parasite_primary == "E_ferrisi",
           relative_weight == min(relative_weight)) 
mean(s$dpi)

# mean dpi at peak weight loss E_ferrisi - challenge
s <- 
    Challenge %>%
    group_by(Mouse_ID, infection) %>% 
    filter(infection == "challenge", Parasite_challenge == "E_ferrisi",
           relative_weight == min(relative_weight)) 
mean(s$dpi)


# How many mice died in the primary infections
primary <- Challenge %>%
    filter(infection == "primary", Parasite_primary == "E_falciformis")
# max weight loss primary
mean(primary$WL_max)

# How many mice died in the primary infections
primary <- Challenge %>%
    filter(infection == "primary", Parasite_primary == "E_ferrisi")
# max weight loss primary
mean(primary$WL_max, na.rm = TRUE)

# How many mice died in the primary infections
primary <- Challenge %>%
    filter(infection == "primary", Parasite_primary == "uninfected")
# max weight loss primary
mean(primary$WL_max)

chale <- Challenge %>%
    filter(infection == "challenge", Parasite_challenge == "E_falciformis")
# max weight loss primary
mean(chale$WL_max)

# How many mice died in the primary infections
chale <- Challenge %>%
    filter(infection == "challenge", Parasite_challenge == "E_ferrisi")
# max weight loss primary
mean(chale$WL_max)

# How many mice died in the primary infections
chale <- Challenge %>%
    filter(infection == "challenge", Parasite_challenge == "uninfected")
# max weight loss primary
mean(chale$WL_max)



challenge <- Challenge %>%
    filter(infection == "challenge")

mean(Challenge$WL_max)

intersection <- intersect(primary$Mouse_ID, challenge$Mouse_ID)


s <- lab %>%
    filter(death == "primary")



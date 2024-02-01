library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggridges)


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

# check the distributions of weight loss for strains

A <- lab %>%
    dplyr::select(Mouse_ID, mouse_strain, relative_weight, WL_max)

# Define colors
colors <- c("TRUE" = "firebrick3", "FALSE" = "steelblue")

# transform mouse strain into factor
lab$mouse_strain <- as.factor(lab$mouse_strain)

lab$mouse_strain <- gsub(pattern = "_", " ", lab$mouse_strain)

# order factor levels
lab$mouse_strain <- factor(lab$mouse_strain, 
                           levels = names(
                               sort(tapply(lab$WL_max, lab$mouse_strain, median))))

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
                   "uninfected" = "blue")

lab$current_infection <- gsub(pattern = "_", replacement = ". ", lab$current_infection)

ggplot(lab %>%
           filter(infection == "challenge"), 
       aes(x = WL_max, y = current_infection, fill = current_infection)) + 
    geom_density_ridges(jittered_points = TRUE, position = position_points_jitter(height = 0), 
                        scale = 0.9, alpha = 0.6, point_shape = 21, point_size = 2, point_alpha = 1) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5, position = position_nudge(x = 0.2)) +
   # coord_flip() +
    theme_minimal() +
    scale_color_manual(values = color_mapping) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.3)) +
    xlab("Maximum relative weight loss") +
    ylab("Parasite strain") -> eimeria_weight

eimeria_weight

ggsave(filename = "figures/eimeria_strains_weight.jpeg", plot = eimeria_weight, 
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


######### parasite vs mouse
# Create the bar plot using ggplot2
ggplot(lab, aes(x = mouse_strain, y = current_infection, fill = current_infection)) +
    geom_segment(aes(xend = GO_Term, yend = 0),  size = 1.5) +
    geom_point(size = 3, shape = 19, color = "mediumvioletred", fill = "white") +
    coord_flip() +
    labs(x = "Enriched GO Terms", y = "-log10(p-value)",
         title = "Gene Ontology Enrichment Analysis") +
    theme_minimal()

Field <- hm %>%
    filter(origin == "Field")

Field_s <- hm %>%
    dplyr::select(c(Mouse_ID, Year, )

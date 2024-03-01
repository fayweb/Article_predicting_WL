
library(tidyr)
library(dplyr)
library(broom)
library(ggplot2)
library(patchwork)



hm <- read.csv("Data/Data_output/imputed_clean_data.csv")


# WOrking with laboratory data only
# Select genes
lab <- hm %>%
    dplyr::filter(origin == "Lab")

# List of dependent variables
dependent_vars <- c("IFNy", "CXCR3", "IL.6", "IL.13",# "IL.10",
                    "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                    "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                    "TICAM1", "TNF")

# Assuming 'lab' is your dataframe
lab$current_infection <- factor(lab$current_infection, levels = c("uninfected", "E_falciformis", "E_ferrisi"))

# Perform regressions
results <- lapply(dependent_vars, function(var) {
    lm_formula <- as.formula(paste(var, "~ current_infection", sep = " "))
    lm(lm_formula, data = lab)
})

names(results) <- dependent_vars  # Name each regression with the name of the dependent variable


# Correctly add 'model' column to each tidied model's dataframe
tidy_models <- do.call(rbind, lapply(names(results), function(name) {
    model_df <- tidy(results[[name]], conf.int = TRUE)
    model_df$model <- name  # Add model name as a new column
    return(model_df)
}))

# Now 'tidy_models' should have a 'model' column correctly populated

# Filter out intercept terms
tidy_models_no_intercept <- tidy_models[!grepl("intercept", tidy_models$term, ignore.case = TRUE),]

tidy_models_no_intercept <- as.data.frame(tidy_models_no_intercept) %>%
    mutate(term = case_when(
        term == "current_infectionE_falciformis" ~ "E. falciformis",
        term == "current_infectionE_ferrisi" ~ "E. ferrisi",
    ))

# Define custom colors for the independent variables
custom_colors <- c("E. falciformis" = "salmon", "E. ferrisi" = "forestgreen")

# Create the coefficient plot
ggplot(tidy_models_no_intercept, aes(x = model, y = estimate, color = term)) +
    coord_flip() +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = custom_colors) +
    theme_classic() +
    theme(axis.text.y = element_text(angle = 45, hjust = 1)) +
    labs(x = "Coefficient Estimate", y = "Gene expression level") +
#, title = "Coefficient Estimates by Gene Expression") +
    theme(legend.title = element_blank(),
          legend.position = "none") -> coef_mmr

print(coef_mmr)

ggsave(filename = "~/GitHub/Article_predicting_WL/figures/coef_plot_lab_genes.jpeg",
       plot = coef_mmr, width = 6, height = 4, dpi = 300)


####################################################
    ####################################################
############################# gene expression distribution
    

lab$current_infection <- gsub(pattern = "_", replacement = ". ", lab$current_infection)

    
color_mapping <- c("E. falciformis" = "salmon", 
                       "E. ferrisi" = "forestgreen", 
                       "uninfected" = "skyblue")
    
    lab %>%
        pivot_longer(cols = all_of(dependent_vars), 
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
        theme(legend.title = element_blank())+
        labs(y = "Density", 
             x = "Gene expression level") -> density_imm
    
    density_imm
    
    ggsave(filename = "~/GitHub/Article_predicting_WL/figures/density_immune_genes.jpeg",
           plot = density_imm, width = 10, height = 8, dpi = 300)
    
    
    #######################
    # combine
    comb <- (density_imm | coef_mmr) +
        plot_layout(guides = 'collect') + # Collect all legends into a single legend
        plot_annotation(tag_levels = 'A') # Add labels (A, B, C, etc.)
    
    # Add a figure title
    comb <- comb + 
        plot_annotation(title = 'Fig. 2', 
                        theme = theme(plot.title = element_text(size = 20, hjust = 0)))
    
    comb <- comb + 
        plot_layout(heights = c(1, 1), 
                    widths = c(2, 1))

    
    # Display the panel figure
    print(comb)
    
    # Save the panel figure
    ggsave('figure_panels/Fig_2_immune_gene_expression_lab.jpeg', 
           comb, width = 12, height = 6, dpi = 300)    
    
  
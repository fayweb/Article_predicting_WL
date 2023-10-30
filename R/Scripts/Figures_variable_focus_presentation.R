
################### IRGM1, SOCS1, MUC2
#focus on some variables for the presentations
# Variables you want to focus on
focus_vars <- c("IRGM1", "SOCS1", "MUC2")

vpg <- vpg %>%
    mutate(focus = ifelse(Variable %in% focus_vars, "Focus", "Fade"))

pca_variables <- ggplot(vpg, aes(x = PC1, y = PC2)) +
    geom_segment(aes(xend = 0, yend = 0), color = "gray50") +
    
    # Plot faded variables first
    geom_point(data = filter(vpg, focus == "Fade"), aes(color = cos2), size = 3, alpha = 0.05) +
    geom_label_repel(data = filter(vpg, focus == "Fade"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf, alpha = 0.3) +
    
    # Plot focus variables on top
    geom_point(data = filter(vpg, focus == "Focus"), aes(color = cos2), size = 3) +
    geom_label_repel(data = filter(vpg, focus == "Focus"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf) +
    
    coord_equal() +
    xlab("PC1 (32.83%)") +
    ylab("PC2 (16.25%)") +
    theme_minimal() + 
    guides(color = guide_colorbar(title = "Squared Distance from Origin")) +
    scale_color_gradientn(colors = gradient_colors, guide = "none")

print(pca_variables)

ggsave(filename = "figures/pca_variables_SOC1_IRGM1_MUC2.jpeg", plot = pca_variables, 
       width = 12, height = 6, dpi = 600)

################### MUC5AC, IL1RN, MPO
#focus on some variables for the presentations
# Variables you want to focus on
focus_vars <- c("MUC5AC", "IL1RN", "MPO")

vpg <- vpg %>%
    mutate(focus = ifelse(Variable %in% focus_vars, "Focus", "Fade"))

pca_variables <- ggplot(vpg, aes(x = PC1, y = PC2)) +
    geom_segment(aes(xend = 0, yend = 0), color = "gray50") +
    
    # Plot faded variables first
    geom_point(data = filter(vpg, focus == "Fade"), aes(color = cos2), size = 3, alpha = 0.05) +
    geom_label_repel(data = filter(vpg, focus == "Fade"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf, alpha = 0.3) +
    
    # Plot focus variables on top
    geom_point(data = filter(vpg, focus == "Focus"), aes(color = cos2), size = 3) +
    geom_label_repel(data = filter(vpg, focus == "Focus"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf) +
    
    coord_equal() +
    xlab("PC1 (32.83%)") +
    ylab("PC2 (16.25%)") +
    theme_minimal() + 
    guides(color = guide_colorbar(title = "Squared Distance from Origin")) +
    scale_color_gradientn(colors = gradient_colors, guide = "none")

print(pca_variables)

ggsave(filename = "figures/pca_variables_MUC5AC_IL1RN_MPO.jpeg", plot = pca_variables, 
       width = 12, height = 6, dpi = 600)

################### IL13
#focus on some variables for the presentations
# Variables you want to focus on
focus_vars <- "IL.13"

vpg <- vpg %>%
    mutate(focus = ifelse(Variable %in% focus_vars, "Focus", "Fade"))

pca_variables <- ggplot(vpg, aes(x = PC1, y = PC2)) +
    geom_segment(aes(xend = 0, yend = 0), color = "gray50") +
    
    # Plot faded variables first
    geom_point(data = filter(vpg, focus == "Fade"), aes(color = cos2), size = 3, alpha = 0.05) +
    geom_label_repel(data = filter(vpg, focus == "Fade"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf, alpha = 0.3) +
    
    # Plot focus variables on top
    geom_point(data = filter(vpg, focus == "Focus"), aes(color = cos2), size = 3) +
    geom_label_repel(data = filter(vpg, focus == "Focus"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf) +
    
    coord_equal() +
    xlab("PC1 (32.83%)") +
    ylab("PC2 (16.25%)") +
    theme_minimal() + 
    guides(color = guide_colorbar(title = "Squared Distance from Origin")) +
    scale_color_gradientn(colors = gradient_colors, guide = "none")

print(pca_variables)

ggsave(filename = "figures/pca_variables_IL.13.jpeg", plot = pca_variables, 
       width = 12, height = 6, dpi = 600)

################## TICAM1, NCR1, PRF1, CXCR3, RETNLB, IL.6, CXCL9, CASP1, MYD88, TNF
#focus on some variables for the presentations
# Variables you want to focus on
focus_vars <- c("TICAM1", "NCR1", "PRF1", "CXCR3", "RETNLB", "IL.6", "CXCL9", 
                "CASP1", "MYD88", "TNF")

vpg <- vpg %>%
    mutate(focus = ifelse(Variable %in% focus_vars, "Focus", "Fade"))

pca_variables <- ggplot(vpg, aes(x = PC1, y = PC2)) +
    geom_segment(aes(xend = 0, yend = 0), color = "gray50") +
    
    # Plot faded variables first
    geom_point(data = filter(vpg, focus == "Fade"), aes(color = cos2), size = 3, alpha = 0.05) +
    geom_label_repel(data = filter(vpg, focus == "Fade"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf, alpha = 0.3) +
    
    # Plot focus variables on top
    geom_point(data = filter(vpg, focus == "Focus"), aes(color = cos2), size = 3) +
    geom_label_repel(data = filter(vpg, focus == "Focus"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf) +
    
    coord_equal() +
    xlab("PC1 (32.83%)") +
    ylab("PC2 (16.25%)") +
    theme_minimal() + 
    guides(color = guide_colorbar(title = "Squared Distance from Origin")) +
    scale_color_gradientn(colors = gradient_colors, guide = "none")

print(pca_variables)

ggsave(filename = "figures/pca_variables_most_genes.jpeg", plot = pca_variables, 
       width = 12, height = 6, dpi = 600)



################################# pca variables focusing on pc1
################## TNF, IDO1, RETNLB, CXCL9, TICAM1, CXCR3, IFNy, MYD88, IL1RN
#focus on some variables for the presentations
# Variables you want to focus on
focus_vars <- c("TNF", "IDO1", "RETNLB", "CXCL9", "TICAM1", "CXCR3", "IFNy", 
                "MYD88", "IL1RN")

vpg <- vpg %>%
    mutate(focus = ifelse(Variable %in% focus_vars, "Focus", "Fade"))

pca_variables <- ggplot(vpg, aes(x = PC1, y = PC2)) +
    geom_segment(aes(xend = 0, yend = 0), color = "gray50") +
    
    # Plot faded variables first
    geom_point(data = filter(vpg, focus == "Fade"), aes(color = cos2), size = 3, alpha = 0.01) +
    geom_label_repel(data = filter(vpg, focus == "Fade"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf, alpha = 0.3) +
    
    # Plot focus variables on top
    geom_point(data = filter(vpg, focus == "Focus"), aes(color = cos2), size = 3) +
    geom_label_repel(data = filter(vpg, focus == "Focus"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf) +
    
    coord_equal() +
    xlab("PC1 (32.83%)") +
    ylab("PC2 (16.25%)") +
    theme_minimal() + 
    guides(color = guide_colorbar(title = "Squared Distance from Origin")) +
    scale_color_gradientn(colors = gradient_colors, guide = "none")

print(pca_variables)

ggsave(filename = "figures/pca_variables_pc1.jpeg", plot = pca_variables, 
       width = 12, height = 6, dpi = 600)



################################# pca variables focusing on pc2
################## IRGM1, SOCS1, MPO, MUC2, IL1RN
#focus on some variables for the presentations
# Variables you want to focus on
focus_vars <- c("IRGM1", "SOCS1", "MPO", "MUC2", "IL1RN")

vpg <- vpg %>%
    mutate(focus = ifelse(Variable %in% focus_vars, "Focus", "Fade"))

pca_variables <- ggplot(vpg, aes(x = PC1, y = PC2)) +
    geom_segment(aes(xend = 0, yend = 0), color = "gray50") +
    
    # Plot faded variables first
    geom_point(data = filter(vpg, focus == "Fade"), aes(color = cos2), size = 3, alpha = 0.01) +
    geom_label_repel(data = filter(vpg, focus == "Fade"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf, alpha = 0.3) +
    
    # Plot focus variables on top
    geom_point(data = filter(vpg, focus == "Focus"), aes(color = cos2), size = 3) +
    geom_label_repel(data = filter(vpg, focus == "Focus"), aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = Inf) +
    
    coord_equal() +
    xlab("PC1 (32.83%)") +
    ylab("PC2 (16.25%)") +
    theme_minimal() + 
    guides(color = guide_colorbar(title = "Squared Distance from Origin")) +
    scale_color_gradientn(colors = gradient_colors, guide = "none")

print(pca_variables)

ggsave(filename = "figures/pca_variables_pc2.jpeg", plot = pca_variables, 
       width = 12, height = 6, dpi = 600)



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
library(tidyverse)
library(randomForest)
library(caret)
library(ggpubr)
library(rfUtilities)

#import data
hm <- read.csv("Data/Data_output/imputed_clean_data.csv")

Gene_v <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", 
            "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", 
            "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

# prepare the lab data
lab <- hm %>% 
    dplyr::filter(origin == "Lab")

#select the imputed gene columns
gene_m <-  lab %>%
    drop_na(delta_ct_cewe_MminusE) %>%
    dplyr::select(c(Mouse_ID, all_of(Gene_v), delta_ct_cewe_MminusE))

genes <- gene_m %>%
    dplyr::select(-Mouse_ID)

gene_W <- lab  %>%
    drop_na(delta_ct_cewe_MminusE) %>%
    dplyr::select(c(all_of(Gene_v), delta_ct_cewe_MminusE))

# Splitting the data
set.seed(333)
training.samples <- createDataPartition(y = gene_W$delta_ct_cewe_MminusE, p = .7, list = FALSE)
train.data <- gene_W[training.samples, ] 
test.data <- gene_W[-training.samples, ] 

# Training the random forest model
set.seed(333)
delta_predict_gene <- randomForest(delta_ct_cewe_MminusE ~., data = train.data, 
                                   proximity = TRUE, ntree = 1000) 

saveRDS(delta_predict_gene, file = "R/Models/delta_predict_gene.RData")

# Cross-validation
predict_delta_cv <- rf.crossValidation(x = delta_predict_gene, xdata = train.data, 
                                       p = 0.10, n = 99, ntree = 501)

par(mfrow=c(2,2))
plot(predict_delta_cv) 


# Making predictions on test data
predictions <- predict(delta_predict_gene, test.data)

result <- test.data
result$predictions <- predictions

# Merging predictions with the original lab dataset
test_lab <- lab %>%
    left_join(result, by = c("delta_ct_cewe_MminusE", 
                             c("IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", 
                                                        "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", 
                                                        "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")))

# Correlation between predicted and actual data
cor_result <- cor(result$delta_ct_cewe_MminusE, result$predictions, 
                  method = c("pearson", "kendall", "spearman"))

cor_result

cor_test <- cor.test(result$delta_ct_cewe_MminusE, result$predictions)

cor_test

spearman_cor <- cor(result$delta_ct_cewe_MminusE, result$predictions, method = "spearman")
spearman_cor

# Checking for ΔCt values
test_lab <- test_lab %>%
    mutate(infected_delta = 
               case_when(
                   delta_ct_cewe_MminusE > -5 ~ "infected",
                   delta_ct_cewe_MminusE <= -5 ~ "uninfected"
               )) %>%
    drop_na(delta_ct_cewe_MminusE)

# Plotting
plotting <- test_lab %>%
    drop_na(delta_ct_cewe_MminusE) %>%
    ggplot(aes(x = predictions, y = delta_ct_cewe_MminusE, color = current_infection)) +
    geom_point(aes(size = delta_ct_cewe_MminusE), alpha = 0.7) +
    labs(
        x = "Predictions: ΔCt Value", 
        y = "Observed: ΔCt Value",
        title = "Relationship between Predicted and Observed ΔCt Value",
        color = "Current Infection",
        size = "Delta Ct value",
        shape = "Delta Ct threshold"
    ) +
    theme_minimal() +
    theme(
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
    ) +
    scale_color_manual(values = c(E_falciformis = "salmon", 
                                  E_ferrisi = "forestgreen", 
                                  uninfected = "deepskyblue")) +
    scale_size_continuous(range = c(2, 10))

plotting

ggsave(filename = "figures/predictions_random_for_lab_delta_ct.jpeg", 
       plot = plotting, width = 8, height = 5, dpi = 1000)


# Calculate linear model and visualize
lm_fit <- lm(delta_ct_cewe_MminusE ~ predictions, data = test_lab)

lm_fit

formula_text <- paste0("delta_ct_cewe_MminusE = ", 
                       round(coef(lm_fit)[1], 2), 
                       ifelse(coef(lm_fit)[2] >= 0, " + ", " - "), 
                       abs(round(coef(lm_fit)[2], 2)), 
                       " * predictions")
# Calculate Spearman's rank correlation (rho)
rho <- cor(test_lab$delta_ct_cewe_MminusE, test_lab$predictions, method = "spearman")
rho_text <- paste0("Rho (Spearman) = ", round(rho, 2))

# Linear regression formula annotation
lm_coef <- coef(lm_fit)
formula_text <- paste0("ΔCt = ", round(lm_coef[1], 2), 
                       ifelse(lm_coef[2] >= 0, " + ", " - "), 
                       abs(round(lm_coef[2], 2)), 
                       " * predictions")

test_lab %>%
    ggplot(aes(x = predictions, y = WL_max)) +
    geom_smooth(method = lm, se = TRUE) +
    labs(x = "Predictions: Maximum weight loss", 
         y = "Observed: Maximum weight loss") +
    geom_point(aes(x = predictions, y = WL_max, size = 0.8, alpha = 0.3)) +
    labs(x = "Predictions: Maximum weight loss", 
         y = "Observed: Maximum weight loss") +
    theme_light() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "none") +
    annotate("text", x = min(test_lab$predictions), y = max(test_lab$WL_max), 
             label = formula_text, hjust = 0, vjust = 4, size = 4, color = "blue") +
    annotate("text", x = min(test_lab$predictions), y = max(test_lab$WL_max), 
             label = cor_text, hjust = 0, vjust = 1.5, size = 4, color = "blue") -> linear_plot

linear_plot


linear_plot

ggsave(filename = "figures/linear_model_of_delta_rf_with_annotations.jpeg", 
       plot = linear_plot, width = 10, height = 6, dpi = 1000)


# Training model on all data (assuming this was your intention with the "toa" comment)
delta_predict_gene_all <- randomForest(delta_ct_cewe_MminusE ~., data = genes, 
                                       proximity = TRUE, ntree = 1000) 
saveRDS(delta_predict_gene_all, "R/Models/predict_delta_all.rds")

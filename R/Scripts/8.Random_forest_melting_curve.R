library(tidyverse)
library(randomForest)
library(caret)
library(ggpubr)
library(rfUtilities)
library(heatmaply)
library(viridis)

# Import data
hm <- read.csv("Data/Data_output/imputed_clean_data.csv")

Gene_v <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", 
            "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", 
            "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

# Prepare the lab data
lab <- hm %>% 
    dplyr::filter(origin == "Lab")

# Create a subset of the data using selected genes and the new target variable
data <- lab %>% select(c(Gene_v, "MC.Eimeria"))

# Convert MC.Eimeria to factor
data$MC.Eimeria <- as.factor(data$MC.Eimeria)

# Split data into training and testing sets
set.seed(289)
train_indices <- sample(1:nrow(data), nrow(data)*0.7)
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]

# Train a Random Forest model
rf_model <- randomForest(MC.Eimeria ~ ., data = train_data, ntree = 500)

# Make predictions on the test set
test_data$predictions <- predict(rf_model, test_data)

# Print summary of the Random Forest model
print(rf_model)

# Confusion Matrix and Model Evaluation
confusion <- confusionMatrix(test_data$predictions, test_data$MC.Eimeria)
print(confusion)

# Convert the confusion matrix to a tidy data frame for ggplot
confusion_melted <- as.data.frame(as.table(confusion$table))
colnames(confusion_melted) <- c("Reference", "Predicted", "n")

# Plot
ggplot(data = confusion_melted, aes(x = Predicted, y = reorder(Reference, desc(Reference)))) +
    geom_tile(aes(fill = n), colour = "white") +
    geom_text(aes(label = sprintf("%d", n)), vjust = 1) +
    scale_fill_gradient(low = "white", high = "Firebrick1") +
    labs(x = 'Predicted', y = 'Actual', fill = 'Number of observations', 
         title = 'Confusion Matrix for MC.Eimeria') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> confusion_plot_MC_Eimeria

confusion_plot_MC_Eimeria

ggsave(filename = "figures/confusion_matrix_predicting_MC_Eimeria.jpeg", 
       plot = confusion_plot_MC_Eimeria, width = 6, height = 4, dpi = 1000)

# Extract Feature Importance
feature_importance <- importance(rf_model)
feature_importance_ordered <- feature_importance[order(-feature_importance[,1]),]

# Train a Random Forest model - all data as training data 
rf_model <- randomForest(MC.Eimeria ~ ., data = data, ntree = 1000)
print(rf_model)

saveRDS(rf_model, "R/Models/predict_MC_Eimeria.rds")

library(tidyverse)
library(randomForest)
library(caret)
library(ggpubr)
library(rfUtilities)
library(heatmaply)


#import data
hm <- read.csv("Data/Data_output/imputed_clean_data.csv")

Gene_v <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", 
            "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", 
            "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

# prepare the lab data
lab <- hm %>% 
    dplyr::filter(origin == "Lab")

# Load data
# Assuming the dataset is stored in a variable named 'lab'
# lab <- read.csv("your_data_path.csv")

# Create a subset of the data using selected genes and the target variable
data <- lab %>% select(c(Gene_v, "current_infection"))

# Split data into training and testing sets
set.seed(123)
train_indices <- sample(1:nrow(data), nrow(data)*0.7)
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]

# Train a Random Forest model
rf_model <- randomForest(current_infection ~ ., data = train_data, ntree = 500)

# Make predictions on the test set
test_data$predictions <- predict(rf_model, test_data)

# Print summary of the Random Forest model
print(rf_model)

# Confusion Matrix and Model Evaluation
confusion <- confusionMatrix(test_data$predictions, test_data$current_infection)
print(confusion)

# Plot confusion matrix

heatmaply::heatmaply(
    as.matrix(confusion$table),
    xlab = "Predicted",
    ylab = "Reference",
    main = "Confusion Matrix",
    notecol = "black"  # change font color of cells
)

# Extract Feature Importance
feature_importance <- importance(rf_model)
feature_importance_ordered <- feature_importance[order(-feature_importance[,1]),]

# Plot Feature Importance
ggplot(data = as.data.frame(feature_importance_ordered), aes(x = reorder(Row.names(feature_importance_ordered), MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Feature Importance", x = "Features", y = "Mean Decrease in Gini") +
    theme_minimal()


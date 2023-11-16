# Load libraries
library(tidyverse)
library(tidyr)
library(dplyr)
library(cowplot)
library(randomForest)
library(varImp)
library(ggplot2)
library(gridGraphics)
library(ggpmisc)
library(caret)
library(ggpubr)
library(ggiraphExtra)
library(ggeffects)
library(rfUtilities) # Implements a permutation test cross-validation for 
# Random Forests models

#import data
hm <- read.csv("Data/Data_output/imputed_clean_data.csv")

#vectors for facs selection
facs_v <- c( "Treg", "CD4", "Treg17", "Th1", "Th17", "CD8",
                "Act_CD8", "IFNy_CD4", "IFNy_CD8")



# prepare the lab data
lab <- hm %>% 
    dplyr::filter(origin == "Lab")


#select the imputed facs columns
facs_m <-  lab %>%
    dplyr::select(c(Mouse_ID, all_of(facs_v), WL_max))

# drop na
facs_m <- facs_m %>%
    drop_na()

# select only the facss
facss <- facs_m %>%
    dplyr::select(-Mouse_ID)

# select the facss and the weight loss
facs_W <- lab  %>%
    dplyr::select(c(all_of(facs_v), WL_max)) %>%
    drop_na()

repeat_cv <- trainControl(method = "repeatedcv", #repeated cross validation
                          number = 5, # 5 fold cross validation
                          repeats = 3)



# split data into training and test
set.seed(333) # this will help us reproduce this random assignment

# in this way we can pick the random numbers
training.samples <- createDataPartition(y = facs_W$WL_max, p = .7, list = FALSE) 

# this is the partiicition! In this case 0.7 = training data and 0.3 = testing
# we don't want to get a list in return
train.data <- facs_W[training.samples, ] 
test.data <- facs_W[-training.samples, ] 


## ----predicting_weight_loss_model---
set.seed(333)


#train the model
WL_predict_facs <- randomForest(WL_max ~., data = train.data, 
                                proximity = TRUE, ntree = 1000) 

# ntree = number of trees     
# save the model 
save(WL_predict_facs, file =  "R/Models/WL_predict_facs.RData")
print(WL_predict_facs)

predict_WL_cv <- rf.crossValidation(x = WL_predict_facs, xdata = train.data, 
                                    p = 0.10, n = 99, ntree = 501)

predict_WL_cv$fit.var.exp

par(mfrow=c(2,2))


##################
##################
########## Plots

root_mean <- plot(predict_WL_cv)

# Root Mean Squared Error (observed vs. predicted) from each Bootstrap 
# iteration (cross-validation)
mean_error <- plot(predict_WL_cv, stat = "mse")

#Percent variance explained from specified fit model
model_var <- plot(predict_WL_cv, stat = "var.exp")

#Mean Absolute Error from each Bootstrapped model
abs_error <- plot(predict_WL_cv, stat = "mae")


#d# ---------------------------------------------------------------------------------------------------
error_random  <- plot(WL_predict_facs)

## ---------------------------------------------------------------------------------------------------
# number of trees with lowest MSE
which.min(WL_predict_facs$mse)

# RMSE of this optimal random forest
sqrt(WL_predict_facs$mse[which.min(WL_predict_facs$mse)])

WL_predict_facs$mtry
oob_error_rate <- WL_predict_facs$mse[WL_predict_facs$ntree]
oob_error_rate <- 1 - sum(diag(WL_predict_facs$confusion)) / sum(WL_predict_facs$confusion)


### Visualize variable importance ---
#Call importance() function on the model model to check how the attributes used 
# as predictors affect our WL_predict_facs
ImpData <- as.data.frame(importance(WL_predict_facs))
ImpData$Var.Names <- row.names(ImpData)
varImp(WL_predict_facs)

#WL_predict_facs$mse

## S3 method for class 'randomForest'
plot(WL_predict_facs, type = "l", main=deparse(substitute(x)))

variable_importance <- varImpPlot(WL_predict_facs)


ggsave(filename = "figures/variable_imporance_random.jpeg", width = 10, height = 8, dpi = 300)


# Get variable importance from the WL_predict_facs fit
ImpData <- as.data.frame(importance(WL_predict_facs))

#The predict() function in R is used to predict the values based on the 
# input data.
predictions <- predict(WL_predict_facs, test.data)

# assign test.data to a new object, so that we can make changes
result <- test.data

#add the new variable of predictions to the result object
result <- cbind(result, predictions)

# what is the correlation between predicted and actual data?
cor(result$WL_max, result$predictions, 
    method = c("pearson", "kendall", "spearman"))

cor.test(result$WL_max, result$predictions)


test_lab <- lab %>%
    left_join(result, by = intersect(colnames(lab), colnames(result)))

test_lab <- test_lab %>%
    drop_na(predictions)


# what is the correlation between predicted and actual data?
cor(result$WL_max, result$predictions, 
    method = c("pearson", "kendall", "spearman"))

model <- lm(predictions ~ WL_max, data = test_lab)
summary(model)

ggpredict(model, terms = c("WL_max")) %>% 
    plot(colors = "blue") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Observed maximum weight loss during infections") +
    ylab("Predicted maximum weight loss during infections") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))-> lm_short

lm_short



model <- lm(predictions ~ WL_max * current_infection, data = test_lab)    
summary(model)

#### Plotting
# Then, define the color for each level of infection
color_mapping <- c("E_falciformis" = "salmon", 
                   "E_ferrisi" = "forestgreen", 
                   "uninfected" = "cornflowerblue")

ggpredict(model, terms = c("WL_max", "current_infection")) %>% 
    plot(colors = "darkorchid") +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Observed maximum weight loss during infections") +
    ylab("Predicted maximum weight loss during infections") +
    theme_minimal() +
    scale_color_manual(values = color_mapping) +
    scale_fill_manual(values = color_mapping) +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)
    ) -> lm_weight_loss_predictions

lm_weight_loss_predictions


model <- lm(predictions ~ WL_max * delta_ct_cewe_MminusE , data = test_lab)
summary(model)

ggpredict(model, terms = c("WL_max", "delta_ct_cewe_MminusE"), interactive=TRUE) %>% 
    plot() +
    labs(title = NULL) +  # This removes the title
    # ggtitle("Effect of PC2 on Predicted Weight Loss") +
    xlab("Observed maximum weight loss during infections") +
    ylab("Predicted maximum weight loss during infections") +
    theme_minimal() +
    scale_color_manual(values = c("darkred", "gold", "violet")) +
    scale_fill_manual(values = c("darkred", "gold", "violet")) +
    theme(
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))-> lm_short

lm_short

### plotting
test_lab %>%
    ggplot(aes(x = predictions, y = WL_max, color = current_infection)) +
    geom_point(aes(size = delta_ct_cewe_MminusE), alpha = 0.7) +
    labs(
        x = "Predictions: Maximum weight loss", 
        y = "Observed: Maximum weight loss",
        #  title = "Relationship between Predicted and Observed Weight Loss",
        #subtitle = "Grouped by Current Infection and Sized by Delta CT Value",
        color = "Treatment group",
        size = "Caecal infection intensities,
      Delta Ct value",
        shape = "Delta Ct treshold"
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
    scale_size_continuous(range = c(2, 10)) -> predictions_random_for_lab

predictions_random_for_lab


ggsave(plot = predictions_random_for_lab, 
       filename = "figures/predictions_random_for_lab.jpeg", width = 8, height = 5,
       dpi = 1000)


# Calculate the linear model
lm_fit <- lm(WL_max ~ predictions, data = test_lab)

# Extract coefficients for the model formula
intercept <- round(coef(lm_fit)[1], 2)
slope <- round(coef(lm_fit)[2], 2)
formula_text <- paste0("WL_max = ", intercept, " ", ifelse(slope >= 0, "+ ", "- "), abs(slope), " * predictions")

# Calculate correlation
cor_value <- round(cor(test_lab$WL_max, test_lab$predictions), 2)
cor_text <- paste0("Rho = ", cor_value)

test_lab   %>%
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

ggsave(filename = "figures/linear_model_of_random_forest.jpeg", plot = linear_plot, 
       width = 10, height = 6,
       dpi = 1000)




#train the model
WL_predict_facs <- randomForest(WL_max ~., data = facs_W, 
                                proximity = TRUE, ntree = 1000) 

WL_predict_facs
# ntree = number of trees     
# save the model 
# toa = trained on all
saveRDS(WL_predict_facs, "R/Models/predict_WL.rds")

print(WL_predict_facs)





figure_panel_2 <- ggarrange(predictions_random_for_lab,
                            ggarrange(lm_weight_loss_predictions, lm_short, 
                                      labels = c("B", "C"), ncol = 2),
                            nrow = 2, labels = "A")



# Adding the title "Figure 1" to the entire arrangement
figure_panel_2 <- annotate_figure(figure_panel_2, 
                                  top = text_grob("Figure 2", size = 14, 
                                                  face = "bold"))

print(figure_panel_2)


ggsave("figure_panels/figure_panel_2.jpeg", figure_panel_2, 
       width = 18, height = 18, dpi = 300)


act_cd8 <- lm(WL_max ~ Act_CD8, data = lab)
summary(act_cd8)

cor(facs_W$Act_CD8, facs_W$WL_max)


x <- hm %>%
    filter(origin == "Field") %>%
    drop_na(Act_CD8)
glimpse()
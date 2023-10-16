#install.packages("optimx", version = "2021-10.12") # this package is required for 
#the parasite load package to work
library(tidyverse)
library(tidyr)
library(dplyr)
library(cowplot)
library(randomForest)
library(ggplot2)
library(caret)
library(ggpubr)
library(rfUtilities) # Implements a permutation test cross-validation for 
# Random Forests models

#import data
hm <- read.csv("Data/Data_output/imputed_clean_data.csv")


Gene_v   <- c("IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF") #"IL.12", "IRG6")


# prepare the lab data
lab <- hm %>% 
  dplyr::filter(origin == "Lab")


#select the imputed gene columns
gene_m <-  lab %>%
  dplyr::select(c(Mouse_ID, all_of(Gene_lab), WL_max))

genes <- gene %>%
  dplyr::select(-Mouse_ID)

gene_W <- lab  %>%
    dplyr::select(c(all_of(Gene_lab), WL_max))

repeat_cv <- trainControl(method = "repeatedcv", #repeated cross validation
                           number = 5, # 5 fold cross validation
                           repeats = 3)

# split data into training and test
set.seed(333) # this will help us reproduce this random assignment

# in this way we can pick the random numbers
training.samples <- createDataPartition(y = gene_W$WL_max, p = .7, list = FALSE) 

# this is the partiicition! In this case 0.7 = training data and 0.3 = testing
# we don't want to get a list in return
train.data <- gene_W[training.samples, ] 
test.data <- gene_W[-training.samples, ] 


## ----predicting_weight_loss_model---
set.seed(333)


#train the model
WL_predict_gene <- randomForest(WL_max ~., data = train.data, 
                                    proximity = TRUE, ntree = 1000) 
# ntree = number of trees     
# save the model 
save(WL_predict_gene, file =  "R/Models/WL_predict_gene.RData")
print(WL_predict_gene)

predict_WL_cv <- rf.crossValidation(x = WL_predict_gene, xdata = train.data, 
                                    p = 0.10, n = 99, ntree = 501)

predict_WL_cv$fit.var.exp

par(mfrow=c(2,2))

plot(predict_WL_cv) 

# Root Mean Squared Error (observed vs. predicted) from each Bootstrap 
# iteration (cross-validation)
plot(predict_WL_cv, stat = "mse")

#Percent variance explained from specified fit model
plot(predict_WL_cv, stat = "var.exp")

#Mean Absolute Error from each Bootstrapped model
plot(predict_WL_cv, stat = "mae")


## ---------------------------------------------------------------------------------------------------
plot(WL_predict_gene)


## ---------------------------------------------------------------------------------------------------
# number of trees with lowest MSE
which.min(WL_predict_gene$mse)

# RMSE of this optimal random forest
sqrt(WL_predict_gene$mse[which.min(WL_predict_gene$mse)])


### Visualize variable importance ---
#Call importance() function on the model model to check how the attributes used 
# as predictors affect our WL_predict_gene
importance(WL_predict_gene)

#WL_predict_gene$mse

## S3 method for class 'randomForest'
plot(WL_predict_gene, type = "l", main=deparse(substitute(x)))
varImpPlot(WL_predict_gene)

# Get variable importance from the WL_predict_gene fit
ImpData <- as.data.frame(importance(WL_predict_gene))
ImpData$Var.Names <- row.names(ImpData)

#The predict() function in R is used to predict the values based on the 
# input data.
predictions <- predict(WL_predict_gene, test.data)

# assign test.data to a new object, so that we can make changes
result <- test.data

#add the new variable of predictions to the result object
result <- cbind(result, predictions)

# what is the correlation between predicted and actual data?
cor(result$WL_max, result$predictions, 
    method = c("pearson", "kendall", "spearman"))

cor.test(result$WL_max, result$predictions)

cor(result$WL_max, result$predictions, 
    method = "spearman")

test_lab <- lab %>%
  left_join(result, by = c("WL_max", "IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF"))

test_lab <- test_lab %>%
  drop_na(predictions)


# what is the correlation between predicted and actual data?
cor(result$WL_max, result$predictions, 
    method = c("pearson", "kendall", "spearman"))

### delta ct We considered ΔCt  = −5 our limit of detection
test_lab <- test_lab %>%
    mutate(infected_delta = 
               case_when(
                   delta_ct_cewe_MminusE > -5 ~ "infected",
                   delta_ct_cewe_MminusE < -5 ~ "uninfected"
               ))


### plotting
test_lab %>%
    drop_na(delta_ct_cewe_MminusE) %>%
    mutate(infected_delta = 
               case_when(
                   delta_ct_cewe_MminusE > -5 ~ "infected",
                   delta_ct_cewe_MminusE < -5 ~ "uninfected"
               )) %>%
    ggplot(aes(x = predictions, y = WL_max, color = current_infection)) +
    # Geom
    geom_point(aes(size = delta_ct_cewe_MminusE, shape = infected_delta), alpha = 0.7) +
    
    # Labels
    labs(
        x = "Predictions: Maximum weight loss", 
        y = "Observed: Maximum weight loss",
        title = "Relationship between Predicted and Observed Weight Loss",
        #subtitle = "Grouped by Current Infection and Sized by Delta CT Value",
        color = "Current Infection",
        size = "Delta Ct value",
        shape = "Delta Ct treshold"
    ) +
    
    # Theme adjustments
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
    
    # Color adjustments as per given values
    scale_color_manual(values = c(E_falciformis = "salmon", 
                                  E_ferrisi = "forestgreen", 
                                  uninfected = "deepskyblue")) +
    
    # Size adjustments
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
WL_predict_gene <- randomForest(WL_max ~., data = gene, 
                                    proximity = TRUE, ntree = 1000) 
# ntree = number of trees     
# save the model 
# toa = trained on all
saveRDS(WL_predict_gene, "R/Models//predict_WL.rds")

print(WL_predict_gene)



####################################################################################
####################################################################################
####################################################################################
####################################################################################
####################################################################################
############################ infection status model

gene_curr <- lab %>%
  dplyr::select(c(Mouse_ID, all_of(Gene_v), current_infection))


gene <- gene_curr %>%
  dplyr::select(-Mouse_ID)

lab$current_infection <- as.factor(lab$current_infection)

# split data into training and test
set.seed(123) # this will help us reproduce this random assignment
# in this way we can pick the random numbers
training.samples <- createDataPartition(y = gene$current_infection, p = .7, list = FALSE)

train.data <- gene[training.samples, ] 
test.data <- gene[-training.samples, ] 


#train the model
model_Parasite <- randomForest(current_infection ~., 
                               data = train.data, proximity = TRUE,
                               ntree = 1500) # number of trees

# save the model 
save(model_Parasite, file =  "r_scripts/models/predict_infecting_parasite.rds")


#The predict() function in R is used to predict the values based on the input 
# data.
predictions_parasite <- predict(model_Parasite, test.data_parasite)
# assign test.data to a new object, so that we can make changes
result_parasite <- test.data_parasite
#add the new variable of predictions to the result object
result_parasite <- cbind(result_parasite, predictions_parasite)
#add the results to a data frame containing test data and the prediction
result_parasite <- cbind(lab2[row.names(result_parasite), ], predictions_parasite)


## ---------------------------------------------------------------------------------------------------

conf_matrix_parasite <- 
  confusionMatrix(
    result_parasite$predictions_parasite,
    reference = result_parasite$infection)

print(conf_matrix_parasite)

conf_matrix_parasite$table

plt <- as.data.frame(conf_matrix_parasite$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

ggplot(plt, aes(x = Prediction, y =  Reference, fill= Freq)) +
        geom_tile() + geom_text(aes(label=Freq)) +
        scale_fill_gradient(low="white", high="darkturquoise") +
        labs(x = "Predictions",y = "Reference") 



## ---------------------------------------------------------------------------------------------------
lab2 <- lab %>% 
  dplyr::filter(infection %in% c("E_ferrisi", "uninfected"))


lab2$infection <- as.factor(lab2$infection)

gene_curr <- lab2 %>%
  dplyr::select(c(Mouse_ID, all_of(Gene_lab), infection))



gene <- gene_curr %>%
  dplyr::select(-Mouse_ID)


# split data into training and test
set.seed(123) # this will help us reproduce this random assignment
# in this way we can pick the random numbers
training.samples <- gene$infection%>%
  createDataPartition(p = .7, list = FALSE) 
train.data_parasite <- gene[training.samples, ] 
test.data_parasite <- gene[-training.samples, ] 




## ---------------------------------------------------------------------------------------------------

#train the model
model_Parasite <- randomForest(infection ~., 
                               data = train.data_parasite, proximity = TRUE,
                      ntree = 1500) # number of trees

# save the model 
save(model_Parasite, file =  "r_scripts/models/predict_infecting_parasite.rds")

print(model_Parasite)



## ---------------------------------------------------------------------------------------------------
model_Parasite_cv <- rf.crossValidation(x = model_Parasite, xdata =  
                                          train.data_parasite, 
                                    p = 0.10, n = 99, ntree = 501)

model_Parasite_cv$fit.var.exp


# Plot cross validation versus model producers accuracy




## ---------------------------------------------------------------------------------------------------
#The predict() function in R is used to predict the values based on the input 
# data.
predictions_parasite <- predict(model_Parasite, test.data_parasite)
# assign test.data to a new object, so that we can make changes
result_parasite <- test.data_parasite
#add the new variable of predictions to the result object
result_parasite <- cbind(result_parasite, predictions_parasite)
#add the results to a data frame containing test data and the prediction
result_parasite <- cbind(lab2[row.names(result_parasite), ], predictions_parasite)


## ---------------------------------------------------------------------------------------------------

conf_matrix_parasite <- 
  confusionMatrix(
    result_parasite$predictions_parasite,
    reference = result_parasite$infection)

print(conf_matrix_parasite)

conf_matrix_parasite$table

plt <- as.data.frame(conf_matrix_parasite$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

ggplot(plt, aes(x = Prediction, y =  Reference, fill= Freq)) +
        geom_tile() + geom_text(aes(label=Freq)) +
        scale_fill_gradient(low="white", high="darkturquoise") +
        labs(x = "Predictions",y = "Reference") 



## ---------------------------------------------------------------------------------------------------
#train the model
model_Parasite <- randomForest(infection ~., data = gene, 
                                    proximity = TRUE, ntree = 1000) 
# ntree = number of trees     
# save the model 
# toa = trained on all
saveRDS(model_Parasite, "r_scripts/models/predict_infecting_parasite.rds")

print(model_Parasite)



## ---------------------------------------------------------------------------------------------------

str(lab$MC.Eimeria)

MC <- lab %>%
  dplyr::select(c(Mouse_ID, all_of(Gene_lab), MC.Eimeria))



MC <- MC %>%
  dplyr::select(-Mouse_ID)


# split data into training and test
set.seed(182) # this will help us reproduce this random assignment
# in this way we can pick the random numbers
training.samples <- MC$MC.Eimeria %>%
  createDataPartition(p = .7, list = FALSE) 

train.data_MC <- MC[training.samples, ] 
test.data_MC <- MC[-training.samples, ] 




## ---------------------------------------------------------------------------------------------------

#train the model
model_MC <- randomForest(MC.Eimeria ~., 
                               data = train.data_MC, proximity = TRUE,
                      ntree = 1500) # number of trees

# save the model 
saveRDS(model_MC, file =  "r_scripts/models/predict_MC.rds")

print(model_MC)



## ---------------------------------------------------------------------------------------------------
model_MC_cv <- rf.crossValidation(x = model_MC, xdata =  
                                          train.data_MC, 
                                    p = 0.10, n = 99, ntree = 501)

model_MC_cv$fit.var.exp


# Plot cross validation versus model producers accuracy




## ---------------------------------------------------------------------------------------------------

#The predict() function in R is used to predict the values based on the input 
# data.
predictions_MC <- predict(model_MC, test.data_MC)

# assign test.data to a new object, so that we can make changes
result_MC <- test.data_MC
#add the new variable of predictions to the result object
result_MC <- cbind(result_MC, predictions_MC)
#add the results to a data frame containing test data and the prediction
result_MC <- cbind(lab[row.names(result_MC), ], predictions_MC)


## ---------------------------------------------------------------------------------------------------

conf_matrix_MC <- 
  confusionMatrix(result_MC$predictions_MC, reference = result_MC$MC.Eimeria)

print(conf_matrix_MC)

conf_matrix_MC$table

plt <- as.data.frame(conf_matrix_MC$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

ggplot(plt, aes(x = Prediction, y =  Reference, fill= Freq)) +
        geom_tile() + geom_text(aes(label=Freq)) +
        scale_fill_gradient(low="white", high="darkturquoise") +
        labs(x = "Predictions",y = "Reference") 



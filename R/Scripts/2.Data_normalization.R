source("R/Scripts/1.Data_prep.R")


## ---------------------------------------------------------------------------------------------------
library(mice)
library(stringr)
library(gridExtra)


## ---------------------------------------------------------------------------------------------------
hm <- data


## Normalizing 
hm$Mouse_ID <- str_replace(hm$Mouse_ID, "_", "")

field <- hm %>%
    dplyr::filter(origin == "Field") 

# select the genes, the mouse identifier and the house keeping gene
gmf <- field[, c("Mouse_ID", Genes_v, "GAPDH")] 

# Remove columns with only NA values
gmf <- gmf %>% select_if(~!all(is.na(.)))

#remove rows with only nas 
gmf <- gmf[!apply(is.na(gmf[-1]), 1, all), ]

gf <- gmf[,-1]

##select same rows in the first table
field <- field %>% 
    filter(Mouse_ID %in% gmf$Mouse_ID)


###############lab
#select the genes and lab mice
lab <- hm %>%
    dplyr::filter(origin == "Lab", Position == "mLN") %>%
    group_by(Mouse_ID, infection) %>%
    filter(dpi == max_dpi)
    

# select the genes, the mouse identifier and the house keeping gene
gml <- lab[, c("Mouse_ID", Genes_v, "PPIB")] 

gml <- unique(gml)

# Remove columns with only NA values
gml <- gml %>% select_if(~!all(is.na(.)))

#remove rows with only nas 
gml <- gml[!apply(is.na(gml[-1]), 1, all), ]

gl <- gml[,-1]

##select same rows in the first table
lab <- lab %>%
    filter(Mouse_ID %in% gml$Mouse_ID)


# looking at patterns of nas)
#pattern_na <-as.data.frame(md.pattern(field_genes))
#field 
sapply(gmf, function(x) sum(is.na(x)))

#lab
sapply(gml, function(x) sum(is.na(x)))

#remove duplicates
lab_prim <- lab %>%
    filter(death == "primary")
lab_chal <- lab %>% 
    filter(death == "challenge", infection == "challenge")

lab <- rbind(lab_prim, lab_chal)
rm(lab_prim, lab_chal)

##############################################################################

# Livak's method of normalisation
# Create a function to calculate the delta delta ct for each gene for each mouse
# 1. Step 1: ΔCt (sample) = Ct(gene of interest) − Ct(reference gene)
# 2. ΔΔct = Δct from step 1 - mean(Δct for the gene)
# 3. 2-ΔΔCt
df <- field

#check the distributions of genes
plot_list <- lapply(names(gf), function(colname){
  ggplot(df, aes_string(x = colname)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(title = paste("Histogram of", colname))
})


#grid.arrange(grobs = plot_list, ncol = 3) 
#placed in hashtags, takes to long to process
calculate_expression <- function(df) {
  # Extract GAPDH column
  GAPDH <- df$GAPDH
  
  # Step 1: Calculate delta ct for all genes except GAPDH
  delta_ct <- df[, sapply(df, is.numeric)]  # Ensure only numeric columns are taken
  delta_ct <- delta_ct[, !colnames(delta_ct) %in% "GAPDH"]  # Exclude GAPDH
  delta_ct <- sapply(delta_ct, function(gene) gene - GAPDH)
  
  # Step 2: Calculate ΔΔct
  mean_dct <- colMeans(delta_ct, na.rm = TRUE)  # Calculate mean delta_ct for each gene
  delta_delta_ct <- t(apply(delta_ct, 1, function(row) row - mean_dct))
  
  # Step 3: Calculate the result
  result <- 2^(-delta_delta_ct)
  result <- round(result, 2)
  return(result)
}

# Use the function
result_df <- calculate_expression(df)

# View the result
print(result_df)

#check the distributions of genes
plot_list <- lapply(names(result_field[,-1]), function(colname){
  ggplot(result_field[,-1], aes_string(x = colname)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(title = paste("Histogram of", colname))
})


#grid.arrange(grobs = plot_list, ncol = 3)

############### same for lab ##################################################
# In the lab we use PPIB
#hist(lab$PPIB)

#check the distributions of genes
plot_list <- lapply(names(gl), function(colname){
  ggplot(gl, aes_string(x = colname)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(title = paste("Histogram of", colname))
})


#grid.arrange(grobs = plot_list, ncol = 3)


calculate_expression <- function(df) {
  # Extract PPIB column
  PPIB <- df$PPIB
  
  # Step 1: Calculate delta ct for all genes except PPIB
  delta_ct <- df[, sapply(df, is.numeric)]  # Ensure only numeric columns are taken
  delta_ct <- delta_ct[, !colnames(delta_ct) %in% "PPIB"]  # Exclude PPIB
  delta_ct <- sapply(delta_ct, function(gene) gene - PPIB)
  
  # Step 2: Calculate ΔΔct
  mean_dct <- colMeans(delta_ct, na.rm = TRUE)  # Calculate mean delta_ct for each gene
  delta_delta_ct <- t(apply(delta_ct, 1, function(row) row - mean_dct))
  
  # Step 3: Calculate the result
  result <- 2^(-delta_delta_ct)
  
  return(as.data.frame(result))
}

# Use the function
result_lab <- calculate_expression(gml)

# View the result
#print(result_lab)
#summary(result_lab)

#check the distributions of genes
plot_list <- lapply(names(result_lab[,-1]), function(colname){
  ggplot(result_lab[,-1], aes_string(x = colname)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(title = paste("Histogram of", colname))
})


grid.arrange(grobs = plot_list, ncol = 3)



#join the normalised data
result <- rbind(result_field, result_lab)


# When working with DDct Values
# A ddct of 0 = no change in fold expression
# negative ddct = upregulation
# postive ddct = downregulation
# Check if each column is numeric and, if so, multiply by -1
#result[] <- lapply(result, function(x) if(is.numeric(x)) -1 * x else x)

# Therefore deciding to transform data by *-1
# then positive = upgregulation - more intuitive for interpretation
rm(gmf,data, gf, gl, gml, hm, result_field, result_lab)

hm <- rbind(field, lab)

#remove non normalised genes
hm <-hm %>%
  dplyr::select(-all_of(Genes_v))

# join the normalised genes
hm <- hm %>%
  left_join(result, by = "Mouse_ID")

rm(lab, field, result)

## ---------------------------------------------------------------------------------------------------
##save the imputed data 
write.csv(hm, "Data/Data_output/2.1.norm_MICE_data_set.csv", row.names = FALSE)


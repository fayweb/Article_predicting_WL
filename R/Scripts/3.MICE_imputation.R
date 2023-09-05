
source("R/Scripts/2.Data_normalization.R")

library(VIM)

hm_genes <- hm[,c("Mouse_ID", Genes_v)]

#pattern_na <-as.data.frame(md.pattern(field_genes))
sapply(hm_genes, function(x) sum(is.na(x)))

genes <- hm_genes[, -1]

# The frequency distribution of the missing cases per variable can be obtained 
# as:
init <- mice(genes, maxit = 0)

#we want to impute only the specific variables
meth <- init$method


## ---------------------------------------------------------------------------------------------------
aggr_plot <- aggr(hm_genes, col=c('navyblue','red'), numbers=TRUE, 
                  sortVars=TRUE, labels=names(hm_genes), cex.axis=.7, gap=3, 
                  ylab=c("Histogram of missing data","Pattern"))


## ---------------------------------------------------------------------------------------------------
marginplot(hm_genes[c(6,8)])


## ---------------------------------------------------------------------------------------------------
# removing il 10
genes <- genes[, !(names(genes) %in% "IL.10")]


# removed(because of large missing numbers)
# m=5 refers to the number of imputed datasets. Five is the default value.
igf <- mice(genes, m = 5, seed = 500) # method = meth,


## ---------------------------------------------------------------------------------------------------
summary(igf)


## ---------------------------------------------------------------------------------------------------
# to check each column with imputed data
## igf$imp$IFNy
#Now we can get back the completed dataset using the complete()
complete_genes <- complete(igf, 1)

#sapply(complete_field, function(x) sum(is.na(x)))
#visualize missingness
vis_dat(complete_genes)


## ---------------------------------------------------------------------------------------------------
# remove IL.10 from the selection vector
Genes_v <- Genes_v[!Genes_v == "IL.10"]

#remove the non imputed genes from our data set
hm <-hm %>%
  dplyr::select(-all_of(Genes_v))

# join the imputed genes
Mouse_ID <- hm_genes$Mouse_ID
result <- data.frame(Mouse_ID, complete_genes)

hm <- hm %>%
  left_join(result, by = "Mouse_ID")

## ---------------------------------------------------------------------------------------------------
plot(igf)


## ---------------------------------------------------------------------------------------------------
xyplot(igf, IFNy ~ IL.13 + IRGM1 + MUC2, pch=18,cex=1)


## ---------------------------------------------------------------------------------------------------
xyplot(igf,IFNy ~ IL.13 + PRF1 + CASP1, pch=18,cex=1)


## ---------------------------------------------------------------------------------------------------
stripplot(igf, pch = c(20,21), cex = 1.2)


## ---------------------------------------------------------------------------------------------------
#bwplot(igf)


## ---------------------------------------------------------------------------------------------------
densityplot(igf)


## ---------------------------------------------------------------------------------------------------
 ##save the imputed data 
write.csv(hm, "Data/Data_output/2.imputed_MICE_data_set.csv", 
          row.names = FALSE)


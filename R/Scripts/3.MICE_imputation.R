
source("R/Scripts/2.Data_normalization.R")

hm <- hm_norm

rm(hm_norm)

hm_genes <- rbind(gmf, gml)

# remove HKG
hm_genes <- hm_genes %>%
  dplyr::select(-c(PPIB, GAPDH))

#dplyr::select(-Mouse_ID)
# looking at patterns of nas)
#pattern_na <-as.data.frame(md.pattern(field_genes))
sapply(hm_genes, function(x) sum(is.na(x)))

genes <- hm_genes[,-1]

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
#genes <- genes %>%
 # dplyr::select(-IL.10)
# removed already at previous step (because of large missing numbers)
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

##transform genes with multiplaction by -1! We want to avoid confusion in the subsequent
#data analysis 
#the higher the numbers represent higher fold expression
complete_genes <- unique(complete_genes)
complete_genes <- complete_genes %>%
    mutate_if(is.numeric, funs(.*-1))

## ---------------------------------------------------------------------------------------------------
#remove the non imputed genes from our data set
hm_selection_g <- hm_selection_g %>%
  dplyr::select(-c("IFNy", "CXCR3", "IL.6", "IL.13", #"IL.10",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF", "origin"))
# add the new imputed genes to the data
hm_selection_g <- hm_selection_g %>%
    left_join(complete_genes, by = "Mouse_ID")

hm_selection_g <- unique(hm_selection_g)

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
write.csv(hm_select, "Data/Data_output/2.imputed_MICE_data_set.csv", 
          row.names = FALSE)


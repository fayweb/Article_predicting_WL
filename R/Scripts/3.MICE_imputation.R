## ---------------------------------------------------------------------------------------------------
library(mice)
library(VIM)



## ---------------------------------------------------------------------------------------------------
hm <- read.csv("Data/Data_output/2.1.norm_MICE_data_set.csv")


## ---------------------------------------------------------------------------------------------------
# Vectors for selecting genes
#Lab genes
# The measurements of IL.12 and IRG6 are done with an other assay and will 
#ignore for now
Gene_lab   <- c("IFNy", "CXCR3", "IL.6", "IL.13", # "IL.10",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF") #"IL.12", "IRG6")

Genes_wild   <- c("IFNy", "CXCR3", "IL.6", "IL.13",# "IL.10", 
                  "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                  "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                  "TICAM1", "TNF") #, "IL.12", "IRG6")

Facs_lab <- c("CD4", "Treg", "Div_Treg", "Treg17", "Th1", 
                    "Div_Th1", "Th17", "Div_Th17", "CD8", "Act_CD8", 
                    "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8","Treg_prop", 
                    "IL17A_CD4")  

Facs_wild <- c( "Treg", "CD4", "Treg17", "Th1", "Th17", "CD8",
                     "Act_CD8", "IFNy_CD4", "IL17A_CD4", "IFNy_CD8")


## ---------------------------------------------------------------------------------------------------
field <- hm %>%
  dplyr::filter(origin == "Field") 

field <- unique(field)
genes_mouse_field <- field %>%
  dplyr::select(c(Mouse_ID, "IFNy", "CXCR3", "IL.6", "IL.13",# "IL.10", 
                  "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                  "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                  "TICAM1", "TNF"))

genes <- genes_mouse_field  %>%
  dplyr::select(-Mouse_ID)
#remove rows with only nas
genes <- genes[,colSums(is.na(genes))<nrow(genes)]
#remove colums with only nas 
genes <- genes[rowSums(is.na(genes)) != ncol(genes), ]
genes_mouse_field <- genes_mouse_field[row.names(genes), ]
##select same rows in the first table
field <- field[row.names(genes), ]


###############lab
#select the genes and lab muce
lab <- hm %>%
  dplyr::filter(origin == "Lab", Position == "mLN") #selecting for mln to avoid
# duplicates
lab <- unique(lab)
gene_lab_mouse <- lab %>%
  dplyr::select(c(Mouse_ID, "IFNy", "CXCR3", "IL.6", "IL.13",# "IL.10", 
                  "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                  "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                  "TICAM1", "TNF"))

gene_lab_mouse <- unique(gene_lab_mouse)

genes_lab <- gene_lab_mouse[, -1]

#remove rows with only nas
genes_lab <- genes_lab[,colSums(is.na(genes_lab))<nrow(genes_lab)]

#remove colums with only nas 
genes_lab <- genes_lab[rowSums(is.na(genes_lab)) != ncol(genes_lab), ]

genes_lab <- unique(genes_lab)

#select same rows in the first table
gene_lab_mouse <- gene_lab_mouse[row.names(genes_lab), ]

##select same rows in the first table
lab <- lab[row.names(genes_lab), ]

hm_genes <- rbind(gene_lab_mouse, genes_mouse_field)

hm_selection_g <- rbind(lab, field)

genes <- hm_genes %>% 
  left_join(hm_selection_g %>%
              dplyr::select(c(Mouse_ID, origin)),
            by = "Mouse_ID")


genes <- genes %>%
  dplyr::select(-Mouse_ID)
  
genes$origin <- as.factor(genes$origin)

#dplyr::select(-Mouse_ID)
# looking at patterns of nas)
#pattern_na <-as.data.frame(md.pattern(field_genes))
sapply(hm_genes, function(x) sum(is.na(x)))


## ---------------------------------------------------------------------------------------------------
# Discarding the origin 
#genes <- genes %>% dplyr::select(-origin)


#had to remove as they were disturbing the imputation: Worms_presence, 
#MC.Eimeria.FEC,  Heligmosomoides_polygurus, Zfy2, Y,  MpiC,
#vis_miss(field)
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


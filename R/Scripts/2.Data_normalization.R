source("R/Scripts/1.Data_prep.R")


## ---------------------------------------------------------------------------------------------------
library(mice)
library(stringr)



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


## ---------------------------------------------------------------------------------------------------


####################### field ##########################
df <- gmf

############### IL.13 
# dct
df <- df %>%
    mutate(IL.13_dct = IL.13 - GAPDH)
# max of dct
dct_max <- max(df$IL.13_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(IL.13_N = 2^ - (IL.13_dct - dct_max)) %>%
    mutate(IL.13_N = round(IL.13_N, digits = 2))

############# I tried writing a nice function but I failed so now I am
# going t repeat this many many times

############### IFNy
# dct
df <- df %>%
    mutate(IFNy_dct = IFNy - GAPDH)

# max of dct
dct_max <- max(df$IFNy_dct, na.rm = TRUE)

#fold gene expression
df <- df %>%
    mutate(IFNy_N = 2^ - (IFNy_dct - dct_max)) %>%
    mutate(IFNy_N = round(IFNy_N, digits = 2))


############### CXCR3
df <- df %>%
    mutate(CXCR3_dct = CXCR3 - GAPDH)

# max of dct
dct_max <- max(df$CXCR3_dct, na.rm = TRUE)

#fold gene expression
df <- df %>%
    mutate(CXCR3_N = 2^ - (CXCR3_dct - dct_max)) %>%
    mutate(CXCR3_N = round(CXCR3_N, digits = 2))

############### IL.6
df <- df %>%
    mutate(IL.6_dct = IL.6 - GAPDH)

# max of dct
dct_max <- max(df$IL.6_dct, na.rm = TRUE)

#fold gene expression
df <- df %>%
    mutate(IL.6_N = 2^ - (IL.6_dct - dct_max)) %>%
    mutate(IL.6_N = round(IL.6_N, digits = 2))

##############  IL1RN
df <- df %>%
    mutate(IL1RN_dct = IL1RN - GAPDH)

# max of dct
dct_max <- max(df$IL1RN_dct, na.rm = TRUE)

#fold gene expression
df <- df %>%
    mutate(IL1RN_N = 2^ - (IL1RN_dct - dct_max)) %>%
    mutate(IL1RN_N = round(IL1RN_N, digits = 2))


##############  CASP1
df <- df %>%
    mutate(CASP1_dct = CASP1 - GAPDH)

# max of dct
dct_max <- max(df$CASP1_dct, na.rm = TRUE)

#fold gene expression
df <- df %>%
    mutate(CASP1_N = 2^ - (CASP1_dct - dct_max)) %>%
    mutate(CASP1_N = round(CASP1_N, digits = 2))


##############  CXCL9
df <- df %>%
    mutate(CXCL9_dct = CXCL9 - GAPDH)

# max of dct
dct_max <- max(df$CXCL9_dct, na.rm = TRUE)

#fold gene expression
df <- df %>%
    mutate(CXCL9_N = 2^ - (CXCL9_dct - dct_max)) %>%
    mutate(CXCL9_N = round(CXCL9_N, digits = 2))

##############  IDO1
df <- df %>%
    mutate(IDO1_dct = IDO1 - GAPDH)
# max of dct
dct_max <- max(df$IDO1_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(IDO1_N = 2^ - (IDO1_dct - dct_max)) %>%
    mutate(IDO1_N = round(IDO1_N, digits = 2))

##############  IRGM1
df <- df %>%
    mutate(IRGM1_dct = IRGM1 - GAPDH)
# max of dct
dct_max <- max(df$IRGM1_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(IRGM1_N = 2^ - (IRGM1_dct - dct_max)) %>%
    mutate(IRGM1_N = round(IRGM1_N, digits = 2))

##############  MPO
df <- df %>%
    mutate(MPO_dct = MPO - GAPDH)
# max of dct
dct_max <- max(df$MPO_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(MPO_N = 2^ - (MPO_dct - dct_max)) %>%
    mutate(MPO_N = round(MPO_N, digits = 2))


##############  MUC2
df <- df %>%
    mutate(MUC2_dct = MUC2 - GAPDH)
# max of dct
dct_max <- max(df$MUC2_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(MUC2_N = 2^ - (MUC2_dct - dct_max)) %>%
    mutate(MUC2_N = round(MUC2_N, digits = 2))

##############  MUC5AC
df <- df %>%
    mutate(MUC5AC_dct = MUC5AC - GAPDH)
# max of dct
dct_max <- max(df$MUC5AC_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(MUC5AC_N = 2^ - (MUC5AC_dct - dct_max)) %>%
    mutate(MUC5AC_N = round(MUC5AC_N, digits = 2))

##############  MYD88
df <- df %>%
    mutate(MYD88_dct = MYD88 - GAPDH)
# max of dct
dct_max <- max(df$MYD88_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(MYD88_N = 2^ - (MYD88_dct - dct_max)) %>%
    mutate(MYD88_N = round(MYD88_N, digits = 2))

##############  NCR1
df <- df %>%
    mutate(NCR1_dct = NCR1 - GAPDH)
# max of dct
dct_max <- max(df$NCR1_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(NCR1_N = 2^ - (NCR1_dct - dct_max)) %>%
    mutate(NCR1_N = round(NCR1_N, digits = 2))


##############  PRF1
df <- df %>%
    mutate(PRF1_dct = PRF1 - GAPDH)
# max of dct
dct_max <- max(df$PRF1_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(PRF1_N = 2^ - (PRF1_dct - dct_max)) %>%
    mutate(PRF1_N = round(PRF1_N, digits = 2))

##############  RETNLB
df <- df %>%
    mutate(RETNLB_dct = RETNLB - GAPDH)
# max of dct
dct_max <- max(df$RETNLB_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(RETNLB_N = 2^ - (RETNLB_dct - dct_max)) %>%
    mutate(RETNLB_N = round(RETNLB_N, digits = 2))



##############  SOCS1
df <- df %>%
    mutate(SOCS1_dct = SOCS1 - GAPDH)
# max of dct
dct_max <- max(df$SOCS1_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(SOCS1_N = 2^ - (SOCS1_dct - dct_max)) %>%
    mutate(SOCS1_N = round(SOCS1_N, digits = 2))

##############  TICAM1
df <- df %>%
    mutate(TICAM1_dct = TICAM1 - GAPDH)
# max of dct
dct_max <- max(df$TICAM1_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(TICAM1_N = 2^ - (TICAM1_dct - dct_max)) %>%
    mutate(TICAM1_N = round(TICAM1_N, digits = 2))

##############  TNF
df <- df %>%
    mutate(TNF_dct = TNF - GAPDH)
# max of dct
dct_max <- max(df$TNF_dct, na.rm = TRUE)
#fold gene expression
df <- df %>%
    mutate(TNF_N = 2^ - (TNF_dct - dct_max)) %>%
    mutate(TNF_N = round(TNF_N, digits = 2))


df -> df_field



## ---------------------------------------------------------------------------------------------------
################################# lab 
# select first the field samples 


df_lab <- gml

############### IFNy
df_lab <- df_lab %>%
    mutate(IFNy_dct = IFNy - PPIB)
# max of dct
dct_max <- max(df_lab$IFNy_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(IFNy_N = 2^ - (IFNy_dct - dct_max)) %>%
    mutate(IFNy_N = round(IFNy_N, digits = 2))



############### CXCR3
df_lab <- df_lab %>%
    mutate(CXCR3_dct = CXCR3 - PPIB)
# max of dct
dct_max <- max(df_lab$CXCR3_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(CXCR3_N = 2^ - (CXCR3_dct - dct_max)) %>%
    mutate(CXCR3_N = round(CXCR3_N, digits = 2))

############### IL.6
df_lab <- df_lab %>%
    mutate(IL.6_dct = IL.6 - PPIB)
# max of dct
dct_max <- max(df_lab$IL.6_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(IL.6_N = 2^ - (IL.6_dct - dct_max)) %>%
    mutate(IL.6_N = round(IL.6_N, digits = 2))


############## IL.13
df_lab <- df_lab %>%
    mutate(IL.13_dct = IL.13 - PPIB)
# max of dct
dct_max <- max(df_lab$IL.13_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(IL.13_N = 2^ - (IL.13_dct - dct_max)) %>%
    mutate(IL.13_N = round(IL.13_N, digits = 2))

##############  IL1RN
df_lab <- df_lab %>%
    mutate(IL1RN_dct = IL1RN - PPIB)
# max of dct
dct_max <- max(df_lab$IL1RN_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(IL1RN_N = 2^ - (IL1RN_dct - dct_max)) %>%
    mutate(IL1RN_N = round(IL1RN_N, digits = 2))

##############  CASP1
df_lab <- df_lab %>%
    mutate(CASP1_dct = CASP1 - PPIB)
# max of dct
dct_max <- max(df_lab$CASP1_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(CASP1_N = 2^ - (CASP1_dct - dct_max)) %>%
    mutate(CASP1_N = round(CASP1_N, digits = 2))


##############  CXCL9
df_lab <- df_lab %>%
    mutate(CXCL9_dct = CXCL9 - PPIB)
# max of dct
dct_max <- max(df_lab$CXCL9_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(CXCL9_N = 2^ - (CXCL9_dct - dct_max)) %>%
    mutate(CXCL9_N = round(CXCL9_N, digits = 2))


##############  IDO1
df_lab <- df_lab %>%
    mutate(IDO1_dct = IDO1 - PPIB)
# max of dct
dct_max <- max(df_lab$IDO1_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(IDO1_N = 2^ - (IDO1_dct - dct_max)) %>%
    mutate(IDO1_N = round(IDO1_N, digits = 2))

##############  IRGM1
df_lab <- df_lab %>%
    mutate(IRGM1_dct = IRGM1 - PPIB)
# max of dct
dct_max <- max(df_lab$IRGM1_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(IRGM1_N = 2^ - (IRGM1_dct - dct_max)) %>%
    mutate(IRGM1_N = round(IRGM1_N, digits = 2))

##############  MPO
df_lab <- df_lab %>%
    mutate(MPO_dct = MPO - PPIB)
# max of dct
dct_max <- max(df_lab$MPO_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(MPO_N = 2^ - (MPO_dct - dct_max)) %>%
    mutate(MPO_N = round(MPO_N, digits = 2))

##############  MUC2
df_lab <- df_lab %>%
    mutate(MUC2_dct = MUC2 - PPIB)
# max of dct
dct_max <- max(df_lab$MUC2_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(MUC2_N = 2^ - (MUC2_dct - dct_max)) %>%
    mutate(MUC2_N = round(MUC2_N, digits = 2))

##############  MUC5AC
df_lab <- df_lab %>%
    mutate(MUC5AC_dct = MUC5AC - PPIB)
# max of dct
dct_max <- max(df_lab$MUC5AC_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(MUC5AC_N = 2^ - (MUC5AC_dct - dct_max)) %>%
    mutate(MUC5AC_N = round(MUC5AC_N, digits = 2))


##############  MYD88
df_lab <- df_lab %>%
    mutate(MYD88_dct = MYD88 - PPIB)
# max of dct
dct_max <- max(df_lab$MYD88_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(MYD88_N = 2^ - (MYD88_dct - dct_max)) %>%
    mutate(MYD88_N = round(MYD88_N, digits = 2))


##############  NCR1
df_lab <- df_lab %>%
    mutate(NCR1_dct = NCR1 - PPIB)
# max of dct
dct_max <- max(df_lab$NCR1_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(NCR1_N = 2^ - (NCR1_dct - dct_max)) %>%
    mutate(NCR1_N = round(NCR1_N, digits = 2))

##############  PRF1
df_lab <- df_lab %>%
    mutate(PRF1_dct = PRF1 - PPIB)
# max of dct
dct_max <- max(df_lab$PRF1_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(PRF1_N = 2^ - (PRF1_dct - dct_max)) %>%
    mutate(PRF1_N = round(PRF1_N, digits = 2))

##############  RETNLB
df_lab <- df_lab %>%
    mutate(RETNLB_dct = RETNLB - PPIB)
# max of dct
dct_max <- max(df_lab$RETNLB_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(RETNLB_N = 2^ - (RETNLB_dct - dct_max)) %>%
    mutate(RETNLB_N = round(RETNLB_N, digits = 2))

##############  SOCS1
df_lab <- df_lab %>%
    mutate(SOCS1_dct = SOCS1 - PPIB)
# max of dct
dct_max <- max(df_lab$SOCS1_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(SOCS1_N = 2^ - (SOCS1_dct - dct_max)) %>%
    mutate(SOCS1_N = round(SOCS1_N, digits = 2))

##############  TICAM1
df_lab <- df_lab %>%
    mutate(TICAM1_dct = TICAM1 - PPIB)
# max of dct
dct_max <- max(df_lab$TICAM1_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(TICAM1_N = 2^ - (TICAM1_dct - dct_max)) %>%
    mutate(TICAM1_N = round(TICAM1_N, digits = 2))

##############  TNF
df_lab <- df_lab %>%
    mutate(TNF_dct = TNF - PPIB)
# max of dct
dct_max <- max(df_lab$TNF_dct, na.rm = TRUE)
#fold gene expression
df_lab <- df_lab %>%
    mutate(TNF_N = 2^ - (TNF_dct - dct_max)) %>%
    mutate(TNF_N = round(TNF_N, digits = 2))

df_lab <- df_lab[,-1]


## ---------------------------------------------------------------------------------------------------
#df_lab <- df_lab %>% 
#  dplyr::select(-c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",                
#"IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO",   "MUC2", "MUC5AC", "MYD88",
#"NCR1", "PRF1", "RETNLB", "SOCS1",   "TICAM1", "TNF", PPIB, contains("_dct")))

# remove ending _N
#df_lab <- df_lab %>%
#  rename_with(~str_remove(.x, "_N"))

#df_field <- df_field %>% 
#  dplyr::select(-c(all_of(Genes_wild), GAPDH, contains("_dct")))

# remove ending _N
#df_field <- df_field %>%
# rename_with(~str_remove(.x, "_N"))



###################### let's try just using DCT
df_lab <- df_lab %>% 
    dplyr::select(-c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10", "IL1RN","CASP1", 
                     "CXCL9", "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", 
                     "NCR1", "PRF1", "RETNLB", "SOCS1",  "TICAM1", "TNF", "PPIB", 
                     contains("_N")))
Mouse_ID <- gml$Mouse_ID

df_lab <- cbind(Mouse_ID, df_lab)

# remove ending _dct
df_lab <- df_lab %>%
    rename_with(~str_remove(.x, "_dct"))

df_field <- df_field %>% 
    dplyr::select(-c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10", "IL1RN","CASP1", 
                     "CXCL9", "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", 
                     "NCR1", "PRF1", "RETNLB", "SOCS1",  "TICAM1", "TNF", "GAPDH", 
                     contains("_N")))

# remove ending
df_field <- df_field %>%
    rename_with(~str_remove(.x, "_dct"))


# add the new genes to the complete data sets 
lab <- lab %>%
    dplyr::select(-all_of(Genes_v)) %>%
    left_join(df_lab, by = "Mouse_ID")


field <- field %>%
    dplyr::select(-c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10", "IL1RN","CASP1", 
                     "CXCL9", "IDO1", "IRGM1", "MPO",  "MUC2", "MUC5AC", "MYD88", 
                     "NCR1", "PRF1", "RETNLB", "SOCS1",  "TICAM1", "TNF")) %>%
    left_join(df_field, by = "Mouse_ID")


hm_norm <- rbind(lab, field)


## ---------------------------------------------------------------------------------------------------
##save the imputed data 
write.csv(hm_norm, "Data/Data_output/2.1.norm_MICE_data_set.csv", row.names = FALSE)


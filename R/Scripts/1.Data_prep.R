#BiocManager::install("limma")

library(mice)
library(stringr)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(tidyr)
library(janitor)
library(visdat)
library(corrplot)
library(RColorBrewer)
library(ggplot2)
library(VIM)
library(limma)


# Import data 
#Challenge = experimental challenge infection data

#SOTA = State of the Art, of wild data


Challenge <- read.csv("Data/Data_input/Challenge_infections.csv")
SOTA <- read.csv("Data/Data_input/SOTA_Data_Product.csv")

# Vectors for selecting genes
#Lab genes
# The measurements of IL.12 and IRG6 are done with an other assay and will 
#ignore for now
Genes_v   <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF") #"IL.12", "IRG6")


# Cleaning

## Challenge

Challenge <- Challenge %>%
    dplyr::mutate(Parasite_primary = case_when(
        primary_infection == "E64" ~ "E_ferrisi",
        primary_infection == "E88" ~ "E_falciformis",
        primary_infection == "Eflab" ~ "E_falciformis",
        primary_infection == "E139" ~ "E_ferrisi",
        primary_infection == "UNI" ~ "uninfected",
        TRUE ~ ""))

Challenge <- Challenge %>%
    dplyr::mutate(Parasite_challenge = case_when(    
        challenge_infection == "E64" ~ "E_ferrisi",
        challenge_infection == "E88" ~ "E_falciformis",
        challenge_infection == "Eflab" ~ "E_falciformis",
        challenge_infection == "E139" ~ "E_ferrisi",
        challenge_infection == "UNI" ~ "uninfected",
        TRUE ~ ""))

Challenge <- Challenge %>%
    dplyr::mutate(infection_history = case_when(
        Parasite_primary == "uninfected" & 
            Parasite_challenge == "uninfected" ~ "uninfected",
        Parasite_primary == "uninfected" & 
            Parasite_challenge == "E_ferrisi" ~ "uninfected_ferrisi",
        Parasite_primary == "uninfected" & 
            Parasite_challenge == "E_falciformis" ~ "uninfected_falciformis",
        Parasite_primary == "E_falciformis" & 
            Parasite_challenge == "E_falciformis" ~ "falciformis_falciformis",
        Parasite_primary == "E_falciformis" & 
            Parasite_challenge == "E_ferrisi" ~ "falciformis_ferrisi",
        Parasite_primary == "E_falciformis" & 
            Parasite_challenge == "uninfected" ~ "falciformis_uninfected",
        Parasite_primary == "E_ferrisi" & 
            Parasite_challenge == "E_falciformis" ~ "ferrisi_falciformis",
        Parasite_primary == "E_ferrisi" & 
            Parasite_challenge == "E_ferrisi" ~ "ferrisi_ferrisi",
        Parasite_primary == "E_ferrisi" &
            Parasite_challenge == "uninfected" ~ "ferrisi_uninfected",
        TRUE ~ ""))


Challenge[sapply(Challenge, is.infinite)] <- NA


# IF a mouse dies before the end of the experiment the weight on the dpis after
# death will be NA. For further analysis I am changing those dpis to NA as well
# so that I can get the max dpi for each mouse

# replace the 0 with NA (incorrect experimental measurement)
Challenge$weight[Challenge$weight == 0] <- NA
Challenge$relative_weight[Challenge$relative_weight == 0] <- NA

# change dpis to NA after death
Challenge$dpi[is.na(Challenge$weight)] <- NA

### Add the variable end weight (relative weight at day of sacrifice)
# start by adding the variable dpi_max which inficates the last day of each mouse
Challenge <- Challenge%>%
    dplyr::filter(!weight == "NA") %>% 
    dplyr::group_by(EH_ID, infection) %>%
    dplyr::mutate(max_dpi = max(dpi), WL_max = (100 - min(relative_weight)))

# Creating death variable
# Identify EH_IDs that exist for both 'primary' and 'challenge'
both_ids <- intersect(
    Challenge %>% filter(infection == "primary") %>% pull(EH_ID),
    Challenge %>% filter(infection == "challenge") %>% pull(EH_ID)
)

# Identify EH_IDs that exist for 'primary' but not for 'challenge'
missing_ids <- setdiff(
    Challenge %>% filter(infection == "primary") %>% pull(EH_ID),
    Challenge %>% filter(infection == "challenge") %>% pull(EH_ID)
)

# Add the 'death' column
Challenge <- Challenge %>%
    mutate(
        death = case_when(
            EH_ID %in% both_ids  ~ "challenge",
            EH_ID %in% missing_ids ~ "primary",
            TRUE ~ ""
        )
    )

#There are two measuremts for CXCR3
# We want to here keep the CXCR3_bio 
Challenge <- Challenge %>% 
    dplyr::select(-CXCR3)

#Now rename the CXCR3_bio to CXCR3
Challenge <- Challenge %>%
    dplyr::rename(CXCR3 = CXCR3_bio)


Challenge <- Challenge %>%
    dplyr::mutate(origin = "Lab")
SOTA <- SOTA %>%
    dplyr::mutate(origin = "Field")

# Adjust the parasite names to fit the lab
SOTA <- SOTA %>%
    dplyr::mutate(eimeriaSpecies = case_when(
        eimeriaSpecies == "Negative" ~ "uninfected",
        eimeriaSpecies == "" ~ "NA",
        TRUE ~ eimeriaSpecies))

# Rename column names to match each other
Challenge <- Challenge %>% 
    dplyr::rename(Mouse_ID = EH_ID, delta_ct_cewe_MminusE = delta, 
                  MC.Eimeria = Eim_MC, Feces_Weight = feces_weight)


# Here I create a new column, where we get the actual infection status
# According to the melting curve for eimeria 
Challenge <- Challenge %>%
  dplyr::mutate(current_infection = case_when(
    Parasite_challenge == "E_ferrisi" & MC.Eimeria == "TRUE" ~ "E_ferrisi",
    Parasite_challenge == "E_ferrisi" & MC.Eimeria == "FALSE" ~ "uninfected",
    Parasite_challenge == "E_falciformis" & MC.Eimeria == "TRUE" ~ "E_falciformis",
    Parasite_challenge == "E_falciformis" & MC.Eimeria == "FALSE" ~ "uninfected",
    Parasite_challenge == "uninfected" & MC.Eimeria == "TRUE" ~ "E_falciformis",
    Parasite_challenge == "uninfected" & MC.Eimeria == "FALSE" ~ "uninfected",
    TRUE ~ Parasite_challenge
  ))


Challenge <- Challenge %>%
    dplyr::mutate(immunization = case_when(
    infection_history == "falciformis_ferrisi" ~ "heterologous",
    infection_history == "ferrisi_falciformis" ~ "heterologous",
    infection_history == "falciformis_uninfected" ~ "uninfected",
    infection_history == "ferrisi_uninfected" ~ "uninfected",
    infection_history == "ferrisi_ferrisi" ~ "homologous",
    infection_history == "falciformis_falciformis" ~ "homologous",
    infection_history == "uninfected_falciformis" ~ "naive",
    infection_history == "uninfected_ferrisi" ~ "naive",
    infection_history == "uninfected" ~ "uninfected",
    TRUE ~ "NA"
))


# Join wild and lab data 
length(intersect(colnames(Challenge), colnames(SOTA)))

#37 intersecting columns

# create a function that is the opposite of intersect
outersect <- function(x, y) {
    sort(c(setdiff(x, y),
           setdiff(y, x)))
}

length(outersect(colnames(Challenge), colnames(SOTA)))

# 148 dissimilar columns

#expected columns:
37 + 148 #185

# now join the two data sets
data <- full_join(Challenge, SOTA, 
                  by = intersect(colnames(SOTA), colnames(Challenge)))

write.csv(x = data, file = "Data/Data_output/full_data_prior_imputation.csv", 
          row.names = FALSE)

data <- data %>%
  dplyr::select(-ends_with("_N"))
  
rm(Challenge, SOTA)

hm <- data


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

hm <- rbind(lab, field)



gene_correlation <- lab %>% 
    filter(infection == "challenge", dpi == max_dpi) %>%
    ungroup() %>%
    dplyr::select(all_of(Genes_v))

# draw correlation between the genes
gene_correlation <- as.matrix(cor(gene_correlation, 
                                  use="pairwise.complete.obs"))

# load the function to calculate the p value for correlations
source("R/Functions/p_value_for_correlations.R")

# matrix of the p-value of the correlatio
p.mat <- cor.mtest(gene_correlation)

corrplot(gene_correlation, 
         method = "circle",  #method of the plot, "color" would show colour gradient
         tl.col = "black", tl.srt=45, #colour of labels and rotation
         col = brewer.pal(n = 8, name ="RdYlBu"), #colour of matrix
         order="hclust", #hclust reordering
         p.mat = p.mat, sig.level = 0.01, insig = "blank",
         addCoef.col = 'black',
         number.cex=0.5, 
         title = "Lab") 
#Add significance level to the correlogram
#remove the values that are insignificant


hm_genes <- hm[,c("Mouse_ID", Genes_v, "GAPDH", "PPIB")]

#pattern_na <-as.data.frame(md.pattern(field_genes))
sapply(hm_genes, function(x) sum(is.na(x)))

genes <- hm_genes[, !colnames(hm_genes) %in% "Mouse_ID"]

##########################################################
#Normalization

genes_matrix <- as.matrix(genes)


genes_matrix <- limma::normalizeQuantiles(genes_matrix)


genes <- as.data.frame(genes_matrix)

# The frequency distribution !f the missing cases per variable can be obtained 
# as:
init <- mice(genes, maxit = 0)

#we want to impute only the specific variables
meth <- init$method

# removing il 10
genes <- genes[, !(names(genes) %in% "IL.10")]

# removed(because of large missing numbers)
# m=5 refers to the number of imputed datasets. Five is the default value.
igf <- mice(genes, m = 5, seed = 500) # method = meth,

summary(igf)

## igf$imp$IFNy
#Now we can get back the completed dataset using the complete()
complete_genes <- complete(igf, 1)

# Transform all columns by multiplying by -1 (higher expression - more positive)
complete_genes[] <- -1 * complete_genes

# join the imputed genes
Mouse_ID <- hm_genes$Mouse_ID
result <- data.frame(Mouse_ID, complete_genes)

hm_imp <- hm %>%
    dplyr::select(-all_of(Genes_v)) %>%
    left_join(result, by = "Mouse_ID")

lab <- hm_imp %>%
    filter(origin== "Lab")

Genes_v <- setdiff(Genes_v, "IL.10")

gene_correlation <- lab[,Genes_v]

# draw correlation between the genes
gene_correlation <- as.matrix(cor(gene_correlation, 
                                  use="pairwise.complete.obs"))

# matrix of the p-value of the correlatio
p.mat <- cor.mtest(gene_correlation)

corrplot(gene_correlation, 
         method = "circle",  #method of the plot, "color" would show colour gradient
         tl.col = "black", tl.srt=45, #colour of labels and rotation
         col = brewer.pal(n = 8, name ="RdYlBu"), #colour of matrix
         order="hclust", #hclust reordering
         p.mat = p.mat, sig.level = 0.01, insig = "blank",
         addCoef.col = 'black',
         number.cex=0.5,
         title = "Lab") 
#Add significance level to the correlogram
#remove the values that are insignificant

hm_imp <- hm_imp %>%
    mutate(GAPDH = GAPDH.y, PPIB = PPIB.y) %>%
    dplyr::select(-c(GAPDH.x, GAPDH.y, PPIB.x, PPIB.y))

# add missing infection intensities
ii <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/CEWE_FECES_infection_intensities")

ii$Mouse_ID <- gsub(pattern = "AA_", replacement = "AA", x = ii$Mouse_ID)


hm_imp <- hm_imp %>%
    left_join(ii, by = intersect(colnames(hm_imp), colnames(ii)))


write.csv(hm_imp, "Data/Data_output/imputed_clean_data.csv", row.names = FALSE)




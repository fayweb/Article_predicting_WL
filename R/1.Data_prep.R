
## ----load_libraries, echo=FALSE, include = FALSE----------------------------------------------------
library(pheatmap)
library(tidyr)
library(dplyr)
library(janitor)
library(visdat)


## ---------------------------------------------------------------------------------------------------
Challenge <- read.csv("Data/Data_input/Challenge_infections.csv")
SOTA <- read.csv("Data/Data_input/SOTA_Data_Product.csv")
# Vectors for selecting genes
#Lab genes
# The measurements of IL.12 and IRG6 are done with an other assay and will 
#ignore for now
Gene_lab   <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF") #"IL.12", "IRG6")
Genes_wild   <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10", 
                  "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                  "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                  "TICAM1", "TNF") #, "IL.12", "IRG6")
Facs_lab <- c("Position", "CD4", "Treg", "Div_Treg", "Treg17", "Th1", 
                    "Div_Th1", "Th17", "Div_Th17", "CD8", "Act_CD8", 
                    "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8","Treg_prop", 
                    "IL17A_CD4")  
Facs_wild <- c( "Treg", "CD4", "Treg17", "Th1", "Th17", "CD8",
                     "Act_CD8", "IFNy_CD4", "IL17A_CD4", "IFNy_CD8")


## ---- echo = FALSE----------------------------------------------------------------------------------
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
### Add the variable end weight (relative weight at day of sacrifice)
# start by adding the variable dpi_max which inficates the last day of each mouse
Challenge <- Challenge %>% 
  dplyr::filter(!weight == "NA") %>%
  dplyr::group_by(EH_ID, infection) %>%
  dplyr::mutate(dpi_max = max(dpi), WL_max = (min(relative_weight) - 100))
#somehow case when dplyr ways didn't work for me and this is the only solution 
#that is functional
#let's filter for the challenge mice
chal <- Challenge %>% filter(infection == "challenge")
#now only select the rows where the dpi is equal to the dpi max for each mouse
chal <- chal[chal$dpi == chal$dpi_max, ] 
#now we can easily add the variable end weight to each mouse (which in now equal
#to the weight on the dpi = dpi_max)
chal <- chal %>% dplyr::mutate(end_rel_weight = (weight/weight_dpi0) * 100)
#let'repeat for the prim 
#let's filter for the challenge mice
prim <- Challenge %>% filter(infection == "primary")
#now only select the rows where the dpi is equal to the dpi max for each mouse
prim <- prim[prim$dpi == prim$dpi_max, ] 
#now we can easily add the variable end weight to each mouse (which in now equal
#to the weight on the dpi = dpi_max)
prim <- prim %>% 
  dplyr::mutate(end_rel_weight = (weight/weight_dpi0) * 100)
c <- rbind(chal, prim)
#now jon it to the challenge infections
c %>% 
  dplyr::select(EH_ID, end_rel_weight) %>%
  right_join(Challenge) -> Challenge
Challenge <- unique(Challenge)
#There are two measuremts for CXCR3
# We want to here keep the CXCR3_bio 
Challenge <- Challenge %>% 
  dplyr::select(-CXCR3)
#Now rename the CXCR3_bio to CXCR3
Challenge <- Challenge %>%
  dplyr::rename(CXCR3 = CXCR3_bio)
rm(c, chal, prim)


## ---- message = FALSE, echo = FALSE-----------------------------------------------------------------
length(intersect(colnames(Challenge), colnames(SOTA)))
#37 intersecting columns
# create a function that is the opposite of intersect
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
length(outersect(colnames(Challenge), colnames(SOTA)))
# 149 dissimilar columns
# add a column that indicates where the samples are from 
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
#expected columns:
37 + 149 #186
# now join the two data sets
data <- full_join(Challenge, SOTA, 
                  by = intersect(colnames(SOTA), colnames(Challenge)))
write.csv(data, "Data/Data_output/1.MICE_cleaned_data.csv", row.names = FALSE)




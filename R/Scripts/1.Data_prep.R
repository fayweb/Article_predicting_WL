library(tidyverse)
library(tidyr)
library(dplyr)
library(janitor)
library(visdat)


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

data <- data %>%
  dplyr::select(-ends_with("_N"))
  
rm(Challenge, SOTA)

write.csv(data, "Data/Data_output/1.MICE_cleaned_data.csv", row.names = FALSE)




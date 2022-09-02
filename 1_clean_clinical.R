library(readxl)
library(tidyverse)
library(naniar)

{
    clinical.pt1 <- read_excel('Duke_Breast_Cancer_MRI--Clinical_min_copy.xlsx')
    birads <- read_excel('Duke_Breast_Cancer_MRI--Clinical_BIRADSdata_copy.xlsx')
    
    # Combine clinical(part1) and birads data
    clinical <- cbind(clinical.pt1, birads)
    clinical_summary <- do.call(cbind, lapply(clinical, summary))
}

#### DEFINE FUNCTIONS ##############################################################################
clean.clinical <- function(clinical.df) {
    # Drop ultrasound data
    clinical.df <- clinical.df %>% select(-contains('US_'))       #removes bolus vol column too
    # Drop MRI technical data
    clinical.df <- clinical.df[, -which(names(clinical.df) %in% c(
        'TR','TE','Slice_Thickness','Rows','Columns','Reconstruction_Diameter'))]
    # Drop mammography data
    clinical.df <- clinical.df %>% select(-contains('Mammo_'))
    clinical.df <- clinical.df[, -which(names(clinical.df) %in% c('Breast_Density','Age_at_mammo'))]
    # Further remove variables
    cols.drop <- c(
        'Metastatic_at_Presentation','Days_to_MRI','Position','Bilateral_Information',
        'Bilateral_diff_rec_status', 'Bilateral_side_annotated','Bilateral_OtherSide_side_of_cancer',
        'Bilateral_OtherSide_oncotype_score','Bilateral_OtherSide_nottingham_grade',
        'Bilateral_OtherSide_ER','Bilateral_OtherSide_PR','Bilateral_OtherSide_HER2',
        'Bilateral_OtherSide_mol_subtype','Multicentric/Multifocal'
        )
    clinical.df <- clinical.df[, -which(names(clinical.df) %in% cols.drop)]
    
    return(clinical.df)
    }
    
clean.race <- function(clinicaldata) {
    clinicaldata$Race_cleaned <- ifelse(
        clinicaldata$Race_and_Ethnicity == 1, 1,
        ifelse(clinicaldata$Race_and_Ethnicity == 2, 2, 
               0)
        )
    clinicaldata <- clinicaldata %>% relocate(Race_cleaned, .after = Menopause)
    clinicaldata$Race_and_Ethnicity <- NULL
    return(clinicaldata)
    }

create.TherapyColumns <- function(cleanedClinical) {
    # Create df of just therapies
    clinical_therapies <- cleanedClinical %>% select(contains('adj'))
    clinical_neoadj <- clinical_therapies[, which(names(clinical_therapies) %in% c(
        'Neoadjuvant_Radiation_Therapy',
        'Neoadjuvant_Chemotherapy',
        'Neoadjuvant_Endocrine_Therapy',
        'Neoadjuvant_Anti-Her2_Neu_Therapy'
        )
        )
        ]
    clinical_neoadj[is.na(clinical_neoadj)] <- 0
    clinical_neoadj$num_neoadjuvant <- rowSums(clinical_neoadj)
    clinical_neoadj$rcd_neoadjuvant <- ifelse(clinical_neoadj$num_neoadjuvant > 0, 1, 0) # 1=yes, 0=no
    clinical_therapies$Received_Neoadjuvant <- clinical_neoadj$rcd_neoadjuvant
    cleanedClinical$Received_Neoadjuvant <- clinical_neoadj$rcd_neoadjuvant
    
    clinical_adj <- clinical_therapies[, which(names(clinical_therapies) %in% c(
        'Adjuvant_Radiation_Therapy',
        'Adjuvant_Chemotherapy',
        'Adjuvant_Endocrine_Therapy',
        'Adjuvant_Anti-Her2_Neu_Therapy'
        )
        )
        ]
    clinical_adj$num_adjuvant <- rowSums(clinical_adj)
    clinical_adj$rcd_adjuvant <- ifelse(clinical_adj$num_adjuvant > 0, 1, 0) # 1=yes, 0=no
    clinical_therapies$Received_Adjuvant <- clinical_adj$rcd_adjuvant
    cleanedClinical$Received_Adjuvant <- clinical_adj$rcd_adjuvant
    
    returns <- list(
        'clinical' = cleanedClinical,
        'neoadj' = clinical_neoadj,
        'adjuvant' = clinical_adj
        )
    return(returns)
    }

# Convert variables to factors
convert.toFactors <- function(combined.clinical) {
    cols.numerics <- c(
        'Patient_ID','age','age_last_contact','Date_of_Birth','Oncotype_score','Days_to_Surgery',
        'Days_to_local_recurrence','Days_to_distant_recurrence','Days_to_death',
        'Days_to_last_local_recurree_free_assessment','Days_to_last_distant_recurree_free_assessment',
        'Days_to_last_contact'
        )
    clinical.factors <- combined.clinical[, -which(names(combined.clinical) %in% cols.numerics)]
    col.factors <- names(clinical.factors)
    # Convert columns to factors
    combined.clinical[col.factors] <- lapply(combined.clinical[col.factors], factor)
    
    return(combined.clinical)
    }
    

calculate.age <- function(clinicalData) {
    clinicalData$age <- clinicalData$Date_of_Birth/(-365.25)
    attach(clinicalData)
    clinicalData <- clinicalData %>% relocate(age, .after = Patient_ID)
    days_since_Dx <- Days_to_last_contact/365.25
    clinicalData$age_last_contact <- age + days_since_Dx
    clinicalData <- clinicalData %>% relocate(age_last_contact, .after = age)
    
    return(clinicalData)
    }

# Standardize continuous variables
std.contvars <- function(clinicalvars) {
    continuous.df <- select_if(clinicalvars, is.numeric)
    continuous.summary <- do.call(cbind, lapply(continuous.df, summary))
    continuous.scaled <- scale(continuous.df, center = TRUE, scale = TRUE)
    continuous.scaled.df <- as.data.frame(continuous.scaled)
    
    returns <- list(
        'continuous' = continuous.df,
        'continuous_summary' = continuous.summary,
        'continuous_scaled' = continuous.scaled.df
        )
    return(returns)
    }


### Clean for basic analyses -----------------------------------------------------------------------
drop.Adjuvant <- function(combined.clinical, adj.df) {
    combined.clinical <- combined.clinical[ , -which(names(combined.clinical) %in% names(adj.df))]
    combined.clinical$Received_Adjuvant <- NULL
    
    return(combined.clinical)
    }

# Drop patients who received more than 1 form of neoadjuvant therapy
drop.multipleNeo <- function(combined.clinical, neoadj.df) {
    combined.clinical$num_neoadj <- neoadj.df$num_neoadjuvant
    combined.clinical <- combined.clinical[!combined.clinical$num_neoadj > 1, ]
    
    return(combined.clinical)
    }

keep.basics <- function(clinical.cleaned) {
    col.select <- c(
        'Patient_ID','age', 'Oncotype_score',
        'Menopause','Race_cleaned','Histologic_type','ER','PR','HER2','Mol_Subtype',
        'Tumor_Location','Staging_Tumor_Size','Staging_Nodes','Staging_Metastasis',
        'Nottingham_grade','Tumor_Grade_(tubule)','Tumor_Grade_(nuclear)','Tumor_Grade_(mitotic)',
        'Neoadjuvant_Chemotherapy','Neoadjuvant_Anti-Her2_Neu_Therapy','Neoadjuvant_Radiation_Therapy',
        'Neoadjuvant_Endocrine_Therapy', 'Known_Ovarian_Status', 'Therapeutic_or_Prophylactic_Oophorectomy',
        'Received_Neoadjuvant'
        )
    clinical.basic <- clinical.cleaned[, names(clinical.cleaned) %in% col.select]
    
    return(clinical.basic)
    }


#### CALL FUNCTIONS ################################################################################

clinical <- clean.clinical(clinical)

clinical <- clean.race(clinical)

returns.combineData <- create.TherapyColumns(clinical)
clinical <- returns.combineData$clinical
neoadjuvant.df <- returns.combineData$neoadj
adjuvant.df <- returns.combineData$adjuvant

clinical <- convert.toFactors(clinical)

clinical <- calculate.age(clinical)

returns.stdCont <- std.contvars(clinical)
continuous <- returns.stdCont$continuous
cont.1Dsummary <- returns.stdCont$continuous_summary
cont.scaled <- returns.stdCont$continuous_scaled


clinical <- drop.Adjuvant(clinical, adjuvant.df)

clinical <- drop.multipleNeo(clinical, neoadjuvant.df)

clinical.full.cleaned <- clinical
write.table(
    clinical.full.cleaned, 
    file='duke_clinical_cleaned.tsv', 
    quote=FALSE, 
    sep='\t', 
    row.names = FALSE
    )

### Save a condensed version!! ---------------------------------------------------------------------
clinical.short <- keep.basics(clinical)
clinical <- clinical.short

write.table(
    clinical, 
    file='duke_clinical_basics.tsv', 
    quote=FALSE, 
    sep='\t', 
    row.names = FALSE
)


#### INDICATE NEOADJUVANT THERAPY TYPE #############################################################
neoadj.df <- clinical.full.cleaned %>% select (contains('neoadj'))
attach(neoadj.df)
neoadj.df$neoadj_type <-
    ifelse(num_neoadj == 1 & Neoadjuvant_Radiation_Therapy == 1, 'R',
           ifelse(num_neoadj == 1 & Neoadjuvant_Endocrine_Therapy == 1, 'E',
                  ifelse(num_neoadj == 1 & Neoadjuvant_Chemotherapy == 1, 'C',
                         ifelse(num_neoadj == 1 & `Neoadjuvant_Anti-Her2_Neu_Therapy` == 1, 'H',
                                'None'

                                # ifelse(num_neoadj == 2 & (Neoadjuvant_Radiation_Therapy == 1 & Neoadjuvant_Chemotherapy == 1), 'C,R',
                                #        ifelse(num_neoadj == 2 & (Neoadjuvant_Endocrine_Therapy == 1 & Neoadjuvant_Chemotherapy == 1), 'C,E',
                                #               ifelse(num_neoadj == 2 & (`Neoadjuvant_Anti-Her2_Neu_Therapy` == 1 & Neoadjuvant_Chemotherapy == 1), 'C,H',
                                #
                                #                      ifelse(num_neoadj == 3 &
                                #                                 (Neoadjuvant_Endocrine_Therapy == 1 & `Neoadjuvant_Anti-Her2_Neu_Therapy` == 1 & Neoadjuvant_Chemotherapy == 1), 'C,E,H',
                                #                             ifelse(num_neoadj == 3 &
                                #                                        (Neoadjuvant_Radiation_Therapy == 1 & `Neoadjuvant_Anti-Her2_Neu_Therapy` == 1 & Neoadjuvant_Chemotherapy == 1), 'C,H,R', 'None'
                                                                        )
                                                                   )
                                                            )
                                                     )
                  #                             )
                  #                      )
                  #               )
                  #        )
                  # )
#            
# # Radiation and endocrine don't occur together
# # Radiation and Anti-HER2 only occur in a 3-combo therapy
# # Endocrine and Anti-HER2 only occur in a 3-combo therapy
# neoadj.df$neoadj_type <- as.factor(neoadj.df$neoadj_type)
# sum(neoadj.df$num_neoadj > 1) 
# sum(neoadj.df$num_neoadj == 1)
# sum(neoadj.df$num_neoadj < 1)
# 
# 
# clinical$Neoadjuvant_Type <- neoadj.df$neoadj_type
# clinical <- clinical %>% relocate(Neoadjuvant_Type, .after = Received_Neoadjuvant)



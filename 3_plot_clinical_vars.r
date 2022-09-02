library(readxl)
library(tidyverse)
library(naniar)
library(ggplot2)
library(DataExplorer)

source('prep_clinical.R')


#### CONTINUOUS ####################################################################################
cont.scaled %>% plot_density()
cont.scaled %>% plot_correlation()

attach(continuous)
par(mfrow=c(2,2))
{
    boxplot(age,
            main = 'Age',
            xlab = 'years',
            horizontal = TRUE
            )
    hist(age,
         main = 'Age',
         xlab = 'years',
         )
    boxplot(age_last_contact,
            main = 'Age of last contact',
            xlab = 'years',
            horizontal = TRUE)
    boxplot(Oncotype_score,
            main = 'Oncotype Dx score',
            xlab = 'score',
            horizontal = TRUE)
    hist(Oncotype_score,
         main = 'Oncotype Score',
         xlab = 'score',
         breaks = 20
        )
    # boxplot(Bilateral_OtherSide_oncotype_score,
    #         main = 'Bilateral Other side Oncotype Dx score',
    #         xlab = 'score',
    #         horizontal = TRUE)
    boxplot(Days_to_Surgery,
            main = 'Days to surgery',
            xlab = 'days',
            horizontal = TRUE)
    boxplot(Days_to_local_recurrence,
            main = 'Days to local recurrence',
            xlab = 'days',
            horizontal = TRUE)
    boxplot(Days_to_distant_recurrence,
            main = 'Days to distant recurrence',
            xlab = 'days',
            horizontal = TRUE)
    boxplot(Days_to_death,
            main = 'Days to death',
            xlab = 'days',
            horizontal = TRUE)
    boxplot(Days_to_last_local_recurree_free_assessment,
            main = 'Days to last local recurrence-free assessment',
            xlab = 'days',
            horizontal = TRUE)
    boxplot(Days_to_last_distant_recurree_free_assessment,
            main = 'Days to last distant recurrence-free assessment',
            xlab = 'days',
            horizontal = TRUE)
}

par(mfrow=c(1,1))
attach(cont.scaled)
boxplot(Days_to_Surgery,
        Days_to_local_recurrence,
        Days_to_distant_recurrence,
        Days_to_death,
        Days_to_last_local_recurree_free_assessment,
        Days_to_last_distant_recurree_free_assessment,
        main = 'Comparison of days to event',
        names = c(
            'surgery',
            'local_recurrence',
            'distant_recurrence',
            'death',
            'last_local_recurree_free_assessment',
            'last_distant_recurree_free_assessment'
            ),
        las = 2,
        ylab = 'days (scaled)'
        )


### MISSINGNESS ####################################################################################
# Missingness graphs: https://cran.r-project.org/web/packages/naniar/vignettes/naniar-visualisation.html
gg_miss_var(clinical, show_pct = TRUE)
gg_miss_fct(clinical[, -1], Received_Neoadjuvant) +
    # scale_fill_gradient2(
    #     '% missing',
    #     mid=muted('purple4'),
    #     high='gold',
    #     midpoint=0,
    #     limits=c(0,  100)
    #     ) +
    labs(title = 'Missing Data in Duke Breast Cancer Cohort') +
    theme(plot.title = element_text(size = 12))


gg_miss_fct(clinical, Neoadjuvant_Type)
gg_miss_fct(clinical, Recurrence)


# How much of Oncotype score is missing?
neoadj.yes <- subset(clinical, Received_Neoadjuvant == 1, select = c(Oncotype_score, Received_Neoadjuvant))
sum(is.na(neoadj.yes$Oncotype_score))
207 / 243 * 100
neoadj.no <- subset(clinical, Received_Neoadjuvant == 0, select = c(Oncotype_score, Received_Neoadjuvant))
sum(is.na(neoadj.no$Oncotype_score))
370 / 592 * 100

# therapy.resp <- c(
#     'Pathologic_stage_T_Neoadjuvant',
#     'Pathologic_stage_N_Neoadjuvant',
#     'Pathologic_stage_M_Neoadjuvant',
#     'Overall_Near-complete_Response_Stricter',
#     'Overall_Near-complete_Response_Looser',
#     'Near-complete_Response_Graded_Measure'
#     )
# clinical.wTherapy.resp <- cbind(clinical, subset(clinical.full.cleaned, select = therapy.resp))
# gg_miss_var(subset(clinical.full.cleaned, select = therapy.resp), show_pct = TRUE)
# gg_miss_fct(clinical.wTherapy.resp, Received_Neoadjuvant)


### FACTOR DISTRIBUTION ############################################################################
attach(clinical)
par(mfrow=c(2,2))

#### Demographic Data ------------------------------------------------------------------------------
{
    table(Menopause)
    menopause_prop <- prop.table(table(Menopause))
    barplot(
        menopause_prop,
        main = 'menopause status',
        ylab = 'proportion',
        names = c('pre-','post-','N/A'),
        ylim = c(0,1)
        )
    
    table(Race_and_Ethnicity)
    race_prop <- prop.table(table(Race_and_Ethnicity))
    barplot(
        race_prop,
        main = 'race/ethnicity status',
        ylab = 'proportion',
        names = c('N/A','white','black','asian','native','hispanic','multi','hawa','amer indian'),
        ylim = c(0,1)
    )

    table(Race_cleaned)
    race_prop <- prop.table(table(Race_cleaned))
    barplot(
        race_prop,
        main = 'race/ethnicity status',
        ylab = 'proportion',
        names = c('other','white','black'),
        ylim = c(0,1)
    )
    
    
    table(Metastatic_at_Presentation)
    metastatic_prop <- prop.table(table(Metastatic_at_Presentation))
    barplot(
        metastatic_prop,
        main = 'metastasis status',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
}

#### Histology Data -------------------------------------------------------------------------------
{
    table(Histologic_type)
    hist_type_prop <- prop.table(table(Histologic_type))
    barplot(
        hist_type_prop,
        main = 'histologic type',
        ylab = 'proportion',
        names = c('DCIS','ductal','lobular','metaplastic',
                  # 'LCIS',
                  'tubular',
                  # 'mixed',
                  # 'micropapillary',
                  'colloid'
                  ),
        ylim = c(0,1)
        )

    table(Mol_Subtype)
    mol_subtype_prop <- prop.table(table(Mol_Subtype))
    barplot(
        mol_subtype_prop,
        main = 'molecular subtype',
        ylab = 'proportion',
        names = c('Luminal A','Luminal B','HER2','TNBC'),
        ylim = c(0,1)
        )
    
    table(ER)
    ER_prop <- prop.table(table(ER))
    barplot(
        ER_prop,
        main = 'ER receptor',
        ylab = 'proportion',
        names = c('negative','positive'),
        ylim = c(0,1)
    )

    table(PR)
    PR_prop <- prop.table(table(PR))
    barplot(
        PR_prop,
        main = 'PR receptor',
        ylab = 'proportion',
        names = c('negative','positive'),
        ylim = c(0,1)
    )
    
    table(HER2)
    HER2_prop <- prop.table(table(HER2))
    barplot(
        HER2_prop,
        main = 'HER2 receptor',
        ylab = 'proportion',
        names = c('negative','positive'),
        ylim = c(0,1)
    )
    
    table(Staging_Tumor_Size)
    stagingT_prop <- prop.table(table(Staging_Tumor_Size))
    barplot(
        stagingT_prop,
        main = 'Staging_Tumor_Size',
        ylab = 'proportion',
        ylim = c(0,1)
    )
    
    table(Staging_Nodes)
    stagingN_prop <- prop.table(table(Staging_Nodes))
    barplot(
        stagingN_prop,
        main = 'Staging_Nodes',
        ylab = 'proportion',
        ylim = c(0,1)
    )
    
    table(Staging_Metastasis)
    stagingM_prop <- prop.table(table(Staging_Metastasis))
    barplot(
        stagingM_prop,
        main = 'Staging_Metastasis',
        ylab = 'proportion',
        ylim = c(0,1)
    )
    
    table(`Tumor_Grade_(tubule)`)
    tumorT_prop <- prop.table(table(`Tumor_Grade_(tubule)`))
    barplot(
        tumorT_prop,
        main = 'Tumor_Grade_(tubule)',
        ylab = 'proportion',
        ylim = c(0,1)
    )
    
    table(`Tumor_Grade_(nuclear)`)
    tumorN_prop <- prop.table(table(`Tumor_Grade_(nuclear)`))
    barplot(
        tumorN_prop,
        main = 'Tumor_Grade_(nuclear)',
        ylab = 'proportion',
        ylim = c(0,1)
    )
    
    table(`Tumor_Grade_(mitotic)`)
    tumorM_prop <- prop.table(table(`Tumor_Grade_(mitotic)`))
    barplot(
        tumorM_prop,
        main = 'Tumor_Grade_(mitotic)',
        ylab = 'proportion',
        ylim = c(0,1)
    )
    
    table(Nottingham_grade)
    nottingham_prop <- prop.table(table(Nottingham_grade))
    barplot(
        nottingham_prop,
        main = 'Nottingham Grade',
        ylab = 'proportion',
        names = c('low','intermediate','high'),
        ylim = c(0,1)
    )
    
    table(Tumor_Location)
    location_prop <- prop.table(table(Tumor_Location))
    barplot(
        location_prop,
        main = 'Tumor Location',
        ylab = 'proportion',
        ylim = c(0,1)
    )
    
    table(`Multicentric/Multifocal`)
    multicentric_prop <- prop.table(table(`Multicentric/Multifocal`))
    barplot(
        multicentric_prop,
        main = 'Multicentric/Multifocal',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
    
}

#### Therapies -------------------------------------------------------------------------------------
###### Neoadjuvant therapy distribution ------------------------------
{
    table(Received_Neoadjuvant)
    neoT_prop <- prop.table(table(Received_Neoadjuvant))
    barplot(
        neoT_prop,
        main = 'received neoadjuvant therapy',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
    
    # Neoadjuvant type
    table(neoadj.df$neoadj_type)
    neoType_prop <- prop.table(table(neoadj.df$neoadj_type))
    barplot(
        neoType_prop,
        main = 'Neoadjuvant Therapy Type',
        ylab = 'proportion',
        ylim = c(0,1)
    )
    types <- c('None','Chemo','Endocrine','HER2','Radiation')#,'C,E','C,H','C,R')
    neoadj.type.df <- neoadj.df[ , 11:12]
    ggplot(neoadj.type.df, aes(neoadj_type)) +
        scale_x_discrete(limits = types) +
        geom_bar() +
        ggtitle('Distribution of Neoadjuvant Therapy Types')
   
   
    # radiation
    table(Neoadjuvant_Radiation_Therapy)
    neoT_rad_prop <- prop.table(table(Neoadjuvant_Radiation_Therapy))
    barplot(
        neoT_rad_prop,
        main = 'neoadjuvant radiation therapy',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
    table(Clinical_Response_Neoadjuvant)
    neoT_rad_clinical_prop <- prop.table(table(Clinical_Response_Neoadjuvant))
    barplot(
        neoT_rad_clinical_prop,
        main = 'clinical response neoadjuvant radiation therapy',
        ylab = 'proportion',
        names = c('complete','not complete','no imaging'),
        las = 2,
        ylim = c(0,1)
    )
    table(Pathologic_Response_Neoadjuvant)
    neoT_rad_path_prop <- prop.table(table(Pathologic_Response_Neoadjuvant))
    barplot(
        neoT_rad_path_prop,
        main = 'pathologic response neoadjuvant radiation therapy',
        ylab = 'proportion',
        names = c('complete','not complete','DCIS only','LCIS only','unavailable'),
        las = 2,
        ylim = c(0,1)
    )
    # chemotherapy
    table(Neoadjuvant_Chemotherapy)
    neoT_chemo_prop <- prop.table(table(Neoadjuvant_Chemotherapy))
    barplot(
        neoT_chemo_prop,
    main = 'neoadjuvant chemotherapy',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
    # endocrine therapy
    table(Neoadjuvant_Endocrine_Therapy)
    neoT_endo_prop <- prop.table(table(Neoadjuvant_Endocrine_Therapy))
    barplot(
        neoT_endo_prop,
        main = 'neoadjuvant endocrine therapy',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
    # anti-Her2 therapy
    table(`Neoadjuvant_Anti-Her2_Neu_Therapy`)
    neoT_her2_prop <- prop.table(table(`Neoadjuvant_Anti-Her2_Neu_Therapy`))
    barplot(
        neoT_her2_prop,
        main = 'neoadjuvant anti-HER2 therapy',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
    
}

###### Adjuvant therapy distribution ------------------------------
# {
#     table(Received_Adjuvant)
#     adjT_prop <- prop.table(table(Received_Adjuvant))
#     barplot(
#         adjT_prop,
#         main = 'received adjuvant therapy',
#         ylab = 'proportion',
#         names = c('no','yes'),
#         ylim = c(0,1)
#     )
#     # radiation
#     table(Adjuvant_Radiation_Therapy)
#     adjT_rad_prop <- prop.table(table(Adjuvant_Radiation_Therapy))
#     barplot(
#         adjT_rad_prop,
#         main = 'adjuvant radiation therapy',
#         ylab = 'proportion',
#         names = c('no','yes'),
#         ylim = c(0,1)
#     )
#     # chemotherapy
#     table(Adjuvant_Chemotherapy)
#     adjT_chemo_prop <- prop.table(table(Adjuvant_Chemotherapy))
#     barplot(
#         adjT_chemo_prop,
#         main = 'adjuvant chemotherapy',
#         ylab = 'proportion',
#         names = c('no','yes'),
#         ylim = c(0,1)
#     )
#     # endocrine therapy
#     table(Adjuvant_Endocrine_Therapy)
#     adjT_endo_prop <- prop.table(table(Adjuvant_Endocrine_Therapy))
#     barplot(
#         adjT_endo_prop,
#         main = 'adjuvant endocrine therapy',
#         ylab = 'proportion',
#         names = c('no','yes'),
#         ylim = c(0,1)
#     )
#     # anti-Her2 therapy
#     table(`Adjuvant_Anti-Her2_Neu_Therapy`)
#     adjT_her2_prop <- prop.table(table(`Adjuvant_Anti-Her2_Neu_Therapy`))
#     barplot(
#         adjT_her2_prop,
#         main = 'neoadjuvant anti-HER2 therapy',
#         ylab = 'proportion',
#         names = c('no','yes'),
#         ylim = c(0,1)
#     ) 
# }

#### Recurrence/Surgery distribution ---------------------------------------------------------------
{
    table(Recurrence)
    recurrence_prop <- prop.table(table(Recurrence))
    barplot(
        recurrence_prop,
        main = 'recurrence',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )  

    table(Surgery)
    surgery_prop <- prop.table(table(Surgery))
    barplot(
        surgery_prop,
        main = 'surgery',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
    
    table(Definitive_Surgery_Type)
    surgery_type_prop <- prop.table(table(Definitive_Surgery_Type))
    barplot(
        surgery_type_prop,
        main = 'surgery type',
        ylab = 'proportion',
        names = c('BCS','mastectomy'),
        ylim = c(0,0.6)
    )
}

# Recurrence and neoadjuvant therapy (0=no, 1=yes)
table(clinical$Received_Neoadjuvant, clinical$Recurrence)


### ADDITIONAL ANALYSES ############################################################################
#### Out of those that received/did not receive neoadjuvant therapy...------------------------------
{
    no_neoadj <- subset(
        clinical,
        Received_Neoadjuvant==0,
        select = c(Received_Adjuvant, Surgery, Definitive_Surgery_Type, Recurrence)
        )
    table(no_neoadj$Received_Adjuvant)
    no_neoadj_adj <- prop.table(table(no_neoadj$Received_Adjuvant))
    # barplot(
    #     no_neoadj_adj,
    #     main = 'Adjuvant therapy following no neoadjuvant therapy',
    #     ylab = 'proportion',
    #     names = c('no','yes'),
    #     ylim = c(0,1)
    # )
    table(no_neoadj$Surgery)
    no_neoadj_surgery <- prop.table(table(no_neoadj$Surgery))
    barplot(
        no_neoadj_surgery,
        main = 'Surgery following no neoadjuvant therapy',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
    table(no_neoadj$Definitive_Surgery_Type)
    no_neoadj_surgery_type <- prop.table(table(no_neoadj$Definitive_Surgery_Type))
    barplot(
        no_neoadj_surgery_type,
        main = 'Surgery type following no neoadjuvant therapy',
        ylab = 'proportion',
        names = c('BCS','mastectomy'),
        ylim = c(0,1)
    )
    table(no_neoadj$Recurrence)
    no_neoadj_rec <- prop.table(table(no_neoadj$Recurrence))
    barplot(
        no_neoadj_rec,
        main = 'Recurrence if no neoadjuvant therapy',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
}

{
    yes_neoadj <- subset(
        clinical,
        Received_Neoadjuvant==1,
        select = c(Received_Adjuvant, Surgery, Definitive_Surgery_Type, Recurrence)
        )
    table(yes_neoadj$Received_Adjuvant)
    # yes_neoadj_adj <- prop.table(table(yes_neoadj$Received_Adjuvant))
    # barplot(
    #     yes_neoadj_adj,
    #     main = 'Adjuvant therapy following neoadjuvant therapy',
    #     ylab = 'proportion',
    #     names = c('no','yes'),
    #     ylim = c(0,1)
    # )
    table(yes_neoadj$Surgery)
    yes_neoadj_surgery <- prop.table(table(yes_neoadj$Surgery))
    barplot(
        yes_neoadj_surgery,
        main = 'Surgery following neoadjuvant therapy',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
    table(yes_neoadj$Definitive_Surgery_Type)
    yes_neoadj_surgery_type <- prop.table(table(yes_neoadj$Definitive_Surgery_Type))
    barplot(
        yes_neoadj_surgery_type,
        main = 'Surgery type following neoadjuvant therapy',
        ylab = 'proportion',
        names = c('BCS','mastectomy'),
        ylim = c(0,1)
    )
    table(yes_neoadj$Recurrence)
    yes_neoadj_rec <- prop.table(table(yes_neoadj$Recurrence))
    barplot(
        yes_neoadj_rec,
        main = 'Recurrence if neoadjuvant therapy',
        ylab = 'proportion',
        names = c('no','yes'),
        ylim = c(0,1)
    )
}
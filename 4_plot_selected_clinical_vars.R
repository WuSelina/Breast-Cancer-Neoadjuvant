library(readxl)
library(tidyverse)
library(naniar)
library(ggplot2)
library(DataExplorer)

source('prep_clinical.R')


### MISSINGNESS ####################################################################################
# Missingness graphs: https://cran.r-project.org/web/packages/naniar/vignettes/naniar-visualisation.html
gg_miss_fct(clinical, Received_Neoadjuvant)


#### CONTINUOUS ####################################################################################
cont.interest <- c('Age','Oncotype_score')

cont.scaled %>% plot_density()
cont.scaled %>% plot_correlation()

attach(continuous)
par(mfrow=c(2,2))

{
    boxplot(age,
            main = 'Age',
            xlab = 'years',
            horizontal = TRUE)
    
    boxplot(Oncotype_score,
            main = 'Oncotype Dx score',
            xlab = 'score',
            horizontal = TRUE)
}


### FACTOR DISTRIBUTION ############################################################################
dis.interest <- c('Menopause','Race_and_Ethnicity','Histologic_type','Mol_Subtype','ER','PR','HER2',
                  '','','','','','',''
                  )


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
    table(Received_Neoadjuvant, Menopause)
    
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
    table(Received_Neoadjuvant, Race_and_Ethnicity)
    
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
    table(Neoadjuvant_Type)
    neoType_prop <- prop.table(table(Neoadjuvant_Type))
    barplot(
        neoType_prop,
        main = 'Neoadjuvant Therapy Type',
        ylab = 'proportion',
        ylim = c(0,1)
    )
    types <- c('None','C','E','H','R')#,'C,E','C,H','C,R')
    neoadj.type.df <- neoadj.df[ , 11:12]
    ggplot(neoadj.type.df, aes(Neoadjuvant_Type)) +
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
    table(Received_Neoadjuvant, Recurrence)
}



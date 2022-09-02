{
    library(readxl)
    library(tidyverse)
    library(mlbench)
    library(caret) 
}

### LOAD DATA ######################################################################################
{
    features.key <- read_excel('Duke_Breast_Cancer_MRI--Imaging_features_key.xlsx')
    features <- read_excel('Duke_Breast_Cancer_MRI--Imaging_Features.xlsx')
    
    # Convert columns to factors from key
    features.key$FeatureGroup <- as.factor(features.key$FeatureGroup)
    levels(features.key$FeatureGroup)
}

# Define Functions ---------------------------------------------------------------------------------
# Function to remove highly correlated features
remove_highCorr <- function(df) {
    correlationMatrix <- cor(df)
    highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.8)
    cat('TBremoved:', length(highlyCorrelated), 'Remaining:', length(df) - length(highlyCorrelated))

    new.df <- subset(df, select = -highlyCorrelated)
    
    returns <- list(
        'correlation_matrix' = correlationMatrix, 
        'high_correlated' = highlyCorrelated,
        'df_dropped_highcorr' = new.df
        )
    return(returns)
}

### FGT ENHANCEMENT FEATURES (284) #################################################################
# Enhancement (78) ---------------------------------------------------------------------------------
{
    fgt_enhance <- subset(features.key, FeatureGroup=='FGT Enhancement', select=c(FeatureNum,FeatureNames))
    fgt_enhance_RFs <- features[, names(features) %in% fgt_enhance$FeatureNames]
    fgt_enhance_RFs <- as.data.frame(fgt_enhance_RFs)
    {
        fgt_enhance_RFs$Grouping_based_mean_of_peak_enhancement_slope_3D_tissue_T1_Group_3 <- as.numeric(
            fgt_enhance_RFs$Grouping_based_mean_of_peak_enhancement_slope_3D_tissue_T1_Group_3    )
        
        fgt_enhance_RFs$Grouping_based_mean_of_washin_slope_3D_tissue_T1_Group_3 <- as.numeric(
            fgt_enhance_RFs$Grouping_based_mean_of_washin_slope_3D_tissue_T1_Group_3)
        
        fgt_enhance_RFs$Grouping_based_mean_of_peak_enhancement_slope_3D_tissue_PostCon_Group_3 <- as.numeric(
            fgt_enhance_RFs$Grouping_based_mean_of_peak_enhancement_slope_3D_tissue_PostCon_Group_3)
        
        fgt_enhance_RFs$Grouping_based_mean_of_washin_slope_3D_tissue_PostCon_Group_3 <- as.numeric(
            fgt_enhance_RFs$Grouping_based_mean_of_washin_slope_3D_tissue_PostCon_Group_3)
    }
    
    
    drop.FGT.enhance <- names(which(colSums(is.na(fgt_enhance_RFs)) > 0))
    fgt_enhance_RFs <- fgt_enhance_RFs[ , -which(names(fgt_enhance_RFs) %in% drop.FGT.enhance)]
    
    # Correlation matrix
    correlation.returns <- remove_highCorr(fgt_enhance_RFs)
    fgt.enhance.df <- correlation.returns$df_dropped_highcorr
    
    # fgt_enhance_RFs_scaled <- scale(fgt_enhance_RFs, center = TRUE, scale = TRUE)
}

# Texture (176) ------------------------------------------------------------------------------------
{
    fgt_texture <- subset(features.key, FeatureGroup=='FGT Enhancement Texture', select=c(FeatureNum,FeatureNames))
    fgt_texture_RFs <- features[, names(features) %in% fgt_texture$FeatureNames]
    fgt_texture_RFs <- as.data.frame(fgt_texture_RFs)
    
    # Correlation matrix
    correlation.returns <- remove_highCorr(fgt_texture_RFs)
    fgt.texture.df <- correlation.returns$df_dropped_highcorr
    
    # fgt_texture_RFs_scaled <- scale(fgt_texture_RFs, center = TRUE, scale = TRUE)
    
}

# Variation (30) -----------------------------------------------------------------------------------
{
    fgt_variation <- subset(features.key, FeatureGroup=='FGT Enhancement Variation', select=c(FeatureNum,FeatureNames))
    fgt_variation_RFs <- features[, names(features) %in% fgt_variation$FeatureNames]
    fgt_variation_RFs <- as.data.frame(fgt_variation_RFs)
    {
        fgt_variation_RFs$Grouping_based_variance_of_peak_enhancement_slope_3D_tissue_PostCon_Group_3 <- as.numeric(
            fgt_variation_RFs$Grouping_based_variance_of_peak_enhancement_slope_3D_tissue_PostCon_Group_3)
        
        fgt_variation_RFs$Grouping_based_variance_of_washin_slope_3D_tissue_PostCon_Group_3 <- as.numeric(
            fgt_variation_RFs$Grouping_based_variance_of_washin_slope_3D_tissue_PostCon_Group_3    )
        
        fgt_variation_RFs$Grouping_based_variance_of_peak_enhancement_slope_3D_tissue_T1_Group_3 <- as.numeric(
            fgt_variation_RFs$Grouping_based_variance_of_peak_enhancement_slope_3D_tissue_T1_Group_3    )
        
        fgt_variation_RFs$Grouping_based_variance_of_washin_slope_3D_tissue_T1_Group_3 <- as.numeric(
            fgt_variation_RFs$Grouping_based_variance_of_washin_slope_3D_tissue_T1_Group_3    )
    }
    
    drop.FGT.var <- names(which(colSums(is.na(fgt_variation_RFs)) > 0))
    fgt_variation_RFs <- fgt_variation_RFs[ , -which(names(fgt_variation_RFs) %in% drop.FGT.var)]
    
    
    # Correlation matrix
    correlation.returns <- remove_highCorr(fgt_variation_RFs)
    fgt.variation.df <- correlation.returns$df_dropped_highcorr
    
    # fgt_variation_RFs_scaled <- scale(fgt_variation_RFs, center = TRUE, scale = TRUE) 
}

# Combine FGT features -----------------------------------------------------------------------------
{
    fgt.combined <- cbind(fgt.enhance.df, fgt.texture.df, fgt.variation.df)
    fgt.scaled.df <- scale(fgt.combined, center = TRUE, scale = TRUE)
    
    # Correlation heatmap
    library(reshape2)
    library(ggplot2)
    corr.mat <- round(cor(fgt.scaled.df),2)
    melted.corr.mat <- melt(corr.mat)
    ggplot(melted.corr.mat, aes(x=Var1, y=Var2, fill=value)) + geom_tile()
    
    # hist(fgt.combined.thenscaled,
    #      main='Distribution of scaled & centered BPE Features',
    #      breaks= 100,
    #      xlab='Radiomic feature value',
    #      xlim=c(-6,6),
    #      freq=FALSE
    #      )
}


#### COORDINATE WITH CLINICAL DATA #################################################################
{
    clinical.OG <- read_excel('Duke_Breast_Cancer_MRI--Clinical_min_copy.xlsx')
    fgt.scaled.df <- as.data.frame(fgt.scaled.df)
    rownames(fgt.scaled.df) <- clinical.OG$Patient_ID
    
    source('prep_clinical.R')
    fgt.cleaned <- subset(fgt.scaled.df, rownames(fgt.scaled.df) %in% clinical$Patient_ID)
    
    # Merge imaging and clinical data
    fgt.cleaned$Patient_ID <- rownames(fgt.cleaned)
    df.merged <- merge(clinical, fgt.cleaned)
}

write.table(
    df.merged, 
    file='duke_clinical_imaging_merged.tsv', 
    quote=FALSE, 
    sep='\t', 
    row.names = FALSE
    )

write.table(
    fgt.cleaned, 
    file='duke_clinical_imaging.tsv', 
    quote=FALSE, 
    sep='\t', 
    row.names = TRUE
)


library(tidyverse)
# library(DataExplorer)

source('prep_imaging.R')    # sources the clinical data as well

fgt.cleaned$Patient_ID <- NULL
fgt.cleaned$rcd_neoadj <- df.merged$Received_Neoadjuvant

### TEST MULTICOLLINEARITY #########################################################################
## Variance Inflation Factors ------------------------------------------------------------------------
# library('olsrr')
# 
# # Create linear regression model
# model.all <- lm(
#     as.numeric(rcd_neoadj) ~ .,
#     data = fgt.cleaned
# )
# summary(model.all)
# 
# testVIF <- ols_vif_tol(model.all)
# testVIF <- testVIF[order(testVIF$VIF), ]
# keep.variables <- subset(testVIF, VIF < 7, select = Variables)      # Need source to say why
# 
# # Remove features with severe multicollinearity
# fgt.selected <- subset(fgt.cleaned, select = names(fgt.cleaned) %in% keep.variables$Variables)
# fgt.selected$rcd_neoadj <- fgt.cleaned$rcd_neoadj


### SELECTION USING ALGORITHMS #####################################################################
# Logit --------------------------------------------------------------------------------------------
library(glmnet)
logit.model <- glm(rcd_neoadj ~., family = 'binomial', data = fgt.selected)
options(scipen=999)
summary(logit.model)

var.importance <- caret::varImp(logit.model)
var.importance$var <- rownames(var.importance)
var.importance <- var.importance[order(var.importance$Overall, decreasing=TRUE), ]
row.names(var.importance) <- NULL


importance <- varImp(logit.model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(var.importance)


### PCA ############################################################################################
# library(devtools)
library(ggbiplot)

img.pca <- prcomp(fgt.selected[,-length(fgt.selected)])
summary(img.pca)
ggbiplot(img.pca, 
         obs.scale = 1,
         var.scale = 1,
         groups = clinical$Received_Neoadjuvant,
         ellipse = TRUE,
         var.axes=FALSE
         )
    # ggtitle('PCA Plot from 37 Features')


### CLUSTERING #####################################################################################
neoadj <- clinical[, c("Received_Neoadjuvant", "Known_Ovarian_Status")]
neoadj$colors_neo <- ifelse(neoadj$Received_Neoadjuvant == 0, 'blue', 'red')
heatmap(
    data.matrix(fgt.selected[,-length(fgt.selected)]),
    margins = c(3,2),
    RowSideColors = neoadj$colors,
    xlab = 'BPE features',
    ylab = 'Sample',
    labRow = c(''),
    labCol = c('')
    )
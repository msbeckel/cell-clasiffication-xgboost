#############################################################################################
#Aim: XGBoost Model for multiclass clasiffication of neuronal cell types
#Author: Maximiliano S. Beckel
#Creación: 2023/12/28
#Última modificación: 2024/01/08
#############################################################################################
# limpio la memoria
rm(list = ls()) # remove all objects
gc() # garbage collection

library(edgeR)
library(ggplot2)
library(matrixStats)
library(data.table)
library(stringr)
library(SingleCellExperiment)
library(scuttle)
library(robustbase)
library(ggplot2)
library(scater)
library(scran)
library(scuttle)
library(bluster)
library(data.table)
library(mlexperiments)
library(mllrnrs)
library(mlbench)
library(MLmetrics)
library(ParamHelpers)
library(caret)
options(bitmapType='cairo')
#set parallel backend
library(parallel)
#library(parallelMap) 
#parallelStartSocket(cpus = 25)


# Aqui se debe poner la carpeta de la computadora local
setwd("/home/maxibeckel/github/cell-classification-xgboost/exp") # Establezco el Working Directory

# cargo elsce object donde voy a entrenar el modelo
sce <- readRDS("../datasets/sce.rds")
cell_type         <- as.factor(colData(sce)[,"cluster"])


#
n = 12
exp = vector("list", length = n)

for(i in 1:n){
  e <- dir()[i]
  load(paste0(getwd(), "/", e, "/resultados.RData"))
  cm_test <- confusionMatrix(xgpred$data$response,xgpred$data$truth)
  #Long table
  tmp <- colnames(cm_test$byClass) %in% c("Specificity", "Pos Pred Value", "Neg Pred Value", "Prevalence", "Detection Rate", "Detection Prevalence")
  #bC <- as.data.frame(stack(cm_test$byClass[,!tmp]))
  bC <- as.data.frame(cm_test$byClass[,!tmp])
  #Change Cell type names
  tmp <- levels(cell_type)
  names(tmp) <- paste0("Class: ", seq(0, length(levels(cell_type))-1))
  rownames(bC) <- tmp[rownames(bC)]
  bC     <- bC[complete.cases(bC), ]
  bC$cell_type <- rownames(bC); rownames(bC) <- NULL
  bC$exp <- e

  exp[[i]] <- bC
  cat(paste0(e, "\tF1 score:\t", F1_Score(xgpred$data$response,xgpred$data$truth), "\n"))
}

experiments = do.call(rbind, exp)

long <- melt(setDT(experiments), id.vars = c("cell_type","exp"), variable.name = "measure")
long <- cbind(long, data.frame(strsplit2(long$exp, "H|T")))
long$X1 <- factor(long$X1, c("NI", "A", "M", "S"))

#all together
ggplot(long, aes(x = exp, y = value, colour = exp)) +
  geom_boxplot() +
  facet_wrap(~measure)+
  labs(title = "", x = "exp", y = "Performance")+
  theme_light()

#Only one measure
ggplot(long[measure == "F1",], aes(x = X1, y = value, colour = X1)) +
  geom_boxplot() +
  facet_wrap(~X2)+
  labs(title = "F1", x = "exp", y = "Performance")+
  theme_light()

#Repeated Measures ANOVA
model <- aov(value~factor(X1)+Error(factor(cell_type)), data = long[measure=="F1" & X2 == 1])
summary(model)


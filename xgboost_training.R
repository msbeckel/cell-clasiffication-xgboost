#############################################################################################
#Aim: XGBoost Model for multiclass clasiffication of neuronal cell types
#Author: Maximiliano S. Beckel
#Creación: 2023/12/28
#Última modificación: 2024/01/18
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
library(mlr)
library(mlexperiments)
library(mllrnrs)
library(mlbench)
library(MLmetrics)
library(ParamHelpers)
library(caret)
library(xgboost)
library(dplyr)
options(bitmapType='cairo')
#set parallel backend
library(parallel)
#library(parallelMap) 
#parallelStartSocket(cpus = 25)

# Main options of the experiment----
PARAM             <- list()
PARAM$imputation  <- "S" 
PARAM$HVG         <- "H3"
PARAM$cell_type   <- "T1"
PARAM$extra       <- ""#extra info of the experiment 
PARAM$experimento <- paste0(PARAM[1:4], collapse = "")
PARAM$seed        <- 44 

PARAM$pb <- list( objective="multi:softprob", eval_metric="auc", nrounds=100L, booster="gbtree", nthread = 25)
PARAM$po <- makeParamSet(makeNumericParam("eta",lower = 0.01,upper = 0.3),
                         makeIntegerParam("max_depth",lower = 3L,upper = 10L), 
                         makeNumericParam("min_child_weight",lower = 1L,upper = 10L), 
                         makeNumericParam("subsample",lower = 0.5,upper = 1), 
                         makeNumericParam("colsample_bytree",lower = 0.5,upper = 1))

PARAM$measures <- list(multiclass.aunp, bac)

#Set directory----

# Aqui se debe poner la carpeta de la computadora local
setwd("/home/maxibeckel/maestria_datos/neuronal-cell-type-classification/") # Establezco el Working Directory

# cargo elsce object donde voy a entrenar el modelo
sce <- readRDS("./datasets/sce.rds")

# creo la carpeta donde va el experimento
dir.create("./exp/", showWarnings = FALSE)
dir.create(paste0("./exp/", PARAM$experimento, "/"), showWarnings = FALSE)

# Establezco el Working Directory DEL EXPERIMENTO
setwd(paste0(getwd(), "/exp/", PARAM$experimento, "/"))


#Make dataset from sce and experiment options
#HVG options
hvgo <- c("HVGs.more", "subset", "HVGs.less")
names(hvgo) <- c("H1", "H2", "H3")
#Cell type options
ctype <- c("cluster", "cluster_agg", "cluster_agg2", "cluster_agg3")
names(ctype) <- c("T1","T2","T3","T4")
#Imputation options
impo <- c("logcounts", "magic", "alra", "saver")
names(impo) <- c("NI", "M", "A", "S")

#
dataset <- as.data.table(t(assay(sce, impo[PARAM$imputation])[rowData(sce)[hvgo[PARAM$HVG]][,1],]))
colnames(dataset) <- rowData(sce)[rowData(sce)[hvgo[PARAM$HVG]][,1], "gene"]
rownames(dataset) <- rownames(colData(sce))
#
cell_type         <- as.factor(colData(sce)[,ctype[PARAM$cell_type]])
dataset$cell_type <- cell_type
#
dataset$clase <- as.integer(dataset$cell_type) - 1L
colnames(dataset) <- make.names(colnames(dataset))
rm(sce);gc()

#Columns names
feature_cols <- setdiff(
  colnames(dataset),
  c("cell_type", "clase")
)
target_col <- "clase"

#Split train and test 
data_split <- splitTools::partition(
  y = dataset[, get(target_col)],
  p = c(train = 0.8, test = 0.2),
  type = "stratified",
  seed = PARAM$seed
)


#Hyperparameter tuning----
datos <- as.data.frame(dataset[,cell_type:=NULL])
traintask <- makeClassifTask (data = datos[data_split$train,],target = target_col)
testtask  <- makeClassifTask (data = datos[data_split$test,],target = target_col)

#Learner
lrn <- makeLearner("classif.xgboost",predict.type = "prob")
lrn$par.vals <- PARAM$pb

#CV
rdesc <- makeResampleDesc("CV",stratify = T,iters=5L)

#search strategy
ctrl <- makeTuneControlRandom(maxit = 10L)

#Parameter tuning
#parameter tuning
mytune <- tuneParams(learner = lrn, 
                     task = traintask, 
                     resampling = rdesc, 
                     measures = list(multiclass.aunp, bac), 
                     par.set = PARAM$po, 
                     control = ctrl, 
                     show.info = T)
paste0("Best performance: ", mytune$y) 

mytune$x

#Train final model

#set hyperparameters
lrn_tune <- setHyperPars(lrn,par.vals = mytune$x)

#train model
xgmodel <- mlr::train(learner = lrn_tune,task = traintask)

#view variable importance plot ----
mat <- xgb.importance (feature_names = feature_cols, model = xgmodel$learner.model)
p1 <- xgb.plot.importance (importance_matrix = mat[1:20])

png("feature_importance.png", width = 6, height = 4, units = "in", res = 600)
p1
dev.off()

#Performance in test dataset----
#predict model
xgpred <- predict(xgmodel,testtask)

#
cm_test <- confusionMatrix(xgpred$data$response,xgpred$data$truth)

# Performance Global
cm_test$overall
F1_Score(xgpred$data$response,xgpred$data$truth)

# Performance por Grupo
apply(cm_test$byClass, 2, summary)

# Plot Performance por Grupo
#Long table
tmp <- colnames(cm_test$byClass) %in% c("Specificity", "Pos Pred Value", "Neg Pred Value", "Prevalence", "Detection Rate", "Detection Prevalence")
bC <- as.data.frame(stack(cm_test$byClass[,!tmp]))

#Change Cell type names
tmp <- levels(cell_type)
names(tmp) <- paste0("Class: ", seq(0, length(levels(cell_type))-1))
bC$row <- tmp[bC$row]
bC     <- bC[complete.cases(bC), ]
#
p2 <- ggplot(bC, aes(x = col, y = value)) +
  geom_boxplot() +
  labs(title = "", x = "Metrica", y = "Performance")+
  theme_light()
png("by_cell_type_performance.png", width = 6, height = 4, units = "in", res = 600)
p2
dev.off()

#
outlier_limits <- bC %>%
  group_by(col) %>%
  summarize(lower_limit = quantile(value, 0.25) - 1.5 * IQR(value),
            upper_limit = quantile(value, 0.75) + 1.5 * IQR(value))

p3 <- ggplot(bC, aes(x = col, y = value)) +
  geom_boxplot() +
  geom_text(data = bC %>%
              left_join(outlier_limits, by = "col") %>%
              filter(value < lower_limit | value > upper_limit),
            aes(label = as.character(row)),
            position = position_jitter(height = 0.03),
            size = 2,
            show.legend = FALSE) +
  labs(title = "", x = "Metrica", y = "Performance")+
  theme_light()

png("by_cell_type_performance_out.png", width = 6, height = 4, units = "in", res = 600)
p3
dev.off()

#Save
save(mytune, xgmodel, xgpred, file = "resultados.RData")

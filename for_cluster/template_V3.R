# DOUBLE CROSS VALIDATION FOR PREDICTION WITH MIXOMICS SPLS OR SPLSDA
# Template for bash script
#

#Clear global environment
rm(list = ls(all.names = TRUE))

#free up memory and report the memory usage
gc()

#set working directory
setwd("/save/amarquissea/analysis/")

#libraries
.libPaths(c(.libPaths(),"/work/amarquissea/R/libraries_4.1.2"))

library(tidyverse) 
library(rsample)
library(mixOmics) 
library(tictoc)

#WARNING purrr::map and mixOmics::map
#WARNING dplyr::select and mixOmics::select

#theme plots
theme_set(theme_minimal())

#paste for log file 
paste("Nb cores :", nbCores)
paste("Outter loop :", repetsOutter, "repetitions -", foldsOutter, "folds")
paste("Inner loop :", repetsInner, "repetitions -", foldsInner, "folds")
paste("PREDICTIONS - PREDICTEURS")

#data path
df_prediction = read_tsv("/save/amarquissea/analysis/Data/V3/DATASET")

#mutate into factor and order levels 
df_prediction = df_prediction %>%
  mutate_at(c("annee","NULOBA_annee","annee_lignee","lignee"), as.factor)

df_prediction = df_prediction %>%
  mutate(annee = fct_relevel(annee, "2018", "2019", "2020")) %>%
  mutate(NULOBA_annee = fct_relevel(NULOBA_annee, "1_2018", "2_2018","3_2018","4_2018","5_2018","6_2018","1_2019","2_2019","3_2019","4_2019","5_2019","1_2020","2_2020","3_2020","4_2020")) %>%
  mutate(annee_lignee = fct_relevel(annee_lignee, "2018_rfi+", "2018_rfi-", "2019_rfi+", "2019_rfi-", "2020_rfi+", "2020_rfi-")) %>%
  mutate(lignee = fct_relevel(lignee, "rfi+", "rfi-"))%>%
  mutate(factRFI = fct_relevel(factRFI, "rfi+", "rfi-"))


#### ------ Hyperparameter selection function ------ ####


# Developped by Q. Le Graverand
# Search for the optimal number of components and the optimal number of variables

select_hyperpar<-function(tuning.object){
  error<-tuning.object[["error.rate"]]
  se<-tuning.object[["error.rate.sd"]]/sqrt(tuning.object[["call"]][["nrepeat"]])
  
  nvar_max<-nrow(error)#number of variables to test
  ncomp_max<-ncol(error)#number of components to test
  
  ncomp_opti<-1#initialise optimal number of components
  nvar_opti<-1#initialise optimal number of components
  moy_ref<-error[nvar_opti,ncomp_opti]#initialise optimal mean
  se_ref<-se[nvar_opti,ncomp_opti]#initialise optimal sd
  
  for (i in 1:ncomp_max){#to check if adding more components improve error rates
    for (j in 1:nvar_max){# to check if selecting more variables improve error rates
      #with a splsda, we want to check if a less parsimonious model has a lower R2, with a spls regression we are looking for lower MAE, MSE
      if (tuning.object[["measure"]]=="BER" | tuning.object[["measure"]]=="MAE" | tuning.object[["measure"]]=="MSE"){
        error_diff<-moy_ref-error[j,i]
      }
      
      #with a spls regression, we want to check if a less parsimonious model has a smaller absolute bias
      if (tuning.object[["measure"]]=="Bias"){
        error_diff<-abs(moy_ref)-abs(error[j,i])
      }
      
      #with a spls regression, we want to check if a less parsimonious model has a higher R2
      if (tuning.object[["measure"]]=="R2"){
        error_diff<-error[j,i]-moy_ref
      }
      
      tol<-se_ref+se[j,i]
      
      #if the gap between the two rates is higher than the difference between the two rates plus their two s.e, then the less parsimonious model is retained
      if (error_diff>tol){
        ncomp_opti<-i
        nvar_opti<-j
        
        moy_ref<-error[j,i]
        se_ref<-se[j,i]
      }
    }
  }
  nvar_opti<-rep(nvar_opti, ncomp_opti)#to repeat ncomp times, we could choose a number per component
  output=list(ncomp_select= ncomp_opti,
              nvar_select=as.numeric(row.names(error)[nvar_opti]))
  
  return(output)
}


#### ------ Sampling ------ ####

set.seed(53) #for repro, erase before running

#Split the dataset into training and testing sets x times, output : all the training and testing sets for the outter loop
cv = rsample::vfold_cv(df_prediction, repeats = repetsOutter, v = foldsOutter, strata = annee_lignee)


#### ------ Inner loop ------ ####

#Inner loop function

inner = function (objet, model, rangeX, rangeY) {
  
  #Predictors
  X.train = training(objet)[rangeX] 
  X.test = testing(objet)[rangeX] 
  
  #Variable to predict
  Y.train = training(objet)[[rangeY]] 
  Y.test = testing(objet)[[rangeY]]
  
  #Tuning step for spls or splsda model on training sets
  if (model == "splsda") {
    
    tuned = tune(method = "splsda", 
                 X= X.train, 
                 Y= Y.train,
                 ncomp = NCOMP,
                 nrepeat = repetsInner,
                 test.keepX = LIST_KEEPX,
                 validation = "Mfold", folds = foldsInner,
                 dist = "max.dist",
                 cpus = nbCores)
  }
  
  if (model == "spls") {
    
    tuned = tune(method = "spls",
                 X= X.train, 
                 Y= Y.train,
                 ncomp = NCOMP,
                 test.keepX = LIST_KEEPX,
                 nrepeat = repetsInner, 
                 validation = "Mfold", folds = foldsInner,
                 dist = "max.dist",
                 measure = "MAE",
                 cpus = nbCores)
  }
  
  gc()
  
  #hyperparameter selection
  hyperpar = select_hyperpar(tuned)
  
  list = list("X.train" = X.train, 
              "X.test" =X.test, 
              "Y.train" =Y.train, 
              "Y.test" =Y.test, 
              "id"=objet$id,
              "hyperparam" = hyperpar)
  
  return(list)
}

#Run inner function on all the sets of the outter loop

paste("Models for PREDICTIONS")
tic()
res = purrr::map(.x = cv$splits, .f = inner, model = "MODEL", rangeX = RANGE_X, rangeY = RANGE_Y)
toc()
names(res) = purrr::map (res, function (x) {paste(x$id[[1]],"_",x$id[[2]], sep = "")})

#### ------ model.fit function ------ ####

# Adapted from Q. Le Graverand
# Fit mixOmics::spls or splsda to a training set 


#function to fit a sPLSDA or sPLS         
model.fit <-function(object, model){
  
  #data partition
  Xtrain<- object$X.train
  Ytrain<- object$Y.train
  
  #to fit a sPLS discriminant analysis model
  if (model=="splsda" | model=="sPLSDA"){
    
    fitted_mod<-splsda(X = Xtrain, Y = Ytrain,
                       ncomp = object$hyperparam$ncomp_select, keepX = object$hyperparam$nvar_select,
                       near.zero.var = F)
  }
  #to fit a sPLS regression model
  if (model=="spls" | model=="sPLS"){
    
    fitted_mod<-spls(X = Xtrain, Y = Ytrain,
                     ncomp = object$hyperparam$ncomp_select, keepX = object$hyperparam$nvar_select,
                     near.zero.var = F)
  }
  return(fitted_mod)
}

#### ------ model.predict function ------ ####

# Adapted from Q. Le Graverand
# Fit mixOmics::spls or splsda to a training set 

#function to carry out predictions from a fitted model
model.predict <- function (object, fitted_model){
  
  #data test
  Xtest<-object$X.test
  
  #sPLSDA
  if (class(fitted_model)[1]== "mixo_splsda"){
    pred_obj<-predict(fitted_model, Xtest)
    Ypredicted<-pred_obj[["class"]][["max.dist"]][,dim(pred_obj[["class"]][["max.dist"]])[2]]
    Ypredicted<-as.character(Ypredicted)
  }
  #sPLSR
  if (class(fitted_model)[1]== "mixo_spls"){
    pred_obj<-predict(fitted_model, Xtest)
    Ypredicted<-pred_obj[["predict"]][,1,dim(pred_obj[["predict"]])[3]]
    Ypredicted<-as.numeric(as.character(Ypredicted))
  }
  return(Ypredicted)
}

#### ------ metrics function ------ ####

# Adapted from Q. Le Graverand
# Evaluate model performances 

#function to compute accuracy metrics, by comparing model predictions to real values
metrics<- function (object, Ypredicted){
  
  #data test
  Ytest<-object$Y.test
  
  #for classification
  if (class(Ypredicted)=="character"){
    #confusion matrix
    conf_mat<-get.confusion_matrix(truth=Ytest, predicted=Ypredicted)
    
    #error rates
    measure<-list(BER=get.BER(conf_mat))#balanced error rate
    for (gr in 1:nrow(conf_mat)){#error rate per class
      measure[[row.names(conf_mat)[gr]]]<-(sum(conf_mat[gr,])-diag(conf_mat)[gr])/sum(conf_mat[gr,])
    }
  }
  
  #for regression
  if (class(Ypredicted)=="numeric"){
    measure<-list(corr_Spearman=cor(Ypredicted, Ytest, method="spearman", use = "complete.obs"),#Spearman correlation
                  R2=cor(Ypredicted, Ytest, method="pearson", use = "complete")**2,#Coefficient of determination (ie square of the Pearson correlation)
                  MSE=mean((Ytest-Ypredicted)**2), #Mean Square Error
                  RMSE=sqrt(mean((Ytest-Ypredicted)**2)),#Root Mean Square Error
                  bias=mean(Ytest-Ypredicted),#Bias
                  MAE=mean(abs(Ytest-Ypredicted)))#Mean Absolute Error
  }
  
  return(measure)
}

#### ------ one.val function ------ ####

# Adapted from Q. Le Graverand
# workflow to acquire all model perfs  

#to predict the outcome Y for one testing set (corresponding to one fold and repetition)
one.val<-function(object,
                  model
){
  #model fitting
  fitted_model<-model.fit (object, model)
  
  #model predictions for the testing set
  models_predictions<-model.predict (object, fitted_model)
  
  #performance metrics of the model (comparing predictions to real values)
  perfs <-list(perfs = metrics (object, models_predictions))
  
  #Variables list stability per comp
  loadings = as.data.frame(t(fitted_model$loadings.star[[1]]))
  loadings = loadings %>% dplyr::mutate("comp" = paste("comp",row_number(), sep="")) %>% relocate(comp)
  loadings = list(loadings = loadings)
  
  #variables of importance for projection 
  vip = list(vip = vip(fitted_model))
  
  appended_object = append(object, perfs)
  appended_object = append(appended_object, vip)
  appended_object = append(appended_object, loadings)
  
  return(appended_object)
  
}

#### ------ Run one.val ------ ####

#run one.val on results from inner loop
paste("Perfs")
tic()
res = purrr::map (.x = res, .f = one.val, model ="MODEL")
toc()


#### ------ Plots and tables ------ ####

model = "MODEL"

#perf 
if (model == "spls") {
         
         perfs_spls = list(
           "id" = names(res), 
           "corr_spearman" = unlist(purrr::map(res, function (x) {x$perfs[[1]]})),
           "R2" = unlist(purrr::map(res, function (x) {x$perfs[[2]]})),
           "MSE" = unlist(purrr::map(res, function (x) {x$perfs[[3]]})),
           "RMSE" = unlist(purrr::map(res, function (x) {x$perfs[[4]]})),
           "bias" = unlist(purrr::map(res, function (x) {x$perfs[[5]]})),
           "MAE" = unlist(purrr::map(res, function (x) {x$perfs[[6]]}))
         )
}

if (model == "splsda") {
         
         perfs_splsda = list(
           "id" = names(res), 
           "BER" = unlist(purrr::map(res, function (x) {x$perfs[[1]]})),
           "rfi+" = unlist(purrr::map(res, function (x) {x$perfs[[2]]})),
           "rfi-" = unlist(purrr::map(res, function (x) {x$perfs[[3]]}))
         )
}

df = as_tibble(t(data.frame(matrix(unlist(perfs_MODEL), nrow=length(perfs_MODEL), byrow = TRUE))))
colnames(df) = names(perfs_MODEL)

df = reshape2::melt(df,id.vars="id")
df$value = as.numeric(df$value)

pdf("/work/amarquissea/output/V3/PREDICTIONS/MODEL_error_PREDICTIONS_PREDICTEURS_CORRECTION_V3.pdf")
df %>%
  ggplot(aes(variable, value)) +
  geom_jitter(color = "light grey", size=1.5, alpha=0.30) +
  geom_boxplot(aes(color=variable))+
  stat_summary(fun = mean, geom="point", shape=20, size=3, color="black", fill="black")+
  labs(y= "Metriques", color = "Mesures" )+
  ggtitle("Mesures pour MODEL PREDICTIONS PREDICTEURS", "CORRECTION")+
  facet_wrap(~variable, scale="free")
dev.off()

write_tsv(df, file = "/work/amarquissea/output/V3/PREDICTIONS/metrics_PREDICTIONS_PREDICTEURS_CORRECTION_V3.tsv")

#summary perfs for log file
df %>%
  group_by(variable) %>%
  summarise(mean = mean(value))

paste("Loading")
loadings = purrr::map(res, function (x) {x$loadings})
loadings = bind_rows(loadings, .id = "rep")

#summary number of components and variables per comp for log file
paste("ncomp")
loadings %>%
  group_by(rep) %>%
  summarise(nbcomp=n())%>%
  ungroup()%>%
  group_by(as.factor(nbcomp))%>%
  summarise(n=n())

pdf("/work/amarquissea/output/V3/PREDICTIONS/MODEL_ncomp_PREDICTIONS_PREDICTEURS_CORRECTION_V3.pdf")
loadings %>%
  group_by(rep)%>%
  summarise(ncomp = n())%>%
  ungroup()%>%
  ggplot(aes(as.factor(ncomp)))+
    geom_bar()+
  ggtitle("ncomp pour MODEL PREDICTIONS PREDICTEURS", "CORRECTION")
dev.off()

long_loadings = reshape2::melt(loadings,id.vars=c("rep", "comp"))

pdf("/work/amarquissea/output/V3/PREDICTIONS/MODEL_keepX_PREDICTIONS_PREDICTEURS_CORRECTION_V3.pdf")
long_loadings %>%
  filter(value != 0) %>%
  group_by(rep, comp)%>%
  summarise(keepX = n()) %>%
  ungroup() %>%
  ggplot(aes(as_factor(keepX))) +
  geom_bar()+
  ggtitle("keepX pour MODEL PREDICTIONS PREDICTEURS", "CORRECTION")
dev.off()

write_tsv(loadings, file = "/work/amarquissea/output/V3/PREDICTIONS/loadings_PREDICTIONS_PREDICTEURS_CORRECTION_V3.tsv")

#vip
vip = purrr::map(res, function (x) {data.frame(x$vip)})
vip = bind_rows(vip, .id = "rep")

write_tsv(vip, file = "/work/amarquissea/output/V3/PREDICTIONS/vip_PREDICTIONS_PREDICTEURS_CORRECTION_V3.tsv")

paste("-------end of script-----------")

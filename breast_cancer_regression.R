library(MASS)
library(FNN)
library(leaps)
library(caTools)
library(glmnet)
library(stats)
library(pls)

## Read Data
## \\ because of Windows
breast_data = as.data.frame(read.csv('.\\DATA\\r_breast_cancer.data'))

## Les donnees comptent differentes sections
## 1) Les moyennes des donnees
## 2) Les Standard Error des donnees (SE)
## 3) Les pires valeurs des donnees (worst)


# On definit la fonction qui nous permettra dextraire les meilleurs subsets
# pour chaque jeu de donnees choisi
getSubsets <- function(formula, train_data, main='', xlab='', ylab=''){
  reg.fit_only_means = regsubsets(formula, data=train_data, method='exhaustive', nvmax = NULL, nbest =1)
  plot(reg.fit_only_means,scale="r2", main=main, xlab=xlab, ylab=ylab + ' r2')
  plot(reg.fit_only_means, scale="adjr2", main=main, xlab=xlab, ylab=ylab + ' adjr2')
  plot(reg.fit_only_means, scale="bic", main=main, xlab=xlab, ylab=ylab + ' bic')
  plot(reg.fit_only_means, scale="Cp", main=main, xlab=xlab, ylab=ylab + ' Cp')
  summary.out <- summary(reg.fit_only_means)
  as.data.frame(summary.out$outmat)
  return(reg.fit_only_means) 
}


#################################
#       CROSS VALIDATION        #
#################################
crossValidate <- function(formula, data, k=5){
  # On divise les donnees
  # en k parties egales
  # Et pour chaque pli on fit un modele
  # sur tous les autres plis
  # On calcule ensuite le r^2 de ce modele
  n = nrow(data)
  folds = sample(1:k,n,replace=TRUE)
  CV = 0
  
  for(k in (1:k)){
    reg = lm(formula, data=data[folds!=k,])
    pred = predict(reg,newdata=data[folds==k,])
    CV = CV + sum((data$Time[folds==k] - pred)^2)
  }
  CV = CV/n
  
  return(CV);
}

#################################
#       SUBSET SELECTION        #
#################################
breast_train = breast_data[1:130,]
breast_test = breast_data[131:194,]

# On calcule les subsets sur les donnees de train
# qui representent 2/3 des donnees totales
wholedata.subsets = getSubsets(breast_train$Time~., breast_train, 'Subset Selection for 2/3rds of whole set', 'Predictors', 'Values')

## =========
## RESULTATS
# Cette analyse nous permet de trouver 4 modeles differents, qui correspondent a r^2, adjusted r^2, Cp et BIC

######
# Cette fonction permet de Cross-Validate
# tous les modeles trouves precedemment avec les donnees de train
batchcrossvalidate <- function(k=5){
  
  # r^2
  r_squared_CV = crossValidate(Time~., breast_train)
  # adj r^2
  adjr_cv = crossValidate(Time ~ radius_mean + perimeter_mean + smoothness_mean + symmetry_mean + 
                            fractal_dimension_mean + radius_se + texture_se + perimeter_se + area_se + 
                            smoothness_se + symmetry_se + smoothness_worst + Tumor_size, breast_train);
  # bic
  bic_cv = crossValidate(Time ~ radius_mean + perimeter_mean + fractal_dimension_mean, breast_train)
  
  # Cp
  Cp_CV = crossValidate(Time ~ radius_mean + perimeter_mean + fractal_dimension_mean + texture_se, breast_train)
  
  # PLot
  barplot(height=c(r_squared_CV, adjr_cv, bic_cv, Cp_CV), names.arg=c('r^2', 'adjusted r^2', 'BIC', 'Cp'))
  
  return (as.data.frame(c(r_squared_CV, adjr_cv, bic_cv, Cp_CV)));
}


## Cette fonction permet de lancer 10 fois
## L'analyse, afin d'obtenir un tableau
## recapitulatif de la cross-validation
## de chaque modele
plot_cross_validations <- function(k=5){
  r = rep(0:10)
  adjr = rep(0:10)
  bic = rep(0:10)
  cp = rep(0:10)
  
  for(i in 1:10){
    res = batchcrossvalidate(k);
    print(res)
    r[i] = res[1,1]
    adjr[i] = res[2, 1]
    bic = res[3, 1]
    cp = res[4, 1]
  }
  data = data.frame(r, adjr, bic, cp)
  names(data) = c('r^2', "adjusted r^2", "bic", "cp")
  return (data)
}

plot_cross_validations(5)




#################################
#        REGULARIZATION         #
#################################

## INIT
x = model.matrix(Time~., breast_data)
y <- as.data.frame(breast_data$Time)

n = nrow(breast_data)
napp = 130
ntst = n - 130

train = sample(1:n, napp)

xapp = x[train, ]
yapp = y[train, ]

xtst = x[-train, ]
ytst = y[-train, ]

# Ridge
cv.out<-cv.glmnet(xapp,yapp,alpha=0)
plot(cv.out)
fit<-glmnet(xapp,yapp,lambda=cv.out$lambda.min,alpha=0)
ridge.pred<-predict(fit,s=cv.out$lambda.min,newx=xtst)
# calcul de r^2
print(mean((ytst-ridge.pred)^2))


# Lasso
cv.out<-cv.glmnet(xapp,yapp,alpha=1)
plot(cv.out)
fit.lasso<-glmnet(xapp,yapp,lambda=cv.out$lambda.min,alpha=0)
lasso.pred<-predict(fit.lasso,s=cv.out$lambda.min,newx=xtst)
# calcul de r^2
print(mean((ytst-lasso.pred)^2))


#################################
# PRINCIPAL COMPONENT ANALYSIS  #
#################################
  ## Princomp
pca = princomp(breast_data, scores=TRUE, cor=TRUE)
summary(pca)
breast_12comp = breast_data[, 1:12]
lm.fit = lm(breast_data$Time~., data=breast_12comp)
lm.pred = predict(lm.fit)
# Calcul de r^2
print(mean((breast_data$Time-lm.pred)^2))

  ## PLS
pcr.fit = pcr(Time~.,data=breast_data, scale=TRUE, validation="CV")
summary(pcr.fit)
validationplot(pcr.fit,val.type = "MSEP",legendpos = "topright")




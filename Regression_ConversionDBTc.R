#linear regression model
# Full data file 
# ==============================================================================
library(tidyverse)
library(skimr)
library(DataExplorer)
library(scales)
library(corrr)
library(foreign)
library(ggplot2)
library(psych)
library(cleandata)
library(dlookr)


# Models
# ==============================================================================
library(glmnet)
library(pls)
library(caret)
library(caTools)
library(Metrics)
#library("glmnet")




#################### Dataset HDS_catalyst #############################

dataor <- read.csv("~/HDS_regression-main/HDS_regression-main/DataregreHDS_DBT_Source.csv", header = TRUE)

#string date
lab <- read.csv("~/HDS_regression-main/HDS_regression-main/ListtextHDS_cat.csv", header = FALSE)
lab <- lab[,-1]

################## Data processing ######################################
datos <- dataor[,6:37]

sel <-replace(datos$Selectivity, datos$Selectivity>3, 0)

datos$Selectivity <-replace(datos$Selectivity, datos$Selectivity>3, max(sel))
C_n <- colnames(datos)

#_________String data_______________
cat_n <- c("No_Metals","Metal_cat1_","Metal_cat2_","Element_1_support_","Element_2_support_",
           "Element_3_support_", "Nano","Mesostructure","Promoters","Impregnation_method","Aditive","Kind_cat","Support_1","Support_2",
           "P_Mo","P_Ni","P_Co","P_W","P_Si","P_Al","P_Ti","Structure_directing_agent")

#__________Numerical data___________
d_N <- c("Temperature","Surface_area","Pore_size","Pressure_H2","Selectivity",
         "Conversion_DBT","Reaction_time","Slab_length","Stacking_grade","Dispersion")

#Data conversion from string to nominals
o=0
for (i in 1:length(cat_n)) {
  o=o+1
  b= i 
  n <- lab[i,]
  a=list()
  d=0
  for (j in 1:8){
    if (n[j]!=""){
      d=d+1
      a[d]=n[j]
  }
  }
  for (k in (1:length(a))){
    f=(k-1)/(length(a)-1)
    datos[,cat_n[i]] <-replace(datos[,cat_n[i]], datos[,cat_n[i]]==a[[k]],round(f,9)) 
  }
}

#Standardization

for (i in 1:length(d_N)){
  Maxc <- max(datos[,d_N[i]])
  Minc  <- min(datos[,d_N[i]])
  A=datos[d_N[i]]
  print(A)
  for (j in 1:length(A) ){
    
    x=(A[j]-Minc)/(Maxc-Minc)
    A[j]=round(x,9)
    datos[,d_N[i]] = A[j]
  }
}

datapros <- datos
skim(datos)

write.csv(datapros, "~/HDS_regression-main/HDS_regression-main/datapros1.csv")
datosp <- read.csv("~/HDS_regression-main/HDS_regression-main/datapros1.csv", header = TRUE)
datosp <- datosp[,-1]

datos <-datosp

# Correlation between variables 
# ==============================================================================
df_correlaciones <- datos %>%
  correlate(method = "pearson") #%>%

#Divide data in trainign and test
# ==============================================================================
set.seed(1235)
id_train <- sample(1:nrow(datos), size = 0.7*nrow(datos), replace = FALSE)

datos_train <- datos[id_train, ]
datos_test  <- datos[-id_train, ]


###########################  RIDGE REGRESSION ##################

# Training and test matrixes
# ==============================================================================

x_train <- model.matrix(Conversion_DBT~., data = datos_train)[, -1]

y_train <- datos_train$Conversion_DBT

x_test <- model.matrix(Conversion_DBT~., data = datos_test)[, -1]

y_test <- datos_test$Conversion_DBT


# Regularization method
model <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  nlambda     = 100,
  standardize = TRUE
)


# Evolution of the lambda  coeffients
# ==============================================================================
regularization <- model$beta %>% 
  as.matrix() %>%
  t() %>% 
  as_tibble() %>%
  mutate(lambda = model$lambda)

regularization <- regularization %>%
  pivot_longer(
    cols = !lambda, 
    names_to = "predictor",
    values_to = "coefficients"
  )

regularization %>%
  ggplot(aes(x = lambda, y = coefficients, color = predictor)) +
  geom_line() +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  labs(title = "Coefficients of regression model") +
  theme_bw() +
  theme(legend.position = "none")


# Evolution of lambda error
# ==============================================================================
set.seed(123)
cv_error <- cv.glmnet(
  x      = x_train,
  y      = y_train,
  alpha  = 0,
  nfolds = 10,
  type.measure = "mse",
  standardize  = TRUE
)

plot(cv_error)

# Improving value of lambda
# ==============================================================================
paste("The best value of lambda:", cv_error$lambda.min)

# ==============================================================================
# The bigger lambda value
paste("Mejor valor de lambda encontrado + 1 desviación estándar:", cv_error$lambda.1se)


# ==============================================================================
model <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 0,
  lambda      = cv_error$lambda.1se,
  standardize = TRUE
)

# Coefficients of model
# ==============================================================================
df_coeficients <- coef(model) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor") %>%
  rename(coeficiente = s0)

df_coeficients %>%
  filter(predictor != "(Intercept)") %>%
  ggplot(aes(x = predictor, y = coeficiente)) +
  geom_col() +
  labs(title = "Coefficients for Ridge model") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90))


# Training predictions
# ==============================================================================
predict_train <- predict(model, newx = x_train)

# Training MSE 
# ==============================================================================
training_mse <- mean((predict_train - y_train)^2)
paste("Training error:", training_mse)


# predict de test
# ==============================================================================
predict_test <- predict(model, newx = x_test)

# test MSE
# ==============================================================================
test_mse_ridge <- mean((predict_test - y_test)^2)
paste("Test error:", test_mse_ridge)




################LASSO REGRESSION########################################


# Training and Test matrixes
# ==============================================================================

x_train <- model.matrix(Conversion_DBT~., data = datos_train)[, -1]
y_train <- datos_train$Conversion_DBT

x_test <- model.matrix(Conversion_DBT ~., data = datos_test)[, -1]
y_test <- datos_test$Conversion_DBT



library("glmnet")
# 
model <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 1,
  nlambda     = 100 ,
  standardize = TRUE
)

#  Evolution of lambda error
# ==============================================================================
regularization <- model$beta %>% 
  as.matrix() %>%
  t() %>% 
  as_tibble() %>%
  mutate(lambda = model$lambda)

regularization <- regularization %>%
  pivot_longer(
    cols = !lambda, 
    names_to = "predictor",
    values_to = "coefficients"
  )

regularization %>%
  ggplot(aes(x = lambda, y = coefficients, color = predictor)) +
  geom_line() +
  scale_x_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  labs(title = "coefficients of model with regularization") +
  theme_bw() +
  theme(legend.position = "none")

# Evolution of lambda error
# ==============================================================================
set.seed(123)
cv_error <- cv.glmnet(
  x      = x_train,
  y      = y_train,
  alpha  = 1,
  nfolds = 10,
  type.measure = "mse" ,
  standardize  = TRUE
)

plot(cv_error)

# The best value found
# ==============================================================================
paste("The best value of lambda:", cv_error$lambda.min)



# ==============================================================================
paste("Mejor valor de lambda encontrado + 1 desviación estándar:", cv_error$lambda.1se)

# Optimum
# ==============================================================================
model <- glmnet(
  x           = x_train,
  y           = y_train,
  alpha       = 1,
  lambda      = cv_error$lambda.1se ,
  standardize = TRUE
)

# Coefficients of model
# ==============================================================================
df_coeficients <- coef(model) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor") %>%
  rename(coeficiente = s0)

df_coeficients %>%
  filter(predictor != "(Intercept)") %>%
  ggplot(aes(x = predictor, y = coeficiente)) +
  geom_col() +
  labs(title = "Coeffients of Lasso model") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90))


df_coeficients %>%
  filter(
    predictor != "(Intercept)",
    coeficiente != 0
  ) 

# Training predict 
# ==============================================================================
predict_train <- predict(model, newx = x_train)

# Training MSE 
# ==============================================================================
training_mse <- mean((predict_train - y_train)^2)
paste("Training Error (mse):", training_mse)

# predict of test
# ==============================================================================
predict_test <- predict(model, newx = x_test)

# Test MSE 
# ==============================================================================
test_mse_lasso <- mean((predict_test - y_test)^2)
paste("Test Error (mse):", test_mse_lasso)
test_mse_lasso

print(test_mse_lasso)


#Ridge and Lasso cross-validation 

set.seed(123)   
library(glmnet) 
library(dplyr)  
library(psych)  

datosp <- read.csv("~/HDS_regression-main/HDS_regression-main/datapros1.csv", header = TRUE)
datosp <- datosp[,-1]


data <-datosp

colnames(data)
datos <- data


y <- datos %>% select(Conversion_DBT) %>%as.matrix()
X <- datos %>% select(-Conversion_DBT) %>% as.matrix()

# Perform 10-fold cross-validation to select lambda ---------------------------
lambdas_to_try <- 10^seq(-3, 5, length.out = 100)

# Setting alpha = 0 implements ridge regression
ridge_cv <- cv.glmnet(X, y, alpha = 0, lambda = lambdas_to_try,
                     standardize = TRUE, nfolds = 10)

# Best cross-validated lambda
lambda_cv <- ridge_cv$lambda.min
# Fit final model, get its sum of squared residuals and multiple R-squared
model_cv <- glmnet(X, y, alpha = 0, lambda = lambda_cv , standardize = TRUE)
y_hat_cv <- predict(model_cv, X)
ssr_cv <- t(y - y_hat_cv) %*% (y - y_hat_cv)
rsq_ridge_cv <- cor(y, y_hat_cv)^2

pred_Ridopt<-predict(model_cv,x_test)
Y_out <- data.frame("Y_ori"= y, "Y_pred"=y_hat_cv)

y1_test <- y_test

Ridlm <- lm(Conversion_DBT~s0,data = Y_out)
dev.new(width=8, height=8, unit="cm")
plot(y,y_hat_cv,col="green", cex.lab=1, cex.axis=1.2, cex.main=1.1, cex.sub=1,main="Ridge regression",xlab="Observed Conversion of DBT",ylab = "Predicted Conversion of DBT",cex=0.8,xlim=c(0,1.1),ylim=c(0,1.14))
points(y1_test,pred_Ridopt, col="blue",cex.lab=1, cex.axis=1, cex.main=1.1, cex.sub=1,cex=1,add=TRUE)
abline(Ridlm,col="black")
print(Ridlm)

Ridtest <- cbind(y1_test,pred_Ridopt)
Ridally <- cbind(y,y_hat_cv)

###############################################################
#Plot Coefficient optimal Ridge 
###############################################################
R_df_coef= coef(model_cv)
print(R_df_coef)

R_df_coef <- coef(model_cv) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor") %>%
  rename(coeficiente = s0)

R_df_coef %>%
  filter(predictor != "(Intercept)") %>%
  ggplot(aes(x = predictor, y = coeficiente)) +
  geom_col() +
  labs(title = "Coefficients for Ridge model") +
  xlab("predictors")+
  ylab("coefficient")+  
  theme_bw() +
  
  theme(axis.text.x = element_text(family="Times",size = 12, angle = 90,color="black"))+
  
  theme(axis.text.y = element_text(family="Times",size = 12,color="black"))

# Use information criteria to select lambda -----------------------------------
#X_scaled <- scale(X)
X_scaled <- X
aic <- c()
bic <- c()
for (lambda in seq(lambdas_to_try)) {
  model <- glmnet(X, y, alpha = 0, lambda = lambdas_to_try[lambda], standardize = TRUE)
  betas <- as.vector((as.matrix(coef(model))[-1, ]))
  resid <- y - (X_scaled %*% betas)
  
  ld <- lambdas_to_try[lambda] * diag(ncol(X_scaled))
  H <- X_scaled %*% solve(t(X_scaled) %*% X_scaled + ld) %*% t(X_scaled)
  df <- tr(H)
  aic[lambda] <- nrow(X_scaled) * log(t(resid) %*% resid) + 2 * df
  bic[lambda] <- nrow(X_scaled) * log(t(resid) %*% resid) + 2 * df * log(nrow(X_scaled))
}


# Optimal lambdas according to both criteria
lambda_aic <- lambdas_to_try[which.min(aic)]
lambda_bic <- lambdas_to_try[which.min(bic)]

# Fit final models, get their sum of squared residuals and multiple R-squared
model_aic <- glmnet(X, y, alpha = 0, lambda = lambda_aic, standardize = TRUE)
y_hat_aic <- predict(model_aic, X)
ssr_aic <- t(y - y_hat_aic) %*% (y - y_hat_aic)
rsq_ridge_aic <- cor(y, y_hat_aic)^2
rmse_ridge <- rmse(y, y_hat_aic)
mae_ridge <- mae(y, y_hat_aic)

model_bic <- glmnet(X, y, alpha = 0, lambda = lambda_bic, standardize = TRUE)
y_hat_bic <- predict(model_bic, X)
ssr_bic <- t(y - y_hat_bic) %*% (y - y_hat_bic)
rsq_ridge_bic <- cor(y, y_hat_bic)^2



# Calculate the weights from univariate regressions
weights <- sapply(seq(ncol(X)), function(predictor) {
  uni_model <- lm(y ~ X[, predictor])
  coeff_variance <- summary(uni_model)$coefficients[2, 2]^2
})


# Heteroskedastic Ridge Regression loss function - to be minimized
hridge_loss <- function(betas) {
  sum((y - X %*% betas)^2) + lambda * sum(weights * betas^2)
}


# Heteroskedastic Ridge Regression function
hridge <- function(y, X, lambda, weights) {
  model_init <- glmnet(X, y, alpha = 0, lambda = lambda, standardize = TRUE)
  betas_init <- as.vector(model_init$beta)
  coef <- optim(betas_init, hridge_loss)$par
  fitted <- X %*% coef
  rsq <- cor(y, fitted)^2
  names(coef) <- colnames(X)
  output <- list("coef" = coef,
                 "fitted" = fitted,
                 "rsq" = rsq)
  return(output)
}


# Fit model to the data for lambda = 0.001
hridge_model <- hridge(y, X, lambda = 0.001, weights = weights)
rsq_hridge_0001 <- hridge_model$rsq


# Perform 10-fold cross-validation to select lambda ---------------------------
lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
# Setting alpha = 1 implements lasso regression
lasso_cv <- cv.glmnet(X, y, alpha = 1, lambda = lambdas_to_try,
                      standardize = TRUE, nfolds = 10)


# Best cross-validated lambda
lambda_cv <- lasso_cv$lambda.min
# Fit final model, get its sum of squared residuals and multiple R-squared
model_cv <- glmnet(X, y, alpha = 1, lambda = lambda_cv, standardize = TRUE)
y_hat_cv <- predict(model_cv, X)
ssr_cv <- t(y - y_hat_cv) %*% (y - y_hat_cv)
rsq_lasso_cv <- cor(y, y_hat_cv)^2
rmse_lasso <- rmse(y, y_hat_cv)
mae_lasso <- mae(y, y_hat_cv)


pred_Lasopt<-predict(model_cv,x_test)
Y_out <- data.frame("Y_ori"= y, "Y_pred"=y_hat_cv)

y1_test <- y_test

Ridlm <- lm(Conversion_DBT~s0,data = Y_out)
plot(y,y_hat_cv,col="green", main="Lasso regression",cex.lab=1, cex.axis=1.2, cex.main=1.1, cex.sub=1,xlab="Observed Conversion of DBT",ylab = "Predicted Conversion of DBT",,cex=0.8,xlim=c(0,1.1),ylim=c(0,1.14))
points(y1_test,pred_Lasopt, col="blue",cex.lab=1, cex.axis=1.2, cex.main=1.1, cex.sub=1,cex=1,add=TRUE)
abline(Ridlm,col="black")
print(Ridlm)

Lastest <- cbind(y1_test,pred_Lasopt)

L_df_coef= coef(model_cv)
print(L_df_coef)

L_df_coef <- coef(model_cv) %>%
  as.matrix() %>%
  as_tibble(rownames = "predictor") %>%
  rename(coeficiente = s0)

L_df_coef %>%
  filter(predictor != "(Intercept)") %>%
  ggplot(aes(x = predictor, y = coeficiente)) +
  geom_col() +
  #labs(title = "Coefficients for Lasso model") +
  xlab("predictors")+
  ylab("coefficient")+
  #axis.title.x = element_text(color="blue", size=14, face="bold")+
  theme_bw() +
  theme(axis.text.x = element_text(family="Times",size = 12, angle = 90,color="black"))+
  
  theme(axis.text.y = element_text(family="Times",size = 12,color="black"))
L_df_coef %>%
  filter(
    predictor != "(Intercept)",
    coeficiente != 0
  ) 

rsq <- cbind("R-squared" = c(rsq_ridge_cv, rsq_lasso_cv))
rownames(rsq) <- c("ridge cross-validated", "lasso cross_validated")
print(rsq)


#####################################Random_Forest###############################
library(randomForest)
set.seed(123)

m1 <- randomForest(formula = Conversion_DBT~., data = datos_train)

m1

plot(m1)

oob.err=double(29)
test.err=double(29)

##mtry is no of Variables randomly chosen at each split
pmtry <-(1:29)
for(mtry in 1:29) 
{
   
  rf=randomForest(formula = Conversion_DBT ~ . , data = datos_train , mtry=mtry) 
  oob.err[mtry] = rf$mse[10] 
  
  pred<-predict(rf,datos_test) #Predictions on Test Set for each Tree

  test.err[mtry]= with(datos_test, mean( (y_test - pred)^2)) #Mean Squared Test Error
  
  cat(mtry," ") #printing the output to the console
}


matplot(1:mtry , cbind(oob.err,test.err), pch=19 , col=c("red","blue"),type="b",ylab="Mean Squared Error",xlab="Number of Predictors Considered at each Split")
#legend("topright",legend=c("Out of Bag Error","Test Error"),pch=19, col=c("red","blue"))

oobt=double(50)
testt=double(50)
pntree <- (1:50)
for (ntree in 1:50)

{
  rf=randomForest(formula = Conversion_DBT ~ . , data = datos_train , mtry=25, ntree=ntree) 

  pred<-predict(rf,datos_test) #Predictions on Test Set for each Tree
  
  testt[ntree]= with(datos_test, mean( (y_test - pred)^2)) #Mean Squared Test Error
  
  cat(ntree," ") #printing the output to the console
}

print(testt)

for ( j in 1:length(testt)){
  if (testt[j] == min(testt))
  {ntree_opt=pntree[j]}
}

print(ntree_opt)



rf_opt=randomForest(formula = Conversion_DBT ~ . , data = datos_train , mtry=25, ntree=ntree_opt)
pred_rfopt<-predict(rf_opt,datos_test)
pred_allrf <-predict(rf_opt,datos)
Y_out <- data.frame("Y_ori"= datos$Conversion_DBT, "Y_pred"=pred_allrf)
rf_r2 <- cor(datos_test$Conversion_DBT, pred_rfopt)^2

rlm <- lm(Y_ori~Y_pred,data = Y_out)

plot(datos$Conversion_DBT,pred_allrf,col="green", cex.lab=1, cex.axis=1.2, cex.main=1.1, cex.sub=1,,cex=1, main="Random Forest regression", ylim=c(0,1), xlim=c(0,1),xlab="Observed Conversion of DBT",ylab = "Predicted Conversion of DBT")
points(datos_test$Conversion_DBT,pred_rfopt, col="blue",cex.lab=1, cex.axis=1.2, cex.main=1.1, cex.sub=1,cex=1, add=TRUE)
abline(rlm,col="black",add=TRUE)

print(rlm)

RFtesty <- cbind(datos_test$Conversion_DBT,pred_rfopt)
RFtestall <- cbind(datos$Conversion_DBT,pred_allrf)

importancia_pred <- rf_opt$importance %>%
  enframe(name = "predictor", value = "importance")

print(importancia_pred)
# plot
ggplot(
  data = importancia_pred,
  aes(x    = reorder(predictor, importance),
      y    = importance,
      fill = importance)
) +
  labs(x = "predictors", title = "Importance predictors") +
  geom_col() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none")+
  theme_bw() +
  theme(axis.text.x = element_text(family="Times",size = 12,color="black"))+
  #theme(axis.text.y = element_text(size = 14, angle = 90))
  #p +
  theme(axis.text.y = element_text(family="Times",size = 12,color="black"))
  

set.seed(1234)



folds <- createFolds(datos$Conversion_DBT, k = 10)

cvRandomForest <- lapply(folds, function(x){
  training_fold <- datos[-x, ]
  test_fold <- datos[x, ]
  RF <- randomForest(formula = Conversion_DBT ~ ., data = training_fold, mtry=25, ntree= ntree_opt)
  
  y_pred <- predict(RF, newdata = test_fold)
  rf_r2 <- cor(test_fold$Conversion_DBT, y_pred)^2
  rf_rmse <- rmse(test_fold$Conversion_DBT, y_pred)
  rf_mae <- mae(test_fold$Conversion_DBT, y_pred)
  
  
  M_RF= list(rf_r2, rf_rmse,rf_mae)
  
  return(M_RF)
  
})


print(cvRandomForest)

RF_rsqall=list()
RF_rmseall = list()
RF_maeall = list()


for (i in 1:10) {
  
  R2 <- as.numeric(cvRandomForest[[i]][1])
  m_rmse <- as.numeric(cvRandomForest[[i]][2])
  M_ae <- as.numeric(cvRandomForest[[i]][3])
  RF_rsqall = c(RF_rsqall,R2)
  RF_rmseall = c(RF_rmseall,m_rmse)
  RF_maeall = c(RF_maeall,M_ae)
  }
RF_rmse=as.numeric(RF_rmseall)
RF_rsq=as.numeric(RF_rsqall)
RF_rsq_mean=mean(RF_rsq)
RF_mae = as.numeric(RF_maeall)
RF_mae_mean = mean(RF_mae)
RF_rmse_mean=mean(RF_rmse)

######################Results####################################

rsq2 <- cbind("R-squared" = c(rsq_ridge_cv, rsq_lasso_cv, RF_rsq_mean),"RMSE"= c(rmse_ridge, rmse_lasso,RF_rmse_mean ),"MAE"=c(mae_ridge,mae_lasso,RF_mae_mean))
rownames(rsq2) <- c("ridge cross-validated", "lasso cross_validated","Random_forest cross_validate")
print(rsq2)



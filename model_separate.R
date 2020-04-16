###################################################################################
#                Separate Model
###################################################################################
library(numDeriv)
library('MASS')
library(nnet)
library(leaps)
library(EvaluationMeasures)
library(glmnet)
library(BeSS)
###################################################################################
#                1 Simulate Data Set
###################################################################################
## Function
fn_1_simulation<-function(par_seed,par_betanum,par_rele_num,par_r,par_samplesize){
  set.seed(par_seed)
  betanum <- par_betanum              ##total number of coefficient
  relevant.b <- par_rele_num          ##total number of relevent coefficient
  unrelavent.b <- betanum-relevant.b  ##total number of unrelevent coefficient
  n <- par_samplesize                 ##sample size
  
  assign("betanum", betanum, envir = .GlobalEnv)
  assign("n", n, envir = .GlobalEnv)
  
  ####X: Relevant Matirx#### 
  fn_generate_x <- function(par_r) {
    r <- par_r
    mm<-diag(1, nrow=betanum) 
    msigma<-matrix(0,nrow = betanum,ncol = betanum)
    for (i in 1:betanum) {
      for (j in 1:betanum) {
        msigma[i,j]<-ifelse(mm[i,j]==1,1,r)
      }
    }
    mean0<-rep(0, times = betanum)
    X<-mvrnorm(n, mu=mean0, Sigma=msigma ,empirical=TRUE)
    ####Return Value####
    assign("X", X, envir = .GlobalEnv)  
  }
  fn_generate_x(par_r) 
  
  ####Y: regression & classification####
  fn_generate_y<-function(par_b0.r, par_betarmin, par_betarmax,par_b0.c, par_betacmin, par_betacmax){
    
    ####Linear regression####
    kr <- 1  ##linear task
    e <- matrix(rnorm(1),nrow = n, ncol = kr) ##Error: e~N(0,1)
    b0.r <- par_b0.r  ##Intercept
    B.r <- matrix(0 ,nrow = kr, ncol = betanum)  ##Coefficient
    yr <- matrix(0 ,nrow = n, ncol = kr)  ##Respond
    for (i in 1:kr) {
      al<-matrix(runif(relevant.b, min = par_betarmin, max = par_betarmax),nrow = 1,ncol = relevant.b)
      bl<-matrix(0,nrow = 1,ncol = unrelavent.b)
      B.r[i,]<-cbind(al,bl)
      yr[,i] <- b0.r + (X%*%B.r[i,])
    }
    
    ####Logistics regression####
    b0.c <- par_b0.c
    a <- matrix(runif(relevant.b, min = par_betacmin, max = par_betacmax),nrow = 1,ncol = relevant.b)
    b <- matrix(0,nrow = 1,ncol = unrelavent.b)
    B.c <- cbind(a,b)                                           
    pi.1 <- (exp(b0.c+X%*%t(B.c))/(1+exp(b0.c+X%*%t(B.c)))) ##P(y=1|X=x)
    pi.2 <- 1/(1+exp(b0.c+X%*%t(B.c)))                        ##P(y=0|X=x)
    
    ####Response Y####  
    ys.r <- yr+e
    #ys.c <- rbinom(n,1,prob=pi.1)  
    ys.c <- ifelse(runif(n) < pi.1, 1,0) ## prob < pi.1 then y=1
    ####Return Value####
    assign("b0.r", b0.r, envir = .GlobalEnv)  
    assign("B.r", B.r, envir = .GlobalEnv)  
    assign("b0.c", b0.c, envir = .GlobalEnv)  
    assign("B.c", B.c, envir = .GlobalEnv)  
    assign("ys.r", ys.r, envir = .GlobalEnv)  
    assign("ys.c", ys.c, envir = .GlobalEnv)
  }
  fn_generate_y(1,5,10,2,5,10) #(par_b0.r, par_betarmin, par_betarmax,par_b0.c, par_betacmin, par_betacmax)
  
  #### DATA ####
  dd<-cbind.data.frame(X,ys.r,ys.c)
  colnames(dd) <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10",
                    "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
                    "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30",
                    "X31","X32","X33","X34","X35","X36","X37","X38","X39","X40",
                    "X41","X42","X43","X44","X45","X46","X47","X48","X49","X50",
                    "ys.r","ys.c")
  return(dd)
}

## Generate
data.sim<-fn_1_simulation(2020,50,10,0.6,2000) 
data.vv <- data.sim[501:1000,]     ##validation data
data.test <- data.sim[1001:2000,]  ##test data
data.ori <- data.sim[1:500,]       ##train data
###################################################################################
#                2 Build Model
###################################################################################
## Function
fn_find_separate_lamda<-function(traindata,testdata,par_lamda_size,par_lamdamin,par_lamdamax) { 
  
  #### Basic Setting 
  par_lamda <- runif(par_lamda_size, min = par_lamdamin, max =par_lamdamax)
  Lossfn_r <- matrix(0,nrow = par_lamda_size,ncol=1) 
  Lossfn_c <- matrix(0,nrow = par_lamda_size,ncol=1) 
  betar <- matrix(0,nrow = par_lamda_size,ncol=betanum+1) 
  betac <- matrix(0,nrow = par_lamda_size,ncol=betanum+1) 
  train_x <- as.matrix(traindata[,1:betanum])   ##assign input data
  train_yr <- as.matrix(traindata[,betanum+1])  ##assign output yr
  train_yc <- as.matrix(traindata[,betanum+2])  ##assign output yc
  
  for (k in 1:par_lamda_size) {
    
    #### Models
    reg_l2<-glmnet(train_x, train_yr, family = "gaussian", alpha = 1,lambda = par_lamda[k])####lasso regression
    logi_l2<-glmnet(train_x, train_yc, family = "binomial", alpha = 1,lambda = par_lamda[k])####logistics regression with l2
    
    #### Calculate Loss: Lr,Lc,Pq 
    x <- as.matrix(testdata[,1:betanum]) 
    y_r <- as.matrix(testdata[,betanum+1])
    y_c <- as.matrix(testdata[,betanum+2])
    r_b0 <- coef(reg_l2)[1]
    r_beta <- coef(reg_l2)[-1]
    c_b0 <- coef(logi_l2)[1]
    c_beta <- coef(logi_l2)[-1]
    penalty_r <- matrix(0,nrow = betanum, ncol=1)
    penalty_c <- matrix(0,nrow = betanum, ncol=1)
    
    Lr <- t(y_r-(r_b0 + x%*%r_beta)) %*% (y_r-(r_b0 + x%*%r_beta))
    Lc <- sum(-((y_c*(c_b0 + x%*%c_beta)) - log(1+exp(c_b0 + x%*%c_beta))))
    for (w in 1:betanum) {penalty_r[w] <- sqrt(r_beta[w]^2)} ###for regression
    for (w in 1:betanum) {penalty_c[w] <- sqrt((r_beta[w]^2)+(c_beta[w]^2))}
    Pg_r <- par_lamda[k]*(sum(penalty_r))
    Pg_c <- par_lamda[k]*(sum(penalty_c))
    
    Lossfn_r[k] <- Lr + Pg_r
    Lossfn_c[k] <- Lc + Pg_c
    
    for (p in 1:betanum+1) { ##p includes bo-b50
      betar[k,p] <- coef(reg_l2)[p]
      betac[k,p] <- coef(logi_l2)[p]
    }
    
  }
  
  #### Results: lamda and its beta (out_reg, out_logi)
  out_reg=cbind(par_lamda,Lossfn_r,betar)
  out_logi=cbind(par_lamda,Lossfn_c,betac)
  assign("out_reg", out_reg, envir = .GlobalEnv)
  assign("out_logi", out_logi, envir = .GlobalEnv)
}

## Generate
result_index_only <- matrix(0,ncol = 9,nrow = 50)
result_esti_beta_only <- matrix(0,ncol =betanum*2+2,nrow = 50)


for (s in 1:50) {
  
  #### Find best tuning parameter  L1/L2
  fn_find_separate_lamda(data.ori,data.vv,1000,0,5)
  
  best_lamda_r <- out_reg[which.min(out_reg[,2]),1] 
  best_lamda_c <- out_logi[which.min(out_logi[,2]),1] 
  
  #### Estimate Coefficients 
  traindata <- data.ori
  train_x <- as.matrix(traindata[,1:betanum])   ##assign input data
  train_yr <- as.matrix(traindata[,betanum+1])  ##assign output yr
  train_yc <- as.matrix(traindata[,betanum+2])  ##assign output yc
  
  onlyreg_l2<-glmnet(train_x, train_yr, family = "gaussian", alpha = 1,lambda = best_lamda_r)
  onlylogi_l2<-glmnet(train_x, train_yc, family = "binomial", alpha = 1,lambda = best_lamda_c)
  
  ###################################################################################
  #                3 Generate Results: Estimate Err & Prediction Err & FNR & FPR
  ###################################################################################
  #### Function
  predict_err <- function(testdata){
    
    #### Basic setting
    testn <- dim(testdata)[1]
    X <- as.matrix(testdata[,1:betanum])
    
    #### Linear regression####
    b0.r_test <- coef(onlyreg_l2)[1]  
    B.r_test <- as.matrix(coef(onlyreg_l2)[-1])
    yr_hat<- b0.r_test + (X%*%B.r_test)
    
    #### Logistics regression####
    b0.c_test<- coef(onlylogi_l2)[1]
    B.c_test <- as.matrix(coef(onlylogi_l2)[-1])
    pi_hat <- (exp(b0.c_test+X%*%B.c_test)/(1+exp(b0.c_test+X%*%B.c_test))) ##P(y=1|X=x)
    yc_hat <- ifelse(runif(testn) < pi_hat, 1,0) ## prob < pi.1 then y=1
    
    err_predict_onlyr <- sum((testdata[,betanum+1]-yr_hat)^2)
    err_predict_onlyc <- mean(testdata[,betanum+2] != yc_hat)
    
    ####Results: estimated y 
    data.est.sep <- cbind.data.frame(X,yr_hat,yc_hat)
    assign("err_predict_onlyr", err_predict_onlyr, envir = .GlobalEnv)
    assign("err_predict_onlyc", err_predict_onlyc, envir = .GlobalEnv)
    assign("data.est.sep", data.est.sep, envir = .GlobalEnv)
  }

  #### Generate: Estimation Error (true beta - estimated beta)
  count1 <- betanum+1
  err_esti_r <- ((b0.r - coef(onlyreg_l2)[1])^2) + sum((B.r - coef(onlyreg_l2)[2:count1])^2) 
  err_esti_c <- ((b0.c - coef(onlylogi_l2)[1])^2) + sum((B.c - coef(onlylogi_l2)[2:count1])^2) 
  err_esti <- err_esti_r + err_esti_c 
  #### Generate: Prediction Error
  predict_err(data.test)
  #### Generate: Variable Selection
  var_selection_onlyr<- regsubsets(yr_hat ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20
                                   +X21+X22+X23+X24+X25+X26+X27+X28+X29+X30+X31+X32+X33+X34+X35+X36+X37+X38+X39+X40
                                   +X41+X42+X43+X44+X45+X46+X47+X48+X49+X50 ,data = data.est.sep, nvmax = 10,really.big=T)
  xx <- data.est.sep[,1:50]
  var_selection_onlyc<- bess.one(xx,data.est.sep$yc_hat ,s=10,family ="binomial")
  
  test_subset_r<-matrix((summary(var_selection_onlyr, matrix.logical = T))$which ,nrow=10,ncol=betanum+1)[10,-1]
  test_subset_c<-matrix(FALSE, nrow=betanum,ncol = 1)
  test_subset_c[which(var_selection_onlyc$beta!=0),1]<-TRUE
  real_subset <- matrix(c(rep(TRUE,10),rep(FALSE,40)))
  
  Rate_FNR_reg<-EvaluationMeasures.FNR(real_subset,test_subset_r)
  Rate_FPR_reg<-EvaluationMeasures.FPR(real_subset,test_subset_r)
  Rate_FNR_logi<-EvaluationMeasures.FNR(real_subset,test_subset_c)
  Rate_FPR_logi<-EvaluationMeasures.FPR(real_subset,test_subset_c)
  
  #### Results: 
  result_index_only[s,] <- cbind(err_esti,err_esti_r,err_esti_c,err_predict_onlyr,err_predict_onlyc,Rate_FNR_reg,Rate_FPR_reg,Rate_FNR_logi,Rate_FPR_logi)
  colnames(result_index_only) <- c("Total Estimate Error","Regression Estimate Error","Logistics Estimate Error",
                              "Regression Prediction Error","Logistic Prediction Error",
                              "Regression FNR","Regression FPR","Logistics FNR","Logistics FPR")
  
  a <- as.matrix(t(coef(onlyreg_l2)))
  b <- as.matrix(t(coef(onlylogi_l2)))
  result_esti_beta_only[s,] <- cbind(a,b)
  print(paste0("s: ", s))
  print(paste0("index: ", result_index_only[s,]))
}


##############################################
colSums(result_index_only)/50

#result_esti_beta_only    #### Estimated beta
#data.est.sep             #### X, yr_hat, yc_hat

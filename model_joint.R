###################################################################################
#                Joint Model
###################################################################################
library(numDeriv)
library('MASS')
library(nnet)
library(leaps)
library(EvaluationMeasures)
library(glmnet)
library(BeSS)
###################################################################################
#                2 Build Model
###################################################################################
## Function
fn_2_find_lamda<-function(traindata,testdata,par_lamda_size,par_lamdamin,par_lamdamax) { 
  
  #### Basic Setting 
  par_lamda <- runif(par_lamda_size, min = par_lamdamin, max =par_lamdamax)
  Lossfn <- matrix(0,nrow = par_lamda_size,ncol=1) 
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
    penalty <- matrix(0,nrow = betanum, ncol=1)
    
    Lr <- t(y_r-(r_b0 + x%*%r_beta)) %*% (y_r-(r_b0 + x%*%r_beta))
    Lc <- sum(-((c_b0 + x%*%c_beta) - log(1+exp(c_b0 + x%*%c_beta))))
    for (w in 1:betanum) {penalty[w] <- sqrt((r_beta[w]^2)+(c_beta[w]^2))}
    Pg <- par_lamda[k]*(sum(penalty))
    
    Lossfn[k] <- Lr + Lc + Pg
    for (p in 1:betanum+1) { ##p includes bo-b50
      betar[k,p] <- coef(reg_l2)[p]
      betac[k,p] <- coef(logi_l2)[p]
    }
    
  }
  
  #### Results: lamda and its beta (out)
  out=cbind(par_lamda,Lossfn,betar,betac)
  assign("out", out, envir = .GlobalEnv)
}
find_Lbarrier <- function(theta) {
  
  #### Basic setting
  count1 <- betanum+1  
  count2 <- betanum+2
  count3 <- betanum+3
  count4 <- betanum*2+2
  r_b0 <- theta[1]
  r_beta <- theta[2:count1]
  c_b0 <- theta[count2]
  c_beta <- theta[count3:count4]
  opti_uj<-matrix(0,nrow=betanum,ncol=1)
  
  for (k in 1:betanum) {opti_uj[k]<-(max(abs(theta[k+1]),abs(theta[k+count2])))+0.1} ## For uj
  uj <- opti_uj
  Lr <- t(y_r-(r_b0 + x%*%r_beta)) %*% (y_r-(r_b0 + x%*%r_beta))  
  
  logterm <- matrix(0,ncol=1,nrow = dim(x)[1])  ## Solving the INF problem.
  for (m in 1:dim(x)[1]) {
    index <- c_b0 + x[m,]%*%c_beta
    if (index > log(10000)) {logterm[m] <- c_b0 + x[m,]%*%c_beta} 
    else if (index < log(0.0001)){logterm[m] <- log(1)} 
    else {logterm[m] <- log(1+exp(c_b0 + x[m,]%*%c_beta))}
  }
  
  #### Calculate LBarrier 
  Lc <- sum(-((c_b0 + x%*%c_beta) - logterm))
  Pg <- (lamda*sum(uj)) + (-(sum(log(r_beta+uj)+log(uj-r_beta)))*(1/tt)) + (-(sum(log(c_beta+uj)+log(uj-c_beta)))*(1/tt))
  return (Lr + Lc + Pg)
}
fn_3_optim<-function(traindata,setlamda, par_t, par_e, par_mu){
  
  #### Basic setting: randomly generate initial theta
  traindata<-as.matrix(traindata)
  assign("lamda", setlamda, envir = .GlobalEnv)
  assign("x", traindata[,1:betanum], envir = .GlobalEnv)
  assign("y_r", traindata[,betanum+1] , envir = .GlobalEnv)
  assign("y_c", traindata[,betanum+2] , envir = .GlobalEnv)
  thetanum <- betanum*3+2
  totalbetanum<-betanum*2+2
  assign("totalbetanum", totalbetanum, envir = .GlobalEnv)
  assign("thetanum", thetanum, envir = .GlobalEnv)
  theta<-matrix(0,ncol = thetanum,nrow=1)
  opti_r_b0<-runif(1,min=0,max=3)   ##theta
  opti_c_b0<-runif(1,min=0,max=3)   ##theta
  opti_r_beta<-as.matrix(runif(betanum,min=range_min,max=range_max),nrow=betanum,ncol=1)   ##theta
  opti_c_beta<-as.matrix(runif(betanum,min=range_min,max=range_max),nrow=betanum,ncol=1)   ##theta
  opti_uj<-matrix(0,nrow=betanum,ncol=1) ##strictly feasible #betanum+1
  
  theta<-c(opti_r_b0,opti_r_beta,opti_c_b0,opti_c_beta)                                           
  for (k in 1:betanum) {opti_uj[k]<-(max(abs(theta[k+1]),abs(theta[k+betanum+2])))+0.1}
  
  #### Basic setting: optim parameter
  tt <- par_t
  m<-2*betanum*2 ##2*betanum*(regression task number+ logistic task number)
  j<-(m/tt)
  
  #### Optim: step1-4 until converge
  repeat{ 
    ####################STEP 1 compute theta by Newton Method########################
    assign("tt", tt, envir = .GlobalEnv) ##"assign" means value tt can be use in function

    step1_opti<-optim(theta, find_Lbarrier, method = "CG",hessian = FALSE, control =list(type=2))#, control = list(abstol=0.1,type=1)
    print(paste0("convergence: ", step1_opti$convergence))
    print(paste0("convergence: ", step1_opti$value))
    
    ####################STEP 2 update theta###########################################
    #print(paste0("old theta: ", theta))
    #print(paste0("tt: ", tt))
    theta<-step1_opti$par
    #print(paste0("new theta: ", theta))
    
    ####################STEP 3 stop criteria##########################################
    j<- m/tt
    print(paste0("m/tt: ", j))
    if (j < par_e) {break}
    ####################STEP 4 increase t#############################################
    tt= (tt*par_mu)
  }
  
  #### Results: beta(theta)
  ans<-list(theta,tt)  
  assign("theta", theta, envir = .GlobalEnv)
  return(ans)
}
## Generate
result_index <- matrix(0,ncol = 9,nrow = 50)
result_esti_beta <- matrix(0,ncol =betanum*2+2,nrow = 50)

for (s in 1:50) {
  
  ####Find best tuning parameter  L1/L2
  fn_2_find_lamda(data.ori,data.vv,1000,0,5)
  
  best_lamda <- out[which.min(out[,2]),1] 
  min_loss <- out[which.min(out[,2]),2]

  #### Optimization 
  data.ori<-as.matrix(data.ori)
  fn_3_optim(data.ori,best_lamda,1,0.01,3) ##setlamda,par_t,par_e,par_mu

  ###################################################################################
  #                3 Generate Results: Estimate Err & Prediction Err & FNR & FPR
  ###################################################################################
  #### Function
  predict_err <- function(testdata){
    
    #### Basic setting
    count1 <- betanum+1
    count2 <- betanum+3
    count3 <- betanum*2+2
    count4 <- betanum+2
    testn <- dim(testdata)[1]
    X <- as.matrix(testdata[,1:betanum])
    
    ####Linear regression####
    b0.r_test <- theta[1]  
    B.r_test <- as.matrix(theta[2:count1])
    yr_hat <- b0.r_test + (X%*%B.r_test)
    
    ####Logistics regression####
    b0.c_test<- theta[count4]
    B.c_test <- as.matrix(theta[count2:count3])
    pi_hat <- (exp(b0.c_test+X%*%B.c_test)/(1+exp(b0.c_test+X%*%B.c_test))) ##P(y=1|X=x)
    yc_hat <- ifelse(runif(testn) < pi_hat, 1,0) ## prob < pi.1 then y=1
    
    err_predict_r <- sum(abs(testdata[,count1]-yr_hat))
    err_predict_c <- mean(testdata[,count4] != yc_hat)
    
    ####Results: estimated y 
    data.est <- cbind.data.frame(X,yr_hat,yc_hat)
    assign("err_predict_r", err_predict_r, envir = .GlobalEnv)
    assign("err_predict_c", err_predict_c, envir = .GlobalEnv)
    assign("data.est", data.est, envir = .GlobalEnv)
  }

  #### Generate: Estimation Error (true beta - estimated beta)
  count1 <- betanum+1
  count2 <- betanum+3
  count3 <- betanum*2+2
  count4 <- betanum+2
  err_esti_r <- (abs(b0.r - theta[1])) + sum(abs(B.r - theta[2:count1])) 
  err_esti_c <- (abs(b0.c - theta[count4])) + sum(abs(B.c - theta[count2:count3])) 
  err_esti <- err_esti_r + err_esti_c 
  #### Generate: Prediction Error
  predict_err(data.test)
  #### Generate: Variable Selection
  var_selection_r<- regsubsets(yr_hat ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20
                             +X21+X22+X23+X24+X25+X26+X27+X28+X29+X30+X31+X32+X33+X34+X35+X36+X37+X38+X39+X40
                             +X41+X42+X43+X44+X45+X46+X47+X48+X49+X50 ,data = data.est, nvmax = 10,really.big=T)
  xx <- data.est[,1:50]
  var_selection_c <- bess.one(xx,data.est$yc_hat ,s=10,family ="binomial")
  
  test_subset_r <- matrix((summary(var_selection_r, matrix.logical = T))$which ,nrow=10,ncol=betanum+1)[10,-1]
  test_subset_c <- matrix(FALSE, nrow=betanum,ncol = 1)
  test_subset_c[which(var_selection_c$beta!=0),1] <- TRUE
  real_subset <- matrix(c(rep(TRUE,10),rep(FALSE,40)))
  
  Rate_FNR_reg<-EvaluationMeasures.FNR(real_subset,test_subset_r)
  Rate_FPR_reg<-EvaluationMeasures.FPR(real_subset,test_subset_r)
  Rate_FNR_logi<-EvaluationMeasures.FNR(real_subset,test_subset_c)
  Rate_FPR_logi<-EvaluationMeasures.FPR(real_subset,test_subset_c)

  #### Results:
  result_index[s,] <- cbind(err_esti,err_esti_r,err_esti_c,err_predict_r,err_predict_c,Rate_FNR_reg,Rate_FPR_reg,Rate_FNR_logi,Rate_FPR_logi)
  colnames(result_index) <- c("Total Estimate Error","Regression Estimate Error","Logistics Estimate Error",
                              "Regression Prediction Error","Logistic Prediction Error",
                              "Regression FNR","Regression FPR","Logistics FNR","Logistics FPR")
  
  result_esti_beta[s,] <- c(theta[1:count3])
  print(paste0("s: ", s))
  print(paste0("index: ", result_index[s,]))
}



#################################################################################
colSums(result_index)/50 ##Average of index


#result_esti_beta    #### Estimated beta
#data.est            #### X, yr_hat, yc_hat



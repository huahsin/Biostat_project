############1. Simulation Dataset ###################################
library(numDeriv)
library('MASS')
library(nnet)
library(leaps)
library(EvaluationMeasures)
library(glmnet)
library(BeSS)
# simulated data2
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

data.sim<-fn_1_simulation(2020,50,10,0,2000)                        ##############################R=0
data.vv <- data.sim[501:1000,]     ##validation data
data.test <- data.sim[1001:2000,]  ##test data
data.ori <- data.sim[1:500,]  ##train data


###########2. find tuning parameter ###################################
result_index <- matrix(0,ncol = 9,nrow = 50)
result_esti_beta<- matrix(0,ncol =betanum*2+2,nrow = 50)
#####Repeat 50 times
for (s in 1:50) {
  
  ####Find best tuning parameter#### ##L1/L2
  fn_find_lamda<-function(traindata,testdata,par_lamda_size,par_lamdamin,par_lamdamax) { 
  
  par_lamda <- runif(par_lamda_size, min = par_lamdamin, max =par_lamdamax)
  Lossfn <- matrix(0,nrow = par_lamda_size,ncol=1) 
  betar <- matrix(0,nrow = par_lamda_size,ncol=betanum+1) 
  betac <- matrix(0,nrow = par_lamda_size,ncol=betanum+1) 
  ####set train & test data  
  train_x <- as.matrix(traindata[,1:betanum]) 
  train_yr <- as.matrix(traindata[,betanum+1])
  train_yc <- as.matrix(traindata[,betanum+2])
  
  for (k in 1:par_lamda_size) {
    ####ridge regression
    reg_l2<-glmnet(train_x, train_yr, family = "gaussian", alpha = 1,lambda = par_lamda[k])
    
    ####logistics regression with l2
    logi_l2<-glmnet(train_x, train_yc, family = "binomial", alpha = 1,lambda = par_lamda[k])
    
    #Generate Loss: Lr,Lc,Pq for kth fold
    #bring test data(x) into new modelthen get yhat  
    x <- as.matrix(testdata[,1:betanum]) 
    y_r <- as.matrix(testdata[,betanum+1])
    y_c <- as.matrix(testdata[,betanum+2])
    r_b0 <- coef(reg_l2)[1]
    r_beta <- coef(reg_l2)[-1]
    c_b0 <- coef(logi_l2)[1]
    c_beta <- coef(logi_l2)[-1]
    penalty <- matrix(0,nrow = betanum, ncol=1)
    
    Lr <- t(y_r-(r_b0 + x%*%r_beta)) %*% (y_r-(r_b0 + x%*%r_beta))
    Lc <- sum(-((y_c*(c_b0 + x%*%c_beta)) - log(1+exp(c_b0 + x%*%c_beta))))
    for (w in 1:betanum) {penalty[w] <- sqrt((r_beta[w]^2)+(c_beta[w]^2))}
    Pg <- par_lamda[k]*(sum(penalty))
    
    Lossfn[k] <- Lr + Lc + Pg
    for (p in 1:betanum+1) { ##p includes bo-b50
      betar[k,p] <- coef(reg_l2)[p]
      betac[k,p] <- coef(logi_l2)[p]
    }
    
  }
  
  out=cbind(par_lamda,Lossfn,betar,betac)
  assign("out", out, envir = .GlobalEnv)
  plot(x=par_lamda,y=Lossfn) 
}
  fn_find_lamda(data.ori,data.vv,1000,0,5)
  
  best_lamda <- out[which.min(out[,2]),1] 
  min_loss <- out[which.min(out[,2]),2]

  ##############3. Optimization ###########################################
  find_Lbarrier <- function(theta) {
    count1 <- betanum+1  
    count2 <- betanum+2
    count3 <- betanum+3
    count4 <- betanum*2+2
    r_b0 <- theta[1]
    r_beta <- theta[2:count1]
    c_b0 <- theta[count2]
    c_beta <- theta[count3:count4]
    opti_uj<-matrix(0,nrow=betanum,ncol=1)
    
    for (k in 1:betanum) {opti_uj[k]<-(max(abs(theta[k+1]),abs(theta[k+count2])))+0.5} ###FOR UJ
    uj <- opti_uj
    Lr <- t(y_r-(r_b0 + x%*%r_beta)) %*% (y_r-(r_b0 + x%*%r_beta))  
    
    logterm <- matrix(0,ncol=1,nrow = dim(x)[1])  ###FOR reduce the dimension of index: solving the INF problem.
    for (m in 1:dim(x)[1]) {
      index <- c_b0 + x[m,]%*%c_beta
      if (index > log(10000)) {logterm[m] <- c_b0 + x[m,]%*%c_beta} 
      else if (index < log(0.0001)){logterm[m] <- log(1)} 
      else {logterm[m] <- log(1+exp(c_b0 + x[m,]%*%c_beta))}
    }
    
    Lc <- sum(-((y_c*(c_b0 + x%*%c_beta)) - logterm))
    Pg <- (lamda*sum(uj)) + (-(sum(log(r_beta+uj)+log(uj-r_beta)))*(1/tt)) + (-(sum(log(c_beta+uj)+log(uj-c_beta)))*(1/tt))
    return (Lr + Lc + Pg)
  }  ##Function
  data.ori<-as.matrix(data.ori)
  fn_4_optim<-function(traindata,setlamda, par_t, par_e, par_mu){
    traindata<-as.matrix(traindata)
    ##Random generate initial theta
    thetanum <- betanum*3+2
    totalbetanum<-betanum*2+2
    assign("totalbetanum", totalbetanum, envir = .GlobalEnv)
    assign("thetanum", thetanum, envir = .GlobalEnv)
    theta<-matrix(0,ncol = thetanum,nrow=1)
    opti_r_b0<-runif(1,min=0,max=3)   ##theta
    opti_c_b0<-runif(1,min=0,max=3)   ##theta
    opti_r_beta<-as.matrix(runif(betanum,min=5,max=10),nrow=betanum,ncol=1)   ##theta
    opti_c_beta<-as.matrix(runif(betanum,min=5,max=10),nrow=betanum,ncol=1)   ##theta
    opti_uj<-matrix(0,nrow=betanum,ncol=1) ##strictly feasible #betanum+1
    
    theta<-c(opti_r_b0,opti_r_beta,opti_c_b0,opti_c_beta)                                           #REMOVE UJ#
    for (k in 1:betanum) {opti_uj[k]<-(max(abs(theta[k+1]),abs(theta[k+betanum+2])))+0.5}
    #for (k in 1:4) {theta[k+10]<-(max(abs(theta[k+1]),abs(theta[k+6])))+0.5}
    tt <- par_t
    m<-2*betanum*2 ##2*betanum*(regression task number+ logistic task number)
    j<-(m/tt)
    
    repeat{ 
      
      ####################STEP 1 compute theta by Newton Method########################
      assign("tt", tt, envir = .GlobalEnv) ##"assign" means value tt can be use in function
      assign("lamda", setlamda, envir = .GlobalEnv)
      assign("x", traindata[,1:betanum], envir = .GlobalEnv)
      assign("y_r", traindata[,betanum+1] , envir = .GlobalEnv)
      assign("y_c", traindata[,betanum+2] , envir = .GlobalEnv)
      
      step1_opti<-optim(theta, find_Lbarrier, method = "CG",hessian = TRUE, control = list(type = 2))#,  control =list(maxit=1000), control = list(abstol=0.1,type=1)
      
      #print(paste0("convergence: ", step1_opti$convergence))
      
      ####################STEP 2 update theta###########################################
      #print(paste0("old theta: ", theta))
      #print(paste0("tt: ", tt))
      theta<-step1_opti$par
      print(paste0("new theta: ", theta))
      
      ####################STEP 3 stop criteria##########################################
      j<- m/tt
      print(paste0("m/tt: ", j))
      
      if (j < par_e) {break}
      ####################STEP 4 increase t#############################################
      tt= (tt*par_mu)
      
    }
    
    ans<-list(theta,tt)  
    assign("theta", theta, envir = .GlobalEnv)
    return(ans)
    
  }
  
  fn_4_optim(data.ori,best_lamda,1,0.001,3) ##setlamda,par_t,par_e,par_mu

  ############## Interest outcome ############################
  count1 <- betanum+1
  count2 <- betanum+3
  count3 <- betanum*2+2
  count4 <- betanum+2
  ##Estimation Error (true beta - estimated beta)
  err_esti_r <- ((b0.r - theta[1])^2) + sum((B.r - theta[2:count1])^2) 
  err_esti_c <- ((b0.c - theta[count4])^2) + sum((B.c - theta[count2:count3])^2) 
  err_esti <- err_esti_r + err_esti_c 
  
  ##Prediction Error
  predict_err <- function(testdata){
    
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
  
    err_predict_r <- sum((testdata[,count1]-yr_hat)^2)
    err_predict_c <- mean(testdata[,count4] != yc_hat)
    err_misclass_table <- table(yc_hat,testdata[,count4])
    
    data.est <- cbind.data.frame(X,yr_hat,yc_hat)
    assign("err_predict_r", err_predict_r, envir = .GlobalEnv)
    assign("err_predict_c", err_predict_c, envir = .GlobalEnv)
    assign("err_misclass_table", err_misclass_table, envir = .GlobalEnv)
    assign("data.est", data.est, envir = .GlobalEnv)
   
  }
  predict_err(data.test)
  
  ##Variable Selection
  var_selection_r<- regsubsets(yr_hat ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20
                             +X21+X22+X23+X24+X25+X26+X27+X28+X29+X30+X31+X32+X33+X34+X35+X36+X37+X38+X39+X40
                             +X41+X42+X43+X44+X45+X46+X47+X48+X49+X50 ,data = data.est, nvmax = 10,really.big=T)
  xx <- data.est[,1:50]
  var_selection_c<- bess.one(xx,data.est$yc_hat ,s=10,family ="binomial")
  
  test_subset_r<-matrix((summary(var_selection_r, matrix.logical = T))$which ,nrow=10,ncol=betanum+1)[10,-1]
  
  test_subset_c<-matrix(FALSE, nrow=betanum,ncol = 1)
  test_subset_c[which(var_selection_c$beta!=0),1]<-TRUE
  real_subset <- matrix(c(rep(TRUE,10),rep(FALSE,40)))
  
  Rate_FNR_reg<-EvaluationMeasures.FNR(real_subset,test_subset_r)
  Rate_FPR_reg<-EvaluationMeasures.FPR(real_subset,test_subset_r)
  Rate_FNR_logi<-EvaluationMeasures.FNR(real_subset,test_subset_c)
  Rate_FPR_logi<-EvaluationMeasures.FPR(real_subset,test_subset_c)

  ############## Output ############################
  ##Result: Estimate Error, Prediction error, False Positive Rate(variable selection)
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






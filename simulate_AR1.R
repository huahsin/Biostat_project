#######Simulation Data - AR(1)#########################################
fn_1_simulation_ar1<-function(par_seed,par_betanum,par_rele_num,par_r,par_samplesize){
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
        msigma[i,j]<-ifelse(mm[i,j]==1,1,r^abs(i-j))
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

data.sim<-fn_1_simulation_ar1(2020,50,10,0.8,2000) 
data.vv <- data.sim[501:1000,]     ##validation data
data.test <- data.sim[1001:2000,]  ##test data
data.ori <- data.sim[1:500,]  ##train data
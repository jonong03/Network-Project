# Clean code: Functions with examples
# This version of Ising model includes self-interaction term (asymmetric Ising)
# For Symetrical Ising, simply set diagonal theta entries as 0
####### Theta's order ########
# t11 | t12 | t13 | t14 | t15
# t21 | t22 | t23 | t24 | t25
# t31 | t32 | t33 | t34 | t35
# t41 | t42 | t43 | t44 | t45
# t51 | t52 | t53 | t54 | t55
##############################

# Load required packages
# gc(); rm(list=ls())
library(dplyr); library(data.table); library(Deriv);library(gee); library(numDeriv)

# Function 1: Density function
dising <- function(vv, p, theta, sign=1, take_log=FALSE){
  if(length(theta) != p + (p*(p-1)/2)){
    cat("wrong dimension of parameter theta\n")
    return()
  }
  #  theta= param
  if(sign==1){a= c(0,1)}
  if(sign==2){a= c(-1,1)}
  
  config0 = expand.grid(rep(list(a),p)) %>% as.matrix
  P <- matrix(0, p, p)
  P[lower.tri(P, diag=TRUE)] <- theta
  P = t(P)
  
  #pmf_inner <- function(y) exp(t(y)%*%P%*%y/2)
  pmf_inner <- function(y) exp(t(y)%*%P%*%y)
  den <- sapply(1:nrow(config0), function(i) pmf_inner(config0[i,])) %>% sum
  
  vv <- as.matrix(vv)
  if(dim(vv)[2] == 1) vv <- t(vv)
  if(take_log){
    num <- (vv%*%P%*%t(vv)) %>% diag   #How to incorporate self-interaction term here?
    return(num - log(den))
  }else{
    num <- exp(vv%*%P%*%t(vv)) %>% diag
    return(num/den)
  }
}

# Function 2: Random number generator
rising <- function(n, p, theta, sign=1){
  
  if(length(theta) != p + p*(p-1)/2){
    cat("wrong dimension of parameter vector\n")
    return()
  }
  
  if(sign==1){a= c(0,1)}
  if(sign==2){a= c(-1,1)}
  
  config0 = expand.grid(rep(list(a),p)) %>% as.matrix
  
  P <- matrix(0, p, p)
  P[lower.tri(P, diag=TRUE)] <- theta
  P = t(P)
  
  #pmf_inner <- function(y) exp(t(y)%*%P%*%y/2)
  pmf_inner <- function(y) exp(t(y)%*%P%*%y)
  pr <- sapply(1:nrow(config0), function(i) pmf_inner(config0[i,]))
  cdf <- (pr/sum(pr)) %>% cumsum
  
  oo <- runif(n)
  oo <- sapply(1:n, function(i) sum(oo[i] > cdf)) + 1
  Y <- matrix(0, n, p)
  for(i in 1:n) Y[i,] <- config0[oo[i],]
  
  return(Y)
}

# Function 3: Max likelihood estimation for thetas. Methods include: newton, logit, gd (gradient descent)
MLE_ising <- function(Y, sign=1, method, lambda=0, maxi=25, start = rep(0, dim(Y)[2]+dim(Y)[2]*(dim(Y)[2]-1)/2), learn_rate=0.5, tol=0.001, corstr= "independence", verbose=FALSE ){
  options(dplyr.summarise.inform = FALSE)
  Y = as.matrix(Y)
  m = dim(Y)[1]
  p = dim(Y)[2]
  colnames(Y)= paste0("x",1:p)
  vars = colnames(Y)
  
  if(sign==1){a= c(0,1)}
  if(sign==2){a= c(-1,1)}
  
  Q = matrix(NA, p, p)
  Q[upper.tri(Q, diag=TRUE)] = 0
  loc = which(Q==0, arr.ind=TRUE,useNames=TRUE) %>% as.data.frame %>% arrange(row)
  t_list = paste0("t_",loc$row, "_",loc$col)
  
  #Method: Maximum likelihood by numerical solution
  if(method=="newton"){
    
    loglik <- function(param) -sum(dising(Y, p, param, sign, take_log=TRUE))
    
    est <- nlm(loglik, rep(0,p+(p*(p-1)/2)), hessian=TRUE, gradtol = 1e-5)
    se <- est$hessian %>% solve %>% diag %>% sqrt
   
    output = cbind(Estimate=est$estimate, se, zv= est$estimate/se, pv=pnorm(-abs(est$estimate/se)))
    rownames(output)= t_list
    
  }
  
  #Method: Maximum likelihood by analytical solution
  if(method=="newton2"){
    
    loglik= function(Y){
      pmf = function(z){
        z= as.matrix(z)
        zz= t(z) %*% z 
        zz= zz[lower.tri(zz, diag= TRUE)]
        return( paste0("exp(",paste0( zz, "*", t_list, collapse="+"),")"))
      }
      config0= expand.grid(rep(list(a),p))
      
      pmf_all = lapply(1:dim(config0)[1], function(i) pmf(config0[i,])) %>% paste(collapse="+")
      part2 = paste0("m*log(",pmf_all,")")
      #part2_f = sapply(1:length(t_list), function(k) Deriv(part2, t_list[k]))
      
      B = (t(Y) %*% Y)
      b = B[lower.tri(B, diag=TRUE)]
      part1 = paste0(t_list,"*",b, collapse="+")
      
      f = parse(text = paste0(part1,"-",part2 ))
      
      return(f)
    }
    ll= loglik(Y)
    
    #Start Newton Method
    x <- start
    x <- setNames(x, t_list)
    llf <- function(k) -eval(ll, envir= as.list(k))
    
    i=0
    diff=rep(99999, p+p*(p-1)/2)
    #Max abs diff < 1e-6
    while (max(abs(diff)) <= 1e-6){
      diff = (hessian(llf, x) %>% solve ) %*% t(jacobian(llf, x)) %>% as.vector
      #?solve  X <- solve(A,B) solve(hessian, jacobian)
      x = x- diff
      
      i = i+1
      if(sum(is.na(diff))>1 | i > 50) { break }
    }
    
    if (sum(is.na(diff))>1){
      output = matrix(NA, ncol=4, nrow= p+p*(p-1)/2) %>% as.data.frame
      rownames(output)<- t_list
    } else {
      se = (solve(hessian(llf, x)) %>% diag %>% sqrt)
      zv = x/ se
      pv = pnorm(-zv)
      output = data.frame(x, se, zv, pv)
      }
    }
  
  #Method: Gradient Descent
  if(method=="gd"){
    
    config0 = expand.grid(rep(list(a),p)) %>% as.matrix
    colnames(config0) = paste0("x",1:p)
    
    #Log-likelihood function
    lik = function(Y){
      
      # Construct Potential equation
      l=0
      for (j in 1: dim(config0)[1]){
        X = as.matrix(config0[j,])
        X_list = X %*% t(X)
        X_list = X_list[upper.tri(X_list, diag=T)]
        
        theta_list2=NULL
        for (i in 1:length(X_list)){
          theta_list2 = c(theta_list2, ifelse(X_list[i]==0, 0 , t_list[i]))
        }
        f = paste0("exp(",paste0(theta_list2,collapse="+"),")")
        l = paste0(l,"+",f)
        
      }
      
      Part2 = paste0("log(",l,")")
      
      A = t(Y) %*% Y
      a = A[upper.tri(A, diag=T)]/dim(Y)[1]
      Part1 = paste0(t_list,"*",a,collapse ="+")
      f = parse(text = paste0(Part1,"-",Part2 ))
      
      return(f)
    }
    likf = lik(Y)
    
    #Inputs for Gradient Descent
    oo = rep(0,p+(p*(p-1)/2))
    names(oo) = t_list
    oo_new = oo
    
    grad = lapply(1:length(t_list), function(i) Deriv(likf, t_list[i]))
    
    #Executing Gradient Descent
    #Three layers: inner - find convergence for thetaj, middle - repeat inner loop for all js, outer - making sure all thetas are coverged
    i=0
    while(i < 1000){
      
      i= i+1
      for(j in 1:length(grad)){
        
        k=1
        while(k < 1000){
          diff = eval(grad[[j]], envir = as.list(oo))
          oo[j] = oo[j]+ (learn_rate/sqrt(k)) * diff
          if((abs(diff)<= tol)) {break} else {k = k+1}
        }
        
      }
      
      if(verbose ==TRUE) { print(oo) } 
      oo_new = rbind(oo_new,oo)
      if (sum(abs(diff(tail(oo_new,n=2))) <= tol) ==length(t_list)) {
        output= data.frame(Estimate= oo)
        break
      }
    }
  }  
  
  #Method: Logistic estimates with variance estimation from data's likelihood 
  if(method=="logit"){
    
    #Design Matrix
    tempm_ = rep(list(matrix(data=0, nrow= p*m, ncol= p)), p)
    Xm_ = NULL
    for (i in 1:p){
      for (j in 1:p){
        tempm_[[i]][(1+(m*(i-1))):(m*i),j] <- Y[,j] # Fill diagonal entries in temp_ with Y[,j]
        tempm_[[i]][(1+(m*(j-1))):(m*j),j] <- Y[,i] # Fill 
      }
      tempm_[[i]][(1+(m*(i-1))):(m*i),i] <- 1 # Fill diagonal entries in temp_ with 1 (Conditional Prob)
      Xm_= cbind(Xm_, tempm_[[i]][,i:p])
      
    }
    colnames(Xm_) <- t_list
    
    model = glm(c(Y)~ 0+Xm_, family= binomial(link="logit"))
    output = summary(model)$coefficients
    
    #VAR_Ising function (incl lambda - ridge penalty)
    VAR_Ising <- function(Y, lambda){
      m <- ncol(Y); n <- nrow(Y); m2 <- m*(m+1)/2
      
      YoY <- t(sapply(1:n, function(i) Y[i,]%x%Y[i,]))
      
      #S <- var(YoY) + lambda * diag(m^2)
      S <- var(YoY)
      
      #eigen(S)$values
      H <- matrix(0, m2, m*m)
      index <- 0
      for(i in 1:m){
        for(j in i:m){
          index <- index + 1
          I <- matrix(0, m, m)
          if(i==j){
            I[i,i] <- 1
          }else{
            I[i,j] <- .5; I[j,i] <- .5
          }
          H[index,] <- c(I)
        }
      }
      V <- solve(H%*%S%*%t(H) + lambda * diag(m2) )/n
      return(V)
    }
    V <- VAR_Ising(Y,lambda=0)
    se_logit <- sqrt(diag(V))
    
    #Inspect which variable is missing (NA)
    na_var<- !(t_list %in% gsub("Xm_","",row.names(output)))
    temp= matrix(NA, ncol=4, nrow= sum(na_var))
    output<- rbind(output, temp)
    row.names(output) <- c(t_list[!na_var], t_list[na_var])
    
    est= as.vector(output[,1])
    output[,2]<- se_logit
    output[,3]<- est/ se_logit
    output[,4]<- pnorm(-abs(est/se_logit))
    
  }   
  
  #Method: Logistic estimates from glm 
  if(method=="logit0"){
    
    #Design Matrix
    tempm_ = rep(list(matrix(data=0, nrow= p*m, ncol= p)), p)
    Xm_ = NULL
    for (i in 1:p){
      for (j in 1:p){
        tempm_[[i]][(1+(m*(i-1))):(m*i),j] <- Y[,j] # Fill diagonal entries in temp_ with Y[,j]
        tempm_[[i]][(1+(m*(j-1))):(m*j),j] <- Y[,i] # Fill 
      }
      tempm_[[i]][(1+(m*(i-1))):(m*i),i] <- 1 # Fill diagonal entries in temp_ with 1 (Conditional Prob)
      Xm_= cbind(Xm_, tempm_[[i]][,i:p])
      
    }
    colnames(Xm_) <- t_list

    model = glm(c(Y)~ 0+Xm_, family= binomial(link="logit"))
    output = summary(model)$coefficients
    
    #Inspect which variable is missing (NA)
    na_var<- !(t_list %in% gsub("Xm_","",row.names(output)))
    temp= matrix(NA, ncol=4, nrow= sum(na_var))
    output<- rbind(output, temp)
    row.names(output) <- c(t_list[!na_var], t_list[na_var])
  }
  
  #Method: Logistic estimates with variance estimation from fisher information score (MLE)
  if(method=="logit1"){
    
    #Design Matrix
    tempm_ = rep(list(matrix(data=0, nrow= p*m, ncol= p)), p)
    Xm_ = NULL
    for (i in 1:p){
      for (j in 1:p){
        tempm_[[i]][(1+(m*(i-1))):(m*i),j] <- Y[,j] # Fill diagonal entries in temp_ with Y[,j]
        tempm_[[i]][(1+(m*(j-1))):(m*j),j] <- Y[,i] # Fill 
      }
      tempm_[[i]][(1+(m*(i-1))):(m*i),i] <- 1 # Fill diagonal entries in temp_ with 1 (Conditional Prob)
      Xm_= cbind(Xm_, tempm_[[i]][,i:p])
      
    }
    colnames(Xm_) <- t_list
    
    model = glm(c(Y)~ 0+Xm_, family= binomial(link="logit"))
    output = summary(model)$coefficients

    #Inspect which variable is missing (NA)
    na_var<- !(t_list %in% gsub("Xm_","",row.names(output)))
    temp= matrix(NA, ncol=4, nrow= sum(na_var))
    output<- rbind(output, temp)
    row.names(output) <- c(t_list[!na_var], t_list[na_var])
    
    #Estimate SE
    loglik= function(Y){
      pmf = function(z){
        z= as.matrix(z)
        zz= t(z) %*% z 
        zz= zz[lower.tri(zz, diag= TRUE)]
        return( paste0("exp(",paste0( zz, "*", t_list, collapse="+"),")"))
      }
      config0= expand.grid(rep(list(a),p))
      
      pmf_all = lapply(1:dim(config0)[1], function(i) pmf(config0[i,])) %>% paste(collapse="+")
      part2 = paste0("m*log(",pmf_all,")")
      #part2_f = sapply(1:length(t_list), function(k) Deriv(part2, t_list[k]))
      
      B = (t(Y) %*% Y)
      b = B[lower.tri(B, diag=TRUE)]
      part1 = paste0(t_list,"*",b, collapse="+")
      
      f = parse(text = paste0(part1,"-",part2 ))
      
      return(f)
    }
    ll= loglik(Y)
    
    est= as.vector(output[,1])
    x <- setNames(est, t_list)
    llf <- function(k) -eval(ll, envir= as.list(k))
    se = (solve(hessian(llf, x)) %>% diag %>% sqrt)

    output[,2]<- se
    output[,3]<- est/ se
    output[,4]<- pnorm(-abs(est/se))
  }  
  
  #Method: Ridge penalty on logistic regression with variance estimation from data's likelihood
  if(method=="ridge"){
    
    #Design Matrix
    tempm_ = rep(list(matrix(data=0, nrow= p*m, ncol= p)), p)
    Xm_ = NULL
    for (i in 1:p){
      for (j in 1:p){
        tempm_[[i]][(1+(m*(i-1))):(m*i),j] <- Y[,j] # Fill diagonal entries in temp_ with Y[,j]
        tempm_[[i]][(1+(m*(j-1))):(m*j),j] <- Y[,i] # Fill 
      }
      tempm_[[i]][(1+(m*(i-1))):(m*i),i] <- 1 # Fill diagonal entries in temp_ with 1 (Conditional Prob)
      Xm_= cbind(Xm_, tempm_[[i]][,i:p])
      
    }
    colnames(Xm_) <- t_list
    
    logistic <- function(y, X, lambda, beta=NULL, lb=10, tol=10^-6, maxite=20){
      p <- ncol(X)
      n <- length(y)
      
      if(is.null(beta)){beta <- rep(0, p); beta[1] <- 1}
      
      for(i in 1:maxite){
        phat <- c(X%*%beta)
        phat <- ifelse(phat > lb, lb, ifelse(phat < -lb, -lb, phat))
        phat <- exp(phat)/(1+exp(phat))
        vhat <- phat*(1-phat)
        b_penalty <- beta
        beta_new <- beta + solve(t(X*vhat)%*%X-lambda*diag(p), t(X)%*%(y-phat)-lambda*beta)
        test <- max(abs(beta-beta_new))
        beta <- beta_new
        if(test < tol) break
      }
      return(beta)
    }
    output = logistic(y= c(Y), X= Xm_, lambda)
    
    VAR_Ising <- function(Y, lambda){
      m <- ncol(Y); n <- nrow(Y); m2 <- m*(m+1)/2
      
      YoY <- t(sapply(1:n, function(i) Y[i,]%x%Y[i,]))
      
      #S <- var(YoY) + lambda * diag(m^2)
      S <- var(YoY)
      
      #eigen(S)$values
      H <- matrix(0, m2, m*m)
      index <- 0
      for(i in 1:m){
        for(j in i:m){
          index <- index + 1
          I <- matrix(0, m, m)
          if(i==j){
            I[i,i] <- 1
          }else{
            I[i,j] <- .5; I[j,i] <- .5
          }
          H[index,] <- c(I)
        }
      }
      V <- solve(H%*%S%*%t(H) + lambda * diag(m2) )/n
      return(V)
    }
    V <- VAR_Ising(Y,lambda=lambda)
    se <- sqrt(diag(V))
    
    output= as.data.frame(output)
    est= as.vector(output[,1])
    output$`Std. Error`<- se
    output$`z value`<- est/ se
    output$`Pr(>|z|)`<- pnorm(-abs(est/se))
    

    
  }
  
  #Method: Log-normal regression
  if(method=="poisson"){
    Y_table = Y %>% as.data.table %>% group_by_all %>% summarize(n=n())
    model = glm(n~.^2, data=Y_table, family=poisson(link="log"))
    output= summary(model)$coefficients
    
    new_name_f = function(oldnames){
      new = gsub(":","",oldnames)
      i = tstrsplit(new,"x")[[2]] %>% as.integer
      j = tstrsplit(new,"x")[[3]] %>% as.integer
      j[which(is.na(j))] = i[which(is.na(j))]
      new = paste0("t","_",i,"_",j)
      new[1]="(Intercept)"
      return(new)
    }
    rownames(output)= new_name_f(rownames(output))

    #Inspect which variable is missing (NA)
    na_var<- !(t_list %in% row.names(output))
    temp= matrix(NA, ncol=4, nrow= sum(na_var))
    output<- rbind(output, temp)
    row.names(output) <- c("(Intercept)",t_list[!na_var], t_list[na_var])
    
  }
  
  #Method: Generalized estimating equations
  if(method=="gee"){
    
    #Design Matrix
    tempm_ = rep(list(matrix(data=0, nrow= p*m, ncol= p)), p)
    Xm_ = NULL
    for (i in 1:p){
      for (j in 1:p){
        tempm_[[i]][(1+(m*(i-1))):(m*i),j] <- Y[,j] # Fill diagonal entries in temp_ with Y[,j]
        tempm_[[i]][(1+(m*(j-1))):(m*j),j] <- Y[,i] # Fill 
      }
      tempm_[[i]][(1+(m*(i-1))):(m*i),i] <- 1 # Fill diagonal entries in temp_ with 1 (Conditional Prob)
      Xm_= cbind(Xm_, tempm_[[i]][,i:p])
      
    }
    colnames(Xm_) <- t_list
    ID= rep(1:dim(Y)[1], dim(Y)[2])
    
    model <- gee(c(Y)~ 0+Xm_, id=ID, family= binomial(link = "logit"), corstr= "independence", maxiter= maxi)
    output = summary(model)$coefficients %>% as.data.frame %>% select(Estimate, se= `Robust S.E.`, zval= `Robust z`) %>% mutate(pv=pnorm(-abs(Estimate/se)))
 
    #model <- geeglm(c(Y)~ 0+Xm_, id=ID, family= binomial(link = "logit"), corstr)
    #output = summary(model)$coefficients %>% as.data.frame %>% select(Estimate, se= `Std.err`, zval= `Wald`, pv= `Pr(>|W|)`)
    
    newnames = gsub("Xm_","",rownames(output))
    rownames(output)<- newnames
  }
  
  #Output reformatting
  output= as.data.frame(output)
  output.names=c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  setnames(output, new= output.names)
  output= output[match(t_list, row.names(output)),]
  
  return(output)
  
}


##############################
# Run Example

# Specify parameters
n= 200; p= 5; param_v= rep(-0.1, p+p*(p-1)/2)

# Random generation
Y= rising(n, p, theta=param_v, sign=1)

# Density function
Yden= dising(Y, p=5, theta=param_v, sign=1, take_log=TRUE)

# Parameter estimation
YMLE= MLE_ising(Y, method="newton2"); YMLE
Ylogit= MLE_ising(Y, method="logit3"); Ylogit

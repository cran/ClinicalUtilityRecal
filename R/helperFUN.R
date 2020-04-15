###### 2: Functions used for optimization
## constraint function
constSNB <- function(alpha,y,p,r,ineqBD){
  lp.recal <- alpha[1] + alpha[2]*logit(p)
  p.recal <- logit.inv(lp.recal)

  sNB <- nb(y = y,p = p.recal,r = r)$snb
  res <-  ineqBD - sNB
  return(res)
}


## objective fcn - negative log liklihood
logLik <- function(alpha,y,p,r,ineqBD){
  z <- logit(p)
  logLik <- sum(-y*(alpha[1] + alpha[2]*z) + log(1 + exp(alpha[1] + alpha[2]*z)))
  return(logLik)
}

## function to optimize in first step
thresh.snb.smooth <- function(t,y,p,r){
  if(length(p)!=length(y)) stop('length of p does not match length y')

  xD <- p[y==1]
  xND <- p[y==0]
  sTPR <- mean(pnorm(q = (xD-t)/0.01))
  sFPR <- mean(pnorm(q = (xND-t)/0.01))
  y.bar <- mean(y)
  sNB <- (sTPR - ((r/(1-r)) * ((1-y.bar)/y.bar) * sFPR))
  return(sNB)
}

#variance needed for lower bound
snbVar.tmax <- function(tVec,y,p,r){
  if(length(p)!=length(y)) stop('length of p does not match length y')
  if(any(tVec <0) | any(tVec >1)) stop('t must be between 0 and 1')

  n <- length(p)
  snbVarVec <- rep(NA,length(tVec))
  for(j in 1:length(tVec)){
    p11 <- mean(y==1 & p >= tVec[j])
    p10 <- mean(y==0 & p >= tVec[j])
    p01 <- mean(y==1 & p < tVec[j])
    p00 <- mean(y==0 & p < tVec[j])

    k <- r/(1-r)
    sum <- (p10*(p10+p01)*k^2) + p11*(p01+(k^2*p10))
    var <- sum/(p11 + p01)^3
    snbVarVec[j] <- sqrt(var/n)
  }
  snbVar.tmax <- snbVarVec
  return(snbVar.tmax)
}


###### Recal Wt fcns ##########
logit.inv <- function(x){1/(1+exp(-x))}
logit <- function(p){log(p/(1-p))}

## Mapping RAW value to tuning parameter
getDelta <- function(delta,lambda,r,rl,ru,rrWt,p,y,rl.raw,ru.raw){
  #[rl,ru] is for weight defintion
  #[rl.raw,ru.raw] is for raw bndwith

  loessCurv <- stats::loess(y~p,span = 2/3,degree=1)
  eventRate <- predict(loessCurv)

  wt <- calWt(r = r,p = p,y = y,lambda = lambda,rl = rl,ru = ru,delta = delta)

  #avg weights
  bnd.wt <- mean(wt[eventRate >= rl.raw & eventRate <= ru.raw])
  ## if band is whole range, make outside band weight 1
  notBnd.wt <-  mean(wt[eventRate < rl.raw | eventRate > ru.raw])

  rr.est.wt <- (notBnd.wt/bnd.wt)
  res <- rr.est.wt - rrWt
  return(res)
}
getLambda <- function(lambda,delta,r,rl,ru,rrWt,p,y,rl.raw,ru.raw){
  #[rl,ru] is for weight defintion
  #[rl.raw,ru.raw] is for raw bndwith

  loessCurv <- stats::loess(y~p,span = 2/3,degree=1)
  eventRate <- predict(loessCurv)

  wt <- calWt(r = r,p = p,y = y,lambda = lambda,rl = rl,ru = ru,delta = delta)

  #avg weights
  bnd.wt <- mean(wt[eventRate >= rl.raw & eventRate <= ru.raw])
  notBnd.wt <-  mean(wt[eventRate < rl.raw | eventRate > ru.raw])

  rr.est.wt <- (notBnd.wt/bnd.wt)
  res <- ifelse(is.nan(rr.est.wt),Inf,rr.est.wt - rrWt)
  return(res)
}


tpr.fun <- function(p,y,r){
  tpr <- mean(p[y==1] >= r)
  return(tpr)
}
fpr.fun <- function(p,y,r){
  fpr <- mean(p[y==0] >= r)
  return(fpr)
}


#### Plotting Function ####

# get sNB curve as function of decesion threshold t
snb.t <- function(par,y,p,r){
  #par is the decision threhsold
  TPR <- mean(p[y==1] > par)
  FPR <- mean(p[y==0] > par)
  y.bar <- mean(y)
  sNB <- (TPR - ((r/(1-r)) * ((1-y.bar)/y.bar) * FPR))
  return(sNB)
}

#####  Constrained Logistic Recalibration ########
##### Functions: Functions to implement constrained calibration approach and compare to standard logistic recalibration


###### 2: Functions used for optimization
## constraint function



## main optimization function
constRecal <- function(y,p,r,int=NULL,alphaLB=c(-10,0),alphaUB=c(10,10),
                        ftol=1e-8,xtol=1e-4,maxeval=1e6){
  #if no intial values specified using logistic recalibrated risks
  if(is.null(int)){
    int <- stdRecal(y = y,p = p)$alpha
  }

  ## Step 1: Get maximum sNB and threshold
  ## get max SNB  and corresponding max value
  ## for one just use which ever one is bigger - probably need to change this
  tMax.optim <- optim(par = r,fn=thresh.snb.smooth,r = r,y = y, p = p, method = "L-BFGS-B",
                control=list(fnscale=-1),lower=0, upper=1)
  tmax <- optimize(f = thresh.snb.smooth,interval = c(0,1),maximum = TRUE,r = r,y = y, p = p)
  sNBMax <- ifelse(tMax.optim$value > tmax$objective, tMax.optim$value, tmax$objective)
  sNBMax.t <- ifelse(tMax.optim$value > tmax$objective, tMax.optim$par, tmax$maximum)

  ## Getting standard error lower bound
  snbMax.se <- snbVar.tmax(tVec = sNBMax.t,
                             y = y,p = p,
                             r = r)
  epsilon <- sNBMax - (snbMax.se)


  constLogLik <- nloptr::nloptr( x0= int,
                         eval_f=logLik,eval_g_ineq = constSNB,eval_g_eq = NULL,
                         lb = alphaLB, ub = alphaUB,
                         opts = list("algorithm" = "NLOPT_GN_ORIG_DIRECT",
                                     "ftol_rel" = ftol,
                                     "print_level" = 0,
                                     "maxeval" = maxeval),
                         y = y, p = p,
                         r = r, ineqBD = epsilon)

  log.phat.const <- constLogLik$solution[1] + (constLogLik$solution[2]*logit(p))
  p.const <- logit.inv(log.phat.const)

  res <- list("alpha"=constLogLik$solution,"y"=y,"p.const"=p.const,
              "snbLB" = epsilon,
              "optimRes" = constLogLik)

  return(res)

}







##### Weighted Recalibration ########
##### Fucntions: Implementation and Evaluation Functions for method #####
##### Date Created: 11/14/18 ######


### Purpose of Code: Functions to implement weighted recalibration approach and compare to
###                  standard logistic recalibration


### Outline: 1) Weight Functions
###          2) Weighted Recalibration Methods Functions
###          3) Tuning parameter selection functions
###          4) Evaluation functions (calibration plots, nb, tpr, fpr)
###          5) Simulation Functions


##### 1: Weight Functions ####
## we use weights based on observed event rates estimate from LOESS regression
## we use exponential decay weights to get smoother weights
## the indicator gives a more general form
## when rl,ru = [0,1] the indicator will be the special case of exponential decay
# weights based on euclidean distance with exponential decay

###### 1: Auxillary Functions ##########

calWt <- function(rl,ru,p,y,r,lambda,delta,returnSmoothedEvent=FALSE){

  if( rl > ru ) stop('lower bound(rl) must be less than upper bound (ru)')
  if( (rl > r) | (ru < r) ) stop('r must be between lower and upper bound')
  if( r > 1 | r < 0 ) stop('r must be between 0 and 1')

  if(length(p)!=length(y)) stop('length of p does not match length y')
  if(lambda < 0 | delta <0 ) stop('tuning parameters must be positive')

  loessCurv <- stats::loess(y~p,span = 2/3,degree=1)
  eventRate <- predict(loessCurv)
  d <- (eventRate - r)^2

  wt <- ifelse(eventRate >= rl & eventRate <= ru,exp(-d/lambda),delta)

  if(returnSmoothedEvent==TRUE){
    res <- list("wt" = wt, "o" = eventRate)
    return(res)
  }
  else if(returnSmoothedEvent==FALSE){
    return(wt)
  }
}

stdRecal <- function(y,p){
  if(length(p)!=length(y)) stop('length of p does not match length y')

  mod <- glm(y ~ logit(p), family="binomial")
  alpha <- mod$coef
  names(alpha) <- c("Intercept","Slope")

  stdRisk <- logit.inv(mod$coef[1] + mod$coef[2]*logit(p))
  res <- list("y"=y,"p.std" = stdRisk, "alpha" = alpha)
  return(res)
}

###### 2: Weighted Recalibration Implementation #####
wtRecal <- function(y,p,r,rl,ru,lambda,delta){
  if( rl > ru ) stop('lower bound(rl) must be less than upper bound (ru)')
  if( (rl > r) | (ru < r) ) stop('r must be between lower and upper bound')
  if( r > 1 | r < 0 ) stop('r must be between 0 and 1')
  if(length(p)!=length(y)) stop('length of p does not match length y')
  if(lambda < 0 | delta <0 ) stop('tuning parameters must be positive')

  log.p <- logit(p)

  ## weighted recalibrated rislambda
  # weights for weighted recalibration
  wt <- calWt(rl = rl,ru = ru,r = r,p = p,y = y,lambda = lambda,delta = delta)

  #coefficients from recalibration method
  mod.wt <- glm(y~log.p, weights = wt, family="quasibinomial")
  alpha.wt <- mod.wt$coef
  wt.conv <- ifelse(mod.wt$converged==TRUE,1,0) #did the model converge

  ## predicted rislambdas from weighted recalibraiton
  log.phat.wt <- alpha.wt[1] + (alpha.wt[2]*log.p)
  phat.wt <- logit.inv(log.phat.wt)

  res <- list("y" = y,
              "p.wt" = phat.wt,
              "alpha.wt" = alpha.wt,
              "wt" = wt,
              "wt.conv" = wt.conv)
  return(res)
}


###### 3: Tuning parameter selection function ########
#### 3a: Tuning Parameter Functions ####

## getting CV grid of lambda/delta's for sequence of RAW
RAWgrid <- function(r,rl,ru,p,y,rawSeq=seq(0.1,0.9,0.1),cvParm,delta=NULL,lambda=NULL,rl.raw,ru.raw){
  #library(lattice)
  if( rl > ru ) stop('lower bound(rl) must be less than upper bound (ru)')
  if( (rl > r) | (ru < r) ) stop('r must be between lower and upper bound')
  if( r > 1 | r < 0 ) stop('r must be between 0 and 1')
  if(any(rawSeq <0)  | any(rawSeq > 1) ) stop("RAW must be in [0,1]")
  if(!cvParm %in% c("lambda","delta")) stop("cvParm must be 'lambda' or 'delta' ")

   ### if searching for k (so delta is fixed)
  if(cvParm=="lambda"){
    #right now delta needs to be specifid (if NULL then just put something)
    deltaPick <- ifelse(is.null(delta),1,delta)
    rrWt.seq <- cbind(rawSeq,NA,deltaPick)
    colnames(rrWt.seq) <- c("RAW","lambda","delta")

    ## getting grid of lambda's for fixed delta to do CV over
    for(i in 1:nrow(rrWt.seq)){
      ## this try catch will run the through the loop without stopping for errors,
      ## but will print at the end if there is any error
      ## The error is because there are no smoothed event rates within the RAW interval
      ## (i.e. R probably too far from prevalence)
      tryCatch({
        rrWt.seq[i,2] <- stats::uniroot(f = getLambda,interval = c(1e-5,1e5),
                                 r = r, rl = rl, ru = ru, rrWt = rawSeq[i],
                                 rl.raw = rl.raw,ru.raw = ru.raw,
                                 p = p,y = y,delta = deltaPick,tol = 1e-06)$root

      }, error=function(e){})
    }
    if(length(rrWt.seq[,2])==sum(is.na(rrWt.seq[,2]))){
      stop("Event rate outside RAW interval, widen RAW interval")
    }
  }

  if(cvParm=="delta"){
    #right now k needs to be specifid (if NULL then just put something)
    lambdaPick <- ifelse(is.null(lambda),1,lambda)
    rrWt.seq <- cbind(rawSeq,lambdaPick,NA)
    colnames(rrWt.seq) <- c("RAW","lambda","delta")

    ## getting grid of k's for fixed delta to do CV over
    for(i in 1:nrow(rrWt.seq)){
      tryCatch({
        rrWt.seq[i,3] <- uniroot(f = getDelta,interval = c(1e-4,1e4),
                                 r = r, rl = rl, ru = ru, rrWt = rawSeq[i],
                                 rl.raw = rl.raw,ru.raw = ru.raw,
                                 p = p,y = y,lambda = lambdaPick)$root

      }, error=function(e){})
    }
    if(length(rrWt.seq[,2])==sum(is.na(rrWt.seq[,2]))){
      stop("Event rate outside RAW interval, widen RAW interval")
    }
  }
  #returns grid with selected tuning parameter sequence
  return(rrWt.seq)
}

#### 3b: Cross-validation Functions using sNB ####
#this depends on k or delta, but this comes from RAW value
cvWtTuning <- function(p,y,r,rl,ru,kFold=5,cvParm,tuneSeq,cv.seed=1111){

  if( rl > ru ) stop('lower bound(rl) must be less than upper bound (ru)')
  if( (rl > r) | (ru < r) ) stop('r must be between lower and upper bound')
  if( r > 1 | r < 0 ) stop('r must be between 0 and 1')
  if(!cvParm %in% c("lambda","delta")) stop("cvParm must be 'lambda' or 'delta' ")

  set.seed(cv.seed)
  flds <- caret::createFolds(y, k = kFold, list = TRUE, returnTrain = TRUE)


  mean.nb <- cbind(tuneSeq,NA)
  colnames(mean.nb) <- c(colnames(tuneSeq),"cv.snb")
  mean.nb <- as.data.frame(mean.nb)
  full <- matrix(nrow=nrow(tuneSeq),ncol=kFold)

  for(i in 1:nrow(tuneSeq)){
    lambda <- tuneSeq[i,"lambda"]
    delta <- tuneSeq[i,"delta"]

    nb.wt <- rep(NA,kFold)
    for(j in 1:kFold){
      #getting training labels
      train.index <- as.numeric(flds[[j]])

      #calculating weights in training set
      pTrain <- p[train.index]
      yTrain <- y[train.index]

      # weights for weighted recalibration
      wt <- calWt(rl = rl, ru = ru,p = pTrain,y = yTrain,r = r,lambda = lambda,delta = delta)
      #calculating weighted recalibrated slope and intercept
      lpTrain <- logit(pTrain)
      mod.wt <- glm(yTrain~lpTrain, weights = wt, family="quasibinomial")
      alpha.wt <- mod.wt$coef

      #estimating recalibration risk score in test fold
      yTest <- y[-train.index]
      pTest <- p[-train.index]
      lpTest <- logit(pTest)

      lpHat.wt <- alpha.wt[1] + (alpha.wt[2]*lpTest)
      pHat.wt <- logit.inv(lpHat.wt)

      nb.wt[j] <- nb(y = yTest,p = pHat.wt,r = r)$snb
    }
    full[i,] <- nb.wt
    mean.nb[i,"cv.snb"] <- mean(nb.wt)
  }
  tune.index <- which.max(mean.nb$cv.snb)
  tune.cv <- mean.nb[tune.index,1]

  res <- list("cv.res" = mean.nb, "cv.param"=tune.cv,"cv.full" = full)
  return(res)
}

# repeat the CV procedure with different sample splittings
cvRepWtTuning <- function(y,p,r,rl,ru,kFold=5,cvRep=25,cvParm,tuneSeq,stdErrRule=TRUE,int.seed=11111){

  if( rl > ru ) stop('lower bound(rl) must be less than upper bound (ru)')
  if( (rl > r) | (ru < r) ) stop('r must be between lower and upper bound')
  if( r > 1 | r < 0 ) stop('r must be between 0 and 1')
  if(!cvParm %in% c("lambda","delta")) stop("cvParm must be 'lambda' or 'delta' ")

  repCV.sNB <- matrix(NA,nrow=nrow(tuneSeq),ncol=cvRep)
  fullList <- vector("list", length = cvRep)

  for (j in 1:cvRep){
    cv.seed <- int.seed + j
    int <- cvWtTuning(y = y,p = p,
                            kFold = kFold,r = r,rl = rl,ru = ru,
                            tuneSeq = tuneSeq, cv.seed = cv.seed,cvParm = cvParm)

    repCV.sNB[,j] <- int$cv.res$cv.snb
    fullList[[j]] <- int$cv.full
  }
  avgCV <- apply(repCV.sNB,1,mean)

  if(stdErrRule==FALSE){
    tune.index <- which.max(avgCV)
    avgCV.res <- cbind(tuneSeq,avgCV)
    max.index <- which.max(avgCV.res[,"avgCV"])

    tune.cv <- avgCV.res[max.index,"avgCV"]
    tune.RAW <- avgCV.res[max.index,"RAW"]
    tune.lambda <- avgCV.res[max.index,"lambda"]
    tune.delta <- avgCV.res[max.index,"delta"]

    res <- list("cv.sNB"=tune.cv, "cv.RAW" = tune.RAW, "cv.lambda" = tune.lambda,
                "cv.delta" = tune.delta,
                "avgCV.res" = avgCV.res, "fullList" = fullList)

  }
  else if(stdErrRule==TRUE){
       sigma.m.lambda <- matrix(unlist(lapply(X = fullList,FUN = function(y){
                                      apply(y,MARGIN = 1,var)
                                      }
                                 )
                          ),ncol = cvRep,byrow = FALSE)
       W <- apply(sigma.m.lambda,1,mean)
       B <- apply(repCV.sNB,MARGIN = 1, var)
       V <- (((kFold - 3)/kFold) * W) +  ((1/cvRep)*B)

       SE <- sqrt(V)
      ## estimates of std err
      avgCV.res <- cbind(tuneSeq,avgCV,SE)

      #bounds for within on SE of maximum
      max.index <- which.max(avgCV.res[,"avgCV"])
      cvMaxUpp <- avgCV.res[max.index,"avgCV"] + avgCV.res[max.index,"SE"]
      cvMaxLow <- avgCV.res[max.index,"avgCV"] - avgCV.res[max.index,"SE"]

      maxRAW <- max(avgCV.res[avgCV.res[,"avgCV"] <= cvMaxUpp & avgCV.res[,"avgCV"] >= cvMaxLow,"RAW"])
      tune.index <- which(avgCV.res[,"RAW"]==maxRAW)

      tune.cv <- avgCV.res[tune.index,"avgCV"]
      names(tune.cv) <- "cv.sNB"
      tune.RAW <- avgCV.res[tune.index,"RAW"]
      tune.lambda <- avgCV.res[tune.index,"lambda"]
      tune.delta <- avgCV.res[tune.index,"delta"]


      colnames(avgCV.res) <- c("RAW","lambda","delta","cvSNB","SE")
      res <- list("cv.sNB"=tune.cv, "cv.RAW" = tune.RAW, "cv.lambda" = tune.lambda,
                "cv.delta" = tune.delta,
                  "avgCV.res" = avgCV.res,
                  "W" = W, "B" = B,"fullList" = fullList)

  }

  return(res)
}



###### 4: Evaluation Functions ########
nb <- function(y,p,r){
  y.bar <- mean(y)
  TPR <- mean(p[y==1]>r)
  FPR <- mean(p[y==0]>r)

  #net benefit
  nb  <- (TPR*y.bar) - ((r/(1-r)) * (1-y.bar) *FPR)
  snb <- ((TPR*y.bar) - ((r/(1-r)) * (1-y.bar) *FPR))/y.bar

  #treat all
  snb.all <- 1 -  ((r/(1-r)) * (1-y.bar) *1)/y.bar

  res <- list("nb"=nb,"snb"=snb,"snb.all"=snb.all)
  return(res)
}





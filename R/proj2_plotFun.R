##### Plotting Functions #######
utils::globalVariables(c("..count.."))


### sNB Recalibration Curve with std error bars


# plotting snb as function of t -- still need to add other potential std error bars
snbRecalPlot <- function(p,p.std,y,r,stdErrThresh=1,ylim=NULL,
                          titlePlot = "Potential sNB Under Recalibration",
                          risk.model.std=TRUE){
  ## sNB recalibration Plot
  t.vec <- seq(0,1,0.0005)
  sNB <- cbind(t.vec,NA)
  for(i in 1:length(t.vec)){
    pick.t <- t.vec[i]
    sNB[i,2] <- snb.t(par= pick.t,y = y,p = p,r = r)
  }

  t.max <- NULL
  t.max$maximum <- sNB[which.max(sNB[,2]),1]
  t.max$objective <- sNB[which.max(sNB[,2]),2]

  sNB.max.se <- snbVar.tmax(tVec = t.max$maximum,y = y,p = p,r = r)
  upp <- t.max$objective + stdErrThresh*sNB.max.se
  low <- t.max$objective - stdErrThresh*sNB.max.se

  ## points marking orig and std recal sNB
  snb.orig <- nb(y = y,p = p,r = r)$snb
  snb.recal <- nb(y = y,p = p.std,r = r)$snb

  if(is.null(ylim)){
    ylim = c(min(c(snb.orig,snb.recal,0)),max(c(snb.orig,snb.recal,0.8)))
  }
  plot(sNB[,1],sNB[,2],type="l",col="black",lwd=2,ylim=ylim,
       xlab="Threshold (t) for Decision Rule",ylab="sNB",
       main=titlePlot)

  #standard error bars for points
  snb.t.orig <- sNB[which.min(abs(snb.orig-sNB[,2])),1]
  sNB.t.orig.se <- snbVar.tmax(tVec = snb.t.orig,y = y,p = p,r = r)

  snb.t.std <- sNB[which.min(abs(snb.recal-sNB[,2])),1]
  sNB.t.std.se <- snbVar.tmax(tVec = snb.t.std,y = y,p = p,r = r)


  points(snb.t.std,
         sNB[which.min(abs(snb.recal-sNB[,2])),2],col="blue",pch=1,cex=1.3,lwd=3)
  points(snb.t.orig,
         sNB[which.min(abs(snb.orig-sNB[,2])),2],col="red",pch=1,cex=1.2,lwd=3)
  #abline(h=upp,lwd=1,col="black",lty=c(2,3,4))
  abline(h=low,lwd=1,col="black",lty=c(2,3,4))
  if(risk.model.std==TRUE){
    arrows(x0 = snb.t.std,y0 = sNB[which.min(abs(snb.recal-sNB[,2])),2] - sNB.t.std.se,
           x1 = snb.t.std,y1 = sNB[which.min(abs(snb.recal-sNB[,2])),2] + sNB.t.std.se, angle = 90,length = 0.1,lwd=1.5,code = 3,col="blue")

    arrows(x0 = snb.t.orig,y0 = sNB[which.min(abs(snb.orig-sNB[,2])),2] - sNB.t.orig.se,
           x1 = snb.t.orig,y1 = sNB[which.min(abs(snb.orig-sNB[,2])),2] + sNB.t.orig.se, angle = 90,length = 0.1,lwd=1.5,code = 3,col="red")
  }
  legend("topleft",paste("Max(sNB) =",round(t.max$objective,3)),bty="n")
  legend("topright",c("Orig Risk Model","Std. Log. Recal. Risk Model",
                      paste(stdErrThresh,"Std Err from Maximum")),
         col=c("red","blue","black"),pch=c(1,1,NA),lwd=c(1.5,1.5),lty = c(NA,NA,2),bty = "n")
}


### calibration plots and histogram
calCurvPlot <- function(y,p,p.std=NULL,p.recal=NULL,stdPlot=FALSE,recalPlot=FALSE,
                         xlim=c(0,1),ylim=c(0,1),
                         label="Original Risk Score",
                         label2 = "Standard Recalibrated Risk Score",
                         label3 = "Weighted/Constrained Recalibrated Risk Score",
                         legendLab = c("Orig.", "Std.", "Wt."),
                         mainTitle="Calibration of Risk Score",
                         hist=TRUE,ylimHist = c(0,0.5),
                         r,rl = -Inf, ru = Inf){

    ### wont work in plot null value so put something outside plotting range
  orig.loess <- data.frame("x"=lowess(x = p,y = y,f = 2/3,iter = 0)$x,"y"=lowess(x = p,y = y,f = 2/3,iter = 0)$y)
  orig.loess$type <- "orig"
  hist.orig <- ggplot2::ggplot(orig.loess,ggplot2::aes(orig.loess$x)) +
    ggplot2::geom_histogram(binwidth = (xlim[2]-xlim[1])/20, ggplot2::aes(y = (..count..)/sum(..count..))) +
    ggplot2::labs(title =NULL, x = label,y="Percentage") +
    ggplot2::geom_vline(xintercept = r,linetype="dotted") +
    ggplot2::geom_vline(xintercept = rl,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
    ggplot2::geom_vline(xintercept = ru,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
    ggplot2::coord_cartesian(xlim=xlim, ylim=ylimHist)

  if(stdPlot==TRUE){
    stdCal.loess <- data.frame("x"=lowess(x = p.std,y = y,f = 2/3,iter = 0)$x,"y"=lowess(x = p.std,y = y,f = 2/3,iter = 0)$y)
    stdCal.loess$type <- "std"

    hist.std <- ggplot2::ggplot(stdCal.loess,ggplot2::aes(stdCal.loess$x)) +
      ggplot2::geom_histogram(binwidth = (xlim[2]-xlim[1])/20, ggplot2::aes(y = (..count..)/sum(..count..))) +
      ggplot2::labs(title =NULL, x = label2,y="Percentage") +
      ggplot2::geom_vline(xintercept = r,linetype="dotted")  +
      ggplot2::geom_vline(xintercept = rl,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::geom_vline(xintercept = ru,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::coord_cartesian(xlim=xlim, ylim=ylimHist)

  }
  else{stdCal.loess <- NULL}

  if(recalPlot==TRUE){
    wtCal.loess <-  data.frame("x"=lowess(x = p.recal,y = y,f = 2/3,iter = 0)$x,"y"=lowess(x = p.recal,y = y,f = 2/3,iter = 0)$y)
    wtCal.loess$type <- "wt"


    hist.wt <- ggplot2::ggplot(wtCal.loess,ggplot2::aes(wtCal.loess$x)) +
      ggplot2::geom_histogram(binwidth = (xlim[2]-xlim[1])/20, ggplot2::aes(y = (..count..)/sum(..count..))) +
      ggplot2::labs(title =NULL, x = label3,y="Percentage") +
      ggplot2::geom_vline(xintercept = r,linetype="dotted")  +
      ggplot2::geom_vline(xintercept = rl,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::geom_vline(xintercept = ru,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::coord_cartesian(xlim=xlim, ylim=ylimHist)

  }
  else(wtCal.loess <- NULL)

  loessDat <- as.data.frame(rbind(orig.loess,stdCal.loess,wtCal.loess))




  if(stdPlot==TRUE & recalPlot==TRUE){
    plot.cal <- ggplot2::ggplot(as.data.frame(loessDat),ggplot2::aes(.data$x,y = .data$y,group=.data$type,col=.data$type)) +
      ggplot2::geom_line() + ggplot2::coord_cartesian(xlim=xlim,ylim=ylim) +
      ggplot2::geom_abline(intercept=0,slope=1,linetype="dashed") +
      ggplot2::geom_vline(xintercept = r,linetype="dotted") +
      ggplot2::geom_hline(yintercept = r,linetype="dotted") +
      ggplot2::geom_vline(xintercept = rl,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::geom_vline(xintercept = ru,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::labs(title = mainTitle, x = "Predicted Risk", y = "Observed Event Rate") +
      ggplot2::scale_colour_discrete(name="Risk Score",
                            breaks=c("orig", "std", "wt"),
                            labels=legendLab)

    if(hist==TRUE){
      suppressWarnings(
        cowplot::plot_grid(plot.cal, NULL,
                hist.orig + ggplot2::geom_line(ggplot2::aes(x = p,y=y,color = "TEST"))
                + ggplot2::scale_color_manual(values = NA) + ggplot2::theme(legend.text = ggplot2::element_blank(), legend.title = ggplot2::element_blank()),
                NULL,
                hist.std + ggplot2::geom_line(ggplot2::aes(x = p,y=y,color = "TEST"))
                + ggplot2::scale_color_manual(values = NA) + ggplot2::theme(legend.text = ggplot2::element_blank(), legend.title = ggplot2::element_blank()),
                NULL,
                hist.wt + ggplot2::geom_line(ggplot2::aes(x = p,y=y,color = "TEST"))
                + ggplot2::scale_color_manual(values = NA) + ggplot2::theme(legend.text = ggplot2::element_blank(), legend.title = ggplot2::element_blank()),
                align = "hv",axis=1, ncol = 1,rel_heights = c(1,-0.2,0.6,-0.2,0.6,-0.2,0.6))
      )
    }
    else{
      suppressWarnings(print(plot.cal))
    }
  }
  else if(stdPlot==TRUE & recalPlot==FALSE){
    plot.cal <- ggplot2::ggplot(as.data.frame(subset(loessDat,subset = loessDat$type!="wt")),
                                ggplot2::aes(.data$x,y = .data$y,group=.data$type,col=.data$type)) +
      ggplot2::geom_line() + ggplot2::coord_cartesian(xlim=xlim,ylim=ylim) +
      ggplot2::geom_abline(intercept=0,slope=1,linetype="dashed") +
      ggplot2::geom_vline(xintercept = r,linetype="dotted") +
      ggplot2::geom_hline(yintercept = r,linetype="dotted") +
      ggplot2::geom_vline(xintercept = rl,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::geom_vline(xintercept = ru,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::labs(title = mainTitle, x = "Predicted Risk", y = "Observed Event Rate") +
      ggplot2::scale_colour_discrete(name="Risk Score",
                            breaks=c("orig", "std"),
                            labels=legendLab[1:2])
    if(hist==TRUE){
      suppressWarnings(
        cowplot::plot_grid(plot.cal,
                NULL,
                hist.orig + ggplot2::geom_line(ggplot2::aes(x = p,y=y,color = "TEST")) +
                  ggplot2::scale_color_manual(values = NA) + ggplot2::theme(legend.text = ggplot2::element_blank(), legend.title = ggplot2::element_blank()),
                NULL,
                hist.std + ggplot2::geom_line(ggplot2::aes(x = p,y=y,color = "TEST")) +
                  ggplot2::scale_color_manual(values = NA) + ggplot2::theme(legend.text = ggplot2::element_blank(), legend.title = ggplot2::element_blank()),
                align = "hv",axis=1, ncol = 1,rel_heights = c(1,-0.2,0.6,-0.2,0.6))
      )
    }
    else{suppressWarnings(print(plot.cal))}
  }

  else if(stdPlot==FALSE & recalPlot==TRUE){
    plot.cal <- ggplot2::ggplot(as.data.frame(subset(loessDat,subset = loessDat$type!="std")),
                                ggplot2::aes(.data$x,y = .data$y,group=.data$type,col=.data$type) ) +
      ggplot2::geom_line() + ggplot2::coord_cartesian(xlim=xlim,ylim=ylim) +
      ggplot2::geom_abline(intercept=0,slope=1,linetype="dashed") +
      ggplot2::geom_vline(xintercept = r,linetype="dotted") +
      ggplot2::geom_hline(yintercept = r,linetype="dotted") +
      ggplot2::geom_vline(xintercept = rl,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::geom_vline(xintercept = ru,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::labs(title = mainTitle, x = "Predicted Risk", y = "Observed Event Rate") +
      ggplot2::scale_colour_discrete(name="Recalibration\nType",
                            breaks=c("orig", "wt"),
                            labels=legendLab[c(1,3)])
    if(hist==TRUE){
      suppressWarnings(
        cowplot::plot_grid(plot.cal,
                NULL,
                hist.orig + ggplot2::geom_line(ggplot2::aes(x = p,y=y,color = "TEST")) +
                  ggplot2::scale_color_manual(values = NA) +theme(legend.text = ggplot2::element_blank(), legend.title = ggplot2::element_blank()),
                NULL,
                hist.wt + ggplot2::geom_line(ggplot2::aes(x = p,y=y,color = "TEST")) +
                  ggplot2::scale_color_manual(values = NA) + ggplot2::theme(legend.text = ggplot2::element_blank(), legend.title = ggplot2::element_blank()),
                align = "hv",axis=1, ncol = 1,rel_heights = c(1,-0.2,0.6,-0.2,0.6))
      )
    }
    else(suppressWarnings(print(plot.cal)))

  }
  else if(stdPlot==FALSE & recalPlot==FALSE){
    plot.cal <- ggplot2::ggplot(as.data.frame(subset(loessDat,subset = loessDat$type=="orig")),
                                ggplot2::aes(.data$x,y = .data$y,group=.data$type,col=.data$type)) +
      ggplot2::geom_line() + ggplot2::coord_cartesian(xlim=xlim,ylim=ylim) +
      ggplot2::geom_abline(intercept=0,slope=1,linetype="dashed") +
      ggplot2::geom_vline(xintercept = r,linetype="dotted") +
      ggplot2::geom_hline(yintercept = r,linetype="dotted") +
      ggplot2::geom_vline(xintercept = rl,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::geom_vline(xintercept = ru,linetype=ifelse(is.infinite(abs(rl)),NA,"dotdash")) +
      ggplot2::labs(title =" Calibration of Risk Score", x = "Predicted Risk", y = "Observed Event Rate") +
      ggplot2::scale_colour_discrete(name="Recalibration\nType",
                            breaks=c("orig"),
                            labels=c(label))
    if(hist==TRUE){
      suppressWarnings(
        cowplot::plot_grid(plot.cal,
                NULL,
                hist.orig + ggplot2::geom_line(ggplot2::aes(x = p,y=y,color = "TEST")) +
                  ggplot2::scale_color_manual(values = NA) + ggplot2::theme(legend.text = ggplot2::element_blank(), legend.title = ggplot2::element_blank()),
                align = "hv",axis=1, ncol = 1,
                rel_heights = c(1,-0.2,0.6)))
    }
    else(suppressWarnings(print(plot.cal)))
  }
}


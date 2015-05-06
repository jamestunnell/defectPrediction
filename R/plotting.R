ts.plot <- function(ts, fname, width = 8, height.each = 2, cex=1.5){
  postscript(file=fname, width=width, height=height.each*ncol(ts),
             onefile=TRUE, horizontal=FALSE)
  par(cex.lab=cex, cex.axis=cex, mfrow=c(ncol(ts),1))
  
  for(i in 1:ncol(ts)){
    y <- as.vector(ts[,i])
    n <- length(y)
    plot(time(ts), y, main=NULL, type="l", xlab = "", ylab = names(ts)[i])
  }
  
  dev.off()
  par(cex.lab=1, cex.axis=1, mfrow=c(1,1))
}

plot.predictions <- function(model, x.range, fname, width = 8, height = 3, cex = 1.5){
  y <- model$data$output
  n <- length(y)
  x <- x.range
  
  postscript(file=fname, width=width, height=height,
             onefile=TRUE, horizontal=FALSE)
  #par(cex.lab=cex, cex.axis=cex, cex=cex)
  plot(x, y,main=NULL, type="l", xlab="week", ylab=colnames(y))
  lines(x,model$estimates$pred, col="red", lty=2)
  #legend("topleft",legend=c("Actual","Fitted"),col=c("black","red"), lty=c(1,2))
  garbage <- dev.off()
  #par(cex.lab=1, cex.axis=1, cex=1)
}
library(dse)

ts.labels <- function(ts.data){
  names.out <- colnames(ts.data$output)
  names.in <- colnames(ts.data$input)
  labs <- list(bugs = names.out[pmatch("Bug",names.out)],
               imps = names.in[pmatch("Imp",names.in)],
               news = names.in[pmatch("Fea",names.in)])
  return(labs)
}

ts.extend.one <- function(ts.data){
  data.ext <- TSdata(
    output = ts.data$output,
    input = rbind(ts.data$input, matrix(c(0,0),nrow=1))
  )
  return(data.ext)
}


forecast.intervals <- function(model, ci){
  alpha <- 1 - ci
  z <- qnorm(1-alpha/2)
  residuals <- (model$estimates$pred - model$data$output)
  mse <- sum(residuals)^2 / length(residuals)
  return(sort(c(-z,z)*sqrt(mse)))
}

# Compute forecasts (w/ confidence intervals) for different
# hypothetical  values of nimps and nnews.
forecast.hypotheticals <- function(
  model.est, data.base, ci, imps.hypoth, news.hypoth){
  data.ext <- ts.extend.one(data.base)
  row.last <- nrow(data.ext$input)
  labs <- ts.labels(data.base)
  
  ci <- rev(sort(ci))
  d <- forecast.intervals(model, ci)
  nd <- length(d)
  nl <- nd + 1
  
  nimps <- length(imps.hypoth)
  nnews <- length(news.hypoth)
  x <- mat.or.vec(nimps*nnews,1)
  y <- mat.or.vec(nimps*nnews,1)
  z <- mat.or.vec(nimps*nnews,nl)
  for(i in 1:nimps){
    for(j in 1:nnews){
      k <- (i-1)*nnews + j
      x[k] <- imps.hypoth[i]
      y[k] <- news.hypoth[j]
      
      data.ext$input[row.last,c(labs$imps,labs$news)] <- c(x[k],y[k])
      fc <- forecast(TSmodel(model.est), data.ext)
      
      z_ <- fc$forecast[[1]][1,]
      lohi <- z_ + d
      for(l in 1:(nd/2)){
        z[k,l] <- lohi[l]
        z[k,nl-l+1] <- lohi[nd-l+1]
      }
      z[k,ceiling(nl/2)] <- z_
    }
  }
  cnames <- rep("",nl)
  for(l in 1:(nd/2)){
    cnames[l] <- paste0(100*ci[l],"% lo")
    cnames[nl-l+1] <- paste0(100*ci[l],"% hi")
  }
  cnames[ceiling(nl/2)] <- "mean"
  colnames(z) <- cnames
  return(list(x=x,y=y,z=z))
}

# forecast.hypotheticals.mean3d <- function(
#   model.est, data.base, imps.hypoth, news.hypoth, fname){
#   
#   data.ext <- ts.extend.one(data.base)
#   row.last <- nrow(data.ext$input)
#   labs <- ts.labels(data.base)
#   
#   x <- NULL; y <- NULL; z <- NULL;
#   for(ir in imps.hypoth){
#     for(nr in news.hypoth){
#       data.ext$input[row.last,c(labs$imps,labs$news)] <- c(ir,nr)
#       fc <- forecast(TSmodel(model.est), data.ext)
#       x <- append(x,ir)
#       y <- append(y,nr)
#       z <- append(z, fc$forecast[[1]][1,])
#     }
#   }
#   library(onion)
#   png(filename = fname, width=800, height=800)
#   p3d(x=x,y=y,z=z,d0 = 1,
#       xlab=labs$imps, ylab=labs$news, zlab=labs$bugs,
#       theta = -60, phi=30, ticktype = "detailed")
#   garbage <- dev.off()
#   return(list(x=x,y=y,z=z))
# }

# # Plot forecast (w/ confidence intervals) for different
# # hypothetical  values of nimps and nnews. Two side-by-side
# # plots are produced, one using hypothetical nimps and actual
# # nnews, and vice-versa.
# forecast.hypotheticals <- function(
#   model.est, data.base, ci, fname, imps.hypoth, news.hypoth){
#   data.ext <- ts.extend.one(data.base)
#   row.last <- nrow(data.ext$input)
#   labs <- ts.labels(data.base)
#   
#   ci <- rev(sort(ci))
#   d <- forecast.intervals(model, ci)
#   nd <- length(d)
#   nl <- nd + 1
#   nnews <- length(news.hypoth)
#   nimps <- length(imps.hypoth)
#   
#   png(filename = fname, width=800, height=400)
#   par(mfrow=c(1,2))
#   
#   x <- imps.hypoth
#   y <- mat.or.vec(nimps,nl)
#   for(i in 1:nimps){
#     data.ext$input[row.last,c(labs$imps,labs$news)] <- c(x[i],0)
#     fc <- forecast(TSmodel(model.est), data.ext)
#     y_ <- fc$forecast[[1]][1,]
#     lohi <- y_ + d
#     for(j in 1:(nd/2)){
#       y[i,j] <- lohi[j]
#       y[i,nl-j+1] <- lohi[nd-j+1]
#     }
#     y[i,ceiling(nl/2)] <- y_
#   }
#   plot(x,y[,ceiling(nl/2)], type="l", ylim=range(y), xlab=labs$imps, ylab=labs$bugs,
#        main=paste(labs$news,"=",0))
#   for(j in 1:(nd/2)){
#     lines(x,y[,j], lty=j+1)
#     lines(x,y[,nl-j+1], lty=j+1)
#   }
#   legend("topleft", legend=c("forecast",paste0(100*ci,"% conf")), lty=1:ceiling(nl/2))
#   
#   x <- news.hypoth
#   y <- mat.or.vec(nnews,nl)
#   for(i in 1:nnews){
#     data.ext$input[row.last,c(labs$imps,labs$news)] <- rev(c(x[i],0))
#     fc <- forecast(TSmodel(model.est), data.ext)
#     y_ <- fc$forecast[[1]][1,]
#     lohi <- y_ + d
#     for(j in 1:(nd/2)){
#       y[i,j] <- lohi[j]
#       y[i,nl-j+1] <- lohi[nd-j+1]
#     }
#     y[i,ceiling(nl/2)] <- y_
#   }
#   plot(x,y[,ceiling(nl/2)], type="l", ylim=range(y), xlab=labs$news, ylab=labs$bugs,
#        main=paste(labs$imps,"=",0))
#   for(j in 1:(nd/2)){
#     lines(x,y[,j], lty=j+1)
#     lines(x,y[,nl-j+1], lty=j+1)
#   }
#   legend("topleft", legend=c("forecast",paste0(100*ci,"% conf")), lty=1:ceiling(nl/2))
#   
#   garbage <- dev.off()
#   par(mfrow=c(1,1))
# }
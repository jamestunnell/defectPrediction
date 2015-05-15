library(dse)
library(xts)
library(car)

#' Pre-modeling: issue loading, sampling, and stationarity testing.
#' 
#' @author James Tunnell
#' 
#' @param issues.file A text file, containing CSV-like table, with software issue data.
#' @param sampling.period Sampling period length, in days.
#' @param ndiff Number of differences to take if needed
#' @param start.date The date to start sampling time series data. If left NULL, then sampling begins as early as possible.
#' @param end.date The date to stop sampling time series data. If left NULL, then sampling ends as late as possible.
#' @param our.dir Path to a directory where stationarity report(s) can be saved as files
#' 
#' @note
#' The issues file that is being used should be a text table, and loadable by
#' \code{read.table(issues.file, header = T)} alone. The table should have 
#' columns for: created, resolved, and type. The type must be one of: "bug", 
#' "improvement", or "newfeature". The created and resolved fields are for 
#' date-time string that can be interpreted using \code{timeDate(mystr)}.
#' 
#' @examples
#' \dontrun{
#' ts.data <- pre.modeling("~/issues.txt", 7, start.date = "2002-01-01",
#'  out.dir="~/testrun")
#' }
#' 
#' @export
pre.modeling <- function(issues.file, sampling.period, ndiff = 1, 
                         start.date = NULL, end.date = NULL, out.dir = NULL){
  issues <- read.table(issues.file, header = T)
  if(is.null(start.date) & is.null(end.date)){
    s <- sample.issues.all(issues, sampling.period)  
  } else if(is.null(end.date)){
    s <- sample.issues.from(issues, sampling.period, start.date)
  } else {
    s <- sample.issues.until(issues, sampling.period, end.date)
  }
  
  ts <- as.xts(data.frame(Bugs=s$bugs, Improvements=s$imps, Features=s$news), s$date)
  
  if(!is.null(out.dir)){
    # cat("Plotting time-series\n")
    fname <- file.path(out.dir, "time_series.eps")
    ts.plot(ts, fname)    
  }
  
  if(!is.null(out.dir)){
    fname <- file.path(out.dir, "stationarity.txt")
    cat("", file = fname, append = F)
  }
  
  ST_TYPE = "constant"
  needs.diffed = F
  for(i in 1:ncol(ts)){
    st <- test.stationarity(ts[,i], type = ST_TYPE, df.level = 1, kpss.level = 10)
    if(!is.null(out.dir)){
      print.stationarity(st, names(ts)[i], fname)
    }
    if(!st$df$stationary | !st$kpss$stationary){
      needs.diffed <- T
    }
  }
  
  ts.diffed <- NULL
  ndiff = if(needs.diffed){ ndiff } else { 0 }
  if(needs.diffed){
    ts.diffed <- diff(ts)
    
    for(i in 1:ncol(ts)){
      names(ts.diffed)[i] <- paste(names(ts.diffed)[i],"(Difference)")
      st <- test.stationarity(ts.diffed[(1+ndiff):nrow(ts.diffed),i], type = ST_TYPE, df.level = 1, kpss.level = 10)
      if(!is.null(out.dir)){
        print.stationarity(st, names(ts.diffed)[i], fname)
      }
      stopifnot(st$df$stationary & st$kpss$stationary)
    }
    
    if(!is.null(out.dir)){
      #   cat("Plotting differenced time-series\n")
      fname <- file.path(out.dir, "time_series_diff.eps")
      ts.plot(ts.diffed[ndiff:nrow(ts),], fname)
      #   cat("\n")
    }
  }
  
  return(list(ts = ts,ndiff = ndiff))
}

#' Time series modeling methodology over a sliding window.
#' 
#' @author James Tunnell
#' 
#' @description
#' A modeling regime for setting up software defect prediction model, built using
#' issue tracking system data. Performs the modeling methodology repeatedly over 
#' a sliding sample window.
#' 
#' @param ts.data An xts time series
#' @param window.size The number of samples to include in a window for modeling.
#' @param ndiff The number of differences to take. Leave at 0 if no differencing is needed.
#' @param conf.levels Confidence level(s) to use in testing model forecast performance. 
#' Should be greater than 0 and less than 100.
#' @param our.dir Path to a directory where plots can be saved as files
#' @param verbose If TRUE, extra info is printed
#' 
#' @return
#' The value returned is a vector with named elements, with each element being
#' the percent of sample windows where model forecasts were within confidence level,
#' and each element name being the confidence level.
#' 
#' @examples
#' \dontrun{
#' ts.data <- pre.modeling("~/issues.txt", 7, start.date = "2002-01-01",
#'  out.dir="~/testrun")
#' result <- model.regime(ts.data, 48, ndiff = 1, out.dir="~/testrun")
#' }
#' @export
model.regime <- function(ts.data, window.size, ndiff=0, 
                         conf.levels=c(75,90), out.dir=NULL, verbose = FALSE){  
  ts <- ts.data
  if(ndiff > 0){
    ts.diffed <- diff(ts, differences = ndiff)
  }
  
  labs <- list(bugs = names(ts)[pmatch("Bug",names(ts))],
               imps = names(ts)[pmatch("Imp",names(ts))],
               news = names(ts)[pmatch("Fea",names(ts))])
  
  ci <- rev(sort(conf.levels)) / 100
  ci.inout <- mat.or.vec(length(ci),2)
  rownames(ci.inout) <- paste0(100*ci,"% conf")
  colnames(ci.inout) <- c("in","out")
  fc.errs <- NULL
  
  for(w.start in (1+ndiff):(nrow(ts)-window.size)){
    s.min <- w.start
    s.max <- w.start + window.size - 1
    s.range <- s.min:s.max
    ts.sub <- ts[s.range]
    
    if(verbose){
      #   cat("=========================================\n")
      #   cat("        Modeling samples", s.min, "to", s.max, "\n")
      #   cat("=========================================\n\n")
      cat(s.min, "to", s.max, ":")      
    }
    
    ts.data <- TSdata(
      output = as.matrix(ts.sub[,labs$bugs]),
      input = as.matrix(ts.sub[,c(labs$imps,labs$news)])
    )
    model <- modeling.methodology(ts.data, verbose = F)
    if(is.null(model)){
      if(verbose){
        cat("No valid models found for this sample range. Skipping.\n")
      }
      next
    }
    
    if(!is.null(out.dir)){
      #   cat("Plotting one-step ahead predictions\n")
      fname <- file.path(out.dir, paste0("one-step_predictions_", s.min, "-", s.max, ".eps"))
      plot.predictions(model, s.range, fname, width = 8, height = 3, cex = 1.35)
    }
    #   cat("=========================================\n")
    #   cat("             Forecasting\n")
    #   cat("=========================================\n\n")
    
    tmp1 <- quantile(ts[s.range,labs$imps], probs = c(0.25,0.75), type=6)
    tmp2 <- quantile(ts[s.range,labs$news], probs = c(0.25,0.75), type=6)
    
    # Used for hypothetical forecasting. Not differenced values!
    # They will be converted if needed
    imps.hypoth <- seq(from = tmp1[['25%']], to = tmp1[['75%']], by = 2)
    news.hypoth <- seq(from = tmp2[['25%']], to = tmp2[['75%']], by = 1)
    
    imps.actual <- as.integer(ts[s.max+1,labs$imps])
    news.actual <- as.integer(ts[s.max+1,labs$news])
    bugs.actual <- as.integer(ts[s.max+1,labs$bugs])
    if(verbose){
      cat("actual imps, news, and bugs:", imps.actual, news.actual, bugs.actual)
    }
    
    if(!(imps.actual %in% imps.hypoth)){
      imps.hypoth <- sort(append(imps.hypoth, imps.actual))
    }
    if(!(news.actual %in% news.hypoth)){
      news.hypoth <- sort(append(news.hypoth, news.actual))
    }
    
    if(ndiff > 0){
      imps.hypoth <- as.integer(imps.hypoth - ts[s.max,labs$imps])
      news.hypoth <- as.integer(news.hypoth - ts[s.max,labs$news])
    }
    if(ndiff > 1){
      imps.hypoth <- as.integer(imps.hypoth - diff(ts)[s.max,labs$imps])
      news.hypoth <- as.integer(news.hypoth - diff(ts)[s.max,labs$news])
    }
    results <- forecast.hypotheticals(model, ts.data, ci=ci,
                                      imps.hypoth=imps.hypoth, news.hypoth=news.hypoth)
    x <- results$x; y <- results$y; z <- results$z
    if(ndiff > 1){
      x <- x + as.integer(diff(ts[,labs$imps])[s.max])
      y <- y + as.integer(diff(ts[,labs$news])[s.max])
      z <- z + as.integer(diff(ts[,labs$bugs])[s.max])
    }
    if(ndiff > 0){
      x <- x + as.integer(ts[s.max,labs$imps])
      y <- y + as.integer(ts[s.max,labs$news])
      z <- z + as.integer(ts[s.max,labs$bugs])
    }
    
    actual.at <- which(x == imps.actual & y == news.actual)
    fc.mean <- z[actual.at, ceiling(ncol(z)/2)]
    fc.errs <- append(fc.errs, fc.mean - bugs.actual)
    
    for(i in 1:floor(ncol(z)/2)){
      ci.lo <- z[actual.at, i]
      ci.hi <- z[actual.at, ncol(z) - i + 1]
      if(bugs.actual >= ci.lo & bugs.actual <= ci.hi){
        ci.inout[i,"in"] <- 1 + ci.inout[i,"in"]
      } else {
        ci.inout[i,"out"] <- 1 + ci.inout[i,"out"]
      }
    }
    
    #   library(onion)
    #   fname <- file.path(out.dir, paste0("forecast_hypotheticals_mean3d_", s.min, "-", s.max, ".png"))
    #   png(filename = fname, width=800, height=800)
    #   p3d(x=x,y=y,z=z[,"mean"],d0 = 1,
    #       xlab="Improvements", ylab="Features", zlab="Bugs",
    #       theta = -120, phi=20, ticktype = "detailed")
    #   garbage <- dev.off()
    
    
    if(!is.null(out.dir)){
      fname <- file.path(out.dir, paste0("forecast_hypotheticals", s.min, "-", s.max, ".csv"))
      forecasts <- data.frame(imps=x, news= y)
      for(cname in colnames(z)){
        forecasts[[cname]] <- z[,cname]
      }
      write.table(forecasts, file = fname, row.names = F, sep = ",")
    }
    if(verbose){
      cat("\n")
    }
  }
  
  if(!is.null(out.dir)){
    fname <- file.path(out.dir, paste0("hist_forecast_errors.eps"))
    postscript(file=fname, width=800, height=600,
               onefile=TRUE, horizontal=FALSE)
    hist(fc.errs, xlab = "Error of forecast mean", main = "", breaks="FD")
    garbage <- dev.off()
    
    fname <- file.path(out.dir, paste0("qq_plot_forecast_errors.eps"))
    postscript(file=fname, width=800, height=600,
               onefile=TRUE, horizontal=FALSE)
    qqPlot(fc.errs)
    garbage <- dev.off()
  }
  retval <- ci.inout[,"in"] / (ci.inout[,"in"] + ci.inout[,"out"])
  return(retval)
}
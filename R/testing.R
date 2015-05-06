require(urca)

test.stationarity <- function(x, type, df.level, kpss.level)
{
  if(type == "constant"){
    df.type = "drift"; kpss.type = "mu"
  } else if (type == "trend"){
    df.type = "trend"; kpss.type = "tau"
  } else {
    cat("Unrecognized type '", type, "'\n")
    return
  }
  
  res.df <- ur.df(x,type=df.type)
  cval.df <- attr(res.df,"cval")
  stat.df <- attr(res.df,"teststat")
  
  df.signif <- list()
  df.stationary <- TRUE
  for(statname in colnames(stat.df)){
    sig.idx <- which(abs(stat.df[,statname]) > abs(cval.df[statname,]))
    if(any(sig.idx)){
      pct <- colnames(cval.df)[min(sig.idx)]
      if(df.level > as.numeric(sub("pct","",pct))){ df.stationary <- FALSE }
      df.signif[statname] <- paste0("<",pct)
    } else {
      pct <- tail(colnames(cval.df),1)
      df.stationary <- FALSE
      df.signif[statname] <- paste0(">",pct)
    }
  }
  
  res.kpss <- ur.kpss(x, type=kpss.type)
  cval.kpss <- attr(res.kpss,"cval")
  stat.kpss <- attr(res.kpss,"teststat")
  sig.idx.kpss <- which(stat.kpss < cval.kpss)
  if(any(sig.idx.kpss)){
    pct <- colnames(cval.kpss)[min(sig.idx.kpss)]
    kpss.stationary <- kpss.level <= as.numeric(sub("pct","",pct))
    kpss.signif <- paste0(">",pct)
  } else{
    kpss.stationary <- F
    kpss.signif <- paste0("<",tail(colnames(cval.kpss),1))
  }
  
  return(list(df=list(statistic = stat.df, stationary = df.stationary, 
                      signif = df.signif, result = res.df),
              kpss=list(statistic = stat.kpss, stationary = kpss.stationary,
                        signif = kpss.signif, result = res.kpss)))
}

print.stationarity <- function(stationarity.results, series.name, fname){
  s <- stationarity.results
  cat("Stationarity of", series.name, "\n", file = fname, append = T)
  cat("----------------------------------------------\n", file = fname, append = T)
  cat("Dickey Fuller test results:", file = fname, append = T)
  for(i in 1:ncol(s$df$statistic)){
    cname <- colnames(s$df$statistic)[i]
    cat(if(i > 1){ "," }, cname,"=", s$df$statistic[1,i],
        paste0("(", s$df$signif[[cname]], ")"), file = fname, append = T)
  }
  cat("\n", file = fname, append = T)
  cat("KPSS:",s$kpss$statistic, paste0("(",s$kpss$signif,")\n"), file = fname, append = T)
  cat("\n", file = fname, append = T)
}
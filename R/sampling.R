require(timeDate)
require(xts)

sample.issues.all <- function(issues, period){
  created <- timeDate(issues$created)
  resolved <- timeDate(issues$resolved)
  minmax <- range(created,resolved)
  s <- sample.issues.dateRange(issues,minmax[1],minmax[2],period)
  return(s)
}

sample.issues.dateRange <- function(issues, startdate, enddate, period){
  startdate <- timeDate(startdate)
  enddate <- timeDate(enddate)
  
  datebreaks <- align(c(startdate,enddate), by=paste0(period,"d"))
  period_days <- as.integer(datebreaks[2]-datebreaks[1], units="days")
  last_datebreak <- timeDate(as.Date(tail(datebreaks,1))+period_days)
  nperiods <- length(datebreaks)
  datebreaks <- append(datebreaks, last_datebreak)
  
  bugs.created <- NULL
  imps.resolved <- NULL
  news.resolved <- NULL
  tsks.created <- NULL
  dates <- NULL
  
  created <- timeDate(issues$created)
  resolved <- timeDate(issues$resolved)
  
  for(i in 1:(nperiods-1)){
    l <- datebreaks[i]
    r <- datebreaks[i+1]
    
    was_created <- created >= l & created < r
    was_resolved <- resolved >= l & resolved < r

    nbugs.created <- length(which(issues$type == "bug" & was_created))
    nimps.resolved <- length(which(issues$type == "improvement" & was_resolved))
    nnews.resolved <- length(which(issues$type == "newfeature" & was_resolved))
    ntsks.created <- length(which(issues$type == "task" & was_created))
    
    bugs.created <- append(bugs.created, nbugs.created)
    imps.resolved <- append(imps.resolved, nimps.resolved)
    news.resolved <- append(news.resolved, nnews.resolved)
    tsks.created <- append(tsks.created, ntsks.created)
    
    dates <- append(dates,as.Date(l))
  }
  
  return(data.frame(bugs=bugs.created, imps=imps.resolved,
                    news=news.resolved, tsks=tsks.created,
                    date=dates))
}
\name{model.regime}
\title{
Time series modeling regime for software defect prediction, using issue tracking system data.
}
\description{
The function starts by extracting time series data from the given issues file. It then evaluates the data for stationarity, differences if needed, and performs the modeling methodology repeatedly over a sliding sample window. The results over the sliding window are reported.
}
\usage{
model.regime(issues.file, sampling.period, window.size, ndiffs = 1, conf.levels = c(75,90) out.dir = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{issues.file}{
A text file, containing CSV-like table, with software issue data
}
  \item{sampling.period}{
Sampling period length, in days
}
  \item{window.size}{
The number of samples to include in a window for modeling
}
  \item{ndiffs}{
The number of differences to take, for non-stationary time series data
}
  \item{conf.levels}{
Confidence level(s) to use in testing model forecast performance
}
  \item{out.dir}{
Path to a directory where plots can be saved as files. Keep NULL if no output files are desired
}
  \item{verbose}{
    If TRUE, print actual future values for all time series variables, at each window.
  }
}
\details{

}
\value{
The value returned is a vector with named elements, with each element being the percent of sample windows where model forecasts were within confidence level, and each element name being the confidence level.}
\author{
James Tunnell
}
\note{
The issues file that is being used should be a text table, and loadable by \code{read.table(issues.file, header = T)} alone. The table should have columns for: created, resolved, and type. The type must be one of: "bug", "improvement", or "newfeature". The created and resolved fields are for date-time string that can be interpreted using \code{timeDate(mystr)}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

\dontrun{
  model.regime(issues.file = "issues.txt",
               sampling.period = 7,
               window.size = 50,
               ndiffs = 1,
               out.dir = "~/modeling"
               )
}
  
## The function is currently defined as
function (issues.file, sampling.period, window.size, ndiffs = 1, 
    out.dir = NULL) 
{
    issues <- read.table(issues.file, header = T)
    s <- sample.issues.all(issues, sampling.period)
    ts <- as.xts(data.frame(Bugs = s$bugs, Improvements = s$imps, 
        Features = s$news), s$date)
    fname <- file.path(out.dir, "time_series.eps")
    ts.plot(ts, fname)
    ST_TYPE = "constant"
    fname <- file.path(out.dir, "stationarity.txt")
    cat("", file = fname, append = F)
    needs.diffed = F
    for (i in 1:ncol(ts)) {
        st <- test.stationarity(ts[, i], type = ST_TYPE, df.level = 1, 
            kpss.level = 10)
        print.stationarity(st, names(ts)[i], fname)
        if (!st$df$stationary | !st$kpss$stationary) {
            needs.diffed <- T
        }
    }
    ndiffs = if (needs.diffed) {
        ndiffs
    }
    else {
        0
    }
    if (needs.diffed) {
        for (i in 1:ncol(ts)) {
            ts[, i] <- diff(ts[, i], differences = ndiffs)
            names(ts)[i] <- paste(names(ts)[i], "(Difference)")
            st <- test.stationarity(ts[(1 + ndiffs):nrow(ts), 
                i], type = ST_TYPE, df.level = 1, kpss.level = 10)
            print.stationarity(st, names(ts)[i], fname)
            stopifnot(st$df$stationary & st$kpss$stationary)
        }
        fname <- file.path(out.dir, "time_series_diff.eps")
        ts.plot(ts[2:nrow(ts)], fname)
    }
    ci <- c(0.9, 0.75)
    labs <- list(bugs = names(ts)[pmatch("Bug", names(ts))], 
        imps = names(ts)[pmatch("Imp", names(ts))], news = names(ts)[pmatch("Fea", 
            names(ts))])
    ci.inout <- mat.or.vec(length(ci), 2)
    rownames(ci.inout) <- paste0(100 * ci, "\% conf")
    colnames(ci.inout) <- c("in", "out")
    fc.errs <- NULL
    for (w.start in (1 + ndiffs):(nrow(ts) - window.size)) {
        s.min <- w.start
        s.max <- w.start + window.size - 1
        s.range <- s.min:s.max
        ts.sub <- ts[s.range]
        cat(s.min, "to", s.max, ":")
        ts.data <- TSdata(output = as.matrix(ts.sub[, labs$bugs]), 
            input = as.matrix(ts.sub[, c(labs$imps, labs$news)]))
        model <- modeling.methodology(ts.data, verbose = F)
        if (is.null(model)) {
            cat("No valid models found for this sample range. Skipping.\n")
            next
        }
        fname <- file.path(out.dir, paste0("one-step_predictions_", 
            s.min, "-", s.max, ".eps"))
        plot.predictions(model, s.range, fname, width = 8, height = 3, 
            cex = 1.35)
        tmp1 <- quantile(s$imps[s.range], probs = c(0.25, 0.75), 
            type = 6)
        tmp2 <- quantile(s$news[s.range], probs = c(0.25, 0.75), 
            type = 6)
        imps.hypoth <- seq(from = tmp1[["25\%"]], to = tmp1[["75\%"]], 
            by = 2)
        news.hypoth <- seq(from = tmp2[["25\%"]], to = tmp2[["75\%"]], 
            by = 1)
        imps.actual <- s$imps[s.max + 1]
        news.actual <- s$news[s.max + 1]
        bugs.actual <- s$bugs[s.max + 1]
        cat("actual imps, news, and bugs:", imps.actual, news.actual, 
            bugs.actual)
        if (!(imps.actual \%in\% imps.hypoth)) {
            imps.hypoth <- sort(append(imps.hypoth, imps.actual))
        }
        if (!(news.actual \%in\% news.hypoth)) {
            news.hypoth <- sort(append(news.hypoth, news.actual))
        }
        if (ndiffs > 0) {
            imps.hypoth <- imps.hypoth - s$imps[s.max]
            news.hypoth <- news.hypoth - s$news[s.max]
        }
        if (ndiffs > 1) {
            imps.hypoth <- imps.hypoth - diff(s$imps)[s.max]
            news.hypoth <- news.hypoth - diff(s$news)[s.max]
        }
        results <- forecast.hypotheticals(model, ts.data, ci = ci, 
            imps.hypoth = imps.hypoth, news.hypoth = news.hypoth)
        x <- results$x
        y <- results$y
        z <- results$z
        if (ndiffs > 1) {
            x <- x + diff(s$imps)[s.max]
            y <- y + diff(s$news)[s.max]
            z <- z + diff(s$bugs)[s.max]
        }
        if (ndiffs > 0) {
            x <- x + s$imps[s.max]
            y <- y + s$news[s.max]
            z <- z + s$bugs[s.max]
        }
        actual.at <- which(x == imps.actual & y == news.actual)
        fc.mean <- z[actual.at, ceiling(ncol(z)/2)]
        fc.errs <- append(fc.errs, fc.mean - bugs.actual)
        for (i in 1:floor(ncol(z)/2)) {
            ci.lo <- z[actual.at, i]
            ci.hi <- z[actual.at, ncol(z) - i + 1]
            if (bugs.actual >= ci.lo & bugs.actual <= ci.hi) {
                ci.inout[i, "in"] <- 1 + ci.inout[i, "in"]
            }
            else {
                ci.inout[i, "out"] <- 1 + ci.inout[i, "out"]
            }
        }
        fname <- file.path(out.dir, paste0("forecast_hypotheticals", 
            s.min, "-", s.max, ".csv"))
        forecasts <- data.frame(imps = x, news = y)
        for (cname in colnames(z)) {
            forecasts[[cname]] <- z[, cname]
        }
        write.table(forecasts, file = fname, row.names = F, sep = ",")
        cat("\n")
    }
    for (i in 1:nrow(ci.inout)) {
        frac <- ci.inout[i, "in"]/sum(ci.inout[i, ])
        cat(paste0(100 * round(frac, 4), "\%"), "in", rownames(ci.inout)[i], 
            "\n")
    }
    fname <- file.path(out.dir, paste0("hist_forecast_errors.eps"))
    postscript(file = fname, width = 800, height = 600, onefile = TRUE, 
        horizontal = FALSE)
    hist(fc.errs, xlab = "Error of forecast mean", main = "", 
        breaks = "FD")
    garbage <- dev.off()
    fname <- file.path(out.dir, paste0("qq_plot_forecast_errors.eps"))
    postscript(file = fname, width = 800, height = 600, onefile = TRUE, 
        horizontal = FALSE)
    qqnorm(fc.errs)
    qqline(fc.errs, distribution = qnorm)
    garbage <- dev.off()
  }
}

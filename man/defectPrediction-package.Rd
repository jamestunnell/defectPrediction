\name{defectPrediction-package}
\alias{defectPrediction-package}
\alias{defectPrediction}
\docType{package}
\title{
Predict software defects using a time series model.
}
\description{
For predicting defects in a software project using a time series modelling approach, given historical data from the project issue tracking system. A vector AR model with exogenous variables (VARX) model is used. The model structure has number of defect (created) as an endogenous variable, and number of new features (resolved) and number of improvements (resolved) as exogenous variables.
}
\details{
\tabular{ll}{
Package: \tab defectPrediction\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-05-05\cr
License: \tab Unlimited\cr
}
The primary function is \code{model.regime}. The function evaluates the time series data for stationarity, differences if needed, and performs the modeling methodology repeatedly over a sliding sample window. The results ofover the sliding window are reported.
}
\author{
James Tunnell

Maintainer: James Tunnell <jamestunnell@gmail.com>
}
\examples{
  \dontrun{
    model.regime(issues.file = "issues.txt",
                 sampling.period = 7,
                 window.size = 50,
                 ndiffs = 1,
                 out.dir = "~/modeling"
                 )
  }
}

\name{spectralGOF-package}
\alias{spectralGOF-package}
\alias{spectralGOF}
\docType{package}
\title{
Spectral goodness of fit for network models

}
\description{
A minimal package, including functions for calculating and visualizing the spectral goodness of fit statistic described in "Spectral Goodness of Fit for Network Models" by Jesse Shore and Benjamin Lubin (Social Networks Vol. 43, October 2015, Pages 16–27).


}
\details{
\tabular{ll}{
Package: \tab spectralGOF\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-08-21\cr
License: \tab GPL-2 \cr
}
The central function in this package is \code{SGOF}, which takes an observed network and a (fitted) model of network structure as arguments, and returns an absolute measure of the model's goodness of fit.  The goodness of fit is calculated on the basis of the how well the spectrum of the graph Laplacian of the observed networks are approximated by the Laplacian of the networks simulated from the fitted model.

Interpretation of the SGOF statistic is analogous in certain ways to interpretation of R^2 from linear regression.  In particular, SGOF=0 indicates no improvement over the null model, and SGOF=1 indicates a perfect model fit.  SGOF < 0 indicates that the null model is a better description of the data than the fitted model.

Additionally, \code{plotSGOFerrors()} is provided for visualizing the output of SGOF analysis.

Until the package is available on CRAN, it is best to check back to people.bu.edu/jccs occasionally for any updates.  
}
\author{
Jesse Shore (jccs@bu.edu)
%[blinded for peer review]
}
\references{
Shore, Jesse and Benjamin Lubin (2014). "Spectral Goodness of Fit for Network Models" Social Networks Vol. 43, October 2015, Pages 16–27
}
\keyword{ package }
\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
%% ~~ simple examples of the most important functions ~~
}

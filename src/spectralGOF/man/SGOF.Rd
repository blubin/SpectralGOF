\name{SGOF}
\alias{SGOF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spectral goodness of fit for network models
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SGOF(fittedERGM = NULL, networkModel = NULL, observedNetwork = NULL, observedSpectrum = NULL, nullDistribution = NULL, nsim = 250)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fittedERGM}{A fitted model of class \code{ergm}.  
%%     ~~Describe \code{fittedERGM} here~~
}
  \item{networkModel}{(if not using ergms) a function representing the fitted network model that can produce a simulated network in the form of an adjacency matrix. At present the function may not require additional arguments to be passed.
%%     ~~Describe \code{networkModel} here~~
}
  \item{observedNetwork}{(if not using ergms) a square, symmetric, adjacency matrix
%%     ~~Describe \code{observedNetwork} here~~
}
  \item{observedSpectrum}{(if not using ergms AND not supplying an \code{observedNetwork}) a vector of eigenvalues of the observed network, sorted into decreasing order.
%%     ~~Describe \code{observedSpectrum} here~~
}
  \item{nullDistribution}{(optionally used to save computation time; can be used with \code{fittedERGM} and \code{networkModel}) An array with some number of rows and \code{n} columns (where \code{n} is the number of nodes in the observed network).  Each row in the array should correspond to the spectrum of a single network generated from the fitted model.  
  
  To save computation time, the \code{nullDistribution} can be extracted from a previous \code{SGOF} or \code{nullSGOF} analysis: \code{mySGOFoutputObject$allNullSpectra}. 
  
%%     ~~Describe \code{nullDistribution} here~~
}
  \item{nsim}{\code{nsim} is relevant for all input types and should be >=1000 for final analyses.  NB if a \code{nullDistribution} is supplied, the number of simulations (rows in \code{nullDistribution}) overrides this argument.
%%     ~~Describe \code{nsim} here~~
}
}
\details{
There are multiple ways to supply the necessary information for this function to operate.

Using ERGMs: If a \code{fittedERGM} is supplied, the only other arguments to specify are \code{nsim} and, optionally, \code{nullDistribution}.

Using other models: If a \code{networkModel} (a function to simulate networks according to a fitted model) is supplied, either \code{observedNetwork} or \code{observedSpectrum} must also be supplied.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A named list with the following elements:
\item{SGOF }{The main result: the mean spectral goodness of fit of the fitted model}
\item{SGOF5 }{The fifth percentile spectral goodness of fit}
\item{SGOF95 }{The ninety-fifth percentile spectral goodness of fit}
\item{observedSpectrum }{The spectrum of the observed network}
\item{allSimSpectra }{an \code{nsim} x \code{n} array containing the complete set of normalized spectra of networks simulated from the fitted model.  Each row is the spectrum from one simulated network.}  
\item{allNullSpectra }{an \code{nsim} x \code{n} array containing the complete set of normalized spectra of networks simulated from the null model.  Each row is the spectrum from one simulated network. OR, if \code{nullDistribution} was supplied, it is simply the same array, returned.} 
\item{nsim  }{the number of simulations that was actually used to calculate the \code{SGOF}. See the description of \code{nsim} in the Arguments section, above. }

There is an S3 \code{print} method defined that simply prints out the mean SGOF with the 5th and 95th percentile results in parentheses.

After running \code{SGOF}, one will usually want to \code{\link{plotSGOFerrors}} on the output to inspect the spectral goodness of fit visually.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Jesse Shore (jccs@bu.edu)
%[blinded for peer review]
}
\references{
Shore, Jesse and Benjamin Lubin (2015). "Spectral Goodness of Fit for Network Models" Social Networks, Vol. 43, pp 16-27
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plotSGOFerrors}}, \code{\link{summarySGOF}},\code{\link{SGOFnull}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

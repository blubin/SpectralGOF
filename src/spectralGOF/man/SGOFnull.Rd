\name{SGOFnull}
\alias{SGOFnull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the null distribution for the SGOF.
%%  ~~function to do ... ~~
}
\description{This function is essentially a more limited version of \code{SGOF}.  \code{SGOFnull} is useful when one wishes to calculate the 5th and 95th percentile values for the null model's SGOF, without having to analyze a fitted model. 


%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SGOFnull(observedNetwork, nsim = 250)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{observedNetwork }{a square adjacency matrix representing the observed network
%%     ~~Describe \code{observedNetwork} here~~
}
  \item{nsim }{the number of networks to simulate from the null model
%%     ~~Describe \code{nsim} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{nullSGOF }{always 0}
\item{allNullSpectra }{an \code{nsim} x \code{n} array containing the complete set of normalized spectra of networks simulated from the null model.  Each row is the spectrum from one simulated network.} 
\item{nullSGOF5 }{The fifth percentile spectral goodness of fit for the null model}
\item{nullSGOF95 }{The ninety-fifth percentile spectral goodness of fit for the null model}
    %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{plotSGOFerrors}
\alias{plotSGOFerrors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visualize errors under null and fitted models for SGOF analysis.
%%  ~~function to do ... ~~
}
\description{Creates visualizations of the spectral goodness of fit of network models.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotSGOFerrors(sgofobject, includeLegend=TRUE, style="flat", insetForFlat=TRUE, insetLoc=c(0.3,0.7,0.025,0.425))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sgofobject}{an SGOF object, generated as the output of the function \code{SGOF}.
%%     ~~Describe \code{sgofobject} here~~
}
  \item{includeLegend}{Whether to generate a legend for the plot.
  

%%     ~~Describe \code{sgofobject} here~~
}
  \item{style}{Choose the style of plot.  Possible values are \code{"original"},\code{"flat"}, and \code{"logflat"}}
  \item{ptcex}{set the point size when plotting spectra.}
  \item{insetForFlat}{whether to include an inset plot giving the original shape of the spectrum}
  \item{insetLoc}{Coordinates of the inset, if using}
}
\details{
The distance between an eigenvalue of the observed network and the corresponding eigenvalue from a simulated network (either from the fitted or null models) is an "error."  This function produces a representative visualization of spectral error based on three spectra: the observed spectrum, the null model spectrum that is closest to the mean Euclidean distance from the observed spectrum, and a fitted model spectrum that is closest to the mean Euclidean distance from the observed spectrum. 

The spectrum of the observed network is plotted, along with colored lines representing 1) error explained by the fitted model (in green) 2) error from the null model that remain as error under the fitted model (in blue) and 3) error introduced by the fitted model (in red).

Additionally, when includeLegend=TRUE, a legend is produced automatically, located either on the top right or bottom left of the plot, depending on the shape of the spectrum.

The default style is "flat" but "original" is also frequently useful.  In the "flat" style, all spectra are plotted after subtracting the spectrum of the observed network.  "original" preserves the original shape of the spectra (ie. does not subtract the spectrum of the observed network).  Since the shape of the spectrum is itself informative, this style is sometimes preferred.  The "logflat" is less frequently called for, but can be useful if one wants to "zoom in" on the small spectrall errors. In general, however, the "flat" or "original" styles should be used.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
Shore, Jesse and Benjamin Lubin (2015). "Spectral Goodness of Fit for Network Models" Social Networks, Volume 43, October 2015, Pages 16-27
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SGOF}}.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

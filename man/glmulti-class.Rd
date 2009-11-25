\name{glmulti-class}
\docType{class}
\alias{glmulti-class}


\title{Class "glmulti"}
\description{Contains the results of a \code{glmulti} analysis.}
\section{Objects from the Class}{
Objects will never be created directly but through calls of \code{glmulti} or by applying
\code{consensus} on a list of \code{glmulti} objects.

}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{Object of class \code{"character"} : the name of the analysis. }
    \item{\code{params}:}{Object of class \code{"list"} : parameter values used when calling \code{glmulti} to produce the object.}
    \item{\code{nbmods}:}{Object of class \code{"integer"} : the number of models that have been found by \code{glmulti}. }
    \item{\code{crits}:}{Object of class \code{"numeric"} : the IC values of the models found, in ascending order. }
    \item{\code{codes}:}{Object of class \code{"data.frame"}: for each model (rows), a description of which terms are present/absent in its formula. }
    \item{\code{K}:}{Object of class \code{"integer"} : for each model, its complexity, from best to worst model. }
    \item{\code{formulas}:}{Object of class \code{"list"} : for each model, its formula, from best to worst. }
    \item{\code{call}:}{Object of class \code{"call"} : the original \code{glmulti} call that produced this object. }
    \item{\code{adi}:}{Object of class \code{"list"} : additional arguments that had been passed to the fitting function through \code{glmulti}. }
  }
}
\references{ }
\author{Vincent Calcagno, McGill University}
\note{}


\seealso{
	 \code{\link{summary.glmulti}}, \code{\link{consensus}}, \code{\link{coef.glmulti}}

}
\examples{
showClass("glmulti")
}
\keyword{classes}

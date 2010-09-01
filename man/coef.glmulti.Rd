\name{coef.glmulti}
\alias{predict.glmulti}
\alias{coef.glmulti}

\title{Multimodel inference with glmulti}
\description{
These functions, applied on a \code{glmulti} object, produce model-averaged coefficients and predictions
from the multiple models in the confidence set (or a subset of them). This allows easy multi-model inference.
}
\usage{
# S3 coef method for class 'glmulti'
coef.glmulti(object, select="all", varweighting="Buckland", ...)

# S3 predict method for class 'glmulti'
predict.glmulti(object, select="all", newdata=NA, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{object}{an object of class \code{glmulti}}
	
	\item{select}{A specification of which models should be used for inference. By default all models are used, see below. }
	
	\item{varweighting}{The method to be used to compute the unconditional variance. "Buckland" implements the approach presented in Buckland et al. 1997. "Johnson" implements a slightly different approach
	 suggested in Johnson \& Omland 2004. The latter results in slightly smaller estimates of the unconditional variance of model coefficients.}
	
	\item{newdata}{New data for which to predict values}
	
   \item{\dots}{ Further arguments.}
}
\details{

\emph{select} can be used to specify which models should be used for inference. By default all are used. If specifying an integer value \emph{x}, only the \emph{x} best models are used.
If a numeric value is provided, if it less than one, models that sum up to \emph{x}\% of evidence weight are used.
If it more than one, models within \emph{x} IC units from the best model are used.

}
\value{


\code{coef} returns a data.frame with model-averaged estimates of the different parameters in the models, as well as their unconditional variance and their importance. 
 

\code{predict} returns a list of three data.frames: the multi-model predictions, their variability across models, and the predictions from each model separetely.


}
\references{ 
Buckland et al. 1997. Model selection: an integral part of inference. Biometrics.
Johnson \& Omland. 2004. Model selection in ecology and evolution. TREE.
 }
\author{Vincent Calcagno, McGill University}

\seealso{ 

\code{\link{glmulti}}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models }
\keyword{regression }


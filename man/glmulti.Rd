\name{glmulti}
\alias{glmulti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ GLM model selection and multimodel inference made easy}
\description{
  This is the main function of package glmulti. Use it instead of glm to automatically compare all possible models.
  Simply indicate the name of your response variable, and the candidate explanatory variables.
  This wrapper for glm() will automatically generate all possible models (under some constraints set by you) and find the best models in terms of some Information Criterion (AIC, AICc or BIC).
  It can handle very large numbers of candidate models, and for this reason	can use a Genetic Algorithm to find the best models when an exhaustive screening of the candidates is not feasible.
}
\usage{
glmulti(y, xr, data,  exclude=c(), intercept=TRUE, marginality=FALSE, level=2, filename="glmulti.output", method="h", crit="aicc", chunk=1, chunks=1, minsize=0, maxsize=-1, minK=0, maxK=-1,  plotty=TRUE, confsetsize=100, popsize=100, mutrate=10^-3,sexrate=0.1,imm=0.3, deltaM=0.05, deltaB=0.05, conseq=5, fitfunc=glm, resumefile = "id",  ...)
}
%- maybe also 'usage' for other objects documented here.


\arguments{
  \item{y}{The name of your dependent variable, as used in your data frame. E.g. "MyResponse" }
  \item{xr}{A vector of the names of your explanatory variables (factors and/or covariates). E.g. c("Age","Height","Sex")}
  \item{data}{ A data frame containing your data (observations as rows, variables as (named) columns. }
  \item{exclude}{This sets a constraint on candidate models. A vector of some terms (main effects and/or interactions) that should be excluded from the candidate models. Use ":" for interactions. E.g. c("Age:Sex","Height:Age") }
  \item{intercept}{Whether to include an intercept in the formulas or not. Default is TRUE, so all formulas will start with "~1"}
  \item{marginality}{Whether to use the general marginality rule or not. Default is FALSE. With TRUE, all predictors found in interaction terms are also included as main effects.}
  \item{level}{The level of interaction between explanatory variables to be considered. Currently only 1 (only main effects) or 2 (main effects plus all pairwise interactions) are supported.}
  \item{filename}{ A file name to be used when exporting the results. Extension ".txt" will be automatically appended. }
  \item{method}{What the function should do: "d" for a simple report of the number of candidate models, "h" for exhaustive screening of the candidates, "g" for a genetic algorithm search, "r" to resume a previous GA simulation.}
  \item{crit}{The Information Criterion to be used ("aic", "aicc" or "bic")}

  \item{chunk}{When splitting an exhaustive screening ("h") approach, this indicates which part (on the \code{chunks} parts) this program should do.}
  \item{chunks}{When splitting an exhaustive screening ("h") approach, this indicates the total number of parts (i.e. into how many pieces the job has been splitted.)}

  \item{minsize}{This sets a constraint on candidate models. Minimal number of TERMS (main effects or interactions) to be included in candidate models (negative = no constraint)}
  \item{maxsize}{This sets a constraint on candidate models. Maximal number of TERMS to be included in candidate models (negative = no constraint)}
  \item{minK}{This sets a constraint on candidate models. Minimal complexity of candidate models (negative = no constraint)}
  \item{maxK}{This sets a constraint on candidate models. Maximal complexity of candidate models (negative = no constraint)}
  
  \item{plotty}{whether to plot the progress of the IC profile when running. Default is TRUE.}
  \item{confsetsize}{ How many models should be returned in the confidence set of models?} 
  
  \item{popsize}{The population size for the genetic algorithm}
  \item{mutrate}{The per locus mutation rate for genetic algorithm, between 0 and 1.}
  \item{sexrate}{The rate of sexual reproduction for the genetic algorithm, between 0 and 1. }
  \item{imm}{The rate of immigration for the genetic algorithm, between 0 and 1.}
  
  \item{deltaM}{The target change in mean IC (defines the stop rules for the GA)}
  \item{deltaB}{The target change in best IC (defines the stop rules for the GA)}
  \item{conseq}{The target successive number of times with no improvement (i.e. target changes have been attained) (defines the stop rule for the GA). The greater it is, the more stringent the stop rule. }
  
  \item{fitfunc}{The fitting function to bee used. Default is glm of course, but you can provide any custom function that matches the parameters/returned values of glm. }
  
  \item{resumefile}{The name of files (without extension) from which to restore the java objects, when using method "r". Default: taken to be filename.}
    
  \item{\dots}{Further options to be passed to the fitting function. E.g. \code{maxit=50} or \code{family=binomial} }
}


\details{
A thorough description of this method and package can be found in Calcagno & de Mazancourt.
As the function step in package MASS, glmulti aims at helping and making automatic the selection of best GLM models, using Information Criteria such as AIC.
Contrary to step, it is not iterative (all models are compared), and is fully automatic (no intervention of the user is required).
It is primarily intended to be used for exploratory analyses with many candidate explanatory variables (say 10 to 20). Calculations will typically be rather long (from one night to three days), but once the calculation is finished you will have a complete picture of the fit of models to your data.
This function uses background java classes contained in the archive glmulti.jar.
Running the function therefore requires a Java Running Environment, and the R package rJava.
The fitting of models is performed by R though, using glm(), or any similar function you may have provided using argument \code{fitfunc}.
}


\value{
 Typically nothing is returned and the results are written in a file (named after \code{filename}). This text file is a data frame (tab separated), that can be read with R, OpenOffice Calc, MS Excel or any similar software. Lines are the models in the candidate set, and columns are properties of these models, most iportantly their GLM formula, their Information Criterion, and their complexity (number of parameters estimated from the data).
 From this file, lots of analyses can be performed without fitting any additional model (e.g. model averaging).
 When using \code{method="d"}, returns the number of models in the candidate set (see Examples). 
 When running a GA, two tiny java files (serialized objects) are also written to the disk at regular intervals. They can be used to resume the calculation (\code{method="r"}) if it was interrupted for any reason. This can also be used to continue a GA with modified parameters (e.g. mutation rate).
 }
 
 
\references{
Burnham & Anderson (2002) Model Selection and Multimodel Inference: an information theoretic approach. Springer.
Calcagno & de Mazancourt (in prep) glmulti: a R package for easy and automatic model selection with Generalized Linear models
Venables & Ripley (1997) Modern applied statistics with S-Plus. Springer.
}


\author{Vincent Calcagno, McGill University, Montreal}

\note{}

\seealso{  \code{\link{step}} and \code{\link{stepAIC}} in package MASS }

\examples{

# For a range of examples see the article by Calcagno & de Mazancourt.
# This simply returns the number of candidate models if one has 2 covariates and 2 factors of three levels each, considers pairwise interactions and does not apply marginality
# First we make a dummy data frame
dummyfactor=as.factor(c(rep("y",5),rep("yn",5),rep("n",10)))
dummycovariate = 1:20
dod = data.frame(2:21,dummyfactor,dummyfactor,dummycovariate,dummycovariate)
# Then we call glmulti
glmulti(names(dod)[1],names(dod)[-1], data=dod, method="d")
# There are 250 models to compare

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{models}

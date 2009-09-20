# Main function glmulti (this is mostly a R front end for background Java classes in package glmulti).
# version 0.5-2

glmulti<-function(y, xr, data,  exclude=c(), intercept=TRUE, marginality=FALSE, level=2, filename="glmulti.output", method="h", crit="aicc", chunk=1, chunks=1, minsize=0, maxsize=-1, minK=0, maxK=-1,  plotty=TRUE, confsetsize=100, popsize=100, mutrate=10^-3,sexrate=0.1,imm=0.3, deltaM=0.05, deltaB=0.05, conseq=5, fitfunc=glm, resumefile = "id",  ...) {

# y : name of dependent variable
# xr: vector of names of explanatory variables (factors and covariates)
# exclude: terms (main effects or interactions) to exclude from the candidate models
# data: the data frame containing your data
# intercept: whether to include an intercept in formulas
# marginality: whether to apply he full marginality rule (a specific kind of constraint on candidate models)
# level: the level of interaction terms to include (1 or 2 only)
# minsize: minimal number of TERMS to be included in candidate models (0 or negative = no constraint)
# maxsize: maximal number of TERMS to be included in candidate models (negative = no constraint)
# minK: minimal complexity of candidate models (0 or negative = no constraint)
# maxK: maximal complexity of candidate models (negative = no constraint)
# method:  h for holistic screening, g for genetic algorithm, r to resume a GA simulation
# confsetsize: how many models to keep in the confidence set of models
# crit: the support criterion to use (aic, aicc or bic)
# popsize: pop size for genetic algorithm
# sexrate: rate of sexual reproduction for genetic algorithm
# mutrate: per locus mutation rate for genetic algorithm
# imm: immigration rate for genetic algorithm
# plotty: whether to draw IC profiles on the fly
# deltaB: target change in best ic (stop rules for GA)
# deltaM: target change in mean ic (stop rules for GA)
# conseq: required successive no improvement (i.e. target change attained) (stop rules for GA)
# fitfunc: the fitting function. Default is glm. You can provide any function that matches the parameters/returned values of glm.
# resumefile: the names of files (without extensions) from which to restore the java objects, when using method "r". Default: taken to be filename.

# ... : further options to be passed to the fitting function


# Functions for calculating the Information Criterion of a model
# Three are provided, accepting logLik, K and n as parameters. You can redefine them if you wish.
faic <- function(lk,k,n) {return(-2*lk+2*k)}
faicc <- function(lk,k,n) {return(-2*lk + 2*k*n/max(n-k-1,0))}
fbic <- function(lk,k,n) {return(-2*lk + k*log(n))}


# some general constants
if (method=="h") DELTAD=100
else DELTAD = 100


# parameters
# remove duplicates
x<-unique(xr)
excluzt<-unique(exclude)
# "qazxcww" is added to avoid zero length arrays that confuse R. It is always ignored.
ssize<-length(data[,1])
write(paste("Sample size:",ssize),file="")
xc<- c("qazxcww")
xq<- c("qazxcww")
nc<-0
nq<-0
nana=names(data)
for (i in x) if (is.factor(data[[which(nana==i)]])) {
nc<-nc+1
xc<-c(xc,i)
} else {
nq<-nq+1
xq<-c(xq,i)
}
excluz = c("qazxcww",excluzt)

if (crit=="aic") support<-faic
if (crit=="aicc") support<-faicc
else support <- fbic

# instanciate java object
if (method=="r") {
# resume a genetic algorithm
# ----------------------------
write("Restoring from files...",file="")
if (resumefile=="id") molly <- .jcall("glmulti/Resumator","Lglmulti/ModelGenerator;","resto",filename)
else molly <- .jcall("glmulti/Resumator","Lglmulti/ModelGenerator;","resto",resumefile)
if (!.jfield(molly,"Z","ok")) {
write("WARNING: Could not restore analysis from files.",file="")
return("Aborted.")
}


} else {
# brand new analysis
molly <- .jnew("glmulti/ModelGenerator",y,.jarray(xc),.jarray(xq),.jarray(excluz),as.integer(level),as.integer(minsize),as.integer(maxsize),as.integer(minK),as.integer(maxK),intercept,marginality)
if (!.jfield(molly,"Z","ok")) {
write("WARNING: oversized candidate set.",file="")
return("Aborted.")
}

# informs it of the number of levels of factors, if any
if (nc>0) {
flevs = integer(nc)
for (i in 1:nc) flevs[i]<-as.integer(nlevels(data[,which(nana==xc[i+1])])) 
.jcall(molly,"V","supplyNbLev",.jarray(flevs))
} else .jcall(molly,"V","supplyNbLev",.jarray(integer(2)))

# informs it about the number of df consumed by the error distrib
options(warn=-1)
.jcall(molly,"V","supplyErrorDF",as.integer(attr(logLik(fitfunc(as.formula(paste(y,"~1")),data=data,...)),"df")-1))
options(warn=0)
}


#
   #
      #


if (method == "d") {
# go for simple diagnostic !
write("\nTASK: Diagnostic of candidate set.",file="")
write(paste(nc,"factor(s)."),file="")
write(paste(nq,"covariate(s)."),file="")
write(paste(.jfield(molly,"I","nbforbxc"),"f exclusion(s)."),file="")
write(paste(.jfield(molly,"I","nbforbxq"),"c exclusion(s)."),file="")
write(paste(.jfield(molly,"I","nbforbxcxc"),"f:f exclusion(s)."),file="")
write(paste(.jfield(molly,"I","nbforbxqxq"),"c:c exclusion(s)."),file="")
write(paste(.jfield(molly,"I","nbforbxcxq"),"f:c exclusion(s)."),file="")
write(paste("Size constraints: min = ",minsize,"max =",maxsize),file="")
write(paste("Complexity constraints: min = ",minK,"max =",maxK),file="")
if (marginality) write("Marginality rule.",file="")
nbcand = .jcall(molly,"I","diagnose")
if (nbcand==-1) write("Your candidate set contains more than 1 billion models.",file="")
else write(paste("Your candidate set contains",nbcand,"models."),file="")

return(nbcand)
}


#
   #
      #



if (method=="h") {
# go for holistic screening !
# ----------------------------
# On initialise le fichier d'output...
write("\nINFO: Preparing output file...",file="")
if (chunks>1) {
filen = paste(filename,"_",chunk,".",chunks,sep="")
} else filen = filename
if (file.exists(paste(filen,".txt",sep=""))) {
write("Caution, such a file already exists.",file="")
dolo=readline("Overwrite? (y/n) ")
if (dolo=="n") return(0)
}
header<-.jcall(molly,"S","makeHeaderH")
if (crit=="aic") header<-paste(header,"K\tAIC\tFormula\n",sep="\t")
else if (crit=="aicc") header<-paste(header,"K\tAICc\tFormula\n",sep="\t")
else header<-paste(header,"K\tBIC\tFormula\n",sep="\t")


write(header,file=paste(filen,".txt",sep=""))


write("INFO: Fitting...",file="")
flush.console()
if(!.jcall(molly,"Z","produceModels",as.integer(chunk),as.integer(chunks))) {
write("!Failed to start Java thread.",file="")
return("Aborted.")
}
lesCrit<-numeric(confsetsize)
lesK<-vector('integer',confsetsize)
lesForms<-vector('character',confsetsize)
lesMods<-vector('character',confsetsize)

# proceed with fitting
curr<-0
sel<-0

while(.jcall(molly,"Z","nextModel")) {
formula<-.jcall(molly,"S","getCurrModel")
#if (!is.null(formula)) {
curr<-curr+1
beber<-fitfunc(as.formula(formula),data=data, ...)
# convergence ?
proceed=TRUE
if (!beber$converged) {
proceed<-FALSE
write("!Fit function failed to converge: model skipped.",file="")
write(paste("skipped:", formula),file="")
}
if(proceed) {
liliac<- logLik(beber)
K<-attr(liliac,"df")
# worthy model?
cricri<-support(liliac[[1]],K,beber$df.null+1)
if (sel<confsetsize) {
sel=sel+1
lesForms[sel]=formula
lesMods[sel]=.jcall(molly,"S","getCurrTerms")
lesCrit[sel]=cricri
lesK[sel]=K
} else {
mini=max(lesCrit)
if (cricri<mini) {
ou=which(lesCrit==mini)[1]
lesForms[ou]=formula
lesMods[ou]=.jcall(molly,"S","getCurrTerms")
lesCrit[ou]=cricri
lesK[ou]=K
}
}
if (curr%%DELTAD == 0) {
write(paste("\nAfter ",curr," models:",sep=""),file="")
bestofsex=min(lesCrit[1:sel])
write(paste("Best model: ",gsub(" ","",lesForms[which(lesCrit == bestofsex)]),sep=""),file="")
write(paste("Crit= ",bestofsex,sep=""),file="")
write(paste("Mean crit= ",mean(lesCrit[1:sel]),sep=""),file="")
flush.console()
if (plotty) {
plot(sort(lesCrit[1:sel]),xlab="Best models",ylab=paste("Support (",crit,")"),main="IC profile")
abline(h=bestofsex+2,col="red")
}
}
}

#} else break
}
# over !
write("\nINFO: Completed! Exporting results...",file="")
for (i in 1:sel) {
write(paste(lesMods[i],lesK[i],lesCrit[i],lesForms[i],sep="\t"),file=paste(filen,".txt",sep=""),append=TRUE)
} 


#
   #
      #



} else {
#Go for genetic algorithm approach !
write("\nTASK: Genetic algorithm in the candidate set.",file="")
# On initialise le fichier d'output...
write("\nINFO: Preparing output file...",file="")
if (chunks>1) {
filen = paste(filename,"_",chunk,".",chunks, sep="")
} else filen = filename
if (file.exists(paste(filen,".txt",sep=""))) {
write("Caution, such a file already exists.",file="")
dolo=readline("Overwrite? (y/n) ")
if (dolo=="n") return(0)
}

write("\nINFO: Initialization...",file="")
currgen = 0
consoude = 0
gogo=TRUE
bestofsex = 10000
minou = 10000
bestofsexN = 1000
minouN = 1000

if (method=="r") popul = .jcall(molly,"[S","initPopAgain", mutrate, imm, sexrate, filen)
else popul = .jcall(molly,"[S","initPop",as.integer(popsize), mutrate, imm, sexrate, as.integer(confsetsize),filen)

header<-.jcall(molly,"S","makeHeaderG")
if (crit=="aic") header<-paste(header,"K\tAIC\tFormula\n",sep="\t")
else if (crit=="aicc") header<-paste(header,"K\tAICc\tFormula\n",sep="\t")
else header<-paste(header,"K\tBIC\tFormula\n",sep="\t")
write(header,file=paste(filen,".txt",sep=""))


write("\nINFO: Algorithm started...",file="")
while (gogo) {
for (i in 1:DELTAD) { 
nbtofit = length(popul)
lesic = numeric(nbtofit)
	if (nbtofit>0) for (m in 1:nbtofit) {
	# fit models
	formula=popul[m]
	beber<-fitfunc(as.formula(formula),data=data, ...)
	liliac<- logLik(beber)
	K<-attr(liliac,"df")
	# convergence ?
	if (beber$converged) lesic[m]=support(liliac[[1]],K,beber$df.null+1)
	else lesic[m]=100000
	}
currgen = currgen+1
popul = .jcall(molly,"[S","nextGeneration",.jarray(lesic))
}
# report about the current state and store things
write(paste("\nAfter ",currgen," generations:",sep=""),file="")
lesCrit = .jcall(molly,"[D","reportConfIC")
bestofsexN=.jcall(molly,"D","reportBestIC")
minouN=.jcall(molly,"D","reportMeanIC")
bestform = .jcall(molly,"S","reportbestModel")
write(paste("Best model: ",gsub(" ","",bestform),sep=""),file="")
write(paste("Crit= ",bestofsexN,sep=""),file="")
write(paste("Mean crit= ",minouN,sep=""),file="")
.jcall(molly,"V","printImage")
flush.console()
if (plotty) {
plot(sort(lesCrit),xlab="Best models",ylab=paste("Support (",crit,")"),main="IC profile")
abline(h=bestofsexN+2,col="red")
}
# Shall we continue ?
if (minouN - minou >= -deltaM && bestofsexN - bestofsex >= -deltaB) consoude = consoude+1 else consoude = 0
if (consoude == conseq) {
write("\nINFO: Improvements in best and average IC in the confidence set have been below the specified goals.",file="")
write("Algorithm is thus declared to have converged.",file="")
gogo=FALSE
} else {
write(paste("Change in best IC:",bestofsexN - bestofsex),file="")
write(paste("Change in mean IC:",minouN - minou),file="")
minou = minouN
bestofsex = bestofsexN
} 
}
# END OF GA !!!!
write("\nINFO: Completed! Exporting results...",file="")
lesForms= .jcall(molly,"[S","reportConfMods")
lesMods = .jcall(molly,"[S","reportConfCodes")
lesKK = .jcall(molly,"[I","reportConfKs")
for (i in 1:length(lesCrit)) {
write(paste(lesMods[i],lesKK[i],lesCrit[i],lesForms[i],sep="\t"),file=paste(filen,".txt",sep=""),append=TRUE)
} 

}


#
   #
      #
      


write("\nBYE: Analysis completed successfully.",file="")

# end of function.
}



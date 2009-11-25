# S4 class definition
setClass("glmulti", representation(name="character", params="list", nbmods="integer", crits="numeric", 
		codes="data.frame", K="integer", formulas="list", call="call", adi="list"))


# global variable
assign("glmultiqaiccvalue",NA,.GlobalEnv) 

	
# S3 function definitions
# utilities
summary.glmulti<-function(object, ...)
{
	who=object@formulas[object@crits==object@crits[1]]
	ww = exp(-(object@crits - object@crits[1])/2)
	ww=ww/sum(ww)
	cucu=function(i) sum(ww[1:i])
	wwc=lapply(1:length(ww),cucu)
	coco = ww%*%as.matrix(object@codes[,!is.na(object@codes[2,])])
	ret=list(name=object@name, method=object@params$method, fitting=object@params$fitfunction, crit=object@params$crit, 
			level=object@params$level, marginality=object@params$marginality,confsetsize=object@params$confsetsize, 
			bestic=object@crits[1], icvalues=object@crits ,	bestmodel=who,modelweights=ww, termweights=coco)
	if (object@params$method == "g") {
		ret=c(ret, generations=object@params$generations, elapsed= (object@params$elapsed)/60)
	}
	ret
}

print.glmulti<-function(x, ...)
{
	write(x@name, file="")
	write(paste("Method: ",x@params$method," / Fitting: ",x@params$fitfunction," / IC used: ", 
			x@params$crit,  sep=""), file="")
	write(paste("Level: ", x@params$level, " / Marginality: ", x@params$marginality,  sep=""), file="")
	write(paste("From ",length(x@crits)," models:",sep=""), file="")
	write(paste("Best IC: ",x@crits[1],sep=""), file="")
	write(paste("Best model:",sep=""), file="")
	who=x@formulas[x@crits==x@crits[1]]
	for (w in who)
		print(w)
	ww = exp(-(x@crits - x@crits[1])/2)
	ww=ww/sum(ww)
	write(paste("Evidence weight: ",ww[1], sep=""), file="")
	write(paste("Worst IC: ",x@crits[length(x@crits)],sep=""), file="")
	write(paste(length(which(x@crits<=x@crits[1]+2)), " models within 2 IC units.", sep=""), file="")
	cucu=function(i) sum(ww[1:i])
	wwc=lapply(1:length(ww),cucu)
	write(paste(length(which(wwc<=0.95)), " models to reach 95% of evidence weight.", sep=""), file="")
	if (x@params$method == "g") {
		write(paste("Convergence after ", x@params$generations, " generations.",  sep=""), file="")
		write(paste("Time elapsed: ",  (x@params$elapsed), " minutes.", sep=""), file="")
	}
}

plot.glmulti<-function(x, type="p", highlight=c(), ...) 
{
	if ( type == "p") {
		plot(x@crits,xlab="Best models",
				ylab=paste("Support (",x@params$crit,")"),pch=19, main="IC profile", ...)
		
		if (length(highlight)!=0) {
			coco = as.matrix(x@codes[,!is.na(x@codes[2,])]) 
			for (tt in highlight) {
				points(y=x@crits[coco[,tt]==1],x=which(coco[,tt]==1) ,pch=19,col="yellow")
				points(y=x@crits[coco[,tt]==1],x=which(coco[,tt]==1) ,pch=21,col="black")
			}
		}
		abline(h=x@crits[1]+2,col="red")
	} else if (type == "w") {
		ww = exp(-(x@crits - x@crits[1])/2)
		ww=ww/sum(ww)
		plot(ww,xlab="Best models",
				ylab=paste("Evidence weight (",x@params$crit,")"),pch=19, main="Profile of model weights", ...)
		cucu=function(i) sum(ww[1:i])
		wwc=lapply(1:length(ww),cucu)
		if (length(highlight)!=0) {
			coco = as.matrix(x@codes[,!is.na(x@codes[2,])]) 
			for (tt in highlight) {
				points(y=ww[coco[,tt]==1],x=which(coco[,tt]==1) ,pch=19,col="yellow")
				points(y=ww[coco[,tt]==1],x=which(coco[,tt]==1) ,pch=21,col="black")
			}
		}
		abline(v=min(which(wwc>=0.95)),col="red")
	} else if (type=="s") {
		ww = exp(-(x@crits - x@crits[1])/2)
		ww=ww/sum(ww)
		coco = as.matrix(x@codes[,!is.na(x@codes[2,])])
		barplot((ww%*%coco), xlab="", ylab=paste("Evidence weight (",x@params$crit,")",sep=""), 
				 main="Importance of model terms",xlim=c(0,1.2*length(x@codes[,!is.na(x@codes[2,])])), las=2,...)
		abline(v=0.1+1.2*which(is.na(x@codes[2,]))-1.2*(1:length(which(is.na(x@codes[2,])))), col="yellow")
	} else {
		warning("plot: Invalid type argument for plotting glmulti objects.")
	}
	
	
}

# model averaging: done through coef
coef.glmulti <- function(object, select="all", varweighting="Buckland", ...) 
{
	ww = exp(-(object@crits - object@crits[1])/2)
	ww=ww/sum(ww)
	cucu=function(i) sum(ww[1:i])
	wwc=lapply(1:length(ww),cucu)
	# what models are to be used
	whom=c()
	if (length(select)>1)
		whom=select
	else if (select=="all")
			whom=1:length(object@crits)
		else if (is.integer(select))
				whom = 1:select
			else if (is.numeric(select)) {
				if (select <=1) 
					whom = which(wwc<=select)
				else 
					whom = which(object@crits <= (object@crits[1]+select))
		}
	mods = object@crits[whom]
	formo = object@formulas[whom]
	# fit selected models
	coffee=list()
	for (i in formo) {
		cak=as.call(list(substitute(match.fun(object@params$fitfunction)), formula=i, data=object@call$data))
		if (length(object@adi)>1)
			for (j in 1:length(object@adi)) {
				cak[[length(names(cak))+1]] = object@adi[[j]]
				names(cak)[length(names(cak))] = names(object@adi)[j]
			}
		modf=eval(cak)
		coffee=c(coffee,list(modf))
	}
	# construct list of coefficients
	coke=lapply(coffee,getfit)
	namou=unique(unlist(lapply(coke,function(x) dimnames(x)[[1]])))

	coconutM=matrix(0,length(formo),length(namou))
	coconutSE=matrix(0,length(formo),length(namou))
	coconutN = numeric(length(namou))
	# get values, deviations, presence for all models
	gettou=function(i) {
		ele=coke[[i]]
		nana = dimnames(ele)[[1]]
		mimi=numeric(3*length(namou))
		if (length(nana) >1) {
			for (k in 1:(length(nana))) {
				mimi[match(nana[k],namou)]=ele[k,1]
				mimi[match(nana[k],namou)+length(namou)]=ele[k,2]
				mimi[match(nana[k],namou)+2*length(namou)]=1
			}
		} else {
			mimi[match(nana[1],namou)]=ele[1]
			mimi[match(nana[1],namou)+length(namou)]=ele[2]
			mimi[match(nana[1],namou)+2*length(namou)]=1
		}
		return(mimi)
	}
	lol=sapply(lapply(1:length(coke),gettou),rbind)
	coconutM = t(lol[1:length(namou),])
	coconutSE = t(lol[(1:length(namou))+length(namou),])
	# NA are set to zero
	coconutM[is.na(coconutM)]=0
	coconutM[is.na(coconutSE)]=0
	coconutN =  t(lol[(1:length(namou))+2*length(namou),])
	nene = matrix(rep(1,length(whom))%*%coconutN, nc=1, dimnames=list( namou, c("Nb models")))
	# construct weight vectors
	waou=ww[whom]/sum(ww[whom])
	waouv= t(matrix(rep(waou,length(namou)),length(namou),length(whom)))*coconutN
	totwaou = waou%*%coconutN
	# weight estimates
	averest = matrix((rep(1,length(whom))%*%(waouv*coconutM))/totwaou, nc=1, dimnames=list( namou, c("Estimate")))
	weighty =  matrix(totwaou, nc=1, dimnames=list( namou, c("Importance")))
	# weight variances
	if (varweighting=="Burnham") {
		squaredevs = rep(1,length(whom))%*%(waouv*((coconutM-t(matrix(rep(averest,length(whom)), length(namou), length(whom))))^2))
		condivars =  rep(1,length(whom))%*%(waouv*(coconutSE^2))
		avervar = matrix(condivars+squaredevs, nc=1, dimnames=list( namou, c("Uncond. variance")))
	} else if (varweighting=="Buckland") {
		squaredevs = ((coconutM-t(matrix(rep(averest,length(whom)), length(namou), length(whom)))))^2
		condivars = coconutSE^2
		avervar = matrix((rep(1,length(whom))%*%(waouv*sqrt(squaredevs+condivars)))^2, nc=1, dimnames=list( namou, c("Uncond. variance")))
	} else 
		avervar = matrix(rep(NA,length(namou)), nc=1, dimnames=list( namou, c("Uncond. variance")))
	averaging = cbind(averest, avervar, nene, weighty)
	
	averaging
}


# model averaged prediction
predict.glmulti <- function(object, select="all", ...) 
{
	ww = exp(-(object@crits - object@crits[1])/2)
	ww=ww/sum(ww)
	cucu=function(i) sum(ww[1:i])
	wwc=lapply(1:length(ww),cucu)
	# what models are to be used
	whom=c()
	if (length(select)>1)
		whom=select
	else if (select=="all")
			whom=1:length(object@crits)
		else if (is.integer(select))
				whom = 1:select
			else if (is.numeric(select)) {
				if (select <=1) 
					whom = which(wwc<=select)
				else 
					whom = which(object@crits <= (object@crits[1]+select))
		}
	mods = object@crits[whom]
	formo = object@formulas[whom]
	# fit selected models
	coffee=list()
	for (i in formo) {
		cak=as.call(list(substitute(match.fun(object@params$fitfunction)), formula=i, data=object@call$data))
		if (length(object@adi)>1)
			for (j in 1:length(object@adi)) {
				cak[[length(names(cak))+1]] = object@adi[[j]]
				names(cak)[length(names(cak))] = names(object@adi)[j]
			}
		modf=eval(cak)
		coffee=c(coffee,list(modf))
	}
	waou=ww[whom]/sum(ww[whom])
	# make predictions
	preds=list()
	for (i in 1:length(formo)) {
		preds = c(preds,list(predict(coffee[[i]])))
	}
	nbpo = length(preds[[1]])
	all = t(matrix(unlist(preds), nr=nbpo))
	dimnames(all)=list(c(), names(preds[[1]]) )
	# handle NA values
	nana = lapply(preds, is.na)
	nbok = numeric(nbpo)
	for (i in 1:length(preds)) {
		nbok = nbok + nana[[i]]
		preds[[i]][is.na(preds[[i]])] = 0
		}
	minou = matrix(waou%*%t(matrix(unlist(preds),nr=nbpo)), dimnames=list(names(preds[[1]]),c() )) 
	vava = matrix(waou%*%t(matrix((unlist(preds)-unlist(rep(minou[1,],length(preds))))^2, nr=nbpo)),dimnames=list(names(preds[[1]]),c()))
	list(averages = minou, variances = vava, all = all, omittedNA = nbok)
	
}

#  S4 function definitions

# generic function
setGeneric("glmulti", function(y, xr, data, exclude=c(), name="glmulti.analysis", intercept=TRUE, marginality=FALSE ,chunk=1, chunks=1, 
		level=2, minsize=0, maxsize=-1, minK=0, maxK=-1, method="h",crit="aicc",confsetsize=100,popsize=100,mutrate=10^-3,
		sexrate=0.1,imm=0.3, plotty=TRUE, report=TRUE, deltaM=0.05, deltaB=0.05, conseq=5, fitfunction="glm", resumefile = "id",  ...) 
		standardGeneric("glmulti"))
		
# write redefinition
setGeneric("write", function(x, file = "data",   ncolumns=if(is.character(x)) 1 else 5, append = FALSE, sep = " ") standardGeneric("write"))
setMethod("write","glmulti", function(x, file , ncolumns, append, sep)
{
		if (!missing(file)&&length(grep(file, "|object"))==1) {
			ir = gsub("\\|object","",file)
			if (ir=="")
				.saveRDS(x, file=x@name)
			else 
				.saveRDS(x, file=ir)
		} else {
			concato=cbind(as.data.frame(x@codes),data.frame(K=x@K),data.frame(IC=x@crits), data.frame(Models=as.character(x@formulas)))
			if (missing(file))
				write.table(concato, file = paste(x@name, ".txt"),  append = append, sep = sep )
			else
				write.table(concato, file = file,  append = append, sep = sep )
		}
})
		
# accessing fitted models for coefficients
setGeneric("getfit", function(object, ...) standardGeneric("getfit"))

setMethod("getfit","glm", function(object, ...)
{
	return(summary(object)$coefficients[,1:2])
})

setMethod("getfit","lm", function(object, ...)
{
	return(summary(object)$coefficients[,1:2])
})

# consensus method
setGeneric("consensus", function (xs, confsetsize=NA, ...)  standardGeneric("consensus"))

setMethod("consensus", signature(xs="list"), function (xs, confsetsize, ...)
{
	lespaul = list()
	for (i in xs) {
	if (class(i)=="glmulti")
		lespaul=c(lespaul, i)
	else if (class(i)=="character") {
			paul = .readRDS(file=i)
			if (class(paul)=="glmulti")
				lespaul=c(lespaul,paul)
		}
	}
	ouou = length(lespaul[[1]]@codes[1,])
	neo = new ("glmulti")
	conca = function (x) {
		cbind(as.data.frame(x@codes),data.frame(K=x@K),data.frame(IC=x@crits), data.frame(formulas=as.character(x@formulas)))
	}
	tot=lapply(lespaul, conca)
	tota = tot[[1]]
	for (h in 2: length(tot))
		tota = rbind(tota, tot[[h]])
	tot = tota
	tot = tot[!duplicated(tot$formulas),]
	tot = tot[order(tot$IC),]
	if (is.na(confsetsize) || length(tot$K)<confsetsize) {
		if (!is.na(confsetsize)&&length(tot$K)<confsetsize)
			warning("Could not gather enough models.")
		neo@K = tot$K
		neo@formulas =  as.list(lapply(unlist(lapply(tot$formulas, function(j) as.character(j))), function(ff) as.formula(c(ff))))
		neo@crits = tot$IC
		neo@codes = tot[,1:ouou]
	} else {
		tot = tot[1:confsetsize, ]
		neo@K = tot$K
		neo@formulas = as.list(lapply(unlist(lapply(tot$formulas, function(j) as.character(j))), function(ff) as.formula(c(ff))))
		neo@crits = tot$IC
		neo@codes = tot[,1:ouou]	
	}
	neo@call = lespaul[[1]]@call
	neo@adi = lespaul[[1]]@adi
	neo@name = paste("consensus of ", length(lespaul), "-", lespaul[[1]]@name, sep="")
	neo@params = lespaul[[1]]@params
	return (neo)
	})

# information criteria
setGeneric("aicc", function(object, ...) standardGeneric("aicc"))
setGeneric("aic", function(object, ...) standardGeneric("aic"))
setGeneric("bic", function(object, ...) standardGeneric("bic"))
setGeneric("qaic", function(object, ...) standardGeneric("qaic"))
setGeneric("qaicc", function(object, ...) standardGeneric("qaicc"))
setMethod("aicc", "lm", function(object, ...)
{
	liliac<- logLik(object)
	k<-attr(liliac,"df")
	n= length(resid(object))
	return(-2*logLik(object) + 2*k*n/max(n-k-1,0))
})

setMethod("bic", signature(object="lm"), function(object, ...)
{
	liliac<- logLik(object)
	k<-attr(liliac,"df")
	n= length(resid(object))
	return(-2*logLik(object) + k*log(n))
})

setMethod("aic", "lm",  function(object, ...)
{
	liliac<- logLik(object)
	k<-attr(liliac,"df")
	return(-2*logLik(object)+2*k)
})
setMethod("qaic", "lm",  function(object, ...)
{
	liliac<- logLik(object)
	k<-attr(liliac,"df")
 	return(liliac / glmultiqaiccvalue + 2 * k)

})

setMethod("qaicc", "lm",  function(object, ...)
{
	liliac<- logLik(object)
	k<-attr(liliac,"df")
	n= length(resid(object))
 	return(liliac / glmultiqaiccvalue + 2*k*n/max(n-k-1,0))

})

	
# interfaces for formulas/models: calls with missing xr argument

setMethod("glmulti","missing",
function(y, xr, data, exclude, name, intercept, marginality ,chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile,  ...) 
{
	write("This is glmulti 0.6-1, november 2009.",file="")
})

setMethod("glmulti",
def=function(y, xr, data, exclude, name, intercept, marginality ,chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile,  ...) 
{
	stop("Improper call of glmulti.")
})


setMethod("glmulti", signature(y="formula", xr= "missing", exclude="missing"),
function(y, xr, data, exclude, name, intercept, marginality ,chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile,  ...) {
	if (missing(data))
		tete = terms(y)
	else
		tete = terms(y, data=data)
	oo = attr(tete,"order")
	dep = as.character(attr(tete,"variables"))[2]
	int = attr(tete,"intercept")
	preds = as.character(attr(tete,"variables"))[-(1:2)]
	if (level==2 && max(oo)>1) {
		# get all possible interactions
		interac = attr(tete,"term.labels")[oo==2]
		neotete = terms(as.formula(paste("h~",paste(preds, collapse="*"))))
		neointerac= attr(neotete,"term.labels")[attr(neotete,"order")==2]
		# get exclusions
		for (i in interac)
			neointerac=neointerac[neointerac!=i]
		# same for main effects
		mama = attr(tete,"term.labels")[oo==1]
		exma = preds
		for (j in mama)
			exma = exma[exma!=j]
		exma = c(exma,neointerac)
	} else {
		preds = attr(tete,"term.labels")[oo==1]
		exma=c(1)
	}
	call = match.call()
	call[[match("y", names(call))]] = dep
	call[[length(names(call))+1]] = preds
	names(call)[length(names(call))] ="xr"
		call[[length(names(call))+1]] = exma
	names(call)[length(names(call))] ="exclude"
	
	if (missing(data)) {
		call[[length(names(call))+1]] = environment(y)
		names(call)[length(names(call))] ="data"
	}
	eval(call)
})

setMethod("glmulti", signature(y="character", xr= "missing", exclude="missing"),
function(y, xr, data, exclude, name, intercept, marginality ,chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile,  ...) {
	call = match.call()
	call[[match("y", names(call))]] = as.formula(y)
	eval(call)
})

setMethod("glmulti", signature(y="glm", xr= "missing", exclude="missing"),
function(y, xr, data, exclude, name, intercept, marginality ,chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile,  ...) {
	call = match.call()
	call[[match("y", names(call))]] = formula(y)
	call[[length(names(call))+1]] = "glm"
	names(call)[length(names(call))] = "fitfunction"
	if (length(y$data)) {
		call[[length(names(call))+1]] = y$data
		names(call)[length(names(call))] = "data"
	}
	eval(call)
})

setMethod("glmulti", signature(y="lm", xr= "missing", exclude="missing"),
function(y, xr, data, exclude, name, intercept, marginality , chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile,  ...) {		
	call = match.call()
	call[[match("y", names(call))]] = formula(y)
	call[[length(names(call))+1]] = "lm"
	names(call)[length(names(call))] = "fitfunction"
	if (length(summary(y)$call$data)) {
		call[[length(names(call))+1]] = summary(y)$call$data
		names(call)[length(names(call))] = "data"
	}
	eval(call)
})



	

# workhorse function (call with building blocks and exclusions, as in earlier versions 0.5-x)
setMethod("glmulti", signature(y="character", xr= "character"),
function(y, xr, data, exclude, name, intercept, marginality ,chunk, chunks, 
		level, minsize, maxsize, minK, maxK, method,crit,confsetsize,popsize,mutrate,
		sexrate,imm, plotty,  report, deltaM, deltaB, conseq, fitfunction, resumefile,  ...) 
{
	if (report) write("Initialization...",file="")

	# some general constants
	if (method=="h") DELTAD=50
	else DELTAD=10


	# handle terms
	# remove duplicates if any
	x<-unique(xr)
	excluzt<-unique(exclude)
	# "qazxcww" is added to avoid zero length arrays that confuse R. It is always ignored.
	databis = model.frame(as.formula(paste(y,"~", paste(x, sep="", collapse="+"),sep="")), data = data)
	ssize<-length(databis[,1])
	xc<- c("qazxcww")
	xq<- c("qazxcww")
	nc<-0
	nq<-0
	nana=names(databis)
	for (i in x) 
		if (is.factor(databis[[which(nana==i)]])) {
			nc<-nc+1
			xc<-c(xc,i)
		} else {
			nq<-nq+1
			xq<-c(xq,i)
		}
	excluz = c("qazxcww",excluzt)
	
	# functions for support and fit
	support = match.fun(crit)
	if (is.function(crit))
		crit = deparse(substitute(crit))
	fitfunc = match.fun(fitfunction)
	if (is.function(fitfunction)) 
  		fitfunction <- deparse(substitute(fitfunction))
	if (crit == "qaic" || crit == "qaicc") {
		# check presence of c estimate
		if (get("glmultiqaiccvalue", pos=globalenv()) == NA) {
			stop("To use QAIC or QAICc you must first provide an estimate of c as the global variable glmultiqaiccvalue.")
		}
	}
	
	
	# instanciate java object
	if (method=="r") {
		# resume a genetic algorithm
		write("Restoring from files...",file="")
		if (resumefile=="id") 
			molly <- .jcall("glmulti/Resumator","Lglmulti/ModelGenerator;","resto",name)
		else molly <- .jcall("glmulti/Resumator","Lglmulti/ModelGenerator;","resto",resumefile)
		if (!.jfield(molly,"Z","ok")) {
			warning("!Could not restore analysis from files.")
			return(-1)
		}

	} else {
		# brand new analysis
		molly <- .jnew("glmulti/ModelGenerator",y,.jarray(xc),.jarray(xq),.jarray(excluz),as.integer(level),as.integer(minsize),
				as.integer(maxsize),as.integer(minK),as.integer(maxK),intercept,marginality)
		if (!.jfield(molly,"Z","ok")) {
			warning("!Oversized candidate set.")
			return(-1)
		}

		# informs java object of the number of levels of factors, if any
		if (nc>0) {
			flevs = integer(nc)
			for (i in 1:nc) 
				flevs[i]<-as.integer(nlevels(databis[,which(nana==xc[i+1])]))
			.jcall(molly,"V","supplyNbLev",.jarray(flevs))
		} else .jcall(molly,"V","supplyNbLev",.jarray(integer(2)))

		# informs it about the number of df consumed by the error distrib
		options(warn=-1)
		.jcall(molly,"V","supplyErrorDF",as.integer(attr(logLik(fitfunc(as.formula(paste(y,"~1")),data=data, ...)),"df")-1))
		options(warn=0)
	}


	#
	   #
		  #


	if (method == "d") {
		# go for simple diagnostic !
		write("TASK: Diagnostic of candidate set.",file="")
		write(paste("Sample size:",ssize),file="")
		write(paste(nc,"factor(s)."),file="")
		write(paste(nq,"covariate(s)."),file="")
		write(paste(.jfield(molly,"I","nbforbxc"),"f exclusion(s)."),file="")
		write(paste(.jfield(molly,"I","nbforbxq"),"c exclusion(s)."),file="")
		write(paste(.jfield(molly,"I","nbforbxcxc"),"f:f exclusion(s)."),file="")
		write(paste(.jfield(molly,"I","nbforbxqxq"),"c:c exclusion(s)."),file="")
		write(paste(.jfield(molly,"I","nbforbxcxq"),"f:c exclusion(s)."),file="")
		write(paste("Size constraints: min = ",minsize,"max =",maxsize),file="")
		write(paste("Complexity constraints: min = ",minK,"max =",maxK),file="")
		if (marginality) 
			write("Marginality rule.",file="")
		nbcand = .jcall(molly,"I","diagnose")
		if (nbcand==-1) 
			write("Your candidate set contains more than 1 billion (1e9) models.",file="")
		else 
			write(paste("Your candidate set contains",nbcand,"models."),file="")

	return(nbcand)
	}


	#
	   #
		  #



	if (method=="h") {
	# go for exhaustive screening 
	if (report) write("TASK: Exhaustive screening of candidate set.",file="")
	# init output object
	resulto <- new("glmulti")
	if (chunks>1) {
		resulto@name = paste(name,"_",chunk,".",chunks,sep="")
	} else 
		resulto@name = name
	resulto@params = list(name=name, intercept=intercept, marginality=marginality ,chunk=chunk, chunks=chunks, 
			level=level, minsize=minsize, maxsize=maxsize, minK=minK, maxK=maxK, method=method,crit=crit, confsetsize=confsetsize, fitfunction=fitfunction)
	resulto@call=match.call()
	resulto@adi=list(...)
	header<-.jcall(molly,"S","makeHeaderH")

	if (report) write("Fitting...",file="")
	flush.console()
	if(!.jcall(molly,"Z","produceModels",as.integer(chunk),as.integer(chunks))) {
		warning("!Failed to start Java thread.")
		return(-1)
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
		curr<-curr+1
		beber<-fitfunc(as.formula(formula), data=data, ...)
		# convergence ?
		proceed=TRUE
		if (fitfunction=="glm" && !beber$converged) {
			proceed<-FALSE
			warning(paste("!fitting function failed to converge. Model skipped: ", formula, sep=""))
		}
		if(proceed) {
			cricri<-support(beber)
			if (sel<confsetsize) {
				sel=sel+1
				lesForms[sel]=formula
				lesMods[sel]=.jcall(molly,"S","getCurrTerms")
				lesCrit[sel]=cricri
				lesK[sel]=attr(logLik(beber),"df")
			} else {
				mini=max(lesCrit)
				if (cricri<mini) {
					ou=which(lesCrit==mini)[1]
					lesForms[ou]=formula
					lesMods[ou]=.jcall(molly,"S","getCurrTerms")
					lesCrit[ou]=cricri
					lesK[ou]=attr(logLik(beber),"df")
				}
			}
			if (curr%%DELTAD == 0) {
				if (report) {
					write(paste("\nAfter ",curr," models:",sep=""),file="")
					bestofsex=min(lesCrit[1:sel])
					write(paste("Best model: ",gsub(" ","",lesForms[which(lesCrit == bestofsex)]),sep=""),file="")
					write(paste("Crit= ",bestofsex,sep=""),file="")
					write(paste("Mean crit= ",mean(lesCrit[1:sel]),sep=""),file="")
					flush.console()
				} 
				if (plotty) {
					plot(sort(lesCrit[1:sel]),xlab="Best models",ylab=paste("Support (",crit,")"),pch=19,main="IC profile")
					abline(h=bestofsex+2,col="red")
				}
			}
		}
	}
	# over !
	if (report) write("Completed.",file="")
	# prepare and return glmulti object
	reglo <- order(lesCrit[1:sel])
	resulto@crits = lesCrit[1:sel][reglo]
	resulto@formulas = lapply(lesForms[1:sel][reglo],as.formula)
	header = strsplit(header, "\t")[[1]]
	kaka = matrix(0, sel, length(header))
	lesMods <- lesMods[1:sel][reglo]
	options(warn=-1)
	for (i in 1:sel)
			kaka[i,] <- as.numeric(strsplit(lesMods[1:sel][i],"\t")[[1]])
	options(warn=0)
	kaka = data.frame(kaka)
	names(kaka)<- header
	resulto@codes = kaka
	resulto@K = as.integer(lesK[1:sel][reglo])
	resulto@nbmods = as.integer(sel)
	return(resulto)
	
	#
	   #
		  #
		  #
		#
	#



	} else {
		#Go for genetic algorithm approach
		write("TASK: Genetic algorithm in the candidate set.",file="")
		# init output object
		resulto <- new("glmulti")
		resulto@name = name
		resulto@params = list(name=name, intercept=intercept, marginality=marginality ,chunk=chunk, chunks=chunks, 
				level=level, minsize=minsize, maxsize=maxsize, minK=minK, maxK=maxK, method=method,crit=crit, confsetsize=confsetsize, fitfunction=fitfunction, popsize=popsize,mutrate=mutrate,
				sexrate=sexrate,imm=imm, plotty=plotty, deltaM=deltaM, deltaB=deltaB, conseq=conseq, resumefile = resumefile)
		resulto@call=match.call()
		resulto@adi=list(...)
		header<-.jcall(molly,"S","makeHeaderG")
		
		write("Initialization...",file="")
		currgen = 0
		consoude = 0
		gogo=TRUE
		bestofsex = 10000
		minou = 10000
		bestofsexN = 1000
		minouN = 1000
		if (method=="r") 
			popul = .jcall(molly,"[S","initPopAgain", mutrate, imm, sexrate, name)
		else popul = .jcall(molly,"[S","initPop",as.integer(popsize), mutrate, imm, sexrate, as.integer(confsetsize), name)
		header<-.jcall(molly,"S","makeHeaderG")

		tini = Sys.time()
		dyniT = numeric(0)
		dyniB = numeric(0)
		dyniM = numeric(0)
		write("Algorithm started...",file="")
		while (gogo) {
			for (i in 1:DELTAD) { 
				nbtofit = length(popul)
				lesic = numeric(nbtofit)
				if (nbtofit>0) 
					for (m in 1:nbtofit) {
						# fit models
						formula=popul[m]
						beber<-glm(as.formula(formula),data=data, ...)
						liliac<- logLik(beber)
						K<-attr(liliac,"df")
						# convergence ?
						if (fitfunction=="glm" && !beber$converged) 
							lesic[m]=10000
						else lesic[m]=support(beber)
						}
				currgen = currgen+1
				popul = .jcall(molly,"[S","nextGeneration",.jarray(lesic))
			}
			# report about the current state and store things
			lesCrit = .jcall(molly,"[D","reportConfIC")
			bestofsexN=.jcall(molly,"D","reportBestIC")
			minouN=.jcall(molly,"D","reportMeanIC")
			bestform = .jcall(molly,"S","reportbestModel")
			dyniT = c(dyniT,as.numeric(Sys.time()-tini))
			dyniB = c(dyniB,bestofsexN)
			dyniM = c(dyniM,minouN)
			if (report) {
				write(paste("\nAfter ",currgen," generations:",sep=""),file="")
				write(paste("Best model: ",gsub(" ","",bestform),sep=""),file="")
				write(paste("Crit= ",bestofsexN,sep=""),file="")
				write(paste("Mean crit= ",minouN,sep=""),file="")
				flush.console()
				}
			.jcall(molly,"V","printImage")
			if (plotty) {
				plot(sort(lesCrit),xlab="Best models",ylab=paste("Support (",crit,")"),pch=19,main="IC profile")
				abline(h=bestofsexN+2,col="red")
			}
			# Shall we continue ?
			if (length(lesCrit)==confsetsize && minouN - minou >= -deltaM && bestofsexN - bestofsex >= -deltaB) 
				consoude = consoude+1
			else consoude = 0
			if (consoude == conseq) {
				write("Improvements in best and average IC have been below the specified goals.",file="")
				write("Algorithm is declared to have converged.",file="")
				gogo=FALSE
			} else {
				if (report) write(paste("Change in best IC:",bestofsexN - bestofsex, "/ Change in mean IC:",minouN - minou),file="")
				minou = minouN
				bestofsex = bestofsexN
			} 
			}
		# END OF GA
		write("Completed.",file="")
		lesForms= .jcall(molly,"[S","reportConfMods")
		lesMods = .jcall(molly,"[S","reportConfCodes")
		lesKK = .jcall(molly,"[I","reportConfKs")
		sel = length(lesCrit)
		# prepare and return glmulti object
		reglo <- order(lesCrit)
		resulto@crits = lesCrit[reglo]
		resulto@formulas = lapply(lesForms[reglo],as.formula)
		header = strsplit(header, "\t")[[1]]
		kaka = matrix(0, sel, length(header))
		lesMods <- lesMods[reglo]
		options(warn=-1)
		for (i in 1:sel)
				kaka[i,] <- as.numeric(strsplit(lesMods[i],"\t")[[1]])
		options(warn=0)
		kaka = data.frame(kaka)
		names(kaka)<- header
		resulto@codes = kaka
		resulto@K = as.integer(lesKK)
		resulto@nbmods = as.integer(sel)
		resulto@params = c(resulto@params, list(generations=currgen, elapsed=as.numeric(Sys.time()-tini), dynat=dyniT, dynab=dyniB, dynam=dyniM))
		return(resulto)

	}


	#
	   #
		  #
		  

	# end of function.
})


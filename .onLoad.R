.onLoad <- function(libname, pkgname) {
.jinit(parameters="-XmX512m",force.init=TRUE)
.jpackage(pkgname)
  
write("-----------------------------------------------------",file="")
write("| This is glmulti v. 1.1                            |",file="")
write("| Vincent Calcagno, march 2009                      |",file="")
write("| Contact: vincent.calcagno[-a-t-]mcgill[-d-o-t-]ca |",file="")
write("-----------------------------------------------------",file="")

}  

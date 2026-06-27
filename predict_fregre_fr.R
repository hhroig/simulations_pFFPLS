predict_fregre_fr<-function(object, new.fdataobj=NULL,...){
  
  if (is.null(object)) stop("No fregre.fd object entered")
  if (is.null(new.fdataobj)) return(object$fitted.values)
  if (object$call[[1]]=="fregre.basis.fr" || object$call[[1]]=="fda.usc::fregre.basis.fr" || object$call[[1]]=="fregre.basis.fr.cv"){
    beta.est<-object$coefficients
    isfdx <- fda::is.fd(new.fdataobj)
    if (isfdx) {
      xcoef<-new.fdataobj$coef
      ncurves<-ncol(xcoef)
    }
    else {
      xfdobj<-fda::Data2fd(argvals =new.fdataobj$argvals, y = t(new.fdataobj$data), basisobj = object$basis.s)
      xcoef<-xfdobj$coef
      ncurves<-ncol(xcoef)
      if (any(new.fdataobj$argvals!=object$x$argvals)) stop("Incorrect argvals")
    }
    H = t(xcoef) %*% object$H 
    beta.xest = beta.est %*% t(H)
    beta.xfd   = fda::fd(beta.xest, object$basis.t)
    if (isfdx) {
      yhat  <-  fda::eval.fd(object$argvals.y,object$alpha.est) %*% matrix(1,1,ncurves) + fda::eval.fd(object$argvals.y, beta.xfd)  
      fitted.values  <- fda::smooth.basis(object$argvals.y, yhat, object$y$basis)$fd 
    }
    else {
      yhat  <-  fda::eval.fd(object$y$argvals,object$alpha.est) %*% matrix(1,1,ncurves) + fda::eval.fd(object$y$argvals, beta.xfd)
      fitted.values <- fda.usc::fdata(t(yhat),new.fdataobj$argvals,new.fdataobj$rangeval,new.fdataobj$names)
    }
    return(fitted.values)
  }
  
  else stop("predict is only implemented for fregre.basis.fr output object")
}

###########################SIMULATION CODE###########################
### Xu Shi (shixu@umich.edu)
### Below is the functions for Section 4 of the manuscript
### Delta_1 is called GEST: biased if MR.wrong | MZ.wrong
### Delta_2 is called IPW: biased if MW.wrong | MZ.wrong
### Delta_3 is called OR: biased if MW.wrong | MR.wrong | MY.wrong
### MLE is called MIAO: biased if MW.wrong | MY.wrong
###########################SIMULATION CODE###########################

# This function is waiting to be rewritten to return a class object.

#' @title Function for estimation
#' @description  The document of this function is yet to be written.
#'
#' @param orig.data original data
#' @param indices indices
#' @param arg a list of arguments
#'
#' @return estimation result
#' @export
#'
#' @examples
#' est_func() // TODO
#'
est_func = function(orig.data,indices,arg=list(MW.wrong=MW.wrong,MR.wrong=MR.wrong,MY.wrong=MY.wrong,MZ.wrong=MZ.wrong,ATE=ATE)){
  ## get data
  data=orig.data[indices,]
  A=data$A
  Z=data$Z
  W=data$W
  Y=data$Y
  X.org = as.matrix(data[,5:ncol(data)])
  MW.wrong=arg$MW.wrong;MR.wrong=arg$MR.wrong;
  MY.wrong=arg$MY.wrong;MZ.wrong=arg$MZ.wrong;
  ATE=arg$ATE
  ## run all methods
  ### IPW
  (IPW.est = IPW(aa=A,zz=Z,ww=W,yy=Y,X.org,MZ.wrong,MW.wrong))
  ### OR
  (OR.est = OR(aa=A,zz=Z,ww=W,yy=Y,X.org,MW.wrong,MR.wrong,MY.wrong))
  ### GEST
  (GEST.est = GEST(aa=A,zz=Z,ww=W,yy=Y,X.org,MZ.wrong,MR.wrong))
  ### MIAO
  (MIAO.est = MIAO(aa=A,zz=Z,ww=W,yy=Y,X.org,MW.wrong,MY.wrong))
  ### MR
  (MR.est = MR(aa=A,zz=Z,ww=W,yy=Y,X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong))
  est = c(
    GEST=GEST.est,
    IPW=IPW.est,
    OR=OR.est,
    MIAO=MIAO.est,
    MR=MR.est,
    ATE.true = ATE
  )
  return(est)
}

###########################SIMULATION CODE###########################
### Xu Shi (shixu@umich.edu)
### Below is the functions for Section 4 of the manuscript
### Delta_1 is called GEST: biased if MR.wrong | MZ.wrong
### Delta_2 is called IPW: biased if MW.wrong | MZ.wrong
### Delta_3 is called OR: biased if MW.wrong | MR.wrong | MY.wrong
### MLE is called MIAO: biased if MW.wrong | MY.wrong
###########################SIMULATION CODE###########################

#' @title Function for estimation
#' @description  The document of this function is yet to be written.
#'
#' @param data The original data. The data must contain columns `A`, `Z`, `W`, `Y`. `A` and `Y` respectively give the treatment arm and the outcome. An auxiliary exposure variable `Z` and an auxiliary outcome variable `W` should be provided as the double negative control variables.
#' @param arg A list of arguments.
#'
#' @return The estimation result, an object of class "proximalr_result", inheriting from the base-type "list". The members include three semiparametric estimators `GEST`\eqn{=\hat\Delta_1}, `IPW`\eqn{=\hat\Delta_2}, and `OR`\eqn{=\hat\Delta_3}, which operate under \eqn{\mathcal{M}_1, \mathcal{M}_2, \mathcal{M}_3}, respectively, the plug-in estimator discussed in Section 2.2.1 which we refer to as the MLE estimator hereafter, and the multiply robust (MR) estimator `MR`=\eqn{\hat\Delta_{\text{mr}}}.
#' @export
#'
#' @examples
#' est_func() # TODO
#'
est_func = function(data,arg=list(MW.wrong,MR.wrong,MY.wrong,MZ.wrong)) {
  A=data$A
  Z=data$Z
  W=data$W
  Y=data$Y
  X.org = as.matrix(data[,5:ncol(data)])
  MW.wrong=arg$MW.wrong;MR.wrong=arg$MR.wrong;
  MY.wrong=arg$MY.wrong;MZ.wrong=arg$MZ.wrong;
  #ATE=arg$ATE
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
    MR=MR.est#,
  #  ATE.true = ATE
  )
  class(est) = c(class(est), "proximalr_result")
  return(est)
}

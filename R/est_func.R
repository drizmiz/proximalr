
#' @importFrom Rdpack reprompt

#' @title Function for estimation
#' @description  The function for calculating semiparametric estimators of \eqn{\Delta} in the case of binary unmeasured confounding and negative control variables.
#'
#' @param data The original data. The data must contain columns `A`, `Z`, `W`, `Y`. `A` and `Y` respectively give the treatment arm and the outcome. An auxiliary exposure variable `Z` and an auxiliary outcome variable `W` should be provided as the double negative control variables. The rest of the columns are known covariates.
#' @param method The chosen estimator. There are six options, one of which, "all", indicates calculating all of the five estimators. The other options are respectively "IPW", "OR", "GEST", "MIAO", and "MR".
#' @param arg A list of arguments indicating the following different model misspecification scenarios.
#' - All models are correctly specified;
#' - \eqn{\mathcal{M}_2} and \eqn{\mathcal{M}_3} are wrong: \eqn{E[W | A,Z,X]} is misspecified by assuming that both \eqn{\xi^W_Z (A,X)} and \eqn{\delta^W_A (Z,X)} are constant;
#' - \eqn{\mathcal{M}_1} and \eqn{\mathcal{M}_3} are wrong: \eqn{R(A,X)} is misspecified by assuming that \eqn{R(A,X)} is a constant;
#' - \eqn{\mathcal{M}_1} and \eqn{\mathcal{M}_2} are wrong: \eqn{f(Z | A,X)} is misspecified by omitting the interaction term \eqn{X_7X_8};
#' - All models are wrong: \eqn{f(Z | A,X)} and \eqn{E[Y|A,Z,X]} are misspecified by omitting the interaction term \eqn{X_7X_8}.
#' @param MW.wrong \eqn{E[W | A,Z,X]} wrong and deltaW in IPW misspecify
#' @param MR.wrong \eqn{R(A,X)} wrong
#' @param MY.wrong \eqn{f(Z | A,X)} wrong
#' @param MZ.wrong \eqn{E[Y|A,Z=0,X]} wrong
#'
#' @return The estimation result, an object of class "proximalr_result", inheriting from the base-type "list". Five estimators of the ATE and their variances are included, among which are three semiparametric estimators `GEST`\eqn{=\hat\Delta_1}, `IPW`\eqn{=\hat\Delta_2}, and `OR`\eqn{=\hat\Delta_3}, which operate under \eqn{\mathcal{M}_1, \mathcal{M}_2, \mathcal{M}_3}, respectively, the plug-in estimator discussed in Section 2.2.1 of (arXiv:1808.04906v3 \[stat.ME\]), which referred to as the MLE estimator hereafter, and the multiply robust (MR) estimator `MR`=\eqn{\hat\Delta_{\text{mr}}}.
#'
#' @references {
#'   \insertRef{Binary}{proximalr}
#' }
#'
#' @export
#'
#' @examples
#' ## The arguments bin_data and bin_args are used as an example.
#' result = est_func(data = bin_data, arg = bin_args)
#' unlist(result)
#' ##     GEST.ATE      GEST.var       IPW.ATE       IPW.var        OR.ATE        OR.var
#' ##   0.27478357    0.01936722    0.64118222   10.85614022    0.39589779    1.65221426
#' ##     MIAO.ATE      MIAO.var       MR.ATE1       MR.ATE2       MR.var1       MR.var2
#' ##   0.35922107    0.24132369   -0.93999111    0.92476932 5980.38682647    5.36019154
#'
est_func = function(data, method="all", arg=list(MW.wrong,MR.wrong,MY.wrong, MZ.wrong)) {
  A=data$A
  Z=data$Z
  W=data$W
  Y=data$Y
  X.org = as.matrix(data[,5:ncol(data)])
  MW.wrong=arg$MW.wrong;MR.wrong=arg$MR.wrong;
  MY.wrong=arg$MY.wrong;MZ.wrong=arg$MZ.wrong;
  #ATE=arg$ATE
  ## run all methods
  if(method == "all"){
    (IPW.est = IPW(aa=A,zz=Z,ww=W,yy=Y,X.org,MZ.wrong,MW.wrong))
    (OR.est = OR(aa=A,zz=Z,ww=W,yy=Y,X.org,MW.wrong,MR.wrong,MY.wrong))
    (GEST.est = GEST(aa=A,zz=Z,ww=W,yy=Y,X.org,MZ.wrong,MR.wrong))
    (MIAO.est = MIAO(aa=A,zz=Z,ww=W,yy=Y,X.org,MW.wrong,MY.wrong))
    (MR.est = MR(aa=A,zz=Z,ww=W,yy=Y,X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong))

    est = c(
      GEST=GEST.est,
      IPW=IPW.est,
      OR=OR.est,
      MIAO=MIAO.est,
      MR=MR.est#,
      #  ATE.true = ATE
    )
  } else if(method == "IPW"){
    ### IPW
    est = (IPW.est = IPW(aa=A,zz=Z,ww=W,yy=Y,X.org,MZ.wrong,MW.wrong))
  } else if(method == "OR"){
    ### OR
    est = (OR.est = OR(aa=A,zz=Z,ww=W,yy=Y,X.org,MW.wrong,MR.wrong,MY.wrong))
  } else if(method == "GEST"){
    ### GEST
    est = (GEST.est = GEST(aa=A,zz=Z,ww=W,yy=Y,X.org,MZ.wrong,MR.wrong))
  } else if(method == "MIAO"){
    ### MIAO
    est = (MIAO.est = MIAO(aa=A,zz=Z,ww=W,yy=Y,X.org,MW.wrong,MY.wrong))
  } else if(method == "MR"){
    ### MR
    est = (MR.est = MR(aa=A,zz=Z,ww=W,yy=Y,X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong))
  } else{
    stop(paste0("Error from est_func. Invaild `method`, value = ", method, "."))
  }

  class(est) = c(class(est), "proximalr_result")
  return(est)
}

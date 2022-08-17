
#' @export
test_once = function(n,p,alpha,beta,MW.wrong,MR.wrong,MY.wrong,MZ.wrong){
  ## generate data
  X.org = rep(1,n)
  for(i in 1:(p-2)){
    #X.org=cbind(X.org,rnorm(n,0,1))
    X.org=cbind(X.org,runif(n))
  }
  X.org = cbind(X.org,X.org[,ncol(X.org)-1]*X.org[,ncol(X.org)])
  A.X = expit(X.org%*%alpha)
  B.X = expit(X.org%*%beta)
  A <- rbinom(n,size=1,prob=A.X)
  Z <- rbinom(n,size=1,prob=expit(tt*A+X.org%*%alpha))
  U <- rbinom(n,size=1,prob=(ma*A+mb*Z+mc*A*Z)/me  )
  W <- rbinom(n,size=1,prob=me*U+B.X)
  Y1 <- rbinom(n,size=1,prob=mr0*me*U+mr1*me*1*U+(mf-ma*(mr0+mr1))*1+B.X)
  Y0 <- rbinom(n,size=1,prob=mr0*me*U+mr1*me*0*U+(mf-ma*(mr0+mr1))*0+B.X)
  Y = A*Y1+(1-A)*Y0
  ## monte carlo ATE
  ATE = mean(Y1)-mean(Y0)
  data = data.frame(A=A,Z=Z,W=W,Y=Y,X.org)
  if(length(which(is.na(W)))>0){
    ##in case P(W|U,X)=me*U+B.X is beyond (0,1)
    data=data[!is.na(W),];n=nrow(data)
  }
  #### run all methods for once
  est=try(
    est_func(
      orig.data=data,indices=1:nrow(data),
      arg=list(MW.wrong=MW.wrong,MR.wrong=MR.wrong,MY.wrong=
                 MY.wrong,MZ.wrong=MZ.wrong,ATE=ATE)
    ),silent=TRUE)
  if(!inherits(est, "try-error")){
    return(list(est=c(est,n=n)))
  }else{
    return(list(est=NULL))
  }
}

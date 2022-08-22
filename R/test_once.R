
generate_data = function(n,p,alpha,beta) {

  ## parameters
  dd = 0.1
  tt = -0.2
  ma=0#0.2
  mb=0.2
  mc=0.2
  mr0=0#0.2
  mr1=0.5 # mr1<1 mr0+mr1>=0
  #mr1=0
  cc=0.5
  (me= ma+mb+mc+dd) ## +0.1 to make sure e-a>b+c>0
  mf = 0#1-(me-ma)*(mr0+mr1)-cc ## -0.1 to make sure 1-(e-a)(r0+r1) > f which also gives 1-(b+c)(r0+r1) > f
  (mt = me-ma-mb-mc)

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

  return(
    data.frame(
      A, Z, W, Y1, Y0, X.org
    )
  )
}

test_once = function(n,p,alpha,beta,MW.wrong,MR.wrong,MY.wrong,MZ.wrong){

  data = generate_data(n,p,alpha,beta)

  #data$Y = A*data$Y1+(1-A)*data$Y0
  ## monte carlo ATE
  ATE = mean(data$Y1)-mean(data$Y0)

  data = data %>%
    mutate(Y = A*Y1 + (1-A)*Y0, .after = "W") %>% select(-c("Y0", "Y1"))

  #data = data.frame(A=A,Z=Z,W=W,Y=Y,X.org)
  if(length(which(is.na(data$W)))>0){
    ##in case P(W|U,X)=me*U+B.X is beyond (0,1)
    data=data[!is.na(data$W),];n=nrow(data)
  }

  est = est_func(data,
             #orig.data=data,#indices=1:nrow(data),
             arg=list(MW.wrong=MW.wrong,MR.wrong=MR.wrong,MY.wrong=
                        MY.wrong,MZ.wrong=MZ.wrong)
             )
  est$ATE.true = ATE

  return(list(est=c(est,n=n)))
}

simulation = function(arg1 = 1, arg2 = 1) {

  wrong.i = as.numeric(arg1)
  myseed =  as.numeric(arg2)
  #my.filepath = ""

  n.rep=4#number of replications per seed
  n=50
  all.wrong = rbind(c(F,F,F,F),c(T,F,F,F),c(F,T,F,F),c(F,F,T,F),c(F,F,T,T))
  all.wrong = all.wrong[wrong.i,]
  MW.wrong = all.wrong[1] ## E[W|AZX] wrong and deltaW in IPW misspecify
  MR.wrong = all.wrong[2] ## R wrong
  MZ.wrong = all.wrong[3] ## f(AZX) wrong
  MY.wrong = all.wrong[4]  ## E[Y|Z=0,AX] wrong #this one never used
  p=10 # 1 intercept, p-2 cov, 1 interaction
  alpha = c(-0.1,rep(-0.1,(p-2)/2),rep(-0.1,(p-2)/2),2)/p
  beta = c(-1,rep(-1/p,length.out=p-1))

  rslt = NULL
  for(i in 1:n.rep+(myseed-1)*n.rep){
    set.seed(i)

    tmp=test_once(n,p,alpha,beta,MW.wrong,MR.wrong,MY.wrong,MZ.wrong)
    rslt=rbind(rslt,unlist(tmp$est))

    # if(i%%10==0){
    #   print(i)
    # }
  }
  colnames(rslt)=names(tmp$est)

  return(list(
    result = rslt,
    arguments = list(MW.wrong, MZ.wrong, MR.wrong, MY.wrong, n, n.rep, myseed)
    ))
}

savesim = function() {

  tmp = simulation()
  ag = tmp$arguments

  my.filepath="./intermediate/"

  #cat(tmp$result)
  res = tmp$result

  save(res,
       file=paste0(my.filepath,
                   "rslt_W_", ag[1] ,"_Z_",ag[2],"_R_",ag[3],"_Y_",ag[4],"_n",ag[5],"_nrep",ag[6],"_seed",ag[7],".RData"))

}

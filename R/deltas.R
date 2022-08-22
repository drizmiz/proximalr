
IPW = function(aa,zz,ww,yy,X.org,MZ.wrong,MW.wrong){
  n=length(aa)
  if(MZ.wrong==T){X.Z = cbind(X.org[,1:(ncol(X.org)-1)])}else{X.Z = X.org}
  Obj.ipw = function(par,A,Z,Y,W,X.org){return(sum(apply(U.ipw(par,A,Z,Y,W,X.org,MZ.wrong,MW.wrong),2,sum)^2))}
  t=nlm(f=Obj.ipw,p=c(rep(0.1,ncol(cbind(aa,X.Z))+ncol(X.org)),0.2,0.2,0.1),
        A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  G.ipw = function(par,A,Z,Y,W,X.org){return(apply(U.ipw(par,A,Z,Y,W,X.org,MZ.wrong,MW.wrong),2,sum))}
  bread=numDeriv::jacobian(func=G.ipw,x=t$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  meat.half=U.ipw(par=t$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org,MZ.wrong,MW.wrong)
  IF = meat.half%*%t(solve(-bread))
  delta = t$estimate[length(t$estimate)]
  delta.var = sum(IF[,ncol(IF)]^2)
  return(

    list(
      ATE=as.numeric(delta),
      var=as.numeric(delta.var)
    )

  )
}

OR = function(aa,zz,ww,yy,X.org,MW.wrong,MR.wrong,MY.wrong){
  n=length(aa)
  if(MY.wrong==T){X.Y = cbind(X.org[,1:(ncol(X.org)-1)])}else{X.Y = X.org}
  Obj.or = function(par,A,Z,Y,W,X.org){return(sum(apply(U.or(par,A,Z,Y,W,X.org,MW.wrong,MR.wrong,MY.wrong),2,sum)^2))}
  if(MR.wrong==F){
    t=nlm(f=Obj.or,p=rep(0,(2+ncol(X.org)+ncol(X.Y)+2)),
          A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  }else{
    t=nlm(f=Obj.or,p=rep(0,(2+ncol(X.org)+ncol(X.Y)+3)),
          A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  }
  G.or = function(par,A,Z,Y,W,X.org){return(apply(U.or(par,A,Z,Y,W,X.org,MW.wrong,MR.wrong,MY.wrong),2,sum))}
  bread=numDeriv::jacobian(func=G.or,x=t$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  meat.half=U.or(par=t$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org,MW.wrong,MR.wrong,MY.wrong)
  IF = meat.half%*%t(solve(-bread))
  delta = t$estimate[length(t$estimate)]
  delta.var = sum(IF[,ncol(IF)]^2)

  return(

    list(
      ATE=as.numeric(delta),
      var=as.numeric(delta.var)
    )

  )
}

MIAO = function(aa,zz,ww,yy,X.org,MW.wrong,MY.wrong){
  n=length(aa)
  if(MY.wrong==T){X.Y = cbind(X.org[,1:(ncol(X.org)-1)])}else{X.Y = X.org}
  Obj.miao = function(par,A,Z,Y,W,X.org){return(sum(apply(U.miao(par,A,Z,Y,W,X.org,MW.wrong,MY.wrong),2,sum)^2))}
  t=nlm(f=Obj.miao,p=c(rep(0.1,2+ncol(X.org)+1+ncol(X.Y)+1)),
        A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  G.miao = function(par,A,Z,Y,W,X.org){return(apply(U.miao(par,A,Z,Y,W,X.org,MW.wrong,MY.wrong),2,sum))}
  bread=numDeriv::jacobian(func=G.miao,x=t$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  meat.half=U.miao(par=t$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org,MW.wrong,MY.wrong)
  IF = meat.half%*%t(solve(-bread))
  delta = t$estimate[length(t$estimate)]
  delta.var = sum(IF[,ncol(IF)]^2)

  return(

    list(
      ATE=as.numeric(delta),
      var=as.numeric(delta.var)
    )

  )
}

GEST = function(aa,zz,ww,yy,X.org,MZ.wrong,MR.wrong){
  n=length(aa)
  if(MZ.wrong==T){X.Z = cbind(X.org[,1:(ncol(X.org)-1)])}else{X.Z = X.org}
  Obj.gest = function(par,A,Z,Y,W,X.org){return(sum(apply(U.gest(par,A,Z,Y,W,X.org,MZ.wrong,MR.wrong),2,sum)^2))}
  if(MR.wrong==F){
    t=nlm(f=Obj.gest,p=c(rep(0,ncol(cbind(aa,X.Z))+ncol(X.org)+2)),
          A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  }else{
    t=nlm(f=Obj.gest,p=c(rep(0,ncol(cbind(aa,X.Z))+ncol(X.org)+3)),
          A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  }
  G.gest = function(par,A,Z,Y,W,X.org){return(apply(U.gest(par,A,Z,Y,W,X.org,MZ.wrong,MR.wrong),2,sum))}
  bread=numDeriv::jacobian(func=G.gest,x=t$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  meat.half=U.gest(par=t$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org,MZ.wrong,MR.wrong)
  IF = meat.half%*%t(solve(-bread))
  delta = t$estimate[length(t$estimate)]
  delta.var = sum(IF[,ncol(IF)]^2)
  return(

    list(
      ATE=as.numeric(delta),
      var=as.numeric(delta.var)
    )

  )
}

MR = function(aa,zz,ww,yy,X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong){
  n=length(aa)
  if(MZ.wrong==T){X.Z = cbind(X.org[,1:(ncol(X.org)-1)])}else{X.Z = X.org}
  if(MY.wrong==T){X.Y = cbind(X.org[,1:(ncol(X.org)-1)])}else{X.Y = X.org}
  Obj.mr = function(par,A,Z,Y,W,X.org){return(sum(apply(U.mr(par,A,Z,Y,W,X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong),2,sum)^2))}
  if(MR.wrong==F){
    t=nlm(f=Obj.mr,p=c(rep(0.1,2+ncol(X.org)+ncol(X.Y)+ncol(X.Z)+1+ncol(X.org)+2)),
          A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  }else{
    t=nlm(f=Obj.mr,p=c(rep(0.1,2+ncol(X.org)+ncol(X.Y)+ncol(X.Z)+1+ncol(X.org)+3)),
          A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  }
  G.mr = function(par,A,Z,Y,W,X.org){return(apply(U.mr(par,A,Z,Y,W,X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong),2,sum))}
  bread=numDeriv::jacobian(func=G.mr,x=t$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  meat.half=U.mr(par=t$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong)
  IF = meat.half%*%t(solve(-bread))
  delta = t$estimate[length(t$estimate)]
  delta.var = sum(IF[,ncol(IF)]^2)
  ##################################################################
  ### MR2 is an attempt to stablize MR estimator and potentially
  ### improve efficiency by fitting the full E[Y|AZX] (E[W|AZX]) model,
  ### but only taking the E[Y|AZ=0X] (E[W|A=0Z=0X]) part
  ### key line of code in U.mr2 is
  ### U3 = c(W-(PWAZX)) *cbind(AZ.W/as.numeric(expit(X.W%*%par.w)*(1-expit(X.W%*%par.w))), X.W)
  ### U4 = c(Y-(PYAZX)) *cbind(AZ.Y/as.numeric(expit(X.Y%*%par.y)*(1-expit(X.Y%*%par.y))), X.Y)
  ##################################################################
  Obj.mr2 = function(par,A,Z,Y,W,X.org){return(sum(apply(U.mr2(par,A,Z,Y,W,X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong),2,sum)^2))}
  if(MR.wrong==F){
    t2=nlm(f=Obj.mr2,p=c(rep(0.1,4+ncol(X.org)+1+ncol(X.Y)+ncol(X.Z)+1+ncol(X.org)+2)),
           A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  }else{
    t2=nlm(f=Obj.mr2,p=c(rep(0.1,4+ncol(X.org)+1+ncol(X.Y)+ncol(X.Z)+1+ncol(X.org)+3)),
           A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  }
  G.mr2 = function(par,A,Z,Y,W,X.org){return(apply(U.mr2(par,A,Z,Y,W,X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong),2,sum))}
  bread2=numDeriv::jacobian(func=G.mr2,x=t2$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org)
  meat.half2=U.mr2(par=t2$estimate,A=aa,Z=zz,Y=yy,W=ww,X.org=X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong)
  IF2 = meat.half2%*%t(solve(-bread2))
  delta2 = t2$estimate[length(t2$estimate)]
  delta2.var = sum(IF2[,ncol(IF2)]^2)

  return(

    list(
      ATE1=as.numeric(delta),
      ATE2=as.numeric(delta2),
      var1=as.numeric(delta.var),
      var2=as.numeric(delta2.var)
    )
  )
}

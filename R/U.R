
#### below are all estimating equations
U.mr = function(par,A,Z,Y,W,X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong){
  if(MZ.wrong==T){X.Z = cbind(X.org[,1:(ncol(X.org)-1)])}else{X.Z = X.org}
  AX.Z=cbind(A,X.Z); n=length(A); X.W = X.org
  if(MY.wrong==T){## E[Y|Z=0,AX] wrong
    X.Y = cbind(X.org[,1:(ncol(X.org)-1)])
  }else{
    X.Y = X.org
  }
  if(MW.wrong==T){
    mycov.delta.Wzx=cbind(rep(1,length(Z)))
    mycov.delta.Wax=cbind(rep(1,length(Z)))
    alphaZ=par[1]; alphaA=par[2]; alphaAZ=0; AZ.W = cbind(Z,A); par.waz=par[1:2]
  }else{
    mycov.delta.Wzx=cbind(rep(1,length(Z)))
    mycov.delta.Wax=cbind(rep(1,length(Z)))
    alphaZ=par[1]; alphaA=0; alphaAZ=par[2]; AZ.W = cbind(Z,A*Z); par.waz=par[1:2]
  }
  par.w = par[(2+1):(2+ncol(X.org))]
  par.y = par[(2+ncol(X.org)+1):(2+ncol(X.org)+ncol(X.Y))]
  par.z = par[(2+ncol(X.org)+ncol(X.Y)+1):(2+ncol(X.org)+ncol(X.Y)+ncol(AX.Z))]
  par.a = par[(2+ncol(X.org)+ncol(X.Y)+ncol(AX.Z)+1):(2+ncol(X.org)+ncol(X.Y)+ncol(AX.Z)+ncol(X.org))]
  t=(2+ncol(X.org)+ncol(X.Y)+ncol(AX.Z)+ncol(X.org))
  if(MR.wrong==F){
    par.r = par[t+1];A.R = cbind(A);A.R2 = cbind(1-A); A1.R=par.r; A0.R = 0
    delta=par[t+2]
  }else{
    par.r = par[(t+1):(t+2)];A.R = cbind(1,A);A.R2 = cbind(1,1-A); A1.R=sum(par.r); A0.R = par.r[1]
    delta=par[t+3]
  }
  PAX = expit(X.org%*%par.a)
  PZAX=expit(AX.Z%*%par.z)
  A2=1-A
  PZA2X=expit(cbind(A=A2,X.Z)%*%par.z)
  PZA1X=expit(cbind(A=1,X.Z)%*%par.z)
  PZA0X=expit(cbind(A=0,X.Z)%*%par.z)
  PZX = PZA1X*PAX+PZA0X*(1-PAX)
  PAZX = (PZA1X*Z+(1-PZA1X)*(1-Z)) * (PAX)  / (PZX*Z+(1-PZX)*(1-Z))
  fAZX = PAZX*A+(1-PAZX)*(1-A)
  fZAX = PZAX*Z+(1-PZAX)*(1-Z)
  fAX = PAX*A+(1-PAX)*(1-A); fA2X = PAX*A2+(1-PAX)*(1-A2)

  Ya1z0x=Ya0z0x=Yaz0x = expit(X.Y%*%par.y)
  R = A.R%*%par.r
  R2 = A.R2%*%par.r
  delta.Wax = alphaZ+alphaAZ*A
  delta.Wzx = alphaA+alphaAZ*Z
  delta.Yax = delta.Wax*R
  Yazx = delta.Yax*Z + Yaz0x
  Wa0z0x = expit(X.W%*%par.w)
  Wazx = alphaA*A+alphaZ*Z+alphaAZ*A*Z+Wa0z0x
  Waz0x = alphaA*A+Wa0z0x
  if(MW.wrong==T){
    condE.delta.Wzx = alphaA
  }else{
    condE.delta.Wzx = (cbind(Z=rep(1,n))%*%alphaAZ)*PZA2X+ (cbind(Z=rep(0,n))%*%alphaAZ)*(1-PZA2X)
  }
  condE.R = A1.R*(1-PAZX)+A0.R*PAZX

  odds.AX = fA2X/fAX
  D_R = (2*Z-1)/fZAX/delta.Wax*(Y-Yaz0x-R*(W-Waz0x))*
    condE.delta.Wzx*odds.AX
  D_delta.Wzx = (2*A-1)/fAZX*(W-Wazx)*condE.R
  bias = R2*delta.Wzx

  U1 = c(Z-expit(AX.Z%*%par.z)) *(AX.Z)
  U2 = c(A-expit(X.org%*%par.a)) *(X.org)
  U3 = c(W-(expit(X.W%*%par.w))) *cbind(X.W) *as.numeric(A==0&Z==0)
  ### note that par.waz.useless in U3 is useless: since we only estimate E[W|Z=0,A=0,X,beta^{W0}]
  U4 = c(Y-(expit(X.Y%*%par.y))) *cbind(X.Y) *as.numeric(Z==0)
  ### note that par.yaz.useless in U4 is useless: since we only estimate E[Y|Z=0,A,X,beta_{mle}]

  U5 = c(  (A-PAX)*(W-(expit(X.W%*%par.w)+alphaZ*Z+alphaA*A+alphaAZ*A*Z))  )*mycov.delta.Wzx
  U6 = c(  (Z-PZAX)*(W-(expit(X.W%*%par.w)+alphaZ*Z+alphaA*A+alphaAZ*A*Z))  )*mycov.delta.Wax
  U7 = c(  Y-Yazx  )*A.R

  U8 = delta - (
    (  (Ya1z0x-Ya0z0x + (alphaZ+alphaAZ*1)*A1.R*Z - (alphaZ+alphaAZ*0)*A0.R*Z ) + (2*A-1)/fAZX*(Y-Yazx)  ) -
      (D_R+D_delta.Wzx+bias)
  )

  U = cbind(U1,U2,U3,U4,U5,U6,U7,U8)
  return(U)
}

U.mr2 = function(par,A,Z,Y,W,X.org,MZ.wrong,MW.wrong,MR.wrong,MY.wrong){
  AZ.Y = cbind(A*Z)
  if(MZ.wrong==T){X.Z = cbind(X.org[,1:(ncol(X.org)-1)])}else{X.Z = X.org}
  AX.Z=cbind(A,X.Z); n=length(A); X.W = X.org
  if(MY.wrong==T){## E[Y|Z=0,AX] wrong
    X.Y = cbind(X.org[,1:(ncol(X.org)-1)])
  }else{
    X.Y = X.org
  }
  if(MW.wrong==T){
    mycov.delta.Wzx=cbind(rep(1,length(Z)))
    mycov.delta.Wax=cbind(rep(1,length(Z)))
    alphaZ=par[1]; alphaA=par[2]; alphaAZ=0; AZ.W = cbind(Z,A);
  }else{
    mycov.delta.Wzx=cbind(rep(1,length(Z)))
    mycov.delta.Wax=cbind(rep(1,length(Z)))
    alphaZ=par[1]; alphaA=0; alphaAZ=par[2]; AZ.W = cbind(Z,A*Z);
  }
  par.waz = par[(ncol(AZ.W)+1):(ncol(AZ.W)+2)] ##note that c(alphaZ,alphaA,alphaAZ) are same role as par.waz but estimated twice; par.waz is useless
  par.w = par[(4+1):(4+ncol(X.org))]
  par.yaz = par[(4+ncol(X.org)+1)]
  par.y = par[(4+ncol(X.org)+1+1):(4+ncol(X.org)+1+ncol(X.Y))]
  par.z = par[(4+ncol(X.org)+1+ncol(X.Y)+1):(4+ncol(X.org)+1+ncol(X.Y)+ncol(AX.Z))]
  par.a = par[(4+ncol(X.org)+1+ncol(X.Y)+ncol(AX.Z)+1):(4+ncol(X.org)+1+ncol(X.Y)+ncol(AX.Z)+ncol(X.org))]
  t=(4+ncol(X.org)+1+ncol(X.Y)+ncol(AX.Z)+ncol(X.org))
  if(MR.wrong==F){
    par.r = par[t+1];A.R = cbind(A);A.R2 = cbind(1-A); A1.R=par.r; A0.R = 0
    delta=par[t+2]
  }else{
    par.r = par[(t+1):(t+2)];A.R = cbind(1,A);A.R2 = cbind(1,1-A); A1.R=sum(par.r); A0.R = par.r[1]
    delta=par[t+3]
  }
  PAX = expit(X.org%*%par.a)
  PZAX=expit(AX.Z%*%par.z)
  A2=1-A
  PZA2X=expit(cbind(A=A2,X.Z)%*%par.z)
  PZA1X=expit(cbind(A=1,X.Z)%*%par.z)
  PZA0X=expit(cbind(A=0,X.Z)%*%par.z)
  PZX = PZA1X*PAX+PZA0X*(1-PAX)
  PAZX = (PZA1X*Z+(1-PZA1X)*(1-Z)) * (PAX)  / (PZX*Z+(1-PZX)*(1-Z))
  fAZX = PAZX*A+(1-PAZX)*(1-A)
  fZAX = PZAX*Z+(1-PZAX)*(1-Z)
  fAX = PAX*A+(1-PAX)*(1-A); fA2X = PAX*A2+(1-PAX)*(1-A2)

  Ya1z0x=Ya0z0x=Yaz0x = expit(X.Y%*%par.y)
  R = A.R%*%par.r
  R2 = A.R2%*%par.r
  delta.Wax = alphaZ+alphaAZ*A
  delta.Wzx = alphaA+alphaAZ*Z
  delta.Yax = delta.Wax*R
  Yazx = delta.Yax*Z + Yaz0x
  Wa0z0x = expit(X.W%*%par.w)
  Wazx = alphaA*A+alphaZ*Z+alphaAZ*A*Z+Wa0z0x
  Waz0x = alphaA*A+Wa0z0x
  if(MW.wrong==T){
    condE.delta.Wzx = alphaA
  }else{
    condE.delta.Wzx = (cbind(Z=rep(1,n))%*%alphaAZ)*PZA2X+ (cbind(Z=rep(0,n))%*%alphaAZ)*(1-PZA2X)
  }
  condE.R = fA2X/fAX
  odds.AX = (1-expit(X.org%*%par.a))/expit(X.org%*%par.a)
  D_R = (2*Z-1)/fZAX/delta.Wax*(Y-Yaz0x-R*(W-Waz0x))*
    condE.delta.Wzx*odds.AX
  D_delta.Wzx = (2*A-1)/fAZX*(W-Wazx)*condE.R
  bias = R2*delta.Wzx


  PWAZX=AZ.W%*%par.waz+expit(X.W%*%par.w)
  PYAZX=AZ.Y%*%par.yaz+expit(X.Y%*%par.y)

  U1 = c(Z-expit(AX.Z%*%par.z)) *(AX.Z)
  U2 = c(A-expit(X.org%*%par.a)) *(X.org)
  U3 = c(W-(PWAZX)) *cbind(AZ.W/as.numeric(expit(X.W%*%par.w)*(1-expit(X.W%*%par.w))), X.W)
  U4 = c(Y-(PYAZX)) *cbind(AZ.Y/as.numeric(expit(X.Y%*%par.y)*(1-expit(X.Y%*%par.y))), X.Y)
  #### rather than the following in U.mr()
  # U3 = c(W-(expit(X.W%*%par.w))) *cbind(X.W) *as.numeric(A==0&Z==0)
  # ### note that par.waz.useless in U3 is useless: since we only take E[W|Z=0,A=0,X,beta^{W0}]
  # U4 = c(Y-(expit(X.Y%*%par.y))) *cbind(X.Y) *as.numeric(Z==0)
  # ### note that par.yaz.useless in U4 is useless: since we only take E[Y|Z=0,A,X,beta_{mle}]

  U5 = c(  (A-PAX)*(W-(expit(X.W%*%par.w)+alphaZ*Z+alphaA*A+alphaAZ*A*Z))  )*mycov.delta.Wzx
  U6 = c(  (Z-PZAX)*(W-(expit(X.W%*%par.w)+alphaZ*Z+alphaA*A+alphaAZ*A*Z))  )*mycov.delta.Wax
  U7 = c(  Y-Yazx  )*A.R

  U8 = delta - (
    (  (Ya1z0x-Ya0z0x + (alphaZ+alphaAZ*1)*A1.R*Z - (alphaZ+alphaAZ*0)*A0.R*Z ) + (2*A-1)/fAZX*(Y-Yazx)  ) -
      (D_R+D_delta.Wzx+bias)
  )

  U = cbind(U1,U2,U3,U4,U5,U6,U7,U8)
  return(U)
}

U.ipw = function(par,A,Z,Y,W,X.org,MZ.wrong,MW.wrong){
  if(MZ.wrong==T){X.Z = cbind(X.org[,1:(ncol(X.org)-1)])}else{X.Z = X.org}
  AX.Z=cbind(A,X.Z)
  par.z = par[1:ncol(AX.Z)]
  par.a = par[(ncol(AX.Z)+1):(ncol(AX.Z)+ncol(X.org))]
  PAX = expit(X.org%*%par.a)
  PZAX=expit(AX.Z%*%par.z)
  A2=1-A
  PZA2X=expit(cbind(A=A2,X.Z)%*%par.z)
  PZA1X=expit(cbind(A=1,X.Z)%*%par.z)
  PZA0X=expit(cbind(A=0,X.Z)%*%par.z)
  PZX = PZA1X*PAX+PZA0X*(1-PAX)
  PAZX = (PZA1X*Z+(1-PZA1X)*(1-Z)) * (PAX)  / (PZX*Z+(1-PZX)*(1-Z))
  fAZX = PAZX*A+(1-PAZX)*(1-A)
  fZAX = PZAX*Z+(1-PZAX)*(1-Z)

  par.w = par[(ncol(AX.Z)+ncol(X.org)+1):(ncol(AX.Z)+ncol(X.org)+2)]
  delta = par[ncol(AX.Z)+ncol(X.org)+3]
  if(MW.wrong==T){
    alphaZ=par.w[1]; alphaA=par.w[2]; alphaAZ=0
  }else{
    alphaZ=par.w[1]; alphaA=0; alphaAZ=par.w[2]
  }
  D.AY = (2*A-1)*Y/(fAZX)
  D.ZY = (2*Z-1)*Y/(fZAX)
  fAX = PAX*A+(1-PAX)*(1-A); fA2X = PAX*A2+(1-PAX)*(1-A2)

  U1 = c(Z-expit(AX.Z%*%par.z)) *(AX.Z)
  U2 = c(A-expit(X.org%*%par.a)) *(X.org)
  U3 = (  (A-expit(X.org%*%par.a))*(W-(alphaA*A+alphaZ*Z+alphaAZ*A*Z))  ) * cbind(rep(1,length(Z)))
  U4 = (  (Z-expit(AX.Z%*%par.z))* (W-(alphaA*A+alphaZ*Z+alphaAZ*A*Z))  ) * cbind(rep(1,length(Z)))
  U5 = rep(delta,length(Y)) - (  D.AY - D.ZY/(alphaZ+alphaAZ*(A))*(alphaA+alphaAZ*(PZA2X))* fA2X/fAX  )
  U = cbind(U1,U2,U3,U4,U5)
  return(U)
}

U.gest = function(par,A,Z,Y,W,X.org,MZ.wrong,MR.wrong){
  if(MZ.wrong==T){X.Z = cbind(X.org[,1:(ncol(X.org)-1)])}else{X.Z = X.org}
  AX.Z=cbind(A,X.Z)
  par.z = par[1:ncol(AX.Z)]
  par.a = par[(ncol(AX.Z)+1):(ncol(AX.Z)+ncol(X.org))]
  t=(ncol(AX.Z)+ncol(X.org))
  if(MR.wrong==F){
    par.r = par[t+1];A.R = cbind(A);A.R2 = cbind(1-A); A1.R=par.r; A0.R = 0
    delta=par[t+2]
  }else{
    par.r = par[(t+1):(t+2)];A.R = cbind(1,A);A.R2 = cbind(1,1-A); A1.R=sum(par.r); A0.R = par.r[1]
    delta=par[t+3]
  }
  PAX = expit(X.org%*%par.a)
  PZAX=expit(AX.Z%*%par.z)
  A2=1-A
  PZA2X=expit(cbind(A=A2,X.Z)%*%par.z)
  PZA1X=expit(cbind(A=1,X.Z)%*%par.z)
  PZA0X=expit(cbind(A=0,X.Z)%*%par.z)
  PZX = PZA1X*PAX+PZA0X*(1-PAX)
  PAZX = (PZA1X*Z+(1-PZA1X)*(1-Z)) * (PAX)  / (PZX*Z+(1-PZX)*(1-Z))
  fAZX = PAZX*A+(1-PAZX)*(1-A)
  fZAX = PZAX*Z+(1-PZAX)*(1-Z)

  U1 = c(Z-expit(AX.Z%*%par.z)) *(AX.Z)
  U2 = c(A-expit(X.org%*%par.a)) *(X.org)
  U3 = c( (Z-expit(AX.Z%*%par.z) )*(Y-W*(A.R%*%par.r))  )*A.R
  D.AY = (2*A-1)*Y/(fAZX)
  D.AW = (2*A-1)*W/(fAZX)
  R2_cond_ZX = (A0.R*PAZX)+(A1.R*(1-PAZX))

  U4 = delta - (  D.AY - (R2_cond_ZX)*D.AW  )
  U = cbind(U1,U2,U3,U4)
  return(U)
}

U.miao = function(par,A,Z,Y,W,X.org,MW.wrong,MY.wrong){
  AZ.Y = cbind(A*Z)
  if(MY.wrong==T){
    X.Y = cbind(X.org[,1:(ncol(X.org)-1)])
  }else{
    X.Y = X.org
  }
  X.W = X.org
  if(MW.wrong==T){
    alphaZ=par[1]; alphaA=par[2]; alphaAZ=0; AZ.W = cbind(Z,A); par.waz=par[1:2]
  }else{
    alphaZ=par[1]; alphaA=0; alphaAZ=par[2]; AZ.W = cbind(Z,A*Z); par.waz=par[1:2]
  }
  par.w = par[(2+1):(2+ncol(X.org))]
  par.yaz=par[(2+ncol(X.org)+1):(2+ncol(X.org)+ncol(AZ.Y))]
  par.y = par[(2+ncol(X.org)+ncol(AZ.Y)+1):(2+ncol(X.org)+ncol(AZ.Y)+ncol(X.Y))]
  delta = par[(2+ncol(X.org)+ncol(AZ.Y)+ncol(X.Y)+1)]

  PWAZX=AZ.W%*%par.waz+expit(X.W%*%par.w)
  PYAZX=AZ.Y%*%par.yaz+expit(X.Y%*%par.y)

  U1 = c(W-(PWAZX)) *cbind(AZ.W/as.numeric(expit(X.W%*%par.w)*(1-expit(X.W%*%par.w))), X.W)
  U2 = c(Y-(PYAZX)) *cbind(AZ.Y/as.numeric(expit(X.Y%*%par.y)*(1-expit(X.Y%*%par.y))), X.Y)
  U3 = delta - (
    Z*par.yaz[1] -
      ((1-A)*1*par.yaz[1]) / (alphaZ+alphaAZ*(1-A)) * (alphaA+alphaAZ*Z)
  )

  U = cbind(U1,U2,U3)
  return(U)
}

U.or = function(par,A,Z,Y,W,X.org,MW.wrong,MR.wrong,MY.wrong){
  AZ.Y = cbind(A*Z)
  if(MY.wrong==T){
    X.Y = cbind(X.org[,1:(ncol(X.org)-1)])
  }else{
    X.Y = X.org
  }
  X.W = X.org
  if(MW.wrong==T){
    alphaZ=par[1]; alphaA=par[2]; alphaAZ=0; AZ.W = cbind(Z,A); par.waz=par[1:2]
  }else{
    alphaZ=par[1]; alphaA=0; alphaAZ=par[2]; AZ.W = cbind(Z,A*Z); par.waz=par[1:2]
  }
  par.w = par[(2+1):(2+ncol(X.org))]
  par.y = par[(2+ncol(X.org)+1):(2+ncol(X.org)+ncol(X.Y))]
  t=(2+ncol(X.org)+ncol(X.Y))
  if(MR.wrong==F){
    par.r = par[t+1];A.R = cbind(A);A.R2 = cbind(1-A); A1.R=par.r*1; A0.R = par.r*0
    delta=par[t+2]
  }else{
    par.r = par[(t+1):(t+2)];A.R = cbind(1,A);A.R2 = cbind(1,1-A); A1.R=sum(par.r); A0.R = par.r[1]
    delta=par[t+3]
  }
  delta.gcomp = delta.Yzx = (alphaZ+alphaAZ*1)*A1.R*Z - (alphaZ+alphaAZ*0)*A0.R*Z
  delta.Wzx = (alphaA+alphaAZ*Z)
  R.or = A.R2%*%par.r
  delta.bias = (R.or*delta.Wzx)
  PWAZX=AZ.W%*%par.waz+expit(X.W%*%par.w)

  U1 = c(W-(PWAZX)) *cbind(AZ.W/as.numeric(expit(X.W%*%par.w)*(1-expit(X.W%*%par.w))), X.W)
  U2 = c(Y-(expit(X.Y%*%par.y))) *cbind(X.Y) *as.numeric(Z==0)
  # ### note that par.yaz is useless so did not show up in U2: since we only take E[Y|Z=0,A,X,beta_{mle}]
  U3 = c(  Y-expit(X.Y%*%par.y)-(A.R%*%par.r)*(W-(alphaA*A+expit(X.W%*%par.w)))  )*(A.R)
  U4 = delta - (
    delta.gcomp-delta.bias
  )
  U = cbind(U1,U2,U3,U4)
  return(U)
}

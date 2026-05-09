R"""
require(mvtnorm)
require(configural) # use of vec2hsfull, vechfull (vectorization and inverse vectorization)


##
##  We assume that all covariance/correlation functions
##  are written as
##  c(x) = ctilde(||x||/phi)  <==> c1(r)=ctilde(r/phi)
##  where phi is called a scale parameter
##

## Just to simplify the call to Bessel functions of order 0,2,4,6,8
J0=function(r) besselJ(r,0);J2=function(r) besselJ(r,2)
J4=function(r) besselJ(r,4);J6=function(r) besselJ(r,6)
J8=function(r) besselJ(r,8)

diff.c2=function(r,phi,order=0,which.cov="Gaussian",d=2,nu=NULL){
  #### ---------------------------------------------------------
  ### Provides the derivatives of the covariance function c2
  ### Arguments
  # r: value at which we want to compute the derivative
  # phi: scale parameter of the covariance function
  # order: order of the derivative to compute
  # which.cov: covariance model
  # d: dimension of the model
  # nu: regularity parameter, for the Matérn covariance 
  ##### Output
  # c2^(order)(r) = (order)-th derivative of c2 given by
  # c2(||x||^2) = E(X(0)X(x))
  ### ----------------------------------------------------------
  if ((which.cov=="Matern") & is.null(nu)) stop('nu must be specified')
  
  switch(which.cov,
    "Gaussian"={
      tmp=2*phi^2
      c1=exp(-r/tmp)
      out=(-1)^(order)*c1/tmp^(order)
      },
    "Matern"={
      z=sqrt(r)*sqrt(2*nu)/phi
      cst=(-1)^order*2/2^order/phi^(2*order) * nu^order / gamma(nu)
      i=which(z==0);j=which(z!=0)
      at.i= gamma(nu-order)/2
      at.j=(z[j]/2)^(nu-order)*besselK(z[j],nu-order)
      out=matrix(0,length(z),length(order))
      out[i,]=cst*at.i
      out[j,]=cst*at.j
      if (length(z)==1) out=c(out)
      },
    "RWM"={
      cst=(-1)^order/2^order/phi^(2*order)*(d/2)^order*gamma(d/2)
      z=sqrt(r)/(phi/sqrt(d))
      i=which(z==0);j=which(z!=0)
      at.i=1/gamma(d/2+order)
      at.j=(z[j]/2)^(-(d/2-1+order)) * besselJ(z[j],d/2-1+order)
      out=matrix(0,length(z),length(order))
      out[i,]=cst*at.i
      out[j,]=cst*at.j
      if (length(z)==1) out=c(out)
    }
  )
  out
}

lambda2p=function(phi,p=1,which.cov="Gaussian",d=2,nu=NULL,r=NULL){
  if (is.null(r)) r=0
  c2p0=diff.c2(r=r,phi,order=p,which.cov=which.cov,d=d,nu=nu)
  out=prod(1:(2*p))/prod(1:p) *(-1)^p*c2p0
  out
}

th.lambda2p=function(phi,p=1,which.cov="Gaussian",d=2,nu=NULL){
  switch(which.cov,
        "Gaussian"={prod(1:(2*p))/prod(1:p)/2^p/phi^(2*p)},
        "Matern"={prod(1:(2*p))/prod(1:p)/2^p/phi^(2*p) * prod(nu/(nu-1:p)) },
        "RWM"={prod(1:(2*p))/prod(1:p)/2^p/phi^(2*p) * prod(d/(d+2*(0:(p-1))))}
  )
}

practical.range=function(rho,phi=NULL,which.crit='max',which.cov="Gaussian",
                         d=2,val=.05,nu=NULL){
  switch(which.cov,
         "Gaussian"={sqrt( -log(val)*(2*phi^2))},
         "Matern"={
           pr=practical.range(rho=rho,phi=phi,which.crit=which.crit,which.cov="Gaussian",
                              d=d,val=val)
           f=function(x){
             (diff.c2(x^2,phi,order=0,which.cov="Matern",d=d,nu=nu)-val)^2
           }
           optimize(f,interval = c(.5*pr,5*pr))$min
         },
         "RWM"={ ## only for d=2, using asymptotic behaviour of J0
           delta=d/2-1
           2*phi*(gamma(delta+1)/(sqrt(pi)*val))^(1/(delta+1/2))
         })
}


rho2phi=function(rho,phi=NULL,d=2,which.cov="Gaussian",which.crit="max",nu=NULL){
### ------------------------------------------------------------------
### Description
### Gives the relation between the parameters rho and phi.
### Arguments:
### rho : real value
### phi : real value (not specified if rho is specified)
### which.cov : used model 
### d : dimension of the model (d=2,3 or 4)
### which.crit :  type of critical points to consider ("all" or "max")
### nu: regularity parameter for the Matern covariance
### Details
### Sets the scale parameter phi for a given rho, 
### or rho for a given phi
### Outputs:  A list with the two values (rho, phi).
### ------------------------------------------------------------------  
  if (is.null(rho) & is.null(phi)) stop("rho or phi must be non-null")
  if (which.cov=="Matern" & is.null(nu)) stop("nu must be specified")
  if (!(d %in% 1:4)) stop("Not implemented for this dimension")
  if (!(which.crit %in% c("all","max"))) stop("Not implemented for this type of points (all or max)")
  if (d==4){ 
    y=rnorm(1e6,0,sqrt(1/3))
    I=mean(pnorm(y)*pnorm(sqrt(2)*y))
  }
  #browser()
  switch(d,
         {#d=1 Constant for max to be checked
           cst<-switch(which.crit,
                       "all"={sqrt(3)/pi},#sqrt(3/(2*pi))*E},  with E=sqrt(2/pi)
                       "max"={sqrt(3)/(2*pi)})#sqrt(3/(2*pi))*E*.5})       
         },
         { #d=2
           cst<-switch(which.crit,
                       "all"={2/(pi*sqrt(3))},
                       "max"={.25*2/(pi*sqrt(3))})       
         },
         { #d=3
           cst<-switch(which.crit,
                       "all"={1/(pi^2*sqrt(2))*116/(12*sqrt(6)) },
                       "max"={1/(pi^2*sqrt(2))*(29-6*sqrt(6))/(12*sqrt(6))})       
         },
         { #d=4
           cst<-switch(which.crit,
                       "all"={1/pi^2*200/(48*sqrt(3))},
                       "max"={1/pi^2*(I*100*pi-57)/(48*sqrt(3)*pi)})        
         })
  if (is.null(phi)){
    f=function(phi){
      C0=diff.c2(r=0,phi,order=1:2,which.cov=which.cov,d=d,nu=nu)
      l2=-2*C0[1];l4=12*C0[2]  ## Second and fourth spectral moments
      rho.phi=(l4/(3*l2))^(d/2)*cst
      (rho.phi-rho)^2
    }
    ff=Vectorize(f)
    phi=optimize(ff,interval=c(1e-6,10),tol=1e-50)$min
  }
  if (is.null(rho)) {
    C0=diff.c2(r=0,phi=phi,order=1:2,which.cov=which.cov,d=d,nu=nu)
    l2=-2*C0[1];l4=12*C0[2]  ## Second and fourth spectral moments
    rho=(l4/(3*l2))^(d/2)*cst
  }
  list(rho=rho,phi=phi)
}

SigmaNablaX0Xre1 <- function(r,phi, d=2,which.cov="Gaussian",nu=NULL){
  ### ----------------------------------------------------------------------------------------------------
  ### Provides Covariance matrix of the upper diagonal and diagonal of gradient X(0) and gradient X(r*e_1)
  ### Arguments
  ### r: value for which the covariance matrix is evaluated
  ### phi: scale parameter
  ### which.cov: used covariance model 
  ### d: dimension of the model (d=1,2,3 or 4)
  ### nu: regularity parameter
  ### Output: A square matrix with d*(d+1) rows.
  ### ----------------------------------------------------------------------------------------------------
  C0=matrix(diff.c2(0,phi=phi,order=,0:4,which.cov=which.cov,d=d,nu=nu))
  Cr=matrix(diff.c2(r^2,phi=phi,order=0:4,which.cov=which.cov,d=d,nu=nu))
  
  if (d==1){
    G1=12*C0[3]
    M= (12*Cr[3] + 8*r^2*Cr[4])^2
    G1= G1 + r^2*C0[2]/(2*(C0[2]^2 - (Cr[2] + 2*r^2*Cr[3] )^2))*M
    G3=12*Cr[3]
    G3=G3+48*r^2*Cr[4] + 16*r^4*Cr[5]
    G3=G3+r^2*(Cr[2] + 2*r^2*Cr[3])/(2*(C0[2]^2 -(Cr[2] + 2*r^2*Cr[3])^2))*M
    out=cbind(c(G1,G3),c(G3,G1))
  } else{
    M <- matrix(0, nrow = d, ncol= d  )
    M[-1,-1] <- 16*Cr[3]^2
    M[1,-1] <- M[-1,1] <- 4*Cr[3]*(12*Cr[3] + 8*r^2*Cr[4])
    M[1,1] <- (12*Cr[3] + 8*r^2*Cr[4])^2
    G1 <- matrix(4*C0[3], nrow=d, ncol=d)
    diag(G1)<- 12*C0[3]
    G1 <- G1 + r^2*C0[2]/(2*(C0[2]^2 - (Cr[2] + 2*r^2*Cr[3] )^2))*M
    G3 <- matrix(4*Cr[3],nrow = d, ncol = d)
    diag(G3) <- 12*Cr[3]
    G3[1,-1] <- G3[1,-1] + 8*r^2*Cr[4]
    G3[-1,1] <- G3[-1,1] + 8*r^2*Cr[4]
    G3[1,1] <- G3[1,1] + 48*r^2*Cr[4] + 16*r^4*Cr[5]
    G3 <- G3 + r^2*(Cr[2] + 2*r^2*Cr[3])/(2*(C0[2]^2 -(Cr[2] + 2*r^2*Cr[3])^2))*M
    G2<-G4<-matrix(0,nr=d*(d-1)/2,d*(d-1)/2)
    diag1=rep(4*C0[3] + 8*r^2*Cr[3]^2*C0[2]/(C0[2]^2 - Cr[2]^2),d-1)
    diag2=rep(4*C0[3],  (d-1)*(d-2)/2)
    diag(G2)=c(diag1,diag2)
    diagt1=rep(4*Cr[3] + 8*r^2*Cr[4] + 8*r^2*Cr[3]^2*Cr[2]/(C0[2]^2 - Cr[2]^2), d-1)
    diagt2=rep(4*Cr[3], (d-1)*(d-2)/2)
    diag(G4)<-c(diagt1,diagt2)
    k <- d*(d-1)/2
    l <- d*(d+1)
    out <- matrix(0, nrow = 2*d+2*k, ncol = 2*d+2*k )
    out[1:d,1:d] <- out[(d+k+1):(2*d+k),(d+k+1):(2*d+k)] <- G1
    out[(d+1):(d+k),(d+1):(d+k)]<-out[(2*d+k+1):(2*d+2*k),(2*d+k+1):(2*d+2*k)]<-G2
    out[(1:d),(d+k+1):(2*d+k)]<-out[(d+k+1):(2*d+k),1:d]<-G3
    out[(d+1):(d+k),(2*d+k+1):(2*d+2*k)]<-out[(2*d+k+1):(2*d+2*k),(d+1):(d+k)]<-G4
    }
  out
}

densGrX0Xre1 <- function(r,phi, d=2,which.cov="Gaussian",nu=NULL){
  ### -----------------------------------------------------------------
  ### Provides the density at 0 of V(r) = (nabla X(0), nabla X(re1))
  ### See  proof of Lemma 4.1, section B.2.1.1 of Azais and Delmas 2022
  ### Arguments
  ### r: value for which the density of V(r) at 0 is evaluated
  ### phi: scale parameter
  ### d: dimension
  ### which.cov: used covariance model
  ### nu: eventually regularity parameter
  ### Output: numerical value
  ### ----------------------------------------------------------------- 
  C0=diff.c2(0,phi=phi,order=1,which.cov=which.cov,d=d,nu=nu) # c2^(1)(0)
  Cr=diff.c2(r^2,phi=phi,order=1:2,which.cov=which.cov,d=d,nu=nu) # c2^(o)(r^2) o=1,2
  
  covGrX0GrXt=matrix(0,nr=2*d,nc=2*d)
  diag(covGrX0GrXt)=-2*C0
  M=matrix(0,nr=d,nc=d)
  diag(M)=-2*Cr[1]
  M[1,1]=-2*Cr[1]-4*r^2*Cr[2]
  covGrX0GrXt[1:d,(d+1):(2*d)]<-covGrX0GrXt[(d+1):(2*d),1:d]<-M
  dmvnorm(x =rep(0,2*d), mean= rep(0,2*d), sigma = covGrX0GrXt)
}

g <- function(vecr,rho,phi=NULL, d=2,which.cov="Gaussian",
              which.crit = "max",B=1e4,do.parallel=FALSE,nu=NULL){
  ### --------------------------------------------------------------------------------------------
  ### Monte-Carlo approximation of the pair correlation function g(r) for different values 
  ### of r
  ### Arguments: 
  ### vecr: vector of values for which we approximate the pcf 
  ### rho: intensity parameter
  ### phi: scale parameter
  ### d: dimension 1,2,3 or 4
  ### which.cov: used covariance model
  ### which.crit: type of critical points
  ### B: number of MC replications 
  ### do.parallel: logical parameter specifying the eventual use of parallel computations
  ### nu: eventually regularity parameter
  ### Output: mainly matrix (2,length(vecr)), 1st row contains mean of Monte-Carlo approximations,
  ### 2nr row contains the standard error of the MC approximation
  ### --------------------------------------------------------------------------------------------
  tmp=rho2phi(rho=rho,phi=phi,d=d,which.cov=which.cov,which.crit=which.crit,nu=nu)
  rho=tmp$rho;phi=tmp$phi
  lr <- length(vecr)
  out <- matrix(0, 2, lr)
  det.minor=function(v,m,d){
    if (m>d) stop('m must be <=d')
    ## det of minor of H_m m=1,2,...,d of the B Hessians H 
    ## where H= is the inverse symmetric vectorization of v
    ## of length d + d(d-1)/2
   
    det2ns=function(v,i){#i given by row
      v[,i[1]]*v[,i[4]]-v[,i[2]]*v[,i[3]]
    } 
    det3ns=function(v,i){#i given by row
      a1= v[,i[1]]*det2ns(v,c(i[5],i[6],i[8],i[9]))
      a2=-v[,i[2]]*det2ns(v,c(i[4],i[6],i[7],i[9]))
      a3= v[,i[3]]*det2ns(v,c(i[4],i[5],i[7],i[8]))
      a1+a2+a3
    }
    det4ns=function(v,i){#i given by row
      a1= v[,i[1]]*det3ns(v,c(i[6],i[7],i[8],i[10],i[11],i[12],i[14],i[15],i[16]))
      a2=-v[,i[2]]*det3ns(v,c(i[5],i[7],i[8],i[9],i[11],i[12],i[13],i[15],i[16]))
      a3= v[,i[3]]*det3ns(v,c(i[5],i[6],i[8],i[9],i[10],i[12],i[13],i[14],i[16]))
      a4=-v[,i[4]]*det3ns(v,c(i[5],i[6],i[7],i[9],i[10],i[11],i[13],i[14],i[15]))
      a1+a2+a3+a4
    }
    if (d==1){out=v[,1]}
    if (d==2){out=switch(m,{v[,1]},{det2ns(v,c(1,3,3,2))})
    }
    if (d==3){out=switch(m,{v[,1]},{det2ns(v,c(1,4,4,2))},
                         {det3ns(v,c(1,4,5,4,2,6,5,6,3))})
    }
    if (d==4){out=switch(m,{v[,1]},{det2ns(v,c(1,5,5,2))},
                         {det3ns(v,c(1,5,6,5,2,8,6,8,3))},
                         {det4ns(v,c(1,5,6,7,5,2,8,9,6,8,3,10,7,9,10,4))})
    }
    out
  }
  
  pb <- txtProgressBar(min = 1, max = lr, style = 3)
  dd=d+d*(d-1)/2 ## Dimension of Hessian matrices
  N01=matrix(rnorm(B*2*dd),nr=2*dd,nc=B) # Generations of N(0,1) identical for all r's 
  if (do.parallel){
    totalCores = detectCores()
    cat('totalCores=',totalCores,'\n')
    cluster <- makeCluster(10)#,type="FORK")
    print(cluster)
  }
  eig.val.mat=NULL
  for (i in 1:lr) {
    setTxtProgressBar(pb, i)
    r=vecr[i]
    S=SigmaNablaX0Xre1(r=r,phi=phi,d=d,which.cov=which.cov,nu=nu)
    eigS=eigen(S);D=matrix(0,nrow(S),ncol(S));
    eigS.val=eigS$values
    if (!(all(eigS.val>0))) {
      eig.val=c(sum(eigS.val<0),mean(eigS.val[eigS.val<0]))
    } else eig.val=c(NA,NA)
    eig.val.mat=rbind(eig.val.mat,eig.val)
    diag(D)=sqrt(pmax(0,eigS.val) );P=eigS$vectors
    V=t(P%*%D%*%N01) # A (B,2*dd) matrix of N(0,S) replications
    #tL=t(chol(S))
    #V=t( tL%*%N01)  
    V0=as.matrix(V[,1:dd])    # Inverse symmetric vectorization of Hessians of X(0)
    V1=as.matrix(V[,-(1:dd)]) # Inverse symmetric vectorization of Hessians of X(re1)
    pGrX0Xre1=densGrX0Xre1(r,phi=phi, d=d,which.cov=which.cov,nu=nu)
    switch(which.crit, 
           "all"={
            adet0=abs(det.minor(V0,d,d))
            adet1=abs(det.minor(V1,d,d))
            pcfs = pGrX0Xre1*adet0*adet1 / rho^2
          },
          "max"={
            pcfs=rep(0,B)
            mat.det0.max<-mat.det1.max<-NULL
            mat.det0.min<-mat.det1.min<-NULL
            for (m in 1:d){
              tmp0=det.minor(V0,m,d);tmp1=det.minor(V1,m,d)
              mat.det0.max<-cbind(mat.det0.max,(-1)^m*tmp0)
              mat.det1.max<-cbind(mat.det1.max,(-1)^m*tmp1)
              mat.det0.min<-cbind(mat.det0.min,tmp0)
              mat.det1.min<-cbind(mat.det1.min,tmp1)
            }
            ## Sylvester's criterion: max rows of mat.det0.max (resp. 1) must be >0
            ##                        min rows of mat.det0.min (resp. 1) must be >0
            M.max=cbind(mat.det0.max,mat.det1.max)
            M.min=cbind(mat.det0.min,mat.det1.min)
            one=matrix(1,nr=2*d)
            ## which replications are s.t. (-1)^m det(Hm)>0 
            imax=which( (2*d-c((M.max>0)%*%one)) ==0)  
            imin=which( (2*d-c((M.min>0)%*%one)) ==0)  
            adet0=abs(mat.det0.max[,d])
            adet1=abs(mat.det1.max[,d])
            iboth=c(imin,imax)
            pcfs[iboth] <- .5* pGrX0Xre1*adet0[iboth]*adet1[iboth]/rho^2
          })
  switch(which.crit,
         "all"={
           m.pcfs=mean(pcfs)
           sd.pcfs=sd(pcfs)
         },
         "max"={
           m.pcfs=mean(pcfs)
           sd.pcfs=sd(pcfs)
         })
  precision <- qnorm(0.975)*sd.pcfs/sqrt(B)
  out[,i]<-c(m.pcfs,precision)
  }
  close(pb)
  if (do.parallel) stopCluster(cluster) 
  test=sum(apply(!is.na(eig.val.mat),2,sum))>0
  if (test){
    eig.val=c(sum(eig.val.mat[,1],na.rm=TRUE),
              min(eig.val.mat[,2],na.rm=TRUE),
              mean(vecr[!is.na(eig.val.mat[,1])]))
  } else eig.val=c(NA,NA,NA)
  list(out=out,eig.val=eig.val)
  
}


#############################################################################
######### unessential/auxiliary functions
#############################################################################

cst.d=function(B=1e4,d=2){
  Jd=matrix(1,d,d)
  Id=matrix(0,d,d);diag(Id)=1
  Idind=matrix(0,d*(d-1)/2,d*(d-1)/2);diag(Idind)=1
  dd=d+d*(d-1)/2
  M=matrix(0,dd,dd);
  M[1:d,1:d]=2*Id+Jd
  if (d>1) M[(d+1):dd,(d+1):dd]=Idind
  z=rmvnorm(n=B,mean=rep(0,dd),sigma=M)
  require(configural)
  fmax=function(l){
    tmp=vechs2full(l[-(1:d)],diagonal=l[1:d])
    abs(det(tmp))*(all(eigen(tmp)$val<0))
  }
  f=function(l){
    tmp=vechs2full(l[-(1:d)],diagonal=l[1:d])
    abs(det(tmp))
  }
  res.all=(apply(z,1,f))
  res.max=(apply(z,1,fmax))
  res.ratio=res.max/res.all
  m=c(mean(res.all),mean(res.max),mean(res.max/res.all))
  s=c(sd(res.all),sd(res.max),sd(res.max/res.all))
  low=m-s*qnorm(.975)/sqrt(B)
  upp=m+s*qnorm(.975)/sqrt(B)
  cbind(m,low,upp)
}

plot.ratio=function(B,dd=1:10){
  out=NULL
  for (d in dd){
    cat('\r d=',d)
    out=rbind(out,cst.d(B=B,d=d)[3,])  
  }
  cat('\n')
  out=log(out);dd=log(dd)
  df=data.frame(d=dd,m=out[,1],low=out[,2],upp=out[,3])
  ggplot(data=df)+
    geom_line(aes(x=d,y=m),col="blue")+
    geom_point(aes(x=d,y=m),col="blue")+
    geom_ribbon(aes(x=d,ymin=low,ymax=upp),fill="blue",alpha=.3)+
    ylab('Ratio of maxima among critical points')#+
  #scale_x_continuous(name="Dimension",breaks=dd)
  #    ylim(c(0,df$upp[1]*1.1))
}


structure.factor <- function(vecr,pcf,sdpcf,veck=NULL,rho,phi=NULL, d=2,
                             which.cov="Gaussian",which.crit = "max",nu=NULL){
 ## Bartlett spectrum or structure factor for isotropic and stationary pcf
  tmp=rho2phi(rho=rho,phi=phi,d=d,which.cov=which.cov,which.crit=which.crit,nu=nu)
  rho=tmp$rho;phi=tmp$phi
  if (is.null(veck)) veck=vecr
  stepr=mean(diff(vecr))
  I=(pcf-1)*vecr^(d/2)*besselJ(vecr%*%t(veck),d/2-1)
  I=apply(I,2,sum)*stepr
  Ilow=(pcf+sdpcf-1)*vecr^(d/2)*besselJ(vecr%*%t(veck),d/2-1)
  Ilow=apply(Ilow,2,sum)*stepr
  Iupp=(pcf-sdpcf-1)*vecr^(d/2)*besselJ(vecr%*%t(veck),d/2-1)
  Iupp=apply(Iupp,2,sum)*stepr
  spec=function(I){
    1+rho*(2*pi)^(d/2)/veck^(d/2-1)*I
  }
  cbind(spec(I),spec(Ilow),spec(Iupp))
}


gammad=function(d,B=1e4){
  out=NULL
  for (b in 1:B){
    L=matrix(0,d,d)
    diag(L)=rnorm(d,0,sqrt(1/3))
    G=vech2full(rnorm(d*(d+1)/2))
    out=c(out,det(G-L)^2)
  }
  mean(out)
}



## TODO
K <- function(vecr, rho, phi=NULL, n=10000, choice = "max"){
  ## right now only for d=2
  tmp <- g(vecr=vecr,rho=rho,phi=phi,d=d,which.cov=which.cov,
           which.crit =which.crit,B=B,nu=nu)
  Y <- tmp$out
  cumsum(2*pi*vecr*Y[1,]*diff(range(Gr))/length(Gr))
}



###### Not used

diff.c1tilde=function(r,order=0,which.cov="RWM",d=2){
  # Description
  # Gives the derivatives of the covariance function c1 (up to the scale parameter $\phi$)
  # Arguments
  # - r : value at which we want to compute the derivative
  # - order : order of the derivative (order = 0 is the function c1)
  # - which.cov : specify the model used (RWM)
  # - d : dimension of the gaussian field (d=2)
  # Details
  #
  # Value
  # - The value of the order-th derivative of c1 evaluated at r
  switch(which.cov,
         "RWM"={
           if (d==2){
             out=switch(order+1,
                        besselJ(r,0), #order=0
                        -1*besselJ(r,1), #order=1
                        1/2 *( besselJ(r,0)-besselJ(r,2) ),
                        -1/4 *(3* besselJ(r,1)-besselJ(r,3) ),
                        1/8 *(3* besselJ(r,0)-4*besselJ(r,2) + besselJ(r,4))
             )}}
  )
  out
}


### obsolete....used to get a continuous version of c2^(o)(r)...but limit was enough
#"RWM"={
#  out=NULL
  ### Only defined for d=2 (to be done for higher d)
#  for (p in order){
    ## The derivatives c_2^p(r^2) have been checked with Wolfram
    ## 1st: 1/(2phi^2) *J0'(r)/r
    ## 2nd: -1/8phi^4*( J0'(r)/r + J2'(r)/r) 
    ## 3rd: 1/(2^6*phi^6)* ( J0'(r)/r + 4/3*J2'(r)/r+1/3*J4'(r)/r) 
    ## 4th: -1/(2^8*3*phi^8)*(J0'(r)/r+3/2J2'(r)/r+3/5*J4'(r)/r+1/10*J6'(r)/r) 
 #   rr=sqrt(r)/phi
#    out=c(out,switch(p+1,
#                     {besselJ(rr,0)},
#                     {-1/(4*phi^2)*(J0(rr)+J2(rr))},
#                     {1/(2^5*phi^4)*(J0(rr)+4/3*J2(rr)+1/3*J4(rr))},
#                     {-1/(2^7*3*phi^6)*(J0(rr)+3/2*J2(rr)+3/5*J4(rr)+1/10*J6(rr))},
#                     {1/(3*2^(11)*phi^8)*(J0(rr)+8/5*J2(rr)+4/5*J4(rr)+8/35*J6(rr)+1/35*J8(rr))}))}
#},

## Matern
#theta=2^(1-nu)/gamma(nu)
#phitilde=phi/sqrt(2*nu)
#rr=sqrt(r)/phitilde
#cst=theta*(-1)^(order)/(2^(order)*phitilde^(2*order))
#i=which(rr==0);j=which(rr!=0)
#out=matrix(0,length(rr),length(order))
#out[i,]=cst*2^(nu-order-1)*gamma(nu-order)
#out[j,]=cst*rr[j]^(nu-order)*besselK(rr[j],nu-order)
#if (length(rr)==1) out=c(out)


## RWM
#    rr=sqrt(r)/(phi/sqrt(d))
#    delta=d/2-1;theta=2^delta*gamma(delta+1)
#    cst=theta*(-1)^(order)/(2*(phi/sqrt(d))^2)^(order)
#    i=which(rr==0);j=which(rr!=0)
#    out=matrix(0,length(rr),length(order))
#    out[i,]=cst/(2^(delta+order)*gamma(delta+order+1))
#    out[j,]=cst*rr[j]^(-delta-order)*besselJ(rr[j],delta+order)
#    if (length(rr)==1) out=c(out) 
"""
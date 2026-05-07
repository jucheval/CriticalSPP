using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..")) # Activate CriticalSPP project environment
using CriticalSPP
using DrWatson

@quickactivate # Activate test project environment
using RCall
using Test

R"""
diff.c2=function(r,phi,order=0,which.cov="Gaussian",d=2,nu=NULL){
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

SigmaNablaX0Xre1 <- function(r,phi, d=2,which.cov="Gaussian",nu=NULL){
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
"""

phirange = range(0.5, 3.0; step=0.5)
rrange = range(0.1, 3.1; step=0.5)
nurange = range(4.1, 6.1; step=0.5)

@testset verbose = true "covariance_hessians_x0_xr" begin
    @testset "Gaussian covariance" begin
        for phi in phirange, d in 1:4, r in rrange
            r_val = rcopy(R"SigmaNablaX0Xre1($r, $phi, which.cov=\"Gaussian\", d=$d)")
            jl_val = CriticalSPP.covariance_hessians_x0_xr(GaussianCovariance(phi, d), r)
            @test isapprox(r_val, jl_val)
        end
    end

    @testset "Matérn covariance" begin
        for phi in phirange, nu in nurange, d in 1:4, r in rrange
            r_val = rcopy(R"SigmaNablaX0Xre1($r, $phi, which.cov=\"Matern\", d=$d, nu=$nu)")
            jl_val = CriticalSPP.covariance_hessians_x0_xr(MaternCovariance(phi, nu, d), r)
            @test isapprox(r_val, jl_val)
        end
    end

    @testset "RWM covariance" begin
        for phi in phirange, d in 2:4, r in rrange
            # d = 1 is excluded because it correspond to the sin-cosine process
            # for which the distribution is degenerated
            r_val = rcopy(R"SigmaNablaX0Xre1($r, $phi, which.cov=\"RWM\", d=$d)")
            jl_val = CriticalSPP.covariance_hessians_x0_xr(RWMCovariance(phi, d), r)
            @test isapprox(r_val, jl_val)
        end
    end
end;
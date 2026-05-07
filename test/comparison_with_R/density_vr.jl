using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..")) # Activate CriticalSPP project environment
using CriticalSPP
using DrWatson

@quickactivate # Activate test project environment
using RCall
using Test

R"""
require(mvtnorm)

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

densGrX0Xre1 <- function(r,phi, d=2,which.cov="Gaussian",nu=NULL){
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
"""

phirange = range(0.5, 3.0; step=0.5)
rrange = range(0.1, 3.1; step=0.5)
nurange = range(4.1, 6.1; step=0.5)

@testset verbose = true "density_vr" begin
    @testset "Gaussian covariance" begin
        for phi in phirange, d in 1:4, r in rrange
            r_val = rcopy(R"densGrX0Xre1($r, $phi, which.cov=\"Gaussian\", d=$d)")
            jl_val = CriticalSPP.density_vr(GaussianCovariance(phi, d), r)
            @test isapprox(r_val, jl_val)
        end
    end

    @testset "Matérn covariance" begin
        for phi in phirange, nu in nurange, d in 1:4, r in rrange
            r_val = rcopy(R"densGrX0Xre1($r, $phi, which.cov=\"Matern\", d=$d, nu=$nu)")
            jl_val = CriticalSPP.density_vr(MaternCovariance(phi, nu, d), r)
            @test isapprox(r_val, jl_val)
        end
    end

    @testset "RWM covariance" begin
        for phi in phirange, d in 2:4, r in rrange
            # d = 1 is excluded because it correspond to the sin-cosine process
            # for which the distribution is degenerated
            r_val = rcopy(R"densGrX0Xre1($r, $phi, which.cov=\"RWM\", d=$d)")
            jl_val = CriticalSPP.density_vr(RWMCovariance(phi, d), r)
            @test isapprox(r_val, jl_val)
        end
    end
end;
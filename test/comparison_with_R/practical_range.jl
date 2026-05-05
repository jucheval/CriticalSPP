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
"""

phirange = range(0.5, 3.0; step=0.5)
valrange = range(0.01, 0.21; step=0.02)
nurange = range(0.5, 3.0; step=0.5)

@testset verbose = true "practical_range" begin
    @testset "Gaussian covariance" begin
        for phi in phirange, d in 1:4, val in range(0.01, 0.21; step=0.02)
            r_val = rcopy(
                R"practical.range(phi=$phi, which.cov=\"Gaussian\", d=$d, val=$val)"
            )
            jl_val = practical_range(GaussianCovariance(phi, d), val)
            @test isapprox(r_val, jl_val)
        end
    end

    @testset "Matérn covariance" begin
        for phi in phirange, nu in nurange, d in 1:4, val in range(0.01, 0.21; step=0.02)
            r_val = rcopy(
                R"practical.range(phi=$phi, which.cov=\"Matern\", d=$d, val=$val, nu=$nu)"
            )
            jl_val = practical_range(MaternCovariance(phi, nu, d), val)
            @test isapprox(r_val, jl_val; rtol=1e-2) # higher tolerance due to numerical optimization
        end
    end
end;
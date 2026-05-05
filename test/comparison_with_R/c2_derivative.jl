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
"""

phirange = range(0.5, 3.0; step=0.5)
srange = range(0.0, 3.0; step=0.5)
nurange = range(0.5, 3.0; step=0.5)

@testset verbose = true "c2_derivative" begin
    @testset "Gaussian covariance" begin
        for phi in phirange, d in 1:4, s in srange, k in 0:10
            r_val = rcopy(R"diff.c2($s, $phi, order=$k, which.cov=\"Gaussian\", d=$d)")
            jl_val = c2_derivative(GaussianCovariance(phi, d), s, k)
            @test isapprox(r_val, jl_val)
        end
    end

    @testset "Matérn covariance" begin
        for phi in phirange, nu in nurange, d in 1:4, s in srange
            for k in 0:(ceil(Int, nu) - 1)
                r_val = rcopy(
                    R"diff.c2($s, $phi, order=$k, which.cov=\"Matern\", d=$d, nu=$nu)"
                )
                jl_val = c2_derivative(MaternCovariance(phi, nu, d), s, k)
                @test isapprox(r_val, jl_val)
            end
        end
    end
end;
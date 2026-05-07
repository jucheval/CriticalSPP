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

rho2phi=function(rho,phi=NULL,d=2,which.cov="Gaussian",which.crit="max",nu=NULL){
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
"""

phirange = range(0.5, 3.0; step=0.5)
type_range = ["all", "max"]
nurange = range(2.5, 4.0; step=0.5)

@testset verbose = true "intensity" begin
    @testset "Gaussian covariance" begin
        for phi in phirange, d in 1:4, type in type_range
            r_val = rcopy(
                R"rho2phi(NULL, $phi, which.cov=\"Gaussian\", d=$d, which.crit=$type)$rho"
            )
            cov = GaussianCovariance(phi, d)
            jl_type = type == "all" ? ALL_CRITICAL : MAX_CRITICAL
            cpp = CriticalPointProcess(cov, jl_type)
            jl_val = intensity(cpp)
            tol = (type == "max" && d == 4) ? 1e-2 : 1e-10 # higher tolerance for max critical points with d=4 due to MC estimation
            @test isapprox(r_val, jl_val, rtol=tol)
        end
    end

    @testset "Matérn covariance" begin
        for phi in phirange, nu in nurange, d in 1:4, type in type_range
            r_val = rcopy(
                R"rho2phi(NULL, $phi, which.cov=\"Matern\", d=$d, nu=$nu, which.crit=$type)$rho",
            )
            cov = MaternCovariance(phi, nu, d)
            jl_type = type == "all" ? ALL_CRITICAL : MAX_CRITICAL
            cpp = CriticalPointProcess(cov, jl_type)
            jl_val = intensity(cpp)
            tol = (type == "max" && d == 4) ? 1e-2 : 1e-10 # higher tolerance for max critical points with d=4 due to MC estimation
            @test isapprox(r_val, jl_val, rtol=tol)
        end
    end

    @testset "RWM covariance" begin
        for phi in phirange, d in 1:4, type in type_range
            r_val = rcopy(
                R"rho2phi(NULL, $phi, which.cov=\"RWM\", d=$d, which.crit=$type)$rho"
            )
            cov = RWMCovariance(phi, d)
            jl_type = type == "all" ? ALL_CRITICAL : MAX_CRITICAL
            cpp = CriticalPointProcess(cov, jl_type)
            jl_val = intensity(cpp)
            tol = (type == "max" && d == 4) ? 1e-2 : 1e-10 # higher tolerance for max critical points with d=4 due to MC estimation
            @test isapprox(r_val, jl_val, rtol=tol)
        end
    end
end;
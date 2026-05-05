using RCall
using Test

R"""
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

# Gaussian covariance
for phi in range(0.5, 3.0; step=0.5), d in 1:4, val in range(0.01, 0.21; step=0.02)
    r_val = rcopy(R"practical.range(phi=$phi, which.cov=\"Gaussian\", d=$d, val=$val)")
    jl_val = practical_range(GaussianCovariance(phi, d), val)
    @test isapprox(r_val, jl_val)
end

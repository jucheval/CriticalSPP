using RCall
using Test

R"""
th.lambda2p=function(phi,p=1,which.cov="Gaussian",d=2,nu=NULL){
  switch(which.cov,
        "Gaussian"={prod(1:(2*p))/prod(1:p)/2^p/phi^(2*p)},
        "Matern"={prod(1:(2*p))/prod(1:p)/2^p/phi^(2*p) * prod(nu/(nu-1:p)) },
        "RWM"={prod(1:(2*p))/prod(1:p)/2^p/phi^(2*p) * prod(d/(d+2*(0:(p-1))))}
  )
}
"""

# Gaussian covariance
for phi in range(0.5, 3.0; step=0.5), d in 1:4, p in 1:15
    r_val = rcopy(R"th.lambda2p($phi, p=$p, which.cov=\"Gaussian\", d=$d)")
    jl_val = spectral_moment(GaussianCovariance(phi, d), p)
    @test isapprox(r_val, jl_val)
end

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..")) # Activate CriticalSPP project environment
using CriticalSPP
using DrWatson

@quickactivate # Activate test project environment
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

phirange = range(0.5, 3.0; step=0.5)
prange = 1:15
nurange = range(0.5, 3.0; step=0.5)

@testset verbose = true "spectral_moment" begin
    @testset "Gaussian covariance" begin
        for phi in phirange, d in 1:4, p in prange
            r_val = rcopy(R"th.lambda2p($phi, p=$p, which.cov=\"Gaussian\", d=$d)")
            jl_val = spectral_moment(GaussianCovariance(phi, d), p)
            @test isapprox(r_val, jl_val)
        end
    end

    @testset "Matérn covariance" begin
        for phi in phirange, nu in nurange, d in 1:4
            for p in 1:(ceil(Int, nu) - 1)
                r_val = rcopy(
                    R"th.lambda2p($phi, p=$p, which.cov=\"Matern\", d=$d, nu=$nu)"
                )
                jl_val = spectral_moment(MaternCovariance(phi, nu, d), p)
                @test isapprox(r_val, jl_val)
            end
        end
    end

    @testset "RWM covariance" begin
        for phi in phirange, d in 1:4, p in prange
            r_val = rcopy(R"th.lambda2p($phi, p=$p, which.cov=\"RWM\", d=$d)")
            jl_val = spectral_moment(RWMCovariance(phi, d), p)
            @test isapprox(r_val, jl_val)
        end
    end
end;
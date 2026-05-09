phirange = range(0.5, 3.0; step=0.5)
prange = 1:15
nurange = range(0.5, 3.0; step=0.5)

@testset "Gaussian covariance" begin
    for phi in phirange, p in prange
        r_val = rcopy(R"th.lambda2p($phi, p=$p, which.cov=\"Gaussian\")")
        jl_val = spectral_moment(GaussianCovariance(phi, 1), p)
        @test isapprox(r_val, jl_val)
    end
end

@testset "Matérn covariance" begin
    for phi in phirange, nu in nurange
        for p in 1:(ceil(Int, nu) - 1)
            r_val = rcopy(R"th.lambda2p($phi, p=$p, which.cov=\"Matern\", nu=$nu)")
            jl_val = spectral_moment(MaternCovariance(phi, nu, 1), p)
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

phirange = range(0.5, 3.0; step=0.5)
valrange = range(0.01, 0.21; step=0.02)
nurange = range(0.5, 3.0; step=0.5)

@testset "Gaussian covariance" begin
    for phi in phirange, val in range(0.01, 0.21; step=0.02)
        r_val = rcopy(
            R"practical.range(rho=NULL, phi=$phi, which.cov=\"Gaussian\", val=$val)"
        )
        jl_val = practical_range(GaussianCovariance(phi, 1), val)
        @test isapprox(r_val, jl_val)
    end
end

@testset "Matérn covariance" begin
    for phi in phirange, nu in nurange, val in range(0.01, 0.21; step=0.02)
        r_val = rcopy(
            R"practical.range(rho=NULL, phi=$phi, which.cov=\"Matern\", val=$val, nu=$nu)"
        )
        jl_val = practical_range(MaternCovariance(phi, nu, 1), val)
        @test isapprox(r_val, jl_val; rtol=1e-2) # higher tolerance due to numerical optimization
    end
end

@testset "RWM covariance" begin
    for phi in phirange, d in [2], val in range(0.01, 0.21; step=0.02)
        r_val = rcopy(
            R"practical.range(rho=NULL, phi=$phi, which.cov=\"RWM\", d=$d, val=$val)"
        )
        jl_val = practical_range(RWMCovariance(phi, d), val)
        @test isapprox(r_val, jl_val)
    end
end
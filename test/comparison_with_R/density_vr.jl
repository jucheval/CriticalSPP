phirange = range(0.5, 3.0; step=0.5)
rrange = range(0.1, 3.1; step=0.5)
nurange = range(4.1, 6.1; step=0.5)

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

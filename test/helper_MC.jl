phirange = range(0.5, 3.0; step=0.5)
rrange = range(0.1, 3.1; step=0.5)
nurange = range(4.1, 6.1; step=0.5)

@testset "Gaussian covariance" begin
    for phi in phirange, d in 1:4, r in rrange
        Σ = CriticalSPP.covariance_hessians_x0_xr(GaussianCovariance(phi, d), r)
        @test isposdef(Σ)
    end
end

@testset "Matérn covariance" begin
    for phi in phirange, nu in nurange, d in 1:4, r in rrange
        Σ = CriticalSPP.covariance_hessians_x0_xr(MaternCovariance(phi, nu, d), r)
        @test isposdef(Σ)
    end
end

@testset "RWM covariance" begin
    for phi in phirange, d in 2:4, r in rrange
        # d = 1 is excluded because it correspond to the sin-cosine process
        # for which the distribution is degenerated
        r == 0.1 && continue
        Σ = CriticalSPP.covariance_hessians_x0_xr(RWMCovariance(phi, d), r)
        @test isposdef(Σ)
    end
end
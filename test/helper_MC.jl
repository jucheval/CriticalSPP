phirange = range(0.5, 2.0; step=0.5)
rrange = range(0.1, 2.1; step=0.5)
nurange = range(4.1, 6.1; step=0.5)

@testset verbose = true "Covariance matrix is positive definite" begin
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
end

@testset verbose = true "Density at zero is positive" begin
    @testset "Gaussian covariance" begin
        for phi in phirange, d in 1:4, r in rrange
            dens = CriticalSPP.density_vr(GaussianCovariance(phi, d), r)
            @test dens > 0
        end
    end

    @testset "Matérn covariance" begin
        for phi in phirange, nu in nurange, d in 1:4, r in rrange
            dens = CriticalSPP.density_vr(MaternCovariance(phi, nu, d), r)
            @test dens > 0
        end
    end

    @testset "RWM covariance" begin
        for phi in phirange, d in 2:4, r in rrange
            # d = 1 is excluded because it correspond to the sin-cosine process
            # for which the distribution is degenerated
            r == 0.1 && continue
            dens = CriticalSPP.density_vr(RWMCovariance(phi, d), r)
            @test dens > 0
        end
    end
end

include("vec_to_mat.jl")
@testset "Determinant of minor matrices" begin
    for d in 1:4, m in 1:d, _ in 1:10
        v = randn(d + d * (d - 1) ÷ 2)
        det_user = CriticalSPP.det_minor(v, m)
        M = vec_to_mat_test(v)
        det_test = det(M[1:m, 1:m])
        @test isapprox(det_user, det_test)
    end
end
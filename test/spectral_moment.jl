phirange = range(0.5, 3.0; step=0.5)
prange = 1:12
nurange = range(0.5, 3.0; step=0.5)

@testset "Gaussian covariance" begin
    for phi in phirange, p in prange
        numerical_val = spectral_moment(GaussianCovariance(phi, 1), p, false)
        closedform_val = spectral_moment(GaussianCovariance(phi, 1), p, true)
        @test isapprox(numerical_val, closedform_val)
    end
end

@testset "Matérn covariance" begin
    for phi in phirange, nu in nurange
        for p in 0:(ceil(Int, nu) - 1)
            numerical_val = spectral_moment(MaternCovariance(phi, nu, 1), p, false)
            closedform_val = spectral_moment(MaternCovariance(phi, nu, 1), p, true)
            @test isapprox(numerical_val, closedform_val)
        end
    end
end

@testset "RWM covariance" begin
    for phi in phirange, d in 1:4, p in prange
        numerical_val = spectral_moment(RWMCovariance(phi, d), p, false)
        closedform_val = spectral_moment(RWMCovariance(phi, d), p, true)
        @test isapprox(numerical_val, closedform_val)
    end
end
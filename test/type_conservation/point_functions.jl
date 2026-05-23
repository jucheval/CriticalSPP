@testset "intensity/scale_from_intensity" begin
    @testset "Gaussian covariance" begin
        for T in _FLOAT_TYPES, D in 1:4, ct in _CRITICAL_TYPES
            phi = T(1.1)

            cov = GaussianCovariance(phi, D)
            cpp = CriticalPointProcess(cov, ct)

            rho = intensity(cpp)
            phi2 = scale_from_intensity(cpp, rho)

            @test rho isa T
            @test phi2 isa T
        end
    end

    @testset "Matern covariance" begin
        for T in _FLOAT_TYPES, D in 1:4, ct in _CRITICAL_TYPES
            phi = T(1.1)
            nu = T(3.0)

            cov = MaternCovariance(phi, nu, D)
            cpp = CriticalPointProcess(cov, ct)

            rho = intensity(cpp)
            phi2 = scale_from_intensity(cpp, rho)

            @test rho isa T
            @test phi2 isa T
        end
    end

    @testset "RWM covariance" begin
        for T in _FLOAT_TYPES, D in 1:4, ct in _CRITICAL_TYPES
            phi = T(1.1)

            cov = RWMCovariance(phi, D)
            cpp = CriticalPointProcess(cov, ct)

            rho = intensity(cpp)
            phi2 = scale_from_intensity(cpp, rho)

            @test rho isa T
            @test phi2 isa T
        end
    end
end

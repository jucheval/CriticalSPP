@testset "Covariance" begin
    @testset "Gaussian covariance" begin
        for T in _FLOAT_TYPES, D in 1:4
            phi = T(1.5)

            g = GaussianCovariance(phi, D)
            @test g isa GaussianCovariance{D,T}
            @test scale(g) isa T
        end
    end

    @testset "Matern covariance" begin
        for T in _FLOAT_TYPES, D in 1:4
            phi = T(1.5)
            nu = T(3.0)

            m = MaternCovariance(phi, nu, D)
            @test m isa MaternCovariance{D,T}
            @test scale(m) isa T
            @test m.nu isa T
        end
    end

    @testset "RWM covariance" begin
        for T in _FLOAT_TYPES, D in 1:4
            phi = T(1.5)

            rwm = RWMCovariance(phi, D)
            @test rwm isa RWMCovariance{D,T}
            @test scale(rwm) isa T
        end
    end
end

@testset "CriticalPointProcess" begin
    for T in _FLOAT_TYPES, D in 1:4, ct in _CRITICAL_TYPES
        cov = GaussianCovariance(T(1.25), D)
        cpp = CriticalPointProcess(cov, ct)

        @test cpp isa CriticalPointProcess{typeof(cov),typeof(ct)}
        @test cpp.cov isa typeof(cov)
        @test critical_type(cpp) isa typeof(ct)
        @test dimension(cpp) == D
    end
end

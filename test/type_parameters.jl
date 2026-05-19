const _CRITICAL_TYPES = (ALL_CRITICAL, MAX_CRITICAL)

@testset "Covariance constructors" begin
    for T in (Float32, Float64, BigFloat), D in 1:4
        phi = T(1.5)

        g = GaussianCovariance(phi, D)
        @test g isa GaussianCovariance{D,T}
        @test scale(g) isa T

        nu = T(3.5)
        m = MaternCovariance(phi, nu, D)
        @test m isa MaternCovariance{D,T}
        @test scale(m) isa T
        @test m.nu isa T

        rwm = RWMCovariance(phi, D)
        @test rwm isa RWMCovariance{D,T}
        @test scale(rwm) isa T
    end
end

@testset "CriticalPointProcess constructor" begin
    for T in (Float32, Float64, BigFloat), D in 1:4, ct in _CRITICAL_TYPES
        cov = GaussianCovariance(T(1.25), D)
        cpp = CriticalPointProcess(cov, ct)

        @test cpp isa CriticalPointProcess{typeof(cov),typeof(ct)}
        @test cpp.cov isa typeof(cov)
        @test critical_type(cpp) isa typeof(ct)
        @test dimension(cpp) == D
    end
end

@testset "intensity/scale_from_intensity" begin
    for T in (Float32, Float64, BigFloat), D in 1:4, ct in _CRITICAL_TYPES
        cov = GaussianCovariance(T(1.1), D)
        cpp = CriticalPointProcess(cov, ct)

        rho = intensity(cpp)
        phi = scale_from_intensity(cpp, rho)

        @test rho isa T
        @test phi isa T

        cov2 = GaussianCovariance(phi, D)
        cpp2 = CriticalPointProcess(cov2, ct)
        @test intensity(cpp2) isa T
    end
end

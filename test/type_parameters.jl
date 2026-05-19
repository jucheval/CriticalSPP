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

@testset "Covariance functors" begin
    for T in (Float32, Float64), D in 1:4
        phi = T(1.2)
        r1 = T(0.8)
        r2 = T(0.0)

        g = GaussianCovariance(phi, D)
        @test g(r1) isa T
        @test g(r2) isa T

        nu = T(2.5)
        m = MaternCovariance(phi, nu, D)
        @test m(r1) isa T
        @test m(r2) isa T

        rwm = RWMCovariance(phi, D)
        @test rwm(r1) isa T
        @test rwm(r2) isa T
    end

    # SpecialFunctions does not currently support BigFloat for Bessel functions
    for D in 1:4
        phi = BigFloat(1.2)
        r = BigFloat(0.8)
        g = GaussianCovariance(phi, D)
        @test g(r) isa BigFloat
    end
end

@testset "c2_derivative" begin
    for T in (Float32, Float64), D in 1:4
        phi = T(1.1)
        s1 = T(0.5)
        s2 = T(0.0)

        g = GaussianCovariance(phi, D)
        for k in 1:3
            @test c2_derivative(g, s1, k) isa T
            @test c2_derivative(g, s2, k) isa T
        end

        nu = T(3.5)
        m = MaternCovariance(phi, nu, D)
        for k in 1:(Int(floor(nu)) - 1)
            @test c2_derivative(m, s1, k) isa T
            @test c2_derivative(m, s2, k) isa T
        end

        rwm = RWMCovariance(phi, D)
        for k in 1:2
            @test c2_derivative(rwm, s1, k) isa T
            @test c2_derivative(rwm, s2, k) isa T
        end
    end

    for D in 1:4
        phi = BigFloat(1.1)
        s1 = BigFloat(0.5)
        s2 = BigFloat(0.0)
        g = GaussianCovariance(phi, D)
        for k in 1:3
            @test c2_derivative(g, s1, k) isa BigFloat
            @test c2_derivative(g, s2, k) isa BigFloat
        end
    end
end

@testset "spectral_moment" begin
    for T in (Float32, Float64), D in 1:4
        phi = T(1.3)

        g = GaussianCovariance(phi, D)
        for p in 1:4
            @test spectral_moment(g, p, true) isa T
            @test spectral_moment(g, p, false) isa T
        end

        nu = T(3.0)
        m = MaternCovariance(phi, nu, D)
        for p in 1:2
            @test spectral_moment(m, p, true) isa T
            @test spectral_moment(m, p, false) isa T
        end

        rwm = RWMCovariance(phi, D)
        for p in 1:3
            @test spectral_moment(rwm, p, true) isa T
            @test spectral_moment(rwm, p, false) isa T
        end
    end

    for D in 1:4
        phi = BigFloat(1.3)
        g = GaussianCovariance(phi, D)
        for p in 1:4
            @test spectral_moment(g, p, true) isa BigFloat
            @test spectral_moment(g, p, false) isa BigFloat
        end
    end
end

@testset "covariance_hessians_x0_xr" begin
    for T in (Float32, Float64), D in 1:2
        phi = T(1.4)
        r = T(0.5)

        cov = GaussianCovariance(phi, D)

        Σ = CriticalSPP.covariance_hessians_x0_xr(cov, r)
        @test eltype(Σ) == T
    end
end

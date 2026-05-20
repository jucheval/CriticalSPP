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
    for T in (Float32, Float64), S in (Float32, Float64), D in 1:4
        phi = T(1.2)
        r1 = S(0.8)
        r2 = S(0.0)
        R = promote_type(T, S)

        g = GaussianCovariance(phi, D)
        @test g(r1) isa R
        @test g(r2) isa R

        nu = T(2.5)
        m = MaternCovariance(phi, nu, D)
        @test m(r1) isa R
        @test m(r2) isa R

        rwm = RWMCovariance(phi, D)
        @test rwm(r1) isa R
        @test rwm(r2) isa R
    end

    # SpecialFunctions does not currently support BigFloat for Bessel functions
    # TODO: add a try-catch to test the BigFloat case when the support is added in SpecialFunctions
    for D in 1:4
        phi = BigFloat(1.2)
        r = BigFloat(0.8)
        g = GaussianCovariance(phi, D)
        @test g(r) isa BigFloat
    end
end

@testset "c2_derivative" begin
    for T in (Float32, Float64), S in (Float32, Float64), D in 1:4
        phi = T(1.1)
        s1 = S(0.5)
        s2 = S(0.0)
        R = promote_type(T, S)

        g = GaussianCovariance(phi, D)
        for k in 1:3
            @test c2_derivative(g, s1, k) isa R
            @test c2_derivative(g, s2, k) isa R
        end

        nu = T(3.5)
        m = MaternCovariance(phi, nu, D)
        for k in 1:(Int(floor(nu)) - 1)
            @test c2_derivative(m, s1, k) isa R
            @test c2_derivative(m, s2, k) isa R
        end

        rwm = RWMCovariance(phi, D)
        for k in 1:2
            @test c2_derivative(rwm, s1, k) isa R
            @test c2_derivative(rwm, s2, k) isa R
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

@testset "practical_range" begin
    for T in (Float32, Float64), S in (Float32, Float64), D in 1:4
        phi = T(1.2)
        val = S(0.5)
        R = promote_type(T, S)

        g = GaussianCovariance(phi, D)
        @test practical_range(g, val) isa R

        nu = T(2.5)
        m = MaternCovariance(phi, nu, D)
        @test practical_range(m, val) isa R

        # TODO: uncomment when practical_range is implemented for RWMCovariance
        # rwm = RWMCovariance(phi, D)
        # @test practical_range(rwm, val) isa R
    end

    for D in 1:4
        phi = BigFloat(1.2)
        val = BigFloat(0.5)
        g = GaussianCovariance(phi, D)
        @test practical_range(g, val) isa BigFloat
    end
end

@testset "covariance_hessians_x0_xr" begin
    for T in (Float32, Float64, BigFloat), S in (Float32, Float64, BigFloat), D in 1:2
        phi = T(1.4)
        r = S(0.5)
        R = promote_type(T, S)

        cov = GaussianCovariance(phi, D)

        Σ = CriticalSPP.covariance_hessians_x0_xr(cov, r)
        @test eltype(Σ) == R
    end
end

@testset "density_vr" begin
    for T in (Float32, Float64, BigFloat), S in (Float32, Float64, BigFloat), D in 1:2
        phi = T(1.4)
        r = S(0.5)
        R = promote_type(T, S)

        cov = GaussianCovariance(phi, D)

        @test CriticalSPP.density_vr(cov, r) isa R
    end
end

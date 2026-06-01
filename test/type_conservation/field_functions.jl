@testset "functor" begin
    @testset "Gaussian covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.2)
            r1 = S(0.8)
            r2 = S(0.0)
            R = promote_type(T, S)

            cov = GaussianCovariance(phi, D)
            @test cov(r1) isa R
            @test cov(r2) isa R
        end
    end

    @testset "Matern covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.2)
            nu = T(2.5)
            r1 = S(0.8)
            r2 = S(0.0)
            R = promote_type(T, S)

            cov = MaternCovariance(phi, nu, D)
            try
                c1_1 = cov(r1)
                c1_2 = cov(r2)
                @test c1_1 isa R
                @test c1_2 isa R
            catch err
                if R == BigFloat && is_bessel_methoderror(err)
                    @test_skip "BigFloat Bessel support missing in SpecialFunctions"
                else
                    rethrow(err)
                end
            end
        end
    end

    @testset "RWM covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.2)
            r1 = S(0.8)
            r2 = S(0.0)
            R = promote_type(T, S)

            cov = RWMCovariance(phi, D)
            try
                c1_1 = cov(r1)
                c1_2 = cov(r2)
                @test c1_1 isa R
                @test c1_2 isa R
            catch err
                if R == BigFloat && is_bessel_methoderror(err)
                    @test_skip "BigFloat Bessel support missing in SpecialFunctions"
                else
                    rethrow(err)
                end
            end
        end
    end
end

@testset "c2_derivative" begin
    @testset "Gaussian covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.1)
            s1 = S(0.5)
            s2 = S(0.0)
            R = promote_type(T, S)

            cov = GaussianCovariance(phi, D)
            for k in 1:3
                @test c2_derivative(cov, s1, k) isa R
                @test c2_derivative(cov, s2, k) isa R
            end
        end
    end

    @testset "Matern covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.1)
            nu = T(3.5)
            s1 = S(0.5)
            s2 = S(0.0)
            R = promote_type(T, S)

            cov = MaternCovariance(phi, nu, D)
            for k in 1:(Int(floor(nu)) - 1)
                try
                    c2_1 = c2_derivative(cov, s1, k)
                    c2_2 = c2_derivative(cov, s2, k)
                    @test c2_1 isa R
                    @test c2_2 isa R
                catch err
                    if R == BigFloat && is_bessel_methoderror(err)
                        @test_skip "BigFloat Bessel support missing in SpecialFunctions"
                    else
                        rethrow(err)
                    end
                end
            end
        end
    end

    @testset "RWM covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.1)
            s1 = S(0.5)
            s2 = S(0.0)
            R = promote_type(T, S)

            cov = RWMCovariance(phi, D)
            for k in 1:3
                try
                    c2_1 = c2_derivative(cov, s1, k)
                    c2_2 = c2_derivative(cov, s2, k)
                    @test c2_1 isa R
                    @test c2_2 isa R
                catch err
                    if R == BigFloat && is_bessel_methoderror(err)
                        @test_skip "BigFloat Bessel support missing in SpecialFunctions"
                    else
                        rethrow(err)
                    end
                end
            end
        end
    end
end

@testset "spectral_moment" begin
    @testset "Gaussian covariance" begin
        for T in _FLOAT_TYPES, D in 1:4
            phi = T(1.3)

            cov = GaussianCovariance(phi, D)
            for p in 1:4
                @test spectral_moment(cov, p) isa T
            end
        end
    end

    @testset "Matern covariance" begin
        for T in _FLOAT_TYPES, D in 1:4
            phi = T(1.3)
            nu = T(3.0)

            cov = MaternCovariance(phi, nu, D)
            for p in 1:4
                @test spectral_moment(cov, p) isa T
            end
        end
    end

    @testset "RWM covariance" begin
        for T in _FLOAT_TYPES, D in 1:4
            phi = T(1.3)

            cov = RWMCovariance(phi, D)
            for p in 1:4
                @test spectral_moment(cov, p) isa T
            end
        end
    end
end

@testset "practical_range" begin
    @testset "Gaussian covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.2)
            val = S(0.5)
            R = promote_type(T, S)

            cov = GaussianCovariance(phi, D)
            @test practical_range(cov, val) isa R
        end
    end

    @testset "Matern covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.2)
            nu = T(2.5)
            val = S(0.5)
            R = promote_type(T, S)

            cov = MaternCovariance(phi, nu, D)
            try
                pr = practical_range(cov, val)
                @test pr isa R
            catch err
                if R == BigFloat && is_bessel_methoderror(err)
                    @test_skip "BigFloat Bessel support missing in SpecialFunctions"
                else
                    rethrow(err)
                end
            end
        end
    end

    @testset "RWM covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 2:4
            phi = T(1.2)
            val = S(0.5)
            R = promote_type(T, S)

            cov = RWMCovariance(phi, D)
            try
                pr = practical_range(cov, val)
                @test pr isa R
            catch err
                if R == BigFloat && is_bessel_methoderror(err)
                    @test_skip "BigFloat Bessel support missing in SpecialFunctions"
                else
                    rethrow(err)
                end
            end
        end
    end
end

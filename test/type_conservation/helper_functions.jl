@testset "covariance_hessians_x0_xr" begin
    @testset "Gaussian covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.4)
            r = S(0.5)
            R = promote_type(T, S)

            cov = GaussianCovariance(phi, D)

            Σ = CriticalSPP.covariance_hessians_x0_xr(cov, r)
            @test eltype(Σ) == R
        end
    end

    @testset "Matern covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.4)
            nu = T(4.5)
            r = S(0.5)
            R = promote_type(T, S)

            cov = MaternCovariance(phi, nu, D)

            try
                Σ = CriticalSPP.covariance_hessians_x0_xr(cov, r)
                @test eltype(Σ) == R
            catch err
                if R == BigFloat && is_bigfloat_bessel_methoderror(err)
                    @test_skip "BigFloat Bessel support missing in SpecialFunctions"
                else
                    rethrow(err)
                end
            end
        end
    end

    @testset "RWM covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.4)
            r = S(0.5)
            R = promote_type(T, S)

            cov = RWMCovariance(phi, D)

            try
                Σ = CriticalSPP.covariance_hessians_x0_xr(cov, r)
                @test eltype(Σ) == R
            catch err
                if R == BigFloat && is_bigfloat_bessel_methoderror(err)
                    @test_skip "BigFloat Bessel support missing in SpecialFunctions"
                else
                    rethrow(err)
                end
            end
        end
    end
end

@testset "density_vr" begin
    @testset "Gaussian covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.4)
            r = S(0.5)
            R = promote_type(T, S)

            cov = GaussianCovariance(phi, D)

            @test CriticalSPP.density_vr(cov, r) isa R
        end
    end

    @testset "Matern covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.4)
            nu = T(2.5)
            r = S(0.5)
            R = promote_type(T, S)

            cov = MaternCovariance(phi, nu, D)

            try
                dens = CriticalSPP.density_vr(cov, r)
                @test dens isa R
            catch err
                if R == BigFloat && is_bigfloat_bessel_methoderror(err)
                    @test_skip "BigFloat Bessel support missing in SpecialFunctions"
                else
                    rethrow(err)
                end
            end
        end
    end

    @testset "RWM covariance" begin
        for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in 1:4
            phi = T(1.4)
            r = S(0.5)
            R = promote_type(T, S)

            cov = RWMCovariance(phi, D)

            try
                dens = CriticalSPP.density_vr(cov, r)
                @test dens isa R
            catch err
                if R == BigFloat && is_bigfloat_bessel_methoderror(err)
                    @test_skip "BigFloat Bessel support missing in SpecialFunctions"
                else
                    rethrow(err)
                end
            end
        end
    end
end

@testset "det_minor" begin
    for T in _FLOAT_TYPES, D in 1:4
        v = rand(T, D + D * (D - 1) ÷ 2)
        @test CriticalSPP.det_minor(v, 1, D) isa T
    end
end

@testset "argument_of_expectation" begin
    for T in _FLOAT_TYPES, D in 1:4, ct in _CRITICAL_TYPES
        ξ0 = rand(T, D + D * (D - 1) ÷ 2)
        ξr = rand(T, D + D * (D - 1) ÷ 2)

        @test CriticalSPP.argument_of_expectation(ct, ξ0, ξr, D) isa T
    end
end

for cov in [GaussianCovariance(2), MaternCovariance(5.5, 2), RWMCovariance(2)]
    @test cov(0.0) == 1.0
    cov == GaussianCovariance(2) && @test cov(1.0) == exp(-0.5)
    @test scale(cov) == 1.0
    @test dimension(cov) == 2

    for type in (ALL_CRITICAL, MAX_CRITICAL)
        cpp = CriticalPointProcess(cov, type)
        @test critical_type(cpp) == type
        @test CriticalSPP.covariance(cpp) == cov
        @test dimension(cpp) == 2

        cpp2 = CriticalPointProcess(cov, type, 2.0)
        @test critical_type(cpp2) == type
        @test dimension(cpp2) == 2
        @test intensity(cpp2) ≈ 2.0
    end
end

@test MaternCovariance(1.5, Float32(1.0), 2) == MaternCovariance(1.5, 1.0, 2)
@test MaternCovariance(Float32(2.0), 2) == MaternCovariance(1.0, 2.0, 2)

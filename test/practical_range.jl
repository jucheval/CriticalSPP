phirange = range(0.5, 3.0; step=0.5)
nurange = range(0.5, 3.0; step=0.5)
valrange = range(0.1, 0.9; step=0.1)

@testset "Gaussian covariance" begin
    for phi in phirange, val in valrange
        cov = GaussianCovariance(phi, 1)
        pr = practical_range(cov, val)
        @test isapprox(cov(pr), val)
    end
end

@testset "Matérn covariance" begin
    for phi in phirange, nu in nurange, val in valrange
        cov = MaternCovariance(phi, nu, 1)
        pr = practical_range(cov, val)
        @test isapprox(cov(pr), val)
    end
end

@testset "RWM covariance" begin
    for phi in phirange, nu in nurange, val in 0.01:0.1:0.41, d in 1:4
        cov = RWMCovariance(phi, d)
        d == 1 && @test_throws ArgumentError(
            "practical range does not make sense for RWM covariance in dimension 1 because the limsup of the covariance is 1 as r goes to infinity",
        ) practical_range(cov, val)
        pr = practical_range(cov, val)
        for _ in 1:5
            pr += rand() * 0.1 * pr
            @test abs(cov(pr)) < val
        end
    end
end

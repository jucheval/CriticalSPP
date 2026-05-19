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

# @testset "RWM covariance" begin
#     for phi in phirange, d in [2], val in valrange # practical range is only implemented for 2D RWM covariance
#         cov = RWMCovariance(phi, d)
#         pr = practical_range(cov, val)
#         @test isapprox(cov(pr), val)
#     end
# end

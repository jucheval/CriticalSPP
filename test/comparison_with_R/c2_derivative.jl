phirange = range(0.5, 3.0; step=0.5)
srange = range(0.0, 3.0; step=0.5)
nurange = range(0.5, 3.0; step=0.5)

@testset "Gaussian covariance" begin
    for phi in phirange, s in srange, k in 0:10
        r_val = rcopy(R"diff.c2($s, $phi, order=$k, which.cov=\"Gaussian\")")
        jl_val = c2_derivative(GaussianCovariance(phi, 1), s, k)
        @test isapprox(r_val, jl_val)
    end
end

@testset "Matérn covariance" begin
    for phi in phirange, nu in nurange, s in srange
        for k in 0:(ceil(Int, nu) - 1)
            r_val = rcopy(R"diff.c2($s, $phi, order=$k, which.cov=\"Matern\", nu=$nu)")
            jl_val = c2_derivative(MaternCovariance(phi, nu, 1), s, k)
            @test isapprox(r_val, jl_val)
        end
    end
end

@testset "RWM covariance" begin
    for phi in phirange, d in 1:4, s in srange, k in 0:10
        r_val = rcopy(R"diff.c2($s, $phi, order=$k, which.cov=\"RWM\", d=$d)")
        jl_val = c2_derivative(RWMCovariance(phi, d), s, k)
        @test isapprox(r_val, jl_val)
    end
end
phirange = range(0.5, 3.0; step=0.5)
type_range = ["all", "max"]
nurange = range(2.5, 4.0; step=0.5)

@testset "Gaussian covariance" begin
    for phi in phirange, d in 1:4, type in type_range
        r_val = rcopy(
            R"rho2phi(NULL, $phi, which.cov=\"Gaussian\", d=$d, which.crit=$type)$rho"
        )
        cov = GaussianCovariance(phi, d)
        jl_type = type == "all" ? ALL_CRITICAL : MAX_CRITICAL
        cpp = CriticalPointProcess(cov, jl_type)
        jl_val = intensity(cpp)
        tol = (type == "max" && d == 4) ? 1e-2 : 1e-10 # higher tolerance for max critical points with d=4 due to MC estimation
        @test isapprox(r_val, jl_val, rtol=tol)
    end
end

@testset "Matérn covariance" begin
    for phi in phirange, nu in nurange, d in 1:4, type in type_range
        r_val = rcopy(
            R"rho2phi(NULL, $phi, which.cov=\"Matern\", d=$d, nu=$nu, which.crit=$type)$rho"
        )
        cov = MaternCovariance(phi, nu, d)
        jl_type = type == "all" ? ALL_CRITICAL : MAX_CRITICAL
        cpp = CriticalPointProcess(cov, jl_type)
        jl_val = intensity(cpp)
        tol = (type == "max" && d == 4) ? 1e-2 : 1e-10 # higher tolerance for max critical points with d=4 due to MC estimation
        @test isapprox(r_val, jl_val, rtol=tol)
    end
end

@testset "RWM covariance" begin
    for phi in phirange, d in 1:4, type in type_range
        r_val = rcopy(R"rho2phi(NULL, $phi, which.cov=\"RWM\", d=$d, which.crit=$type)$rho")
        cov = RWMCovariance(phi, d)
        jl_type = type == "all" ? ALL_CRITICAL : MAX_CRITICAL
        cpp = CriticalPointProcess(cov, jl_type)
        jl_val = intensity(cpp)
        tol = (type == "max" && d == 4) ? 1e-2 : 1e-10 # higher tolerance for max critical points with d=4 due to MC estimation
        @test isapprox(r_val, jl_val, rtol=tol)
    end
end

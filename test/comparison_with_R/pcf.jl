phirange = range(0.5, 3.0; step=0.5)
type_range = ["all", "max"]
nurange = range(4.5, 6.5; step=0.5)
n_MC = 5
rrange = 0.5:0.5:3.0

nlag = length(rrange)
pcf = Vector{Float64}(undef, nlag)
stderr = Vector{Float64}(undef, nlag)

@testset "Gaussian covariance" begin
    for phi in phirange, d in 1:4, type in type_range
        for r in rrange
            rs = [r]
            dd = d + d * (d - 1) ÷ 2
            N01 = randn(2dd, n_MC)

            cov = GaussianCovariance(phi, d)
            jl_type = type == "all" ? ALL_CRITICAL : MAX_CRITICAL
            cpp = CriticalPointProcess(cov, jl_type)
            pcf, stderr, ξ = CriticalSPP._pair_correlation_single(cpp, rs[1], N01)
            jl_val = [pcf'; stderr']

            V = ξ'
            r_val = rcopy(R"g.test($rs,rho=NULL,phi=$phi, d=$d,which.cov=\"Gaussian\",
            which.crit=$type,B=$n_MC,V=$V)$out")
            tol = (type == "max" && d == 4) ? 1e-2 : 1e-5 # higher tolerance for max critical points with d=4 due to MC estimation
            @test isapprox(r_val, jl_val, rtol=tol)
        end
    end
end

@testset "Matérn covariance" begin
    for phi in phirange, nu in nurange, d in 1:4, type in type_range
        for r in rrange
            rs = [r]
            dd = d + d * (d - 1) ÷ 2
            N01 = randn(2dd, n_MC)

            cov = MaternCovariance(phi, nu, d)
            jl_type = type == "all" ? ALL_CRITICAL : MAX_CRITICAL
            cpp = CriticalPointProcess(cov, jl_type)
            pcf, stderr, ξ = CriticalSPP._pair_correlation_single(cpp, rs[1], N01)
            jl_val = [pcf'; stderr']

            V = ξ'
            r_val = rcopy(R"g.test($rs,rho=NULL,phi=$phi, d=$d,which.cov=\"Matern\",
            which.crit=$type,nu=$nu,B=$n_MC,V=$V)$out")
            tol = (type == "max" && d == 4) ? 1e-2 : 1e-5 # higher tolerance for max critical points with d=4 due to MC estimation
            @test isapprox(r_val, jl_val, rtol=tol)
        end
    end
end

@testset "RWM covariance" begin
    for phi in phirange, d in 1:4, type in type_range
        for r in rrange
            rs = [r]
            dd = d + d * (d - 1) ÷ 2
            N01 = randn(2dd, n_MC)

            cov = RWMCovariance(phi, d)
            jl_type = type == "all" ? ALL_CRITICAL : MAX_CRITICAL
            cpp = CriticalPointProcess(cov, jl_type)
            pcf, stderr, ξ = CriticalSPP._pair_correlation_single(cpp, rs[1], N01)
            jl_val = [pcf'; stderr']

            V = ξ'
            r_val = rcopy(R"g.test($rs,rho=NULL,phi=$phi, d=$d,which.cov=\"RWM\",
            which.crit=$type,B=$n_MC,V=$V)$out")
            tol = (type == "max" && d == 4) ? 1e-2 : 1e-5 # higher tolerance for max critical points with d=4 due to MC estimation
            @test isapprox(r_val, jl_val, rtol=tol)
        end
    end
end

phirange = range(0.5, 3.0; step=0.5)
rhorange = range(1.0, 100.0; step=5.0)
nurange = range(2.5, 4.0; step=0.5)
valrange = range(0.1, 0.9; step=0.1)
typerange = [ALL_CRITICAL, MAX_CRITICAL]

@testset "Gaussian covariance" begin
    for phi in phirange, rho in rhorange, type in typerange
        cov = GaussianCovariance(1.0, 1)
        cpp = CriticalPointProcess(cov, type)
        phi = scale_from_intensity(cpp, rho)
        cov = GaussianCovariance(phi, 1)
        cpp = CriticalPointProcess(cov, type)
        @test isapprox(intensity(cpp), rho)
    end
end

@testset "Matérn covariance" begin
    for phi in phirange, nu in nurange, rho in rhorange, type in typerange
        cov = MaternCovariance(1.0, nu, 1)
        cpp = CriticalPointProcess(cov, type)
        phi = scale_from_intensity(cpp, rho)
        cov = MaternCovariance(phi, nu, 1)
        cpp = CriticalPointProcess(cov, type)
        @test isapprox(intensity(cpp), rho)
    end
end

@testset "RWM covariance" begin
    for phi in phirange, d in 1:4, rho in rhorange, type in typerange
        cov = RWMCovariance(1.0, d)
        cpp = CriticalPointProcess(cov, type)
        phi = scale_from_intensity(cpp, rho)
        cov = RWMCovariance(phi, d)
        cpp = CriticalPointProcess(cov, type)
        @test isapprox(intensity(cpp), rho)
    end
end

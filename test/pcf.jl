@testset "pcf_single" begin
    phi = 1.0
    nu = 5.0
    d = 2
    r = 2.0

    cov = GaussianCovariance(phi, d)
    cpp = CriticalPointProcess(cov, MAX_CRITICAL)
    N01 = reshape([1.05, -1.34, 0.3, -0.4, -2.19, 0.36], 6, 1)
    pcf, stderr = CriticalSPP._pair_correlation_single(cpp, r, N01)
    @test pcf == 2.4085090393054207

    cov = MaternCovariance(phi, nu, d)
    cpp = CriticalPointProcess(cov, ALL_CRITICAL)
    N01 = reshape([1.69, -0.89, -0.07, 0.66, 0.07, 0.2], 6, 1)
    pcf, stderr = CriticalSPP._pair_correlation_single(cpp, r, N01)
    @test pcf == 0.20518007023706936

    cov = RWMCovariance(phi, d)
    cpp = CriticalPointProcess(cov, MAX_CRITICAL)
    N01 = reshape([-0.05, 0.55, -0.53, -0.27, -0.09, 2.09], 6, 1)
    pcf, stderr = CriticalSPP._pair_correlation_single(cpp, r, N01)
    @test pcf == 1.792705212503323
end

@testset "pcf_serial" begin
    phi = 2.0
    nu = 6.0
    d = 3
    r = [1.0, 2.0]

    cov = GaussianCovariance(phi, d)
    cpp = CriticalPointProcess(cov, MAX_CRITICAL)
    output = pair_correlation_function(MersenneTwister(1), cpp, r; parallel=:serial)
    @test output == ((
        rs=[1.0, 2.0],
        pcf=[0.03639499323379063, 0.16750513503909623],
        stderr=[1.3601721582588697, 5.9095070705187975],
        n_MC=100000,
    ))

    cov = MaternCovariance(phi, nu, d)
    cpp = CriticalPointProcess(cov, ALL_CRITICAL)
    N01 = reshape([1.69, -0.89, -0.07, 0.66, 0.07, 0.2], 6, 1)
    output = pair_correlation_function(MersenneTwister(1), cpp, r; parallel=:serial)
    @test output == (
        rs=[1.0, 2.0],
        pcf=[1.186940633944181, 0.955240250562449],
        stderr=[5.193355701163643, 3.4092897463239487],
        n_MC=100000,
    )

    cov = RWMCovariance(phi, d)
    cpp = CriticalPointProcess(cov, MAX_CRITICAL)
    N01 = reshape([-0.19, 1.29, -1.59, -0.22, -0.29, -1.22], 6, 1)
    output = pair_correlation_function(MersenneTwister(1), cpp, r; parallel=:serial)
    @test output == (
        rs=[1.0, 2.0],
        pcf=[0.0043908073876960255, 0.018212279321504647],
        stderr=[0.20744815003457243, 0.7006183367647199],
        n_MC=100000,
    )
end

@testset "pcf_threads" begin
    phi = 2.0
    nu = 6.0
    d = 3
    r = [1.0, 2.0]

    cov = GaussianCovariance(phi, d)
    cpp = CriticalPointProcess(cov, MAX_CRITICAL)
    output = pair_correlation_function(
        MersenneTwister(1), cpp, r; parallel=:threads, n_MC=50000, show_progress=false
    )
    # I suspect that eigen decomposition is not stable when using threads
    # That is why I used isapprox here
    @test isapprox(output.pcf, [0.03090374446194397, 0.18199479091006585], rtol=1e-2)
    @test isapprox(output.stderr, [0.9619811259882077, 5.966405406771248], rtol=1e-2)

    cov = MaternCovariance(phi, nu, d)
    cpp = CriticalPointProcess(cov, ALL_CRITICAL)
    N01 = reshape([1.69, -0.89, -0.07, 0.66, 0.07, 0.2], 6, 1)
    output = pair_correlation_function(
        MersenneTwister(1), cpp, r; parallel=:threads, n_MC=50000, show_progress=false
    )
    @test isapprox(output.pcf, [1.2112271927892637, 0.9522771225814033], rtol=1e-2)
    @test isapprox(output.stderr, [5.537312015544417, 3.525118765066465], rtol=1e-2)

    cov = RWMCovariance(phi, d)
    cpp = CriticalPointProcess(cov, MAX_CRITICAL)
    N01 = reshape([-0.19, 1.29, -1.59, -0.22, -0.29, -1.22], 6, 1)
    output = pair_correlation_function(
        MersenneTwister(1), cpp, r; parallel=:threads, n_MC=50000, show_progress=false
    )
    @test isapprox(output.pcf, [0.004544825814447406, 0.01913452648289884], rtol=1e-2)
    @test isapprox(output.stderr, [0.23776984264483733, 0.83360883610212], rtol=1e-2)
end

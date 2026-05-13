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
    @test pcf == 0.21688480435720767

    cov = RWMCovariance(phi, d)
    cpp = CriticalPointProcess(cov, MAX_CRITICAL)
    N01 = reshape([-0.19, 1.29, -1.59, -0.22, -0.29, -1.22], 6, 1)
    pcf, stderr = CriticalSPP._pair_correlation_single(cpp, r, N01)
    @test pcf == 0.0973558013558179
end

@testset "pcf_serial" begin
    phi = 2.0
    nu = 6.0
    d = 3
    r = [1.0, 2.0]

    cov = GaussianCovariance(phi, d)
    cpp = CriticalPointProcess(cov, MAX_CRITICAL)
    output = pair_correlation_function(MersenneTwister(1), cpp, r; parallel=:serial)
    @test output == (
        rs=[1.0, 2.0],
        pcf=(0.03639499323379063, 1.3601721582588697),
        stderr=(0.16750513503909623, 5.9095070705187975),
        n_MC=100000,
    )

    cov = MaternCovariance(phi, nu, d)
    cpp = CriticalPointProcess(cov, ALL_CRITICAL)
    N01 = reshape([1.69, -0.89, -0.07, 0.66, 0.07, 0.2], 6, 1)
    output = pair_correlation_function(MersenneTwister(1), cpp, r; parallel=:serial)
    @test output == (
        rs=[1.0, 2.0],
        pcf=(1.185825438681694, 5.18352815817258),
        stderr=(0.9491433166161567, 3.1353830626716523),
        n_MC=100000,
    )

    cov = RWMCovariance(phi, d)
    cpp = CriticalPointProcess(cov, MAX_CRITICAL)
    N01 = reshape([-0.19, 1.29, -1.59, -0.22, -0.29, -1.22], 6, 1)
    output = pair_correlation_function(MersenneTwister(1), cpp, r; parallel=:serial)
    @test output == (
        rs=[1.0, 2.0],
        pcf=(0.004403745642155812, 0.21541086193483894),
        stderr=(0.018468257258804625, 0.7090098848892278),
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
        MersenneTwister(1), cpp, r; parallel=:threads, n_MC=50000, progressbar=false
    )
    @test output == (
        rs=[1.0, 2.0],
        pcf=[0.03090374446194397, 0.18199479091006585],
        stderr=[0.9619811259882077, 5.966405406771248],
        n_MC=50000,
    )

    cov = MaternCovariance(phi, nu, d)
    cpp = CriticalPointProcess(cov, ALL_CRITICAL)
    N01 = reshape([1.69, -0.89, -0.07, 0.66, 0.07, 0.2], 6, 1)
    output = pair_correlation_function(
        MersenneTwister(1), cpp, r; parallel=:threads, n_MC=50000, progressbar=false
    )
    @test output == (
        rs=[1.0, 2.0],
        pcf=[1.2112271927892637, 0.9522686923515971],
        stderr=[5.50894290323668, 3.107075807749902],
        n_MC=50000,
    )

    cov = RWMCovariance(phi, d)
    cpp = CriticalPointProcess(cov, MAX_CRITICAL)
    N01 = reshape([-0.19, 1.29, -1.59, -0.22, -0.29, -1.22], 6, 1)
    output = pair_correlation_function(
        MersenneTwister(1), cpp, r; parallel=:threads, n_MC=50000, progressbar=false
    )
    @test output == (
        rs=[1.0, 2.0],
        pcf=[0.004667458716599647, 0.019432243383040823],
        stderr=[0.2548776852703169, 0.8368493788235667],
        n_MC=50000,
    )
end

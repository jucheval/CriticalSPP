using BenchmarkTools, CriticalSPP

const SUITE = BenchmarkGroup()

# Create hierarchy of benchmarks:
SUITE["pcf"] = BenchmarkGroup()

cpps = [
    CriticalPointProcess(GaussianCovariance(1.0, 4), MAX_CRITICAL),
    CriticalPointProcess(MaternCovariance(1.0, 5.5, 3), ALL_CRITICAL),
    CriticalPointProcess(RWMCovariance(1.0, 2), ALL_CRITICAL),
]
names = ["Gaussian/d=4/MAX", "Matern/d=3/ALL", "RWM/d=2/ALL"]

for r in [0.1, 2.0], (cpp, name) in zip(cpps, names)
    d = dimension(cpp)
    dd = 2 * (d + d * (d - 1) ÷ 2)
    SUITE["pcf"][name]["r=$r"] = @benchmarkable(
        CriticalSPP._pair_correlation_single($cpp, $r, N01),
        evals = 10,
        samples = 1000,
        setup = (N01 = randn(Float64, $dd, Int(1e3)))
    )
end

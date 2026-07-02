using BenchmarkTools
using CriticalSPP

using Logging
global_logger(NullLogger()) # used to remove the logging info messages

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
    rs = [r]
    SUITE["pcf"][name]["r=$r"] = @benchmarkable(
        pair_correlation_function($cpp, $rs; n_MC=10000, show_progress=false),
        evals = 1,
        samples = 1000
    )
end

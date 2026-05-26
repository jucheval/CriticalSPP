cov = GaussianCovariance(Float32(1.0), 2)
type = MAX_CRITICAL
cpp = CriticalPointProcess(cov, type)
@test CriticalSPP.innertype(cov) == Float32
@test CriticalSPP.innertype(cpp) == Float32
cov2 = GaussianCovariance(1.0, 2)
@test CriticalSPP.convert_innertype(Float64, cov) == cov2
@test CriticalSPP.convert_innertype(Float64, cpp) == CriticalPointProcess(cov2, type)

cov = MaternCovariance(1.0, 5.5, 2)
type = ALL_CRITICAL
cpp = CriticalPointProcess(cov, type)
@test CriticalSPP.innertype(cov) == Float64
@test CriticalSPP.innertype(cpp) == Float64
cov2 = MaternCovariance(Float32(1.0), Float32(5.5), 2)
@test CriticalSPP.convert_innertype(Float32, cov) == cov2
@test CriticalSPP.convert_innertype(Float32, cpp) == CriticalPointProcess(cov2, type)

cov = RWMCovariance(2)
type = MAX_CRITICAL
cpp = CriticalPointProcess(cov, type)
@test CriticalSPP.innertype(cov) == Float64
@test CriticalSPP.innertype(cpp) == Float64
cov2 = RWMCovariance(Float32(1.0), 2)
@test CriticalSPP.convert_innertype(Float32, cov) == cov2
@test CriticalSPP.convert_innertype(Float32, cpp) == CriticalPointProcess(cov2, type)

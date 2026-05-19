module CriticalSPP

using Bessels: gamma, besselk, besselj
using Roots: find_zero
using LinearAlgebra: Symmetric, Diagonal, diagind, diagm, eigen
using Distributions: MvNormal, pdf
using Random: AbstractRNG, default_rng, randn
using ProgressMeter: Progress, next!
using StatsBase: mean

# Covariances
export CovarianceSPP
export GaussianCovariance, MaternCovariance, RWMCovariance
export dimension, scale, practical_range, c2_derivative, spectral_moment

# Critical point processes
export AbstractCriticalType
export MaxCritical, AllCritical, MAX_CRITICAL, ALL_CRITICAL
export CriticalPointProcess
export critical_type, intensity, scale_from_intensity, pair_correlation_function

include("field.jl")
include("point.jl")

end

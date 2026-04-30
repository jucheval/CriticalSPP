"""
    CovarianceSPP{D,T<:Real}

Common interface for all covariance models used to get a critical spatial point process.
"""
abstract type CovarianceSPP{D,T<:Real} end

# Parameters
"""
    dimension(cov)

Return the dimension `D` of the covariance model `cov`.
"""
dimension(::CovarianceSPP{D}) where {D} = D

"""
    scale(cov)

Return the scale parameter `phi` of the covariance model `cov`.
"""
scale(cov::CovarianceSPP) = cov.phi

# Helper function to check that the dimension of the covariance model is less than 4
@inline function _check_dimension(::Val{D}) where {D}
    (1 <= D <= 4) || throw(ArgumentError("dimension D must satisfy 1 <= D <= 4"))
    return D
end

# Quantities derived from the covariance function
"""
    (cov::CovarianceSPP)(r)

Return the value of the covariance `cov` at lag `r`, i.e. ``c_1(r)``.
"""
function (cov::CovarianceSPP) end

"""
    practical_range(cov, val)

Return the practical range of the covariance model `cov` corresponding to the value `val`. More precisely, the practical range is defined as the distance `r` such that cov(`r`) = `val`.
"""
function practical_range end

"""
    c2_derivative(cov, s, k)

Return the `k`-th derivative of the auxiliary covariance function c₂ evaluated at `s = r^2`. More precisely, ``c_2(s) = c_1(\\sqrt{s})`` and `c2_derivative(cov, s, k)` returns ``c_2^{(k)}(s)``.
"""
function c2_derivative end

"""
    spectral_moment(cov, p)

Return the `2p`-th spectral moment of the covariance model `cov`.
"""
function spectral_moment end

# Helper functions for Monte Carlo estimation
"""
    covariance_gradient_x0_xr(cov, r)

Return the covariance matrix of the upper diagonal and diagonal of ∇X(0) and ∇X(`r`*e_1), where X is a Gaussian field with covariance `cov`.
"""
function covariance_gradient_x0_xr end

"""
    density_vr(cov, r)

Return the density at 0 of V(`r`) = (∇X(0), ∇X(`r`*e_1)), where X is a Gaussian field with covariance `cov`.
"""
function density_vr end

#-----------------
# IMPLEMENTATIONS
#-----------------
include("field/gaussian.jl")
include("field/matern.jl")
include("field/rwm.jl")
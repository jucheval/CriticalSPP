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
    (1 <= D <= 4) || throw(DomainError(D, "dimension D must satisfy 1 <= D <= 4"))
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

Return the practical range of the covariance model `cov` corresponding to the value `val`. More precisely, when `cov` is a decreasing covariance function, the practical range is defined as the distance `r₀` such that `cov(r₀) = val`. Otherwise, the practical range is such that `abs(cov(r)) <= val == true` for all `r >= r₀`.
"""
function practical_range end

"""
    c2_derivative(cov, s, k)

Return the `k`-th derivative of the auxiliary covariance function c₂ evaluated at `s = r^2`. More precisely, ``c_2(s) = c_1(\\sqrt{s})`` and `c2_derivative(cov, s, k)` returns ``c_2^{(k)}(s)``.
"""
function c2_derivative end

"""
    spectral_moment(cov, p, [closedform])

Return the `2p`-th spectral moment of the covariance model `cov`. If `closedform=true` (the default), it uses a closed-form expression. Otherwise, it is computed using `c2_derivative(cov, 0, p)`.
"""
# FIXME: numerical computation goes out of bounds for large p.
# To see it, fix prange = 1:15 in test/spectral_moment.jl.
function spectral_moment(cov::CovarianceSPP, p::Integer, closedform::Bool)
    if closedform
        return spectral_moment(cov, p)
    else
        c2p0 = c2_derivative(cov, 0, p)
        return prod((p + 1):(2 * p)) * (-1)^p * c2p0
    end
end

#-----------------
# IMPLEMENTATIONS
#-----------------
include("field/gaussian.jl")
include("field/matern.jl")
include("field/rwm.jl")

#-----------------
# Helper functions for Monte Carlo estimation
#-----------------
include("field/helper_MC.jl")

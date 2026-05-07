"""
    GaussianCovariance(phi, d)

Gaussian covariance with scale parameter `phi` used for a field of dimension `d`. It is defined as

```math
c(r) = \\exp(- r^2 / (2\\phi^2))
```
"""
struct GaussianCovariance{D,T<:Real} <: CovarianceSPP{D,T}
    phi::T
end

# Constructor
function GaussianCovariance(phi::T, d::Integer) where {T<:Real}
    _check_dimension(Val(d))
    phi > zero(T) || throw(DomainError(phi, "the scale parameter must be positive"))
    return GaussianCovariance{d,T}(phi)
end

# Functor
function (cov::GaussianCovariance)(r)
    return exp(-r^2 / (2 * scale(cov)^2))
end

# Practical range
function practical_range(cov::GaussianCovariance, val)
    (0 < val < 1) || throw(DomainError(val, "the value must satisfy 0 < val < 1"))
    return sqrt(-2 * scale(cov)^2 * log(val))
end

# c₂ derivative
function c2_derivative(cov::GaussianCovariance, s, k::Integer)
    tmp = 2 * scale(cov)^2
    return (-1 / tmp)^k * exp(-s / tmp)
end

# Spectral moment
function spectral_moment(cov::GaussianCovariance, p::Integer)
    return _spectral_moment_gaussian(p, scale(cov))
end

function _spectral_moment_gaussian(p::Integer, ϕ)
    return prod((1:2:(2 * p))) / ϕ^(2 * p)
end

# Internals
_string(::GaussianCovariance) = "gaussian"
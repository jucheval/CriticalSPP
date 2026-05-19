"""
    GaussianCovariance([phi], d)

Gaussian covariance with scale parameter `phi` used for a field of dimension `d`. It is defined as

```math
c(r) = \\exp(- r^2 / (2\\phi^2))
```

The scale parameter `phi` is optional and defaults to 1.0.
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

GaussianCovariance(d::Integer) = GaussianCovariance(1.0, d)

# Functor
function (cov::GaussianCovariance)(r)
    ϕ = scale(cov)

    return exp(-r^2 / (2ϕ^2))
end

# Practical range
function practical_range(cov::GaussianCovariance, val)
    (0 < val < 1) || throw(DomainError(val, "the value must satisfy 0 < val < 1"))
    ϕ = scale(cov)

    return sqrt(-2 * ϕ^2 * log(val))
end

# c₂ derivative
function c2_derivative(cov::GaussianCovariance, s, k::Integer)
    ϕ = scale(cov)
    tmp = 2 * ϕ^2

    return (-1 / tmp)^k * exp(-s / tmp)
end

# Spectral moment
function spectral_moment(cov::GaussianCovariance, p::Integer)
    return prod((1:2:(2 * p))) / scale(cov)^(2 * p)
end

# Internals
function constant_λ₄_over_3λ₂(::GaussianCovariance)
    return 1
end

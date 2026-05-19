"""
    RWMCovariance([phi], d)

Random Wave Model covariance with scale parameter `phi` used for a field of dimension `d`. It is defined as

```math
c(r) = \\Gamma(d/2) (\\frac{r}{2\\phi'})^{-(d/2-1)} J_{d/2-1}(r / \\phi'), \\quad \\phi' = \\phi/\\sqrt{d},
```
where ``J_\\nu`` is the Bessel function of the first kind.

The scale parameter `phi` is optional and defaults to 1.0.
"""
struct RWMCovariance{D,T<:Real} <: CovarianceSPP{D,T}
    phi::T
end

# Constructor
function RWMCovariance(phi::T, d::Integer) where {T<:Real}
    _check_dimension(Val(d))
    phi > zero(T) || throw(DomainError(phi, "phi must be positive"))
    return RWMCovariance{d,T}(phi)
end

RWMCovariance(d::Integer) = RWMCovariance(1.0, d)

# Functor
function (cov::RWMCovariance)(r)
    ϕ = scale(cov)
    D = dimension(cov)

    # guard r close to 0 to avoid numerical explosion
    r ≈ 0 && return 1.0

    δ = r * √D / ϕ
    J = besselj(D / 2 - 1, δ)
    Γ = gamma(D / 2)

    return Γ * (δ / 2)^(-(D / 2 - 1)) * J
end

# Practical range (only for d=2, using asymptotic behaviour of J0)
# FIXME: currently it does not satisfy cov(practical_range(cov, val)) ≈ val
# implement numerical root finding?
function practical_range(cov::RWMCovariance, val)
    ϕ = scale(cov)
    D = dimension(cov)

    (0 < val < 1) || throw(DomainError(val, "the value must satisfy 0 < val < 1"))
    D == 2 || throw(DomainError(D, "The dimension D of RWMCovariance must be 2"))

    return 2ϕ * (val * sqrt(pi))^(-2)
    # expression from the R code
    # delta = D / 2 - 1
    # 2 * ϕ * (gamma(delta + 1) / (sqrt(pi) * val))^(1 / (delta + 1 / 2))
end

# c₂ derivative
function c2_derivative(cov::RWMCovariance, s, k::Integer)
    D = dimension(cov)
    ϕ = scale(cov)

    cst = (-1)^k / 2^k / ϕ^(2k) * (D / 2)^k * gamma(D / 2)

    # guard s close to 0 to avoid numerical explosion
    s ≈ 0 && return cst / gamma(D / 2 + k)

    α = D / 2 - 1 + k
    δ = √s * √D / ϕ

    return cst * (δ / 2)^(-α) * besselj(α, δ)
end

# Spectral moment
function spectral_moment(cov::RWMCovariance, p::Integer)
    D = dimension(cov)
    ϕ = scale(cov)
    return prod((1:2:(2p))) / ϕ^(2p) * D^p / prod(D .+ 2 * (0:(p - 1)))
end

# Internals
function constant_λ₄_over_3λ₂(cov::RWMCovariance)
    D = dimension(cov)
    return D//(D + 2)
end

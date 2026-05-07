"""
    MaternCovariance(phi, nu, d)

Matérn covariance with scale parameter `phi` and smoothness parameter `nu` used for a field of dimension `d`. It is defined as

```math
c(r) = \\frac{2^{1-\\nu}}{\\Gamma(\\nu)} \\left(\\frac{r \\sqrt{2\\nu}}{\\phi}\\right)^\\nu K_\\nu\\left(\\frac{r \\sqrt{2\\nu}}{\\phi}\\right),
```
where ``K_\\nu`` is the modified Bessel function of the second kind.
"""
struct MaternCovariance{D,T<:Real} <: CovarianceSPP{D,T}
    phi::T
    nu::T
end

# Constructor
function MaternCovariance(phi::T, nu::T, d::Integer) where {T<:Real}
    _check_dimension(Val(d))
    phi > zero(T) || throw(DomainError(phi, "phi must be positive"))
    nu > zero(T) || throw(DomainError(nu, "nu must be positive"))
    return MaternCovariance{d,T}(phi, nu)
end

# Functor
### adapted from MaternVariogram in GeoStats.jl
function (cov::MaternCovariance)(r)
    ν = cov.nu
    ϕ = scale(cov)

    # guard r close to 0 to avoid numerical explosion
    r ≈ 0 && return 1.0

    δ = r * √(2ν) / ϕ
    Β = besselk(ν, δ)
    Γ = gamma(ν)

    return 2^(1 - ν) / Γ * δ^ν * Β
end

# Practical range
function practical_range(cov::MaternCovariance, val)
    (0 < val < 1) || throw(DomainError(val, "the value must satisfy 0 < val < 1"))

    initial_guess = practical_range(GaussianCovariance(scale(cov), dimension(cov)), val)
    f(x) = cov(x) - val

    return find_zero(f, initial_guess)
end

# c₂ derivative
function c2_derivative(cov::MaternCovariance, s, k::Integer)
    ν = cov.nu
    ϕ = scale(cov)
    k >= ν &&
        throw(ArgumentError("the order k must be less than the smoothness parameter nu"))

    cst = (-1)^k * 2 / 2^k / ϕ^(2k) * ν^k / gamma(ν)

    # guard s close to 0 to avoid numerical explosion
    s ≈ 0 && return cst * gamma(ν - k) / 2

    δ = √(s) * √(2ν) / ϕ
    return cst * (δ / 2)^(ν - k) * besselk(ν - k, δ)
end

# Spectral moment
function spectral_moment(cov::MaternCovariance, p::Integer)
    ν = cov.nu
    ϕ = scale(cov)
    p >= ν && throw(
        ArgumentError("the half-order p must be less than the smoothness parameter nu")
    )

    return prod((1:2:(2p))) / ϕ^(2p) * ν^p / prod(ν .- (1:p))
end
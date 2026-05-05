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
    r ≈ 0 && return gamma(ν) / 2
    # r += eps(typeof(r))

    δ = r * √(2ν) / ϕ
    Β = besselk(ν, δ)
    Γ = gamma(ν)

    return 2^(1 - ν) / Γ * δ^ν * Β
end

# Practical range
function practical_range(cov::MaternCovariance, val)
    initial_guess = practical_range(GaussianCovariance(scale(cov), dimension(cov)), val)
    f(x) = (c2_derivative(cov, x^2, 0) - val)^2
    result = optimize(x -> f(first(x)), [initial_guess])
    return first(result.minimizer)
end

# c₂ derivative
function c2_derivative(cov::MaternCovariance, s, k::Integer)
    ν = cov.nu
    ϕ = scale(cov)
    k >= ν &&
        throw(ArgumentError("the order k must be less than the smoothness parameter nu"))

    # guard s close to 0 to avoid numerical explosion
    s ≈ 0 && return gamma(ν - k) / 2

    δ = √(s) * √(2ν) / ϕ

    cst = (-1)^k * 2 / 2^ν / ϕ^(2k) * ν^k / gamma(ν)
    return cst * δ^(ν - k) * besselk(ν - k, δ)
end

# Spectral moment
function spectral_moment(cov::MaternCovariance, p::Integer)
    ν = cov.nu
    ϕ = scale(cov)
    return prod((1:2:(2p))) / ϕ^(2p) * ν^p / prod(ν .- (1:p))
end
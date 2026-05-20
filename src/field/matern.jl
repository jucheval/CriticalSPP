"""
    MaternCovariance([phi], nu, d)

Matérn covariance with scale parameter `phi` and smoothness parameter `nu` used for a field of dimension `d`. It is defined as

```math
c(r) = \\frac{2^{1-\\nu}}{\\Gamma(\\nu)} \\left(\\frac{r \\sqrt{2\\nu}}{\\phi}\\right)^\\nu K_\\nu\\left(\\frac{r \\sqrt{2\\nu}}{\\phi}\\right),
```
where ``K_\\nu`` is the modified Bessel function of the second kind.

The scale parameter `phi` is optional and defaults to 1.0.
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

MaternCovariance(nu::T, d::Integer) where {T<:Real} = MaternCovariance(1.0, nu, d)

# Functor
### adapted from MaternVariogram in GeoStats.jl
# TODO: make convert functions so that this one is an interface
function (cov::MaternCovariance{D,T})(r::S) where {D,T,S}
    R = promote_type(T, S)
    return MaternCovariance{D,R}(cov.phi, cov.nu)(R(r))
end
function (cov::MaternCovariance{D,T})(r::T) where {D,T}
    ν = cov.nu
    ϕ = scale(cov)

    # guard r close to 0 to avoid numerical explosion
    r ≈ 0 && return one(T)

    δ = r * √(2 * ν) / ϕ
    Β = besselk(ν, δ)
    Γ = gamma(ν)

    return 2^(1 - ν) / Γ * δ^ν * Β
end

# Practical range
function practical_range(cov::MaternCovariance{D,T}, val::S) where {D,T,S}
    R = promote_type(T, S)
    return practical_range(MaternCovariance{D,R}(cov.phi, cov.nu), R(val))
end
function practical_range(cov::MaternCovariance{D,T}, val::T) where {D,T}
    (0 < val < 1) || throw(DomainError(val, "the value must satisfy 0 < val < 1"))

    initial_guess = practical_range(GaussianCovariance(scale(cov), dimension(cov)), val)
    f(x) = cov(x) - val

    return find_zero(f, initial_guess)
end

# c₂ derivative
# TODO first !!!!!!!!!!     COMMIT      !!!!!!!
# TODO: make convert functions so that this one is an interface
function c2_derivative(cov::MaternCovariance{D,T}, s::S, k::Integer) where {D,T,S}
    R = promote_type(T, S)
    return c2_derivative(MaternCovariance{D,R}(cov.phi, cov.nu), R(s), k)
end
function c2_derivative(cov::MaternCovariance{D,T}, s::T, k::Integer) where {D,T}
    ν = cov.nu
    ϕ = scale(cov)
    k < ν ||
        throw(ArgumentError("the order k must be less than the smoothness parameter nu"))

    cst = (-one(T))^k * 2 / 2^k / ϕ^(2k) * ν^k / gamma(ν)

    # guard s close to 0 to avoid numerical explosion
    s ≈ 0 && return cst * gamma(ν - T(k)) / 2

    δ = √s * √(2 * ν) / ϕ
    return cst * (δ / 2)^(ν - k) * besselk(ν - k, δ)
end

# Spectral moment
function spectral_moment(cov::MaternCovariance, p::Integer)
    ν = cov.nu
    ϕ = scale(cov)
    return prod((1:2:(2p))) / ϕ^(2p) * ν^p / prod(ν .- (1:p))
end

# Internals
function constant_λ₄_over_3λ₂(cov::MaternCovariance)
    ν = cov.nu

    ν > 2 || throw(ArgumentError("the smoothness parameter nu must be greater than 2"))
    return ν / (ν - 2)
end

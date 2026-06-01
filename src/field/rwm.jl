"""
    RWMCovariance([phi], d)

Random Wave Model covariance with scale parameter `phi` used for a field of dimension `d`. It is defined as

```math
c(r) = \\Gamma(d/2) (\\frac{r}{2\\phi'})^{-(d/2-1)} J_{d/2-1}(r / \\phi'), \\quad \\phi' = \\phi/\\sqrt{d},
```
where ``J_\\nu`` is the Bessel function of the first kind.

### Arguments
- `phi::Real=1.0`: scale parameter.
- `d::Int`: spatial dimension, must satisfy `1 <= d <= 4`.

### Returns
- `RWMCovariance`: covariance model instance.

### Examples
```jldoctest
julia> cov = RWMCovariance(1.0, 2);

julia> scale(cov)
1.0

julia> dimension(cov)
2

julia> cov(0.0)
1.0
```
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

# Convert
function _convert_innertype(::Type{S}, cov::RWMCovariance{D,T}) where {D,T,S<:Real}
    return RWMCovariance{D,S}(S(cov.phi))
end

# Functor
function (cov::RWMCovariance{D,T})(r::T) where {D,T}
    ϕ = scale(cov)
    d = T(D)

    # guard r close to 0 to avoid numerical explosion
    r ≈ 0 && return one(T)

    δ = r * sqrt(d) / ϕ
    α = d / 2 - one(T)
    J = besselj(α, δ)
    Γ = gamma(d / 2)

    return Γ * (δ / 2)^(-α) * J
end

# Practical range
function practical_range(cov::RWMCovariance{1,T}, val) where {T}
    return throw(
        ArgumentError(
            "practical range does not make sense for RWM covariance in dimension 1 because the limsup of the covariance is 1 as r goes to infinity",
        ),
    )
end
function practical_range(cov::RWMCovariance{D,T}, val) where {D,T}
    (0 < val < 1) || throw(DomainError(val, "the value must satisfy 0 < val < 1"))

    d = T(D)
    β = (d + 1) / 2 - one(T)
    Γ = gamma(d / 2)

    return scale(cov) * 2 * T(pi)^(-1 / 2β) * d^(-1 / 2) * Γ^(1 / β) / val^(1 / β)
end

# c₂ derivative
function c2_derivative(cov::RWMCovariance{D,T}, s::T, k::Integer) where {D,T}
    ϕ = scale(cov)
    d = T(D)

    cst = (-one(T))^k / 2^k / ϕ^(2k) * (d / 2)^k * gamma(d / 2)

    # guard s close to 0 to avoid numerical explosion
    s ≈ 0 && return cst / gamma(d / 2 + k)

    α = d / 2 - one(T) + k
    δ = √s * √d / ϕ

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

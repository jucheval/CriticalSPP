"""
    CovarianceSPP{D,T<:Real}

Common interface for covariance models used to get a critical spatial point process

### Notes
- `cov::CovarianceSPP` is a functor (see below)
- `D` is the field dimension.
- `T` is the scalar type used for covariance parameters and evaluations.

### Methods (`cov::CovarianceSPP`):
```julia
dimension(cov)             # spatial dimension
scale(cov)                 # scale parameter
practical_range(cov, val)  # practical range corresponding to target level val
c2_derivative(cov, s, k)   # k-th derivative of auxiliary covariance c₂ at s
spectral_moment(cov, p)    # 2p-th spectral moment
```
"""
abstract type CovarianceSPP{D,T<:Real} end

# Parameters
"""
    dimension(cov)

Return the spatial dimension `D` of covariance model `cov`.

### Arguments
- `cov::CovarianceSPP`: covariance model.

### Returns
- `Int`: model dimension.
"""
dimension(::CovarianceSPP{D}) where {D} = D

"""
    scale(cov::CovarianceSPP{D,T}) where {D,T} -> T

Return the scale parameter `phi` of covariance model `cov`.

### Arguments
- `cov::CovarianceSPP`: covariance model with field `phi`.

### Returns
- `T`: scale parameter.
"""
scale(cov::CovarianceSPP) = cov.phi

"""
    otherparameters(cov::CovarianceSPP{D,T}) where {D,T} -> NTuple{N,T}

Return the parameters of covariance model `cov` other than the scale parameter `phi`.

### Arguments
- `cov::CovarianceSPP`: covariance model with field `phi`.

### Returns
- `NTuple{N,T}`: `N` parameters other than the scale parameter.
"""
otherparameters(cov::CovarianceSPP) = ()

# Helper functions
@inline function _check_dimension(::Val{D}) where {D}
    D isa Int || throw(DomainError(typeof(D), "type of D must be Int"))
    (1 <= D <= 4) || throw(DomainError(D, "dimension D must satisfy 1 <= D <= 4"))
    return D
end

"""
    innertype(cov::CovarianceSPP)

Return the type `T` of the parameters of covariance model `cov`.
"""
innertype(::CovarianceSPP{D,T}) where {D,T} = T

function convert_innertype(::Type{S}, cov::CovarianceSPP{D,T}) where {D,T,S<:Real}
    return _convert_innertype(S, cov)
end

# Quantities derived from the covariance function
"""
    (cov::CovarianceSPP)(r)

Evaluate covariance function c₁ at lag `r`.

### Arguments
- `cov::CovarianceSPP`: covariance model.
- `r::Real`: lag.

### Returns
- `Real`: covariance value `c_1(r)`.

### Notes
- The return type is the promoted type between the covariance parameters and `r`.
"""
function (cov::CovarianceSPP{D,T})(r::S) where {D,T,S}
    # TODO: For small lags, the condition number of the covariance matrix is large
    # leading to numerical explosion. It concerns all 3 models.
    # See https://github.com/JuliaEarth/GeoStatsFunctions.jl/issues/80
    R = promote_type(T, S)
    covR = convert_innertype(R, cov)
    return covR(R(r))
end

"""
    practical_range(cov, val)

Return a practical range `r₀` such that covariance values beyond `r₀` are below
target level `val`.

### Arguments
- `cov::CovarianceSPP`: covariance model.
- `val::Real`: target covariance level, typically in `(0, 1)`.

### Returns
- `Real`: practical range.

### Notes
- For models with monotone decreasing covariance, `r₀` solves `cov(r₀) = val`.
- For Random Wave models, the covariance is oscillatory with explicit decreasing asymptotic envelope (except for d=1). In this case, `r₀` solves `envelope(r₀) = val`.
- The return type is the promoted type between the covariance parameters and `val`.
"""
function practical_range(cov::CovarianceSPP{D,T}, val::S) where {D,T,S}
    R = promote_type(T, S)
    covR = convert_innertype(R, cov)
    return practical_range(covR, R(val))
end

"""
    c2_derivative(cov, s, k)

Return the `k`-th derivative of auxiliary covariance c₂ at `s`, i.e. c₂(`s`) = c₁(√`s`), where c₁ is the standard covariance function.

### Arguments
- `cov::CovarianceSPP`: covariance model.
- `s::Real`: squared lag.
- `k::Integer`: derivative order.

### Returns
- `Real`: value of c₂⁽ᵏ⁾(s).

### Notes
- The return type is the promoted type between the covariance parameters and `s`.
"""
function c2_derivative(cov::CovarianceSPP{D,T}, s::S, k::Integer) where {D,T,S}
    R = promote_type(T, S)
    covR = convert_innertype(R, cov)
    return c2_derivative(covR, R(s), k)
end

"""
    spectral_moment(cov, p, [closedform])

Return the `2p`-th spectral moment of covariance model `cov`.

### Arguments
- `cov::CovarianceSPP`: covariance model.
- `p::Integer`: moment index.
- `closedform::Bool`: if `true`, use model-specific closed form; if `false`, compute
    from `c2_derivative(cov, 0, p)`.

### Returns
- `Real`: spectral moment of order `2p`.

### Notes
- The derivative-based path can be numerically unstable for large `p`.
"""
function spectral_moment(cov::CovarianceSPP, p::Integer, closedform::Bool)
    # FIXME: numerical computation goes out of bounds for large p.
    # To see it, fix prange = 1:15 in test/spectral_moment.jl.
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

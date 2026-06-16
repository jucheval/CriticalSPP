"""
    intensity(cpp)

Return the intensity of the critical spatial point process `cpp`.

### Arguments
- `cpp::CriticalPointProcess`: critical point process model.

### Returns
- `Real`: process intensity.

### Examples
```jldoctest
julia> cpp = CriticalPointProcess(GaussianCovariance(1.0, 2), MAX_CRITICAL);

julia> intensity(cpp)
0.09188814923696535
```
"""
function intensity(cpp::CriticalPointProcess)
    cov = covariance(cpp)
    d = dimension(cov)
    type = critical_type(cpp)
    ϕ = scale(cov)
    T = innertype(cpp)

    K = convert(T, constant_λ₄_over_3λ₂(cov))
    cst = convert(T, INTENSITY_CONSTANT_DICT[(d, type)])
    return cst * K^(T(d) / T(2)) * ϕ^(-T(d))
end

"""
    scale_from_intensity(cpp, rho)

Return covariance scale `phi` that matches target intensity `rho` for model `cpp`.

### Arguments
- `cpp::CriticalPointProcess`: model template defining covariance family and critical type.
- `rho::Real`: target intensity, must satisfy `rho > 0`.

### Returns
- `Real`: inferred scale parameter `phi`.

### Notes
- The scale of the covariance embedded in `cpp` is not used.
- Other parameters (dimension, critical type, and covariance shape parameters) are reused.
- The return type is the promoted type between the covariance parameters and `rho`.

### Examples
```jldoctest
julia> cpp0 = CriticalPointProcess(GaussianCovariance(1.0, 2), MAX_CRITICAL);

julia> phi = scale_from_intensity(cpp0, 0.02);

julia> scale_from_intensity(cpp0, 1.0)
0.3031305811642325
```
"""
function scale_from_intensity(cpp::CriticalPointProcess, rho::Real)
    rho > 0 || throw(DomainError(rho, "intensity rho must be positive"))

    cov = covariance(cpp)
    T = typeof(scale(cov))
    ρ = convert(T, rho)
    K = convert(T, constant_λ₄_over_3λ₂(cov))

    type = critical_type(cpp)
    D = dimension(cov)

    cst = convert(T, INTENSITY_CONSTANT_DICT[(D, type)])

    return √K * (cst / ρ)^(inv(T(D)))
end

# Internals

_I = 0.301208 # Quasi Monte Carlo estimation of E( Φ(Y) * Φ(√2Y) ) with Y ~ N(0,1/3)
# See the file computation_I.jl

const INTENSITY_CONSTANT_DICT = Dict(
    (1, ALL_CRITICAL) => √3 / π,
    (1, MAX_CRITICAL) => √3 / (2π),
    (2, ALL_CRITICAL) => 2 / (π√3),
    (2, MAX_CRITICAL) => 0.25 * 2 / (π√3),
    (3, ALL_CRITICAL) => 1 / (π^2 * √2) * 116 / (12√6),
    (3, MAX_CRITICAL) => 1 / (π^2 * √2) * (29 - 6√6) / (12√6),
    (4, ALL_CRITICAL) => 1 / π^2 * 100 / (24√3),
    (4, MAX_CRITICAL) => 1 / π^2 * (_I * 100π - 57) / (48√3 * π),
)

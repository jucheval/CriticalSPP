"""
    intensity(cpp)

Return the intensity of the critical spatial point process `cpp`.
"""
function intensity(cpp::CriticalPointProcess)
    cov = covariance(cpp)
    d = dimension(cov)
    type = critical_type(cpp)

    λ₂ = spectral_moment(cov, 1)
    λ₄ = spectral_moment(cov, 2)
    cst = INTENSITY_CONSTANT_DICT[(d, type)]
    return cst * (λ₄ / (3 * λ₂))^(d / 2)
end

"""
    scale_from_intensity(cpp, rho)

Return the scale parameter `phi` corresponding to the intensity `rho` for the critical spatial point process model inferred from `cpp`. 

*Remark:*
The scale parameter of the input underlying covariance `cpp.cov` is not used, but other necessary informations (e.g. type of critical points, dimension and smoothness parameter) are inferred from `cpp`.
"""
function scale_from_intensity(cpp::CriticalPointProcess, rho::Real)
    rho > 0 || throw(DomainError(rho, "intensity rho must be positive"))

    cov = covariance(cpp)
    K = constant_λ₄_over_3λ₂(cov)

    type = critical_type(cpp)
    D = dimension(cov)

    cst = INTENSITY_CONSTANT_DICT[(D, type)]

    return √K * (cst / rho)^(1 / D)
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
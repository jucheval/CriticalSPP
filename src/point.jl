# Type of critical points

"""
    AbstractCriticalType

Common interface for the type of critical points to consider in a critical spatial point process.
"""
abstract type AbstractCriticalType end

"""
    MaxCritical

Type for considering only local maxima, i.e. with index equal to the dimension of the field `D`, as critical points.
"""
struct MaxCritical <: AbstractCriticalType end

"""
    AllCritical

Type for considering all critical points.
"""
struct AllCritical <: AbstractCriticalType end

const MAX_CRITICAL = MaxCritical()
const ALL_CRITICAL = AllCritical()

"""
    CriticalPointProcess(cov, type)

A critical spatial point process with latent GRF with covariance `cov` (with known scale parameter φ) and type of critical points `type` (possible choices are `MAX_CRITICAL` and `ALL_CRITICAL`).
"""
struct CriticalPointProcess{C<:CovarianceSPP,CT<:AbstractCriticalType}
    cov::C
end

function CriticalPointProcess(
    cov::C, ::CT=MAX_CRITICAL
) where {C<:CovarianceSPP,CT<:AbstractCriticalType}
    return CriticalPointProcess{C,CT}(cov)
end

"""
    CriticalPointProcess(init_cov, type, rho)

When the scale parameter φ of the covariance model is not known, it can be inferred from the target intensity `rho` of the critical point process. 

*Remark:*
The scale parameter of the input covariance `init_cov` is not used, but other necessary informations (e.g. dimension and smoothness parameter) are inferred from `init_cov`.
"""
function CriticalPointProcess(
    init_cov::C, type::CT, rho::Real
) where {C<:CovarianceSPP,CT<:AbstractCriticalType}
    phi = scale_from_intensity(CriticalPointProcess(init_cov, type), rho)
    cov_scaled = C(phi)
    return CriticalPointProcess(cov_scaled, type)
end

"""
    critical_type(cpp)

Return the type of critical points considered in the critical spatial point process `cpp`.
"""
critical_type(::CriticalPointProcess{<:CovarianceSPP,CT}) where {CT<:AbstractCriticalType} =
    CT()

"""
    covariance(cpp)

Return the covariance model underlying the critical spatial point process `cpp`.
"""
covariance(cpp::CriticalPointProcess) = cpp.cov

#-----------------
# IMPLEMENTATIONS
#-----------------
include("point/intensity.jl")
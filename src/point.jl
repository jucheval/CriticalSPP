# Type of critical points
"""
    AbstractCriticalType

Common interface for the type of critical points to consider in a critical spatial point process.
"""
abstract type AbstractCriticalType end

"""
    MaxCritical

Type marker for considering only local maxima as critical points.

# Notes
- Use singleton `MAX_CRITICAL` in constructors.
"""
struct MaxCritical <: AbstractCriticalType end

"""
    AllCritical

Type marker for considering all critical points.

# Notes
- Use singleton `ALL_CRITICAL` in constructors.
"""
struct AllCritical <: AbstractCriticalType end

const MAX_CRITICAL = MaxCritical()
const ALL_CRITICAL = AllCritical()

# Critical spatial point process model
"""
    CriticalPointProcess(cov, type)

A critical spatial point process from covariance model `cov` and critical type `type`.

# Arguments
- `cov::CovarianceSPP`: covariance model with known scale parameter `phi`.
- `type::AbstractCriticalType=MAX_CRITICAL`: critical-point selector.

# Returns
- `CriticalPointProcess`: model instance.

# Notes
- Supported type singletons are `MAX_CRITICAL` and `ALL_CRITICAL`.
- The selected critical type is encoded in the parametric type `CT`.
"""
struct CriticalPointProcess{C<:CovarianceSPP,CT<:AbstractCriticalType}
    cov::C
end

function CriticalPointProcess(
    cov::C, (::CT)=MAX_CRITICAL
) where {C<:CovarianceSPP,CT<:AbstractCriticalType}
    return CriticalPointProcess{C,CT}(cov)
end

"""
    CriticalPointProcess(init_cov, type, rho)

Build a critical point process targeting intensity `rho` by inferring covariance scale.

# Arguments
- `init_cov::CovarianceSPP{D,T}`: covariance model used to define dimension and shape parameters.
- `type::AbstractCriticalType`: critical-point selector.
- `rho::Real`: target intensity, must be positive.

# Returns
- `CriticalPointProcess`: model with scale `phi::T` inferred from `rho`.

# Notes
- The scale of `init_cov` is ignored.
- Model-specific parameters other than scale (for example `nu` for Matérn) are reused.
- The scale of the underlying covariance model is of type `T` whatever the type of `rho`.
"""
function CriticalPointProcess(
    init_cov::C, type::CT, rho::Real
) where {C<:CovarianceSPP,CT<:AbstractCriticalType}
    phi = scale_from_intensity(CriticalPointProcess(init_cov, type), rho)
    cov_scaled = C(phi)
    return CriticalPointProcess(cov_scaled, type)
end

# Parameters
"""
    critical_type(cpp)

Return the type of critical points considered in the critical spatial point process `cpp`.

# Arguments
- `cpp::CriticalPointProcess`: critical point process model.

# Returns
- `AbstractCriticalType`: one of `MaxCritical()` or `AllCritical()`.
"""
function critical_type(
    ::CriticalPointProcess{<:CovarianceSPP,CT}
) where {CT<:AbstractCriticalType}
    return CT()
end

"""
    covariance(cpp)

Return the covariance model underlying the critical spatial point process `cpp`.

# Arguments
- `cpp::CriticalPointProcess`: critical point process model.

# Returns
- `CovarianceSPP`: underlying covariance model.
"""
covariance(cpp::CriticalPointProcess) = cpp.cov

"""
    dimension(cpp)

Return the dimension of the critical Point Process `cpp`, which is the same as the dimension of its covariance model.

# Arguments
- `cpp::CriticalPointProcess`: critical point process model.

# Returns
- `Int`: spatial dimension.
"""
dimension(cpp::CriticalPointProcess) = dimension(covariance(cpp))

# Helper functions
"""
    innertype(cpp::CriticalPointProcess)

Return the type `T` of the parameters of the underlying covariance model `cpp.cov`.
"""
innertype(cpp::CriticalPointProcess) = innertype(covariance(cpp))

function convert_innertype(::Type{S}, cpp::CriticalPointProcess) where {S<:Real}
    cov = covariance(cpp)
    type = critical_type(cpp)
    return CriticalPointProcess(convert_innertype(S, cov), type)
end

#-----------------
# IMPLEMENTATIONS
#-----------------
include("point/intensity.jl")
include("point/pcf.jl")

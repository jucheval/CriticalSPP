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
    CriticalPointProcess{C,CT}

A critical spatial point process defined by a covariance model `C` and a type of critical points `CT`.
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
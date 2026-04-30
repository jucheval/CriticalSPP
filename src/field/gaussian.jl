struct GaussianCovariance{D,T<:Real} <: CovarianceSPP{D,T}
    phi::T
end

function GaussianCovariance{D}(phi::T) where {D,T<:Real}
    _check_dimension(Val(D))
    phi > zero(T) || throw(DomainError(phi, "phi must be positive"))
    return GaussianCovariance{D,T}(phi)
end
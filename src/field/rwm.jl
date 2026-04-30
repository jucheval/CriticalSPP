struct RWMCovariance{D,T<:Real} <: AbstractCovarianceModel{D,T}
    phi::T
end

function RWMCovariance{D}(phi::T) where {D,T<:Real}
    _check_dimension(Val(D))
    phi > zero(T) || throw(DomainError(phi, "phi must be positive"))
    return RWMCovariance{D,T}(phi)
end
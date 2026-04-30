struct MaternCovariance{D,T<:Real} <: AbstractCovarianceModel{D,T}
    phi::T
    nu::T
end

function MaternCovariance{D}(phi::T, nu::T) where {D,T<:Real}
    _check_dimension(Val(D))
    phi > zero(T) || throw(DomainError(phi, "phi must be positive"))
    nu > zero(T) || throw(DomainError(nu, "nu must be positive"))
    return MaternCovariance{D,T}(phi, nu)
end
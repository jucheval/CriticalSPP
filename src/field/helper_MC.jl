# Helper functions for Monte Carlo estimation
"""
    covariance_hessians_x0_xr(cov, r)

Return the covariance matrix of the upper diagonal and diagonal of the Hessians ∇²X(0) and ∇²X(`r`e₁) given that ∇X(0) = ∇X(`r`e₁) = 0, where X is a Gaussian field with covariance `cov`. See Lemma 8 in Azaïs & Delmas (2022).

### Arguments
- `cov::CovarianceSPP{D,T}`: covariance model.
- `r::T`: lag along first coordinate direction.

### Returns
- `Symmetric{T}`: covariance matrix of the concatenated vectorized Hessians.

### Notes
- Type `T` is common to the two arguments which is not restrictive in a regular use case.
"""
function covariance_hessians_x0_xr(cov::CovarianceSPP{D,T}, r::T) where {D,T}
    # -----
    # Notations from Lemma 8 in Azaïs & Delmas (2022)
    d = D
    ρ = r

    𝐫′0 = c2_derivative(cov, zero(T), 1)
    𝐫′′0 = c2_derivative(cov, zero(T), 2)

    𝐫′ = c2_derivative(cov, ρ^2, 1)
    𝐫′′ = c2_derivative(cov, ρ^2, 2)
    𝐫′′′ = c2_derivative(cov, ρ^2, 3)
    𝐫′′′′ = c2_derivative(cov, ρ^2, 4)
    # -----

    # Convention: set at least the upper diagonal and diagonal, 
    # then symmetrize with `Symmetric(..., :U)`.
    if d == 1
        Γ₁ = 12 * 𝐫′′0
        M = (12 * 𝐫′′ + 8ρ^2 * 𝐫′′′)^2
        Γ₁ += ρ^2 * 𝐫′0 / (2 * (𝐫′0^2 - (𝐫′ + 2ρ^2 * 𝐫′′)^2)) * M
        Γ₃ = 12 * 𝐫′′
        Γ₃ += 48ρ^2 * 𝐫′′′ + 16ρ^4 * 𝐫′′′′
        Γ₃ += ρ^2 * (𝐫′ + 2ρ^2 * 𝐫′′) / (2 * (𝐫′0^2 - (𝐫′ + 2ρ^2 * 𝐫′′)^2)) * M
        upper_diag = [Γ₁ Γ₃; Γ₃ Γ₁]
    else
        M = zeros(T, d, d)
        M[2:end, 2:end] .= 16 * 𝐫′′^2
        M[1, 2:end] .= 4 * 𝐫′′ * (12 * 𝐫′′ + 8ρ^2 * 𝐫′′′)
        M[2:end, 1] .= 4 * 𝐫′′ * (12 * 𝐫′′ + 8ρ^2 * 𝐫′′′)
        M[1, 1] = (12 * 𝐫′′ + 8ρ^2 * 𝐫′′′)^2

        Γ₁ = 4 * 𝐫′′0 * ones(T, d, d)
        Γ₁[diagind(Γ₁)] .= 12 * 𝐫′′0
        Γ₁ += ρ^2 * 𝐫′0 / (2 * (𝐫′0^2 - (𝐫′ + 2ρ^2 * 𝐫′′)^2)) * M

        Γ₃ = 4 * 𝐫′′ * ones(T, d, d)
        Γ₃[diagind(Γ₃)] .= 12 * 𝐫′′
        Γ₃[1, 2:end] .+= 8ρ^2 * 𝐫′′′
        Γ₃[2:end, 1] .+= 8ρ^2 * 𝐫′′′
        Γ₃[1, 1] += 48ρ^2 * 𝐫′′′ + 16ρ^4 * 𝐫′′′′
        Γ₃ += ρ^2 * (𝐫′ + 2ρ^2 * 𝐫′′) / (2 * (𝐫′0^2 - (𝐫′ + 2ρ^2 * 𝐫′′)^2)) * M

        D₁ = fill(4 * 𝐫′′0 + 8ρ^2 * 𝐫′′^2 * 𝐫′0 / (𝐫′0^2 - 𝐫′^2), d - 1)
        D₂ = fill(4 * 𝐫′′0, (d - 1) * (d - 2) ÷ 2)
        Γ₂ = diagm([D₁; D₂])

        D̃₁ = fill(4 * 𝐫′′ + 8ρ^2 * (𝐫′′′ + 𝐫′′^2 * 𝐫′ / (𝐫′0^2 - 𝐫′^2)), d - 1)
        D̃₂ = fill(4 * 𝐫′′, (d - 1) * (d - 2) ÷ 2)
        Γ₄ = diagm([D̃₁; D̃₂])

        k = d * (d - 1) ÷ 2
        upper_diag = zeros(T, 2(d + k), 2(d + k))

        upper_diag[1:d, 1:d] = Γ₁
        upper_diag[(d + k + 1):(2d + k), (d + k + 1):(2d + k)] = Γ₁
        upper_diag[(d + 1):(d + k), (d + 1):(d + k)] = Γ₂
        upper_diag[(2d + k + 1):(2d + 2k), (2d + k + 1):(2d + 2k)] = Γ₂
        upper_diag[(1:d), (d + k + 1):(2 * d + k)] = Γ₃
        upper_diag[(d + 1):(d + k), (2 * d + k + 1):(2 * d + 2 * k)] = Γ₄
    end

    return Symmetric(upper_diag, :U)
end

"""
    density_vr(cov, r)

Return the density at 0 of V(`r`) = (∇X(0), ∇X(`r`*e₁)), where X is a Gaussian field with covariance `cov`. See proof of Lemma 8 (beginning of section B.2.1) in Azaïs & Delmas (2022).

### Arguments
- `cov::CovarianceSPP{D,T}`: covariance model.
- `r::T`: lag along first coordinate direction.

### Returns
- `T`: Gaussian density value at zero for `V(r)`.

### Notes
- Type `T` is common to the two arguments which is not restrictive in a regular use case.
"""
function density_vr(cov::CovarianceSPP{D,T}, r::T) where {D,T}
    # -----
    # Notations from Lemma 8 in Azaïs & Delmas (2022)
    d = D
    ρ = r

    𝐫′0 = c2_derivative(cov, zero(T), 1)
    𝐫′ = c2_derivative(cov, ρ^2, 1)
    𝐫′′ = c2_derivative(cov, ρ^2, 2)
    # -----

    # Convention: set at least the upper diagonal and diagonal, 
    # then symmetrize with `Symmetric(..., :U)`.
    upper_diag = diagm(fill(-2𝐫′0, 2d))

    upper_diag[1, d + 1] = -2𝐫′ - 4ρ^2 * 𝐫′′
    for i in 2:d
        upper_diag[i, d + i] = -2𝐫′
    end

    Σ = Symmetric(upper_diag, :U)
    return pdf(MvNormal(Σ), zeros(T, 2d))
end

"""
    det_minor(v, m, d)

Compute the determinant of the m x m minor of a d x d symmetric matrix defined by the entries of `v`, where `v` is a vector of length d + d*(d-1)/2 containing the diagonal and upper diagonal entries of the symmetric matrix. 

### Arguments
- `v::AbstractVector`: vectorized symmetric matrix with specific layout (see notes below)
- `m::Integer`: size of leading principal minor.
- `d::Integer`: full matrix size.

### Returns
- `Real`: determinant of the leading `m x m` principal minor.

### Notes
- The layout of `v` is specific to the Hessian representation used in
    Azaïs and Delmas (2022), not the usual packed symmetric layout.
- The minor is always taken on the first `m` rows and columns.

### Examples
```jldoctest
julia> v = [2.0, 3.0, 1.0]; # corresponds to matrix [2.0 1.0; 1.0 3.0]

julia> CriticalSPP.det_minor(v, 2, 2)
5.0
```
"""
function det_minor(v::AbstractVector, m::Integer, d::Integer)
    if m == 1
        return v[1]
    elseif m == 2
        return v[1] * v[2] - v[d + 1] * v[d + 1]
    elseif m == 3
        return v[1] * (v[2] * v[3] - v[2d] * v[2d]) -
               v[d + 1] * (v[d + 1] * v[3] - v[d + 2] * v[2d]) +
               v[d + 2] * (v[d + 1] * v[2d] - v[d + 2] * v[2])
    else
        d1 =
            v[2] * v[3] * v[4] - v[2] * v[10] * v[10] - v[8] * v[8] * v[4] +
            2 * v[8] * v[9] * v[10] - v[9] * v[9] * v[3]
        d2 =
            v[5] * v[3] * v[4] - v[5] * v[10] * v[10] - v[8] * v[6] * v[4] +
            v[8] * v[7] * v[10] +
            v[9] * v[6] * v[10] - v[9] * v[7] * v[3]
        d3 =
            v[5] * v[8] * v[4] - v[5] * v[9] * v[10] - v[2] * v[6] * v[4] +
            v[2] * v[7] * v[10] +
            v[9] * v[6] * v[9] - v[8] * v[7] * v[9]
        d4 =
            v[5] * v[8] * v[10] - v[5] * v[9] * v[3] - v[2] * v[6] * v[10] +
            v[2] * v[7] * v[3] +
            v[9] * v[6] * v[8] - v[8] * v[7] * v[8]
        return v[1] * d1 - v[5] * d2 + v[6] * d3 - v[7] * d4
    end
end

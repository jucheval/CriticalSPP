# Helper functions for Monte Carlo estimation
"""
    covariance_hessians_x0_xr(cov, r)

Return the covariance matrix of the upper diagonal and diagonal of the Hessians ∇²X(0) and ∇²X(`r`e₁) given that ∇X(0) = ∇X(`r`e₁) = 0, where X is a Gaussian field with covariance `cov`. See Lemma 8 in Azaïs & Delmas (2022).
"""
function covariance_hessians_x0_xr(cov::CovarianceSPP, r)
    # -----
    # Notations from Lemma 8 in Azaïs & Delmas (2022)
    ρ = r

    𝐫′0 = c2_derivative(cov, 0.0, 1)
    𝐫′′0 = c2_derivative(cov, 0.0, 2)

    𝐫′ = c2_derivative(cov, ρ^2, 1)
    𝐫′′ = c2_derivative(cov, ρ^2, 2)
    𝐫′′′ = c2_derivative(cov, ρ^2, 3)
    𝐫′′′′ = c2_derivative(cov, ρ^2, 4)
    # -----

    d = dimension(cov)

    # Convention: set at least the upper diagonal and diagonal, 
    # then symmetrize with `Symmetric(..., :U)`.
    if d == 1
        Γ₁ = 12 * 𝐫′′0
        M = (12 * 𝐫′′ + 8ρ^2 * 𝐫′′′)^2
        Γ₁ = Γ₁ + ρ^2 * 𝐫′0 / (2 * (𝐫′0^2 - (𝐫′ + 2ρ^2 * 𝐫′′)^2)) * M
        Γ₃ = 12 * 𝐫′′
        Γ₃ = Γ₃ + 48ρ^2 * 𝐫′′′ + 16ρ^4 * 𝐫′′′′
        Γ₃ = Γ₃ + ρ^2 * (𝐫′ + 2ρ^2 * 𝐫′′) / (2 * (𝐫′0^2 - (𝐫′ + 2ρ^2 * 𝐫′′)^2)) * M
        upper_diag = [Γ₁ Γ₃; Γ₃ Γ₁]
    else
        M = zeros(d, d)
        M[2:end, 2:end] .= 16 * 𝐫′′^2
        M[1, 2:end] .= 4 * 𝐫′′ * (12 * 𝐫′′ + 8ρ^2 * 𝐫′′′)
        M[2:end, 1] .= 4 * 𝐫′′ * (12 * 𝐫′′ + 8ρ^2 * 𝐫′′′)
        M[1, 1] = (12 * 𝐫′′ + 8ρ^2 * 𝐫′′′)^2

        Γ₁ = 4 * 𝐫′′0 * ones(d, d)
        Γ₁[diagind(Γ₁)] .= 12 * 𝐫′′0
        Γ₁ += ρ^2 * 𝐫′0 / (2 * (𝐫′0^2 - (𝐫′ + 2ρ^2 * 𝐫′′)^2)) * M

        Γ₃ = 4 * 𝐫′′ * ones(d, d)
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
        upper_diag = zeros(2(d + k), 2(d + k))

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

Return the density at 0 of V(`r`) = (∇X(0), ∇X(`r`*e_1)), where X is a Gaussian field with covariance `cov`. See proof of Lemma 8 (beginning of section B.2.1) in Azaïs & Delmas (2022).
"""
function density_vr(cov::CovarianceSPP, r)
    # -----
    # Notations from Lemma 8 in Azaïs & Delmas (2022)
    ρ = r

    𝐫′0 = c2_derivative(cov, 0.0, 1)
    𝐫′ = c2_derivative(cov, ρ^2, 1)
    𝐫′′ = c2_derivative(cov, ρ^2, 2)
    # -----

    d = dimension(cov)

    # Convention: set at least the upper diagonal and diagonal, 
    # then symmetrize with `Symmetric(..., :U)`.
    upper_diag = diagm(fill(-2𝐫′0, 2d))

    upper_diag[1, d + 1] = -2𝐫′ - 4ρ^2 * 𝐫′′
    for i in 2:d
        upper_diag[i, d + i] = -2𝐫′
    end

    Σ = Symmetric(upper_diag, :U)
    return pdf(MvNormal(Σ), zeros(2d))
end

"""
    det_minor(v, m)

Compute the determinant of the m x m minor of a d x d symmetric matrix defined by the entries of `v`, where `v` is a vector of length d + d*(d-1)/2 containing the diagonal and upper diagonal entries of the symmetric matrix. 

*Remark:* `v` does not correspond to the usual vectorization of a symmetric matrix, but rather to the specific layout used in Lemma 8 of Azaïs & Delmas (2022) for the Hessians. The order of entries in `v` is as follows: first the diagonal entries (d of them), then the upper diagonal entries in row-major order (d*(d-1)/2 of them). The minor is defined by the first m rows and columns of the original symmetric matrix.
"""
function det_minor(v::Vector, m::Integer)
    d = Int((-1 + isqrt(1 + 8 * length(v))) ÷ 2)
    if d + d * (d - 1) / 2 != length(v)
        throw(
            DomainError(
                v,
                "The vector length must correspond to the number of unique entries in a symmetric d x d matrix, i.e. d + d*(d-1)/2.",
            ),
        )
    end
    if m > d
        throw(
            DomainError(
                m,
                "The minor order m must be less than or equal to the dimension d inferred from the vector length.",
            ),
        )
    end

    inds = DET_MINOR_INDICES[(d, m)]
    return det_at_positions(v, inds)
end

const DET_MINOR_INDICES = Dict(
    (1, 1) => (1,),
    (2, 1) => (1,),
    (2, 2) => (1, 3, 3, 2),
    (3, 1) => (1,),
    (3, 2) => (1, 4, 4, 2),
    (3, 3) => (1, 4, 5, 4, 2, 6, 5, 6, 3),
    (4, 1) => (1,),
    (4, 2) => (1, 5, 5, 2),
    (4, 3) => (1, 5, 6, 5, 2, 8, 6, 8, 3),
    (4, 4) => (1, 5, 6, 7, 5, 2, 8, 9, 6, 8, 3, 10, 7, 9, 10, 4),
)

"""
    det_at_positions(v::Vector, i::NTuple{N,Int})

Compute the determinant of a d x d matrix defined by the entries of `v` at indices `i`, where `i` is a tuple of N = d^2 integers corresponding to the positions of the entries in the vectorized symmetric matrix. The order of indices in `i` should correspond to the layout of the d x d minor in the original symmetric matrix. It is implemented for `N` equal to 1, 4, 9, and 16, corresponding to 1x1, 2x2, 3x3, and 4x4 minors, respectively.
"""
det_at_positions(v::Vector, i::NTuple{1,Int}) = v[i[1]]

det_at_positions(v::Vector, i::NTuple{4,Int}) = v[i[1]] * v[i[4]] - v[i[2]] * v[i[3]]

function det_at_positions(v::Vector, i::NTuple{9,Int})
    return v[i[1]] * det_at_positions(v, (i[5], i[6], i[8], i[9])) -
           v[i[2]] * det_at_positions(v, (i[4], i[6], i[7], i[9])) +
           v[i[3]] * det_at_positions(v, (i[4], i[5], i[7], i[8]))
end

function det_at_positions(v::Vector, i::NTuple{16,Int})
    return v[i[1]] * det_at_positions(
        v, (i[6], i[7], i[8], i[10], i[11], i[12], i[14], i[15], i[16])
    ) -
           v[i[2]] * det_at_positions(
        v, (i[5], i[7], i[8], i[9], i[11], i[12], i[13], i[15], i[16])
    ) +
           v[i[3]] * det_at_positions(
        v, (i[5], i[6], i[8], i[9], i[10], i[12], i[13], i[14], i[16])
    ) -
           v[i[4]] *
           det_at_positions(v, (i[5], i[6], i[7], i[9], i[10], i[11], i[13], i[14], i[15]))
end
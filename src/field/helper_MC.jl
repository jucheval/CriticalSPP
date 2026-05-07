# Helper functions for Monte Carlo estimation
"""
    covariance_hessians_x0_xr(cov, r)

Return the covariance matrix of the upper diagonal and diagonal of the Hessians ∇²X(0) and ∇²X(`r`e₁) given that ∇X(0) = ∇X(`r`e₁) = 0, where X is a Gaussian field with covariance `cov`. See Lemma 8 in Azaïs & Delmas (2022).
"""
function covariance_hessians_x0_xr(cov::CovarianceSPP, r)
    C0 = map(k -> c2_derivative(cov, 0.0, k), 0:4)
    Cr = map(k -> c2_derivative(cov, r^2, k), 0:4)

    d = dimension(cov)

    # Convention: set at least the upper diagonal and diagonal, 
    # then symmetrize with `Symmetric(..., :U)`.
    if d == 1
        Γ₁ = 12 * C0[3]
        M = (12 * Cr[3] + 8 * r^2 * Cr[4])^2
        Γ₁ = Γ₁ + r^2 * C0[2] / (2 * (C0[2]^2 - (Cr[2] + 2 * r^2 * Cr[3])^2)) * M
        Γ₃ = 12 * Cr[3]
        Γ₃ = Γ₃ + 48 * r^2 * Cr[4] + 16 * r^4 * Cr[5]
        Γ₃ =
            Γ₃ +
            r^2 * (Cr[2] + 2 * r^2 * Cr[3]) /
            (2 * (C0[2]^2 - (Cr[2] + 2 * r^2 * Cr[3])^2)) * M
        upper_diag = [Γ₁ Γ₃; Γ₃ Γ₁]
    else
        M = zeros(d, d)
        M[2:end, 2:end] .= 16 * Cr[3]^2
        M[1, 2:end] .= 4 * Cr[3] * (12 * Cr[3] + 8r^2 * Cr[4])
        M[2:end, 1] .= 4 * Cr[3] * (12 * Cr[3] + 8r^2 * Cr[4])
        M[1, 1] = (12 * Cr[3] + 8r^2 * Cr[4])^2

        Γ₁ = 4 * C0[3] * ones(d, d)
        Γ₁[diagind(Γ₁)] .= 12 * C0[3]
        Γ₁ += r^2 * C0[2] / (2 * (C0[2]^2 - (Cr[2] + 2r^2 * Cr[3])^2)) * M

        Γ₃ = 4 * Cr[3] * ones(d, d)
        Γ₃[diagind(Γ₃)] .= 12 * Cr[3]
        Γ₃[1, 2:end] .+= 8r^2 * Cr[4]
        Γ₃[2:end, 1] .+= 8r^2 * Cr[4]
        Γ₃[1, 1] += 48r^2 * Cr[4] + 16r^4 * Cr[5]
        Γ₃ += r^2 * (Cr[2] + 2r^2 * Cr[3]) / (2 * (C0[2]^2 - (Cr[2] + 2r^2 * Cr[3])^2)) * M

        D₁ = fill(4 * C0[3] + 8r^2 * Cr[3]^2 * C0[2] / (C0[2]^2 - Cr[2]^2), d - 1)
        D₂ = fill(4 * C0[3], (d - 1) * (d - 2) ÷ 2)
        Γ₂ = diagm([D₁; D₂])

        D̃₁ = fill(
            4 * Cr[3] + 8r^2 * (Cr[4] + Cr[3]^2 * Cr[2] / (C0[2]^2 - Cr[2]^2)), d - 1
        )
        D̃₂ = fill(4 * Cr[3], (d - 1) * (d - 2) ÷ 2)
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

Return the density at 0 of V(`r`) = (∇X(0), ∇X(`r`*e_1)), where X is a Gaussian field with covariance `cov`.
"""
function density_vr(cov::CovarianceSPP, r) end
"""
    pair_correlation_function(cpp, rs; kwargs...)

Estimate the pair correlation function of the critical point process `cpp` at
lags `rs` with Monte Carlo.

# Arguments
- `cpp::CriticalPointProcess`: critical point process model.
- `rs::AbstractVector{<:Real}`: vector of lags.

# Keyword Arguments
- `n_MC::Integer=100_000`: number of Monte Carlo replications per lag.
- `parallel::Symbol=:auto`: execution policy, one of `:auto`, `:serial`, `:threads`.
- `rng::Random.AbstractRNG=Random.default_rng()`: random number generator.

# Returns
- `NamedTuple`: fields `rs`, `pcf`, `stderr`, `eig_val`, `n_MC`.
"""
function pair_correlation_function(
    cpp::CriticalPointProcess,
    rs::AbstractVector{<:Real};
    n_MC::Integer=100_000,
    parallel::Symbol=:auto,
    rng::AbstractRNG=default_rng(),
)
    _validate_nMC(n_MC)
    _validate_parallel(parallel)

    nlag = length(rs)
    if nlag == 0
        return (
            rs=Float64[],
            pcf=Float64[],
            stderr=Float64[],
            eig_val=Tuple{Int,Float64}[],
            n_MC=n_MC,
        )
    end

    use_threads = _thread_policy(parallel, nlag, n_MC)
    pcf, stderr, eig_val = if use_threads
        _pair_correlation_threaded(cpp, rs, n_MC, rng)
    else
        _pair_correlation_serial(cpp, rs, n_MC, rng)
    end

    return (rs=collect(Float64, rs), pcf=pcf, stderr=stderr, eig_val=eig_val, n_MC=n_MC)
end

@inline function _validate_parallel(parallel::Symbol)
    (parallel == :auto || parallel == :serial || parallel == :threads) ||
        throw(ArgumentError("parallel must be :auto, :serial, or :threads"))
    return parallel
end

@inline function _validate_nMC(n_MC::Integer)
    n_MC > 0 || throw(DomainError(n_MC, "n_MC must be positive"))
    return n_MC
end

function _thread_policy(parallel::Symbol, nlag::Integer, n_MC::Integer)
    nworkers = Threads.nthreads()
    if parallel == :serial
        @info "Starting pcf estimation with 1 thread for $nlag lag(s) and $n_MC Monte Carlo replications per lag."
        return false
    elseif parallel == :threads
        @info "Starting pcf estimation with $nworkers threads for $nlag lag(s) and $n_MC Monte Carlo replications per lag."
        return true
    else
        # Auto policy: thread only if total work is large enough and threads exist.
        work = nlag * n_MC
        if work >= Int(1e7)
            @info "Starting pcf estimation with $nworkers threads for $nlag lag(s) and $n_MC Monte Carlo replications per lag."
            return true
        else
            @info "Starting pcf estimation with 1 thread for $nlag lag(s) and $n_MC Monte Carlo replications per lag."
            return false
        end
    end
end

function _pair_correlation_serial(
    cpp::CriticalPointProcess, rs::AbstractVector{<:Real}, n_MC::Integer, rng::AbstractRNG
)
    D = dimension(cpp)
    nlag = length(rs)
    dd = D + D * (D - 1) ÷ 2

    pcf = Vector{Float64}(undef, nlag)
    stderr = Vector{Float64}(undef, nlag)
    eig_val = Vector{Tuple{Int,Float64}}(undef, nlag)
    N01 = randn(rng, 2dd, n_MC)

    for i in ProgressBar(eachindex(rs))
        pcf[i], stderr[i], eig_val[i] = _pair_correlation_single(cpp, rs[i], N01)
    end

    return pcf, stderr, eig_val
end

function _pair_correlation_threaded(
    cpp::CriticalPointProcess, rs::AbstractVector{<:Real}, n_MC::Integer, rng::AbstractRNG
)
    D = dimension(cpp)
    nlag = length(rs)
    dd = D + D * (D - 1) ÷ 2

    pcf = Vector{Float64}(undef, nlag)
    stderr = Vector{Float64}(undef, nlag)
    eig_val = Vector{Tuple{Int,Float64}}(undef, nlag)
    N01 = randn(rng, 2dd, n_MC)

    Threads.@threads for i in ProgressBar(eachindex(rs))
        pcf[i], stderr[i], eig_val[i] = _pair_correlation_single(cpp, rs[i], N01)
    end

    return pcf, stderr, eig_val
end

"""
    _pair_correlation_single(cpp, r, N01)

Return g(`r`) estimated by Monte Carlo (and its standard error) for a single lag `r`. `N01` is a pre-generated matrix with i.i.d. N(0,1) entries used for the Monte Carlo estimation. It has dimensions `(2D', n_MC)` where `D' = D + D*(D-1)/2` is the dimension of the vectorised Hessian.
"""
function _pair_correlation_single(cpp::CriticalPointProcess, r::Real, N01::AbstractMatrix)
    D = dimension(cpp)
    cov = covariance(cpp)
    dd = D + D * (D - 1) ÷ 2 # dimension of Hessian matrices
    size(N01, 1) == 2dd ||
        throw(DimensionMismatch("N01 must have $(2dd) rows for dimension D=$D"))
    n_MC = size(N01, 2)

    Σ = covariance_hessians_x0_xr(cov, r)
    λ, P = eigen(Σ)
    if !all(λ .> 0)
        eig_val = (sum(λ .< 0), mean(λ[λ .< 0]))
    else
        eig_val = (0, NaN)
    end

    diag = Diagonal(sqrt.(max.(λ, 0)))
    ξ = P * diag * N01  # ξ is 2dd x n_MC

    # Welford's online algorithm to compute mean and standard error in a single pass without storing all values.
    ξ0 = @view ξ[1:dd, 1]
    ξr = @view ξ[(dd + 1):(2 * dd), 1]
    value = argument_of_expectation(critical_type(cpp), ξ0, ξr)
    M = value # running mean
    S = 0.0 # running sum of squared deviations
    count = 1
    for i in 2:n_MC
        ξ0 = @view ξ[1:dd, i]
        ξr = @view ξ[(dd + 1):(2 * dd), i]
        value = argument_of_expectation(critical_type(cpp), ξ0, ξr)

        count += 1
        new_M = M + (value - M) / count
        S = S + (value - M) * (value - new_M)
        M = new_M
    end

    p∇X0∇Xr = density_vr(cov, r)
    ρ = intensity(cpp)
    factor = p∇X0∇Xr / (ρ^2)
    M *= factor
    S *= factor^2

    pcf = M
    stderr = sqrt(S / (n_MC - 1))
    # stderr *= 1.959964 / sqrt(n_MC) # uncomment if testing compliance with R code

    return (pcf, stderr, eig_val)
    # return (pcf, stderr, ξ) # uncomment if testing compliance with R code
end

"""
    argument_of_expectation(<:AbstractCriticalType, ξ0, ξr)

Compute the argument of the expectation in Eq (35) of Azaïs & Delmas (2022) for the given type of critical points, where `ξ0` and `ξr` represent the Hessians ∇²X(0) and ∇²X(t) in the vectorised form used in Lemma 8 of Azaïs & Delmas (2022).
"""
function argument_of_expectation(::AllCritical, ξ0::AbstractVector, ξr::AbstractVector)
    length(ξ0) == length(ξr) ||
        throw(DimensionMismatch("ξ0 and ξr must have the same length"))

    d = d_from_vec(ξ0)
    det0 = det_minor(ξ0, d)
    detr = det_minor(ξr, d)
    return abs(det0) * abs(detr)
end

function argument_of_expectation(::MaxCritical, ξ0::AbstractVector, ξr::AbstractVector)
    length(ξ0) == length(ξr) ||
        throw(DimensionMismatch("ξ0 and ξr must have the same length"))

    d = d_from_vec(ξ0)

    m = 1
    tmp0 = det_minor(ξ0, m)
    tmpr = det_minor(ξr, m)
    ismax = (tmp0 < 0) && (tmpr < 0) # if true, the point is a local maximum candidate
    ismin = (tmp0 > 0) && (tmpr > 0) # if true, the point is a local minimum candidate
    # a local minimum corresponds to a positive definite Hessian
    # i.e. all leading principal minors are positive
    while ismin && m < d
        m += 1
        tmp0 = det_minor(ξ0, m)
        tmpr = det_minor(ξr, m)
        ismin = (tmp0 > 0) && (tmpr > 0)
    end
    # a local maximum corresponds to minus Hessian being positive definite 
    # i.e. all leading principal minors of -H are positive
    while ismax && m < d
        m += 1
        tmp0 = det_minor(ξ0, m)
        tmpr = det_minor(ξr, m)
        ismax = ((-1)^m * tmp0 > 0) && ((-1)^m * tmpr > 0)
    end

    if ismax || ismin
        # By symmetry, the quantity for local maxima is the half of 
        # the quantity for all extrema
        return 0.5 * abs(det_minor(ξ0, d)) * abs(det_minor(ξr, d))
    else
        return 0.0
    end
end

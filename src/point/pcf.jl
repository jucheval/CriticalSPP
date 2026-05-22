"""
    pair_correlation_function([rng], cpp, rs; kwargs...)

Estimate the pair correlation function of the critical point process `cpp` at
lags `rs` with Monte Carlo using the random number generator `rng`.

### Arguments
- `rng::AbstractRNG`: random number generator.
- `cpp::CriticalPointProcess`: critical point process model.
- `rs::AbstractVector{<:Real}`: lags where the function is estimated.

### Keywords
- `n_MC::Integer=100_000`: number of Monte Carlo replications per lag.
- `parallel::Symbol=:auto`: execution policy, one of `:auto`, `:serial`, `:threads`.
- `show_progress::Bool=true`: show progress bar.

### Returns
- `NamedTuple`: fields `rs`, `pcf`, `stderr`, `n_MC`.

### Notes
- `pcf` is the Monte Carlo estimation.
- `stderr` is the Monte Carlo standard deviation estimate of the estimator (not a confidence interval width).
- The number of threads is determined by the `JULIA_NUM_THREADS` environment variable (see https://docs.julialang.org/en/v1/manual/multi-threading/)
- Threaded execution can lead to small numerical differences versus serial mode.

### Examples
```jldoctest
julia> using Random

julia> rng = MersenneTwister(42);

julia> cpp = CriticalPointProcess(GaussianCovariance(1.0, 2), MAX_CRITICAL);

julia> out = pair_correlation_function(rng, cpp, [1.0, 2.0]; n_MC=1000, parallel=:serial, show_progress=false);

julia> out.pcf
2-element Vector{Float64}:
 0.07345406896204186
 0.845105195553286
```
"""
function pair_correlation_function(
    rng::AbstractRNG,
    cpp::CriticalPointProcess,
    rs::AbstractVector{<:Real};
    n_MC::Integer=100_000,
    parallel::Symbol=:auto,
    show_progress::Bool=true,
)
    _check_nMC(n_MC)
    _check_parallel(parallel)

    nlag = length(rs)
    if nlag == 0
        return (rs=Float64[], pcf=Float64[], stderr=Float64[], n_MC=n_MC)
    end

    use_threads = _thread_policy(parallel, nlag, n_MC)
    pcf, stderr = if use_threads
        _pair_correlation_threaded(rng, cpp, rs, n_MC, show_progress)
    else
        _pair_correlation_serial(rng, cpp, rs, n_MC, show_progress)
    end

    return (rs=collect(Float64, rs), pcf=pcf, stderr=stderr, n_MC=n_MC)
end

function pair_correlation_function(
    cpp::CriticalPointProcess, rs::AbstractVector{<:Real}; kwargs...
)
    return pair_correlation_function(default_rng(), cpp, rs; kwargs...)
end

@inline function _check_parallel(parallel::Symbol)
    (parallel == :auto || parallel == :serial || parallel == :threads) ||
        throw(ArgumentError("parallel must be :auto, :serial, or :threads"))
    return parallel
end

@inline function _check_nMC(n_MC::Integer)
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
    rng::AbstractRNG,
    cpp::CriticalPointProcess,
    rs::AbstractVector{<:Real},
    n_MC::Integer,
    show_progress::Bool,
)
    D = dimension(cpp)
    nlag = length(rs)
    dd = D + D * (D - 1) ÷ 2

    pcf = similar(rs, Float64)
    stderr = similar(rs, Float64)
    N01 = randn(rng, 2dd, n_MC)

    p = Progress(nlag; enabled=show_progress)
    for i in eachindex(rs)
        pcf[i], stderr[i] = _pair_correlation_single(cpp, rs[i], N01)
        next!(p)
    end

    return pcf, stderr
end

function _pair_correlation_threaded(
    rng::AbstractRNG,
    cpp::CriticalPointProcess,
    rs::AbstractVector{<:Real},
    n_MC::Integer,
    show_progress::Bool,
)
    D = dimension(cpp)
    nlag = length(rs)
    dd = D + D * (D - 1) ÷ 2

    pcf = similar(rs, Float64)
    stderr = similar(rs, Float64)
    N01 = randn(rng, 2dd, n_MC)

    p = Progress(nlag; enabled=show_progress)
    Threads.@threads for i in eachindex(rs)
        pcf[i], stderr[i] = _pair_correlation_single(cpp, rs[i], N01)
        next!(p)
    end

    return pcf, stderr
end

"""
    _pair_correlation_single(cpp, r, N01)

Return g(`r`) estimated by Monte Carlo (and its standard error) for a single lag `r`. `N01` is a pre-generated matrix with i.i.d. N(0,1) entries used for the Monte Carlo estimation. It has dimensions `(2D', n_MC)` where `D' = D + D*(D-1)/2` is the dimension of the vectorised Hessian.

### Arguments
- `cpp::CriticalPointProcess`: critical point process model.
- `r::Real`: single lag.
- `N01::AbstractMatrix`: standard normal samples of shape `(2D', n_MC)`.

### Returns
- `(pcf, stderr)`: Monte Carlo estimate and associated standard deviation estimate.
"""
function _pair_correlation_single(cpp::CriticalPointProcess, r::Real, N01::AbstractMatrix)
    D = dimension(cpp)
    cov = covariance(cpp)
    dd = D + D * (D - 1) ÷ 2 # dimension of Hessian matrices
    n_MC = size(N01, 2)

    Σ = covariance_hessians_x0_xr(cov, r)
    λ, P = eigen(Σ)
    # # The following is used for debugging only
    # # it checks whether the covariance matrix is positive definite
    # if !all(λ .> 0) 
    #     eig_val = (sum(λ .< 0), mean(λ[λ .< 0]))
    # else
    #     eig_val = (0, NaN)
    # end

    diag = Diagonal(sqrt.(max.(λ, 0)))
    ξ = P * diag * N01  # ξ is 2dd x n_MC

    # Welford's online algorithm to compute mean and standard error in a single pass without storing all values.
    ξ0 = @view ξ[1:dd, 1]
    ξr = @view ξ[(dd + 1):(2 * dd), 1]
    value = argument_of_expectation(critical_type(cpp), ξ0, ξr, D)
    M = value # running mean
    S = 0.0 # running sum of squared deviations
    count = 1
    for i in 2:n_MC
        ξ0 = @view ξ[1:dd, i]
        ξr = @view ξ[(dd + 1):(2 * dd), i]
        value = argument_of_expectation(critical_type(cpp), ξ0, ξr, D)

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

    return pcf, stderr
    # return (pcf, stderr, ξ) # uncomment (and comment previous line) if testing compliance with R code
end

"""
    argument_of_expectation(<:AbstractCriticalType, ξ0, ξr, d)

Compute the argument of the expectation in Eq (35) of Azaïs & Delmas (2022) for the given type of critical points, where `ξ0` and `ξr` represent the Hessians ∇²X(0) and ∇²X(t) in the vectorised form used in Lemma 8 of Azaïs & Delmas (2022). The dimension `d` is passed to avoid recomputing it from the vector length.

### Arguments
- `<:AbstractCriticalType`: type of critical points (MAX_CRITICAL or ALL_CRITICAL).
- `ξ0::AbstractVector`, `ξr::AbstractVector`: vectorized Hessians at two locations.
- `d::Integer`: field dimension.

### Returns
- `Real`: integrand value used in Monte Carlo estimation.
"""
function argument_of_expectation(
    ::AllCritical, ξ0::AbstractVector, ξr::AbstractVector, d::Integer
)
    det0 = det_minor(ξ0, d, d)
    detr = det_minor(ξr, d, d)
    return abs(det0) * abs(detr)
end

function argument_of_expectation(
    ::MaxCritical, ξ0::AbstractVector, ξr::AbstractVector, d::Integer
)
    m = 1
    tmp0 = det_minor(ξ0, m, d)
    tmpr = det_minor(ξr, m, d)
    ismax = (tmp0 < 0) && (tmpr < 0) # if true, the point is a local maximum candidate
    ismin = (tmp0 > 0) && (tmpr > 0) # if true, the point is a local minimum candidate
    # a local minimum corresponds to a positive definite Hessian
    # i.e. all leading principal minors are positive (Sylvester's criterion)
    while ismin && m < d
        m += 1
        tmp0 = det_minor(ξ0, m, d)
        tmpr = det_minor(ξr, m, d)
        ismin = (tmp0 > 0) && (tmpr > 0)
    end
    # a local maximum corresponds to minus Hessian being positive definite 
    # i.e. all leading principal minors of -H are positive
    while ismax && m < d
        m += 1
        tmp0 = det_minor(ξ0, m, d)
        tmpr = det_minor(ξr, m, d)
        ismax = ((-1)^m * tmp0 > 0) && ((-1)^m * tmpr > 0)
    end

    if ismax || ismin
        # In that branch, tmp0 = det_minor(ξ0, d, d) and tmpr = det_minor(ξr, d, d) 
        return 0.5 * abs(tmp0) * abs(tmpr)
        # By symmetry, the quantity for local maxima is the half of 
        # the quantity for all extrema
    else
        return 0.0
    end
end

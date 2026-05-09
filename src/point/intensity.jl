"""
    intensity(cpp)

Return the intensity of the critical spatial point process `cpp`.
"""
function intensity(cpp::CriticalPointProcess)
    cov = covariance(cpp)
    d = dimension(cov)
    type = critical_type(cpp)

    λ₂ = spectral_moment(cov, 1)
    λ₄ = spectral_moment(cov, 2)
    cst = _intensity_constant_dict[(d, type)]
    return cst * (λ₄ / (3 * λ₂))^(d / 2)
end

"""
    scale_from_intensity(cpp, rho)

Return the scale parameter `phi` corresponding to the intensity `rho` for the critical spatial point process model inferred from `cpp`. More precisely, the scale parameter of the input underlying covariance `cpp.cov` is not used, but other necessary informations (e.g. type of critical points, dimension and smoothness parameter) are inferred from `cpp`.
"""
function scale_from_intensity(cpp::CriticalPointProcess, rho::Real)
    rho > 0 || throw(DomainError(rho, "intensity rho must be positive"))

    cov = covariance(cpp)
    if cov isa MaternCovariance
        cov.nu > 2 || throw(
            DomainError(
                cpp,
                "the smoothness parameter ν of the Matern covariance must be greater than 2",
            ),
        )
    end

    type = critical_type(cpp)
    D = dimension(cov)

    # re-implementation of the intensity function 
    # to avoid CovarianceSPP and CriticalPointProcess constructors 
    # in the root-finding loop
    cst = _intensity_constant_dict[(D, type)]

    if cov isa GaussianCovariance
        sp_moment = _spectral_moment_gaussian
        args = ()
    elseif cov isa MaternCovariance
        sp_moment = _spectral_moment_matern
        args = (cov.nu,)
    elseif cov isa RWMCovariance
        sp_moment = _spectral_moment_RWM
        args = (D,)
    else
        throw(ArgumentError("Unsupported covariance type: $(typeof(cov))"))
    end

    function f(ϕ)
        λ₂ = sp_moment(1, ϕ, args...)
        λ₄ = sp_moment(2, ϕ, args...)
        rho′ = cst * (λ₄ / (3 * λ₂))^(D / 2)
        return rho′ - rho
    end

    return find_zero(f, (1e-6, 1e6))
end

# """
#     pair_correlation_function(cpp, rs; kwargs...)

# Estimate the pair correlation function of the critical point process `cpp` at
# lags `rs` with Monte Carlo.

# # Arguments
# - `cpp::CriticalPointProcess`: critical point process model.
# - `rs::AbstractVector{<:Real}`: vector of lags.

# # Keyword Arguments
# - `n_MC::Integer=10_000`: number of Monte Carlo replications per lag.
# - `parallel::Symbol=:auto`: execution policy, one of `:auto`, `:serial`, `:threads`.
# - `num_threads::Union{Nothing,Integer}=nothing`: requested thread cap. This does not change Julia's global thread count.
# - `rng::Random.AbstractRNG=Random.default_rng()`: random number generator.

# # Returns
# - `NamedTuple`: fields `r`, `pcf`, `se`, `n_MC`.
# """
# function pair_correlation_function(
#     cpp::CriticalPointProcess,
#     rs::AbstractVector{<:Real};
#     n_MC::Integer=10_000,
#     parallel::Symbol=:auto,
#     num_threads::Union{Nothing,Integer}=nothing,
#     rng::Random.AbstractRNG=Random.default_rng(),
# )
#     _validate_nMC(n_MC)
#     _validate_parallel(parallel)
#     _validate_num_threads(num_threads)

#     nlag = length(rs)
#     if nlag == 0
#         return (r=Float64[], pcf=Float64[], se=Float64[], n_MC=n_MC)
#     end

#     use_threads, nworkers = _thread_policy(parallel, num_threads, nlag, n_MC)
#     pcf, se = if use_threads
#         _pair_correlation_threaded(cpp, rs, n_MC, nworkers, rng)
#     else
#         _pair_correlation_serial(cpp, rs, n_MC, rng)
#     end

#     return (r=collect(Float64, rs), pcf=pcf, se=se, n_MC=n_MC)
# end

# @inline function _validate_parallel(parallel::Symbol)
#     (parallel == :auto || parallel == :serial || parallel == :threads) ||
#         throw(ArgumentError("parallel must be :auto, :serial, or :threads"))
#     return parallel
# end

# @inline function _validate_num_threads(num_threads::Union{Nothing,Integer})
#     isnothing(num_threads) && return nothing
#     num_threads > 0 || throw(DomainError(num_threads, "num_threads must be positive"))
#     return num_threads
# end

# @inline function _validate_nMC(n_MC::Integer)
#     n_MC > 0 || throw(DomainError(n_MC, "n_MC must be positive"))
#     return n_MC
# end

# function _thread_policy(
#     parallel::Symbol, num_threads::Union{Nothing,Integer}, nlag::Integer, n_MC::Integer
# )
#     available = Base.Threads.nthreads()
#     requested = isnothing(num_threads) ? available : min(Int(num_threads), available)

#     if parallel == :serial
#         return false, 1
#     elseif parallel == :threads
#         return requested > 1, max(requested, 1)
#     else
#         # Auto policy: thread only if total work is large enough and threads exist.
#         work = nlag * n_MC
#         use_threads = (requested > 1) && (work >= 50_000)
#         return use_threads, (use_threads ? requested : 1)
#     end
# end

# function _pair_correlation_serial(
#     cpp::CriticalPointProcess,
#     rs::AbstractVector{<:Real},
#     n_MC::Integer,
#     rng::Random.AbstractRNG,
# )
#     nlag = length(rs)
#     pcf = Vector{Float64}(undef, nlag)
#     se = Vector{Float64}(undef, nlag)
#     for i in eachindex(rs)
#         pcf[i], se[i] = _pair_correlation_single(cpp, rs[i], n_MC, rng)
#     end
#     return pcf, se
# end

# function _pair_correlation_threaded(
#     cpp::CriticalPointProcess,
#     rs::AbstractVector{<:Real},
#     n_MC::Integer,
#     nworkers::Integer,
#     rng::Random.AbstractRNG,
# )
#     nlag = length(rs)
#     g = Vector{Float64}(undef, nlag)
#     se = Vector{Float64}(undef, nlag)

#     chunks = _chunk_ranges(nlag, nworkers)
#     seeds = rand(rng, UInt, length(chunks))
#     tasks = Vector{Task}(undef, length(chunks))

#     for j in eachindex(chunks)
#         rchunk = chunks[j]
#         seed = seeds[j]
#         tasks[j] = Base.Threads.@spawn begin
#             local_rng = Random.Xoshiro(seed)
#             for i in rchunk
#                 g[i], se[i] = _pair_correlation_single(cpp, rs[i], n_MC, local_rng)
#             end
#             nothing
#         end
#     end

#     foreach(fetch, tasks)
#     return g, se
# end

# function _chunk_ranges(n::Integer, k::Integer)
#     k_eff = max(1, min(Int(k), Int(n)))
#     q, r = divrem(n, k_eff)
#     ranges = Vector{UnitRange{Int}}(undef, k_eff)
#     start = 1
#     for j in 1:k_eff
#         len = q + (j <= r ? 1 : 0)
#         stop = start + len - 1
#         ranges[j] = start:stop
#         start = stop + 1
#     end
#     return ranges
# end

# """
#     _pair_correlation_single(cpp, r, n_MC, rng)

# Return `(g_hat(r), se_hat(r))` for a single lag `r`.

# Implement this method with the actual Monte Carlo estimator logic.
# """
# function _pair_correlation_single(
#     ::CriticalPointProcess, ::Real, ::Integer, ::Random.AbstractRNG
# )
#     throw(ArgumentError("Implement _pair_correlation_single for your Monte Carlo kernel"))
# end

_I = 0.30120734 # MC estimation of E( Φ(Y) * Φ(√2Y)) with Y ~ N(0,1/3) on 1e9 samples

_intensity_constant_dict = Dict(
    (1, ALL_CRITICAL) => √3 / π,
    (1, MAX_CRITICAL) => √3 / (2π),
    (2, ALL_CRITICAL) => 2 / (π√3),
    (2, MAX_CRITICAL) => 0.25 * 2 / (π√3),
    (3, ALL_CRITICAL) => 1 / (π^2 * √2) * 116 / (12√6),
    (3, MAX_CRITICAL) => 1 / (π^2 * √2) * (29 - 6√6) / (12√6),
    (4, ALL_CRITICAL) => 1 / π^2 * 100 / (24√3),
    (4, MAX_CRITICAL) => 1 / π^2 * (_I * 100π - 57) / (48√3 * π),
)

sp_moment_dict = Dict(
    "gaussian" => _spectral_moment_gaussian,
    "matern" => _spectral_moment_matern,
    "RWM" => _spectral_moment_RWM,
)
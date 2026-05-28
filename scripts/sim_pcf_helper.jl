# Two dictionaries to map string values used in `params`
const TYPE = Dict("all" => ALL_CRITICAL, "max" => MAX_CRITICAL)
const COVARIANCE = Dict(
    "Gaussian" => d -> GaussianCovariance(d),
    ["Matern(" .* string.(nu) .* ")" .=> (d -> MaternCovariance(nu, d)) for nu in NUs]...,
    "RWM" => d -> RWMCovariance(d),
)

# Take a dictionary of parameters, run the simulation and return a dictionary with the results
function make_pcf(dict::Dict)
    @unpack rho, type, cov, d, nMC = dict
    cpp = CriticalPointProcess(COVARIANCE[cov](d), TYPE[type], rho)

    fmax = type == "all" ? 1.5 : 2.0
    fmin = 0.005
    val = 0.1
    # val = which.cov == "RWM" ? 0.2 : 0.1
    pr = practical_range(cpp.cov, val)
    rmax = pr * fmax
    rmin = pr * fmin

    rs = range(rmin, rmax; length=Lr)
    pcf = pair_correlation_function(cpp, rs; n_MC=nMC)

    output = copy(dict)
    output["rs"] = rs
    output["pcf"] = pcf.pcf
    output["stderr"] = pcf.stderr
    return output
end

using DrWatson
@quickactivate
using CriticalSPP

# relative path to the directory to store results
savedir = joinpath(@__DIR__, "..", "data", "sims", "pcf")

# Length of the range of lag values (used in `make_pcf`)
const Lr = 200

# Parameters:
const NUs = [2.5, 3.5, 4.5]
params = Dict(
    "rho" => 100.0,
    "type" => ["all", "max"],
    "cov" => ["Gaussian"; "Matern(" .* string.(NUs) .* ")"; "RWM"],
    "d" => [1, 2, 3],
    "nMC" => Int(1e4),
)
dicts = dict_list(params)
### remove RWM with d=1 (pcf is not well-defined)
popat!(dicts, 28) # type="max"
popat!(dicts, 25) # type="all"

include("sim_pcf_helper.jl") # define `make_pcf`

# Run simulations if not already done
for config in dicts
    @produce_or_load(make_pcf, config, savedir)
end

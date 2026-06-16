## Preamble
using DrWatson
@quickactivate
using CriticalSPP
using Distributions
using DataFrames, DataFramesMeta
using CairoMakie, AlgebraOfGraphics

# relative paths to the directories
loaddir = joinpath(@__DIR__, "..", "data", "sims", "pcf")
plotdir = joinpath(@__DIR__, "..", "plots")

# collect each file as a row
df_wide = collect_results(loaddir)

## Data wrangling
# reshape `df_wide` because each row contains a vector of 
# lag values, pcf values and standard errors
df = DataFrame()
for row in eachrow(df_wide)
    @unpack d, rs, nMC, pcf, stderr, type, rho, cov = row
    for (r, p, s) in zip(rs, pcf, stderr)
        push!(df, (d=d, r=r, nMC=nMC, pcf=p, stderr=s, type=type, rho=rho, cov=cov))
    end
end

# compute the confidence interval for the pcf estimation
@chain df begin
    @rtransform! :a = quantile(Normal(), 0.975) ./ sqrt.(:nMC) # multiplicative factor for the confidence interval of the exact recovery rate
    @rtransform! :pcf_lower = :pcf - :a * :stderr
    @rtransform! :pcf_upper = :pcf + :a * :stderr
    @select! Not([:a, :nMC, :stderr]) # remove useless columns
    @rsubset!(:pcf < 3.0) # remove exploding values
end

## Plot
using SpecialFunctions

function funI(rs, pcfs, d, rho)
    sB1 = 2 * pi^(d / 2) / gamma(d / 2)
    step = mean(diff(rs))
    return 1 .+ rho * step * cumsum(sB1 * rs .^ (d - 1) .* (pcfs .- 1))
end

df_I = @chain df begin
    @groupby :type :d :cov
    @transform @astable begin
        d, rho = first(:d), first(:rho)
        :I = funI(:r, :pcf, d, rho)
        :I_lower = funI(:r, :pcf_lower, d, rho)
        :I_upper = funI(:r, :pcf_upper, d, rho)
    end
end

# labels for the facets columns
const D_LABELS = Dict(1 => L"d=1", 2 => L"d=2", 3 => L"d=3")
@rtransform! df_I :d = D_LABELS[:d]
# labels for the facets rows
const TYPE_LABELS = Dict("all" => L"\mathcal{L}=\{0,...,d\}", "max" => L"\mathcal{L}=\{d\}")
@rtransform! df_I :type = TYPE_LABELS[:type]

# solid lines for the pcf
lines = visual(Lines) * mapping(:r, :I; color=:cov, row=:type, col=:d)
# bands for the confidence intervals
band =
    visual(Band; alpha=0.3) * mapping(:r, :I_lower, :I_upper; color=:cov, row=:type, col=:d)
# dashed line for the reference value of 1
href = visual(HLines; color=:black, linestyle=:dash) * mapping(1.0)
# combine the three layers into a single plot object
plt = data(df_I) * (lines + band) + href

# plot
set_theme!(theme_ggplot2())
fig, grid = draw(
    plt,
    scales(; X=(; label=L"r"), Y=(; label=L"I_{\mathcal{L}}(r)"));
    figure=(; size=(800, 500)),
    facet=(; linkxaxes=:minimal, linkyaxes=:none),
    legend=(; position=:top, titlesize=0),
)

## Save
safesave(joinpath(plotdir, "I.pdf"), fig)

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

# helper dataframe to compute the minimal r range for each type and dimension
rminmax_df = @chain df begin
    @groupby :type :d :cov
    @combine :rmax = maximum(:r)
    @groupby :type :d
    @combine :rminmax = minimum(:rmax)
end
# filter dataframe to get same r ranges on each facet

## Plot
### define a new dataframe for the plot so that `df` keeps all data
df_plot = @chain df begin
    leftjoin(rminmax_df; on=[:type, :d])
    @rsubset! :r <= :rminmax
    @select! Not(:rminmax)
end

# labels for the facets columns
const D_LABELS = Dict(1 => L"d=1", 2 => L"d=2", 3 => L"d=3")
@rtransform! df_plot :d = D_LABELS[:d]
# labels for the facets rows
const TYPE_LABELS = Dict("all" => L"\mathcal{L}=\{0,...,d\}", "max" => L"\mathcal{L}=\{d\}")
@rtransform! df_plot :type = TYPE_LABELS[:type]

# solid lines for the pcf
lines = visual(Lines) * mapping(:r, :pcf; color=:cov, row=:type, col=:d)
# bands for the confidence intervals
band =
    visual(Band; alpha=0.3) *
    mapping(:r, :pcf_lower, :pcf_upper; color=:cov, row=:type, col=:d)
# dashed line for the reference value of 1
href = visual(HLines; color=:black, linestyle=:dash) * mapping(1.0)
# combine the three layers into a single plot object
plt = data(df_plot) * (lines + band) + href

# plot
set_theme!(theme_ggplot2())
fig, grid = draw(
    plt,
    scales(; X=(; label=L"r"), Y=(; label=L"g_{\mathcal{L}}(r)"));
    figure=(; size=(800, 500)),
    facet=(; linkxaxes=:minimal, linkyaxes=:none),
    legend=(; position=:top, titlesize=0),
)

## Save
safesave(joinpath(plotdir, "pcf.pdf"), fig)

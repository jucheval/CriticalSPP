## Preamble
using DrWatson
@quickactivate
using CriticalSPP
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

## Plot
# keep only small values for asymptotic behavior
### define a new dataframe so that `df` keeps all data
df_plot = @rsubset df :r <= 0.02

# labels for the facets columns
const D_LABELS = Dict(1 => L"d=1", 2 => L"d=2", 3 => L"d=3")
@rtransform! df_plot :d = D_LABELS[:d]
# labels for the facets rows
const TYPE_LABELS = Dict("all" => L"\mathcal{L}=\{0,...,d\}", "max" => L"\mathcal{L}=\{d\}")
@rtransform! df_plot :type = TYPE_LABELS[:type]

# solid lines and markers for the pcf
lines_with_markers =
    (visual(Lines) + visual(Scatter) * mapping(; marker=:cov)) *
    mapping(:r => log, :pcf => log; color=:cov, row=:type, col=:d)
# apply the layer to `df_plot`
plt = data(df_plot) * lines_with_markers

# plot
set_theme!(theme_ggplot2())
fig, grid = draw(
    plt,
    scales(; X=(; label=L"\log(r)"), Y=(; label=L"\log(g_{\mathcal{L}}(r))"));
    figure=(; size=(800, 500)),
    facet=(; linkxaxes=:all, linkyaxes=:none),
    legend=(; position=:top, titlesize=0),
)

## Save
safesave(joinpath(plotdir, "pcflog.pdf"), fig)

## Slopes of log-log curves as a table
using GLM
using Latexify
df_table = @rsubset df :r <= 0.02
@rsubset! df_table (:d > 1 || :r <= 0.002) :pcf > 0.0

df_table = DataFrame([
    (
        cov=gdf[1, :cov],
        d=gdf[1, :d],
        type=gdf[1, :type],
        slope=coef(lm(@formula(log(pcf) ~ log(r)), gdf))[2],
    ) for gdf in groupby(df_table, [:type, :d, :cov])
])
@rtransform! df_table :slope = round(:slope; digits=2)

@chain df_table begin
    sort(:type)
    @rtransform!(:col_name = string(:type, "_", :d))
    unstack(:cov, :col_name, :slope)
    latexify(; env=:table, booktabs=true, latex=true)
    print()
end

## Preamble
using DrWatson
@quickactivate
using CriticalSPP
using DataFrames, DataFramesMeta
using CairoMakie, AlgebraOfGraphics

## Generate dataframe
# Parameters
const NUs = [2.5, 3.5, 4.5]
const RHO = 100.0
const RMAX = Dict(1 => 0.03, 2 => 0.3, 3 => 0.5)
const Lr = 1000
const FACTOR = Dict(1 => 1e-9, 2 => 1e-4, 3 => 1e-3)
const COVS = ["Gaussian"; "Matern(" .* string.(NUs) .* ")"; "RWM"]
const COVARIANCE = Dict(
    "Gaussian" => d -> GaussianCovariance(d),
    ["Matern(" .* string.(nu) .* ")" .=> (d -> MaternCovariance(nu, d)) for nu in NUs]...,
    "RWM" => d -> RWMCovariance(d),
)

# Dataframe
df = DataFrame()
for d in 1:3
    covs = [GaussianCovariance(d), MaternCovariance.(NUs, d)..., RWMCovariance(d)]
    for cov_str in COVS
        cpp = CriticalPointProcess(COVARIANCE[cov_str](d), ALL_CRITICAL, RHO)
        cov = CriticalSPP.covariance(cpp)
        rmax = RMAX[d]
        rs = range(0.0, rmax; length=Lr)
        c₂′r² = map(r -> c2_derivative(cov, r^2, 1), rs)
        c₂′′r² = map(r -> c2_derivative(cov, r^2, 2), rs)
        c₁′′r = 2 * c₂′r² + 4 * rs .^ 2 .* c₂′′r²
        diff = FACTOR[d] * (c₁′′r[1]^2 .- c₁′′r .^ 2)
        append!(df, DataFrame(; d=d, r=rs, diff=diff, diff_type="c1", cov=cov_str))
        diff = FACTOR[d] * (c₂′r²[1]^2 .- c₂′r² .^ 2)
        append!(df, DataFrame(; d=d, r=rs, diff=diff, diff_type="c2", cov=cov_str))
    end
end

## Plot
df_ymax = @chain df begin
    @rsubset :cov != "RWM"
    @groupby :d :diff_type
    @combine :ymax = 1.1 * maximum(:diff)
end

df_plot = @chain df begin
    leftjoin(df_ymax; on=[:d, :diff_type])
    @rtransform! :diff = :diff <= :ymax ? :diff : NaN
    @select! Not(:ymax)
end

# labels for the facets columns
const D_LABELS = Dict(1 => L"d=1", 2 => L"d=2", 3 => L"d=3")
@rtransform! df_plot :d = D_LABELS[:d]
# labels for the facets rows
const ROW_LABELS = Dict(
    "c1" => L"c_1\prime\prime(0)^2 - c_1\prime\prime(r)^2",
    "c2" => L"c_2\prime(0)^2 - c_2\prime(r)^2",
)
@rtransform! df_plot :diff_type = ROW_LABELS[:diff_type]

# solid lines for the curves
lines = visual(Lines) * mapping(:r, :diff; color=:cov, row=:diff_type, col=:d)
# dashed line for the reference value of 0
href = visual(HLines; color=:black, linestyle=:dash) * mapping(0.0)
# combine the three layers into a single plot object
plt = data(df_plot) * lines + href

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
safesave(joinpath(plotdir, "c1c2.pdf"), fig)

using DrWatson, Test
using CriticalSPP
using LinearAlgebra

@quickactivate # Activate test project environment
using RCall

# Execute the following command in the terminal:
# julia --project=. --color=yes test/comparison_with_R/runtests.jl

# Run test suite
ti = time()
include("call_Rcode.jl")
println("Starting tests - Compliance with R code")
println("----")
@testset verbose = true "c2_derivative" begin
    include("c2_derivative.jl")
end
println("----")
@testset verbose = true "spectral_moment" begin
    include("spectral_moment.jl")
end
println("----")
@testset verbose = true "practical_range" begin
    include("practical_range.jl")
end
println("----")
@testset verbose = true "intensity" begin
    include("intensity.jl")
end
println("----")
@testset verbose = true "density_vr" begin
    include("density_vr.jl")
end
println("----")
@testset verbose = true "covariance_hessians" begin
    include("covariance_hessians.jl")
end
println("----")
ti = time() - ti
println("Test took total time of:")
println(round(ti / 60; digits=3), " minutes")

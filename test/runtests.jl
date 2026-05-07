using DrWatson, Test
using CriticalSPP
using LinearAlgebra

# The Pkg manager ]test macro does not work for DrWatson projects
# We need to include the test files manually
# For instance, execute the following command in the terminal:
# julia --project=. --color=yes test/runtests.jl

# Run test suite
ti = time()
println("Starting tests")
println("----")
@testset verbose = true "Spectral moments (numerical and closed-form compliance)" begin
    include("spectral_moment.jl")
end
println("----")
@testset verbose = true "Practical range (inverse of covariance)" begin
    include("practical_range.jl")
end
println("----")
@testset verbose = true "Scale from intensity (inverse of intensity)" begin
    include("scale_intensity.jl")
end
println("----")
@testset verbose = true "Helper functions for MC estimation" begin
    include("helper_MC.jl")
end
println("----")
ti = time() - ti
println("Test took total time of:")
println(round(ti / 60; digits=3), " minutes")

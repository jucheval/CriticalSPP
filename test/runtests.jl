using DrWatson, Test
using CriticalSPP

# The Pkg manager ]test macro does not work for DrWatson projects
# We need to include the test files manually
# For instance, execute the following command in the terminal:
# julia --project=. --color=yes test/runtests.jl

# Run test suite
println("Starting tests\n")
ti = time()

@testset verbose = true "Spectral moments (numerical and closed-form comply)" begin
    include("spectral_moment.jl")
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti / 60; digits=3), " minutes")

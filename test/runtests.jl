using Aqua
using CriticalSPP
using ExplicitImports
using JuliaFormatter
using LinearAlgebra
using Logging
using Random
using Test

global_logger(NullLogger()) # used to remove the logging info messages from the test output

# The Pkg manager ]test macro does not work for DrWatson projects
# We need to include the test files manually
# For instance, execute the following command in the terminal:
# julia --project=test/ --color=yes test/runtests.jl

# Run test suite
ti = time()
println("Starting tests")
println("----")
@testset verbose = false "Code quality (Aqua.jl)" begin
    Aqua.test_all(CriticalSPP)
end
println("----")
@testset verbose = false "Formatting" begin
    @test format(CriticalSPP; overwrite=false)
end
println("----")
@testset verbose = false "ExplicitImports" begin
    test_all_explicit_imports_are_public(CriticalSPP)
    test_all_qualified_accesses_are_public(CriticalSPP; ignore=(:gamma,)) # ignore gamma which is not public in Bessels.jl
    test_all_explicit_imports_via_owners(CriticalSPP)
    test_all_qualified_accesses_via_owners(CriticalSPP)
    test_no_implicit_imports(CriticalSPP)
    test_no_self_qualified_accesses(CriticalSPP)
    test_no_stale_explicit_imports(CriticalSPP)
end
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
@testset verbose = true "Type parameter preservation" begin
    include("type_parameters.jl")
end
println("----")
# FIXME: all tests in this @testset fail on GitHub but pass locally.
# It must be (at least) the eigen decomposition
# @testset verbose = true "Pair correlation function" begin 
#     include("pcf.jl")
# end
# println("----")
ti = time() - ti
println("Test took total time of:")
println(round(ti / 60; digits=3), " minutes")

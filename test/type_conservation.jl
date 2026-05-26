const _CRITICAL_TYPES = (ALL_CRITICAL, MAX_CRITICAL)
const _FLOAT_TYPES = (Float32, Float64, BigFloat)

# Helper functions used in try/catch blocks
### Missing BigFloat support for Bessel functions in SpecialFunctions.jl
function is_bessel_methoderror(err)
    return err isa MethodError &&
           (err.f === CriticalSPP.besselk || err.f === CriticalSPP.besselj)
end
### Missing eigen! method for Symmetric{BigFloat} in LinearAlgebra.jl
function is_eigen!_methoderror(err)
    return err isa MethodError && err.args[2] === eigen!
end

@testset "Constructors" begin
    include("type_conservation/constructors.jl")
end

@testset "Inner type and conversion" begin
    include("type_conservation/innertype.jl")
end

@testset "Random field functions" begin
    include("type_conservation/field_functions.jl")
end

@testset "Point process functions" begin
    include("type_conservation/point_functions.jl")
end

@testset "Helper functions" begin
    include("type_conservation/helper_functions.jl")
end

@testset "PCF function" begin
    # BigFloat is removed because there is no method for eigen decomposition
    # of Symmetric{BigFloat} matrices in LinearAlgebra
    for T in _FLOAT_TYPES, S in _FLOAT_TYPES, D in [2], ct in _CRITICAL_TYPES
        cov = GaussianCovariance(T(1.5), D)
        cpp = CriticalPointProcess(cov, ct)
        rs = S.(1:2)
        R = promote_type(T, S)

        try
            output = pair_correlation_function(cpp, rs; show_progress=false, n_MC=10)
            @test eltype(output.rs) == R
            @test eltype(output.pcf) == R
            @test eltype(output.stderr) == R
        catch err
            if R == BigFloat && is_eigen!_methoderror(err)
                @test_skip "BigFloat eigen! support missing in LinearAlgebra"
            else
                rethrow(err)
            end
        end
    end
end

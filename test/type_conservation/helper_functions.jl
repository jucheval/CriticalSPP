@testset "det_minor" begin
    for T in _FLOAT_TYPES, D in 1:4
        v = rand(T, D + D * (D - 1) ÷ 2)
        @test CriticalSPP.det_minor(v, 1, D) isa T
    end
end

@testset "argument_of_expectation" begin
    for T in _FLOAT_TYPES, D in 1:4, ct in _CRITICAL_TYPES
        ξ0 = rand(T, D + D * (D - 1) ÷ 2)
        ξr = rand(T, D + D * (D - 1) ÷ 2)

        @test CriticalSPP.argument_of_expectation(ct, ξ0, ξr, D) isa T
    end
end

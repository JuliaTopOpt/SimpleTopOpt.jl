using MAT
using SimpleTopOpt
using Test

"""
Unit test runner
"""

@testset "Topflow Suite" begin
    include("unit_tests/ut_topflow.jl")
end




@testset "Top88 Suite" begin

    @testset "Unit tests" begin
        #include("unit_tests/ut_top88.jl")
    end

    @testset "Integration" begin
        #include("unit_tests/int_top88.jl")
    end

end

@testset "TopH Suite" begin

    @testset "Unit tests" begin
        #include("unit_tests/ut_toph.jl")
    end

    @testset "Integration" begin
        #include("unit_tests/int_toph.jl")
    end

end

@testset "Benchmarking" begin

    @testset "top88" begin
        #include("benchmarks/bm_top88.jl")
    end

    @testset "toph" begin
        #include("benchmarks/bm_toph.jl")
    end

end

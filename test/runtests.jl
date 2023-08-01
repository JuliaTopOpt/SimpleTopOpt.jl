using MAT
using SimpleTopOpt
using Test

"""
Unit test runner
"""

all_tests = (isempty(ARGS) || "all" in ARGS)


if all_tests || "topflow" in ARGS
    println("Beginning topflow suite...")

    @testset "Topflow Suite" begin
        # include("unit_tests/ut_topflow.jl")
    end

    @testset "Topflow integration test suite" begin
        # include("unit_tests/symbolics/ut_integration.jl")
    end

    @testset "Topflow symbolic suite" begin
        include("unit_tests/symbolics/ut_symbolics.jl")
    end

    println("Finished Topflow Suite...")
end

if all_tests || "top88" in ARGS
    println("Beginning Top88 suite...")

    @testset "Top88 Suite" begin

        @testset "Unit tests" begin
            include("unit_tests/ut_top88.jl")
        end

        @testset "Integration" begin
            include("unit_tests/int_top88.jl")
        end
    end

    println("Finished Top88 suite...")
end


if all_tests || "toph" in ARGS
    println("Beginning TopH suite...")

    @testset "TopH Suite" begin

        @testset "Unit tests" begin
            include("unit_tests/ut_toph.jl")
        end

        @testset "Integration" begin
            include("unit_tests/int_toph.jl")
        end
    end

    println("Finished TopH suite...")
end


if all_tests || "benchmark" in ARGS

    println("Beginning benchmarks...")

    @testset "Benchmarking" begin

        @testset "top88" begin
            include("benchmarks/bm_top88.jl")
        end

        @testset "toph" begin
            include("benchmarks/bm_toph.jl")
        end
    end

    println("Finished benchmarks...")
end

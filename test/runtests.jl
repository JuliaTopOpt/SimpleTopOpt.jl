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
        include("unit_tests/topflow/ut_topflow.jl")
    end

    @testset "Topflow intermediates Problem type 1" begin
        include("unit_tests/topflow/ut_intermediates_type1.jl")
    end

    @testset "Topflow intermediates Problem type 2" begin
        # include("unit_tests/topflow/ut_intermediates_type2.jl")
    end

    @testset "Topflow end-to-end tests" begin
        include("unit_tests/topflow/int_topflow.jl")
    end

    @testset "Topflow symbolic integration tests" begin
        # include("unit_tests/topflow/ut_integration.jl")
    end
    @testset "MATLAB ffi tests" begin
        # include("unit_tests/topflow/matlab_ffi_tests/ffi_tests.jl")
    end
    @testset "Topflow symbolic suite" begin
        # include("unit_tests/topflow/symbolics/ut_symbolics.jl")
    end


    println("Finished Topflow Suite...")
end

if all_tests || "top88" in ARGS
    println("Beginning Top88 suite...")
    @testset "Top88 Suite" begin
        @testset "Unit tests" begin
            include("unit_tests/top88/ut_top88.jl")
        end
        @testset "Integration" begin
            include("unit_tests/top88/int_top88.jl")
        end
    end
    println("Finished Top88 suite...")
end


if all_tests || "toph" in ARGS
    println("Beginning TopH suite...")
    @testset "TopH Suite" begin
        @testset "Unit tests" begin
            include("unit_tests/toph/ut_toph.jl")
        end
        @testset "Integration" begin
            include("unit_tests/toph/int_toph.jl")
        end
    end
    println("Finished TopH suite...")
end


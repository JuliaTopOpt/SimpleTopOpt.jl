using BenchmarkTools
using LinearAlgebra
using MAT
using SimpleTopOpt
using Test


@testset "Top88 Suite" begin
    # First comparison on 60x40
    vars = matread("mat_cases/top88_60_40_04_3_2_1.mat")
    x1 = vars["x"]
    x1h,_,_ = top88(60, 40, 0.4, 3.0, 2.0, true, false)

    @test size(x1h) == size(x1)
    ss = size(x1h)
    num_elements = ss[1] * ss[2]

    @test (mean(x1h) - mean(x1)) ≈ 0 atol=0.001
    @test (norm(x1h - x1))/num_elements ≈ 0 atol=0.01

    # Second comparison on 30x30
    vars = matread("mat_cases/top88_30_30_04_3_2_1.mat")
    x2 = vars["ans"]
    x2h,_,_ = top88(30, 30, 0.4, 3.0, 2.0, true, false)

    @test size(x2h) == size(x2)
    ss = size(x2h)
    num_elements = ss[1] * ss[2]

    @test (mean(x2h) - mean(x2)) ≈ 0 atol=0.001
    @test (norm(x2h-x2))/num_elements ≈ 0 atol=0.01

    # Benchmarking
    io = IOContext(stdout)

    println("\n -- Benchmarking for 4x4")
    b1 = @benchmark top88(4, 4, 0.4, 3.0, 2.0, true, false);
    show(io, MIME("text/plain"), b1)
    # MATLAB takes approx 0.4 seconds

    println("\n -- Benchmarking for 10x10")
    b2 = @benchmark top88(10, 10, 0.4, 3.0, 2.0, true, false)
    show(io, MIME("text/plain"), b2)
    # MATLAB takes approx 1.6 seconds

    println("\n -- Benchmarking for 16x16")
    b3 = @benchmark top88(16, 16, 0.4, 3.0, 2.0, true, false)
    show(io, MIME("text/plain"), b3)
    # MATLAB takes approx 0.0353 seconds
end

@testset "TopH comparison" begin

    # First comparison on 40x40
    vars = matread("mat_cases/toph_40_40_04_3_12.mat")
    x3 = vars["x"]

    x3h,_,_ = toph(40, 40, 0.4, 3.0, 1.2)

    @test size(x3h) == size(x3)

    ss = size(x3h)
    num_elements = ss[1] * ss[2]

    residuals = abs.(x3h - x3)
    @test (mean(residuals)) ≈ 0 atol=1e-6
    @test (norm(residuals))/num_elements ≈ 0 atol=1e-6

    # Second comparison on 20x20
    vars = matread("mat_cases/toph_20_20_04_3_12.mat")
    x4 = vars["ans"]

    x4h,_,_ = toph(20, 20, 0.4, 3.0, 1.2)

    @test size(x4h) == size(x4)
    
    ss = size(x4h)
    num_elements = ss[1] * ss[2]

    residuals = abs.(x4h - x4)
    @test (mean(residuals)) ≈ 0 atol=1e-6
    @test (norm(residuals))/num_elements ≈ 0 atol=1e-6

    # Benchmarking
    io = IOContext(stdout)

    # NOTE: multiples of 20 for compatibility without needing to impose floor/ ceil restrictions

    println("\n -- Benchmarking for 20x20")
    b1 = @benchmark toph(20, 20, 0.4, 3.0, 2.0)
    show(io, MIME("text/plain"), b1)
    # Takes matlab about t=2.61 to reach change=0.01

    println("\n -- Benchmarking for 16x16")
    b2 = @benchmark toph(40, 40, 0.4, 3.0, 2.0)
    show(io, MIME("text/plain"), b2)
    # Takes MATLAB about t=7.50 to reach change=0.009
end
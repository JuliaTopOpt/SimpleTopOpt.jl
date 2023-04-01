using BenchmarkTools
using LinearAlgebra
using MAT
using SimpleTopOpt
using Test

@testset "Top88 comparison" begin

    println("Starting Top88 comparison suite...")

    vars = matread("mat_cases/top88_60_40_04_3_2_1.mat")
    x1 = vars["x"]

    x1h,_,_ = top88(60, 40, 0.4, 3.0, 2.0, true, false)

    @test size(x1h) == size(x1)

    ss = size(x1h)
    num_elements = ss[1] * ss[2]

    @test (mean(x1h) - mean(x1)) ≈ 0 atol=0.001
    @test (norm(x1h - x1))/num_elements ≈ 0 atol=0.01
    println("Finished first case...")

    vars = matread("mat_cases/top88_30_30_04_3_2_1.mat")
    x1 = vars["ans"]

    x1h,_,_ = top88(30, 30, 0.4, 3.0, 2.0, true, false)

    @test size(x1h) == size(x1)
    
    ss = size(x1h)
    num_elements = ss[1] * ss[2]

    @test (mean(x1h) - mean(x1)) ≈ 0 atol=0.001
    @test (norm(x1h-x1))/num_elements ≈ 0 atol=0.01
    println("Finished second case...")
end

@testset "Top88 benchmarking" begin
    println("Starting Top88 benchmarking suite")

    println("Finished the benchmarking suite")
end


@testset "TopH comparison" begin
    println("Starting TopH comparison suite...")

    vars = matread("mat_cases/toph_40_40_04_3_12.mat")
    x1 = vars["x"]

    x1h,_,_ = toph(40, 40, 0.4, 3.0, 1.2)

    @test size(x1h) == size(x1)

    ss = size(x1h)
    num_elements = ss[1] * ss[2]

    @test (mean(x1h) - mean(x1)) ≈ 0 atol=0.001
    @test (norm(x1h - x1))/num_elements ≈ 0 atol=0.01
    println("Finished first case...")

    vars = matread("mat_cases/toph_80_80_04_3_12.mat")
    x1 = vars["ans"]

    x1h,_,_ = toph(80, 80, 0.4, 3.0, 1.2)

    @test size(x1h) == size(x1)
    
    ss = size(x1h)
    num_elements = ss[1] * ss[2]

    @test (mean(x1h) - mean(x1)) ≈ 0 atol=0.001
    @test (norm(x1h-x1))/num_elements ≈ 0 atol=0.01
    println("Finished the comparison suite")
end


@testset "TopH benchmarking" begin
    println("Starting TopH benchmarking suite")


end






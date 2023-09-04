using MAT
using SimpleTopOpt
using Test
using LinearAlgebra
using Statistics

"""
Integration test suite for Top88
"""

domain_1 = Top88Domain(60, 40)
domain_2 = Top88Domain(30, 30)

SIMP = ModifiedSIMPParameters(penal=3.0)
optimizer = OptimalityCriteria()
sensitivity_filter = SensitivityFilter()
density_filter = DensityFilter()

sensitivity_problem_1 = Top88Problem(domain_1, SIMP, optimizer, sensitivity_filter)
sensitivity_problem_2 = Top88Problem(domain_2, SIMP, optimizer, sensitivity_filter)

density_problem_1 = Top88Problem(domain_1, SIMP, optimizer, density_filter)
density_problem_2 = Top88Problem(domain_2, SIMP, optimizer, density_filter)

@testset "60x40 Integration Suite" begin
    # First comparison on 60x40
    vars = matread("mat_cases/top88_60_40_04_3_2_1.mat")
    t_design_field = vars["x"]

    df_sol = SimpleTopOpt.Top88.optimize(sensitivity_problem_1)
    design_field = df_sol.design

    @test df_sol.converged == true

    @test size(t_design_field) == size(design_field)
    ss = size(design_field)
    num_elements = ss[1] * ss[2]

    @test_broken (mean(design_field) - mean(t_design_field)) ≈ 0 atol = 1e-10
    @test_broken sum(abs.(design_field - t_design_field)) / num_elements ≈ 0 atol=1e-10
    @test_broken (norm(design_field - t_design_field)) / num_elements ≈ 0 atol = 1e-10
end

@testset "30x30 Integration Suite" begin
    # Second comparison on 30x30
    vars = matread("mat_cases/top88_30_30_04_3_2_1.mat")
    t_design_field = vars["ans"]

    df_sol = SimpleTopOpt.Top88.optimize(sensitivity_problem_2)
    design_field = df_sol.design

    @test size(t_design_field) == size(design_field)

    ss = size(design_field)
    num_elements = ss[1] * ss[2]

    @test_broken (mean(design_field) - mean(t_design_field)) ≈ 0 atol = 1e-10
    @test_broken sum(abs.(design_field - t_design_field)) / num_elements ≈ 0 atol=1e-10
    @test_broken (norm(design_field - t_design_field)) / num_elements ≈ 0 atol = 1e-10
end

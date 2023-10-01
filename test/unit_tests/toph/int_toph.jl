using MAT
using SimpleTopOpt
using Test
using LinearAlgebra
using Statistics

"""
Integration test suite for TopH
"""

domain_1 = TophDomain(40, 40)
domain_2 = TophDomain(20, 20)

SIMP = ModifiedSIMPParameters(penal=3.0)
optimizer = OptimalityCriteria()
sensitivity_filter = SensitivityFilter()
density_filter = DensityFilter()

sensitivity_problem_1 = TophProblem(domain_1, SIMP, optimizer, sensitivity_filter)
sensitivity_problem_2 = TophProblem(domain_2, SIMP, optimizer, sensitivity_filter)

density_problem_1 = TophProblem(domain_1, SIMP, optimizer, density_filter)
density_problem_2 = TophProblem(domain_2, SIMP, optimizer, density_filter)

@testset "40x40 Comparison" begin
    # First comparison on 40x40
    vars = matread("mat_cases/toph_40_40_04_3_12.mat")
    t_design_field = vars["x"]

    sol = SimpleTopOpt.TopH.optimize(sensitivity_problem_1)
    design_field = sol.design

    @test sol.converged == true

    @test size(design_field) == size(t_design_field)

    ss = size(t_design_field)
    num_elements = ss[1] * ss[2]

    @test_broken mean(design_field - t_design_field) ≈ 0 atol = 1e-10
    @test_broken norm(design_field - t_design_field) / num_elements ≈ 0 atol = 1e-10
    @test_broken sum(abs.(design_field - t_design_field)) / num_elements ≈ 0 atol=1e-10
end

@testset "20x20 Comparison" begin
    # Second comparison on 20x20
    vars = matread("mat_cases/toph_20_20_04_3_12.mat")
    t_design_field = vars["ans"]

    sol = SimpleTopOpt.TopH.optimize(sensitivity_problem_2)
    design_field = sol.design

    @test size(t_design_field) == size(design_field)

    ss = size(t_design_field)
    num_elements = ss[1] * ss[2]

    @test_broken mean(design_field - t_design_field) ≈ 0 atol=1e-10
    @test_broken norm(design_field - t_design_field) / num_elements ≈ 0 atol=1e-10
    @test_broken sum(abs.(design_field - t_design_field)) / num_elements ≈ 0 atol=1e-10
end

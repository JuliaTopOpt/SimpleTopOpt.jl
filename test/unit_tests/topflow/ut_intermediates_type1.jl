"""
This test suite is focused on the intermediates for the
topflow optimization loop
"""

using Test
using SimpleTopOpt
using SimpleTopOpt.TopFlow: OCUpdate, newton, compute_adjoint_solution, compute_sensitivities
using MAT
using LinearAlgebra

const Lx = 1.0
const Ly = 1.0
const nely = 30
const volfrac = 1/3
const Uin = 1e0
const rho = 1e0
const mu = 1e0
const conit = 50

domain = TopflowDomain(Lx, Ly, nely)
bkman = BrinkmanPenalizationParameters(mu)
cont = SimpleTopOpt.TopflowContinuation(bkman, volfrac, conit)
fea = SimpleTopOpt.TopflowFEA(domain)
optimizer = OptimalityCriteria()
bc = SimpleTopOpt.DoublePipeBC(domain, fea, Uin)
problem_container = DoublePipeProblem(domain, volfrac, optimizer)

### TODO -- the next segment should not be copy-pasted
EN = Diagonal(I, fea.doftot)
ND = copy(EN)
ND[bc.fixedDofs, bc.fixedDofs] .= 0.0
EN -= ND

# Vectors for free DOFs
alldofs = 1:fea.doftot
freedofs = setdiff(alldofs, bc.fixedDofs)

### Initialization
# Solution vector
S = zeros(fea.doftot)
L = copy(S)
S[bc.fixedDofs] = bc.DIR[bc.fixedDofs]

# Design Field
xPhys = volfrac * ones(domain.nely, domain.nelx)

# Counters
loop = 0
loopcont = 0
nlittot = 0
chcnt = 0

# Change
change = Inf
objOld = Inf

# Continuation
qastep = 1
qa = cont.qavec[1]

# Vectorized constants 
dxv = domain.dx * ones(1, fea.neltot)
dyv = domain.dy * ones(1, fea.neltot)
muv = problem_container.physics.mu * ones(1, fea.neltot)
rhov = problem_container.physics.rho * ones(1, fea.neltot)

alpha =
    bkman.alphamin .+
    (bkman.alphamax - bkman.alphamin) * (1 .- xPhys[:]) ./ (1 .+ qa * xPhys[:])
dalpha =
    (qa * (bkman.alphamax - bkman.alphamin) * (xPhys[:] .- 1)) ./
    (xPhys[:] * qa .+ 1) .^ 2 .-
    (bkman.alphamax - bkman.alphamin) ./ (xPhys[:] * qa .+ 1)

alpha_T = collect(alpha')
dalpha_T = collect(dalpha')


@testset "Double checking initializations" begin

    @test size(alpha_T) == (1, 900)
    @test size(dalpha_T) == (1, 900)
    @test norm(alpha_T - (250.0 * ones(1,900))) ≈ 0 atol=1e-10
    @test sum(abs.(alpha_T - (250.0 * ones(1,900)))) ≈ 0 atol=1e-10

    @test size(xPhys) == (30, 30)
    @test norm(xPhys - (1/3) * ones(30,30)) ≈ 0 atol=1e-10
    @test sum(abs.(xPhys - (1/3) * ones(30,30))) ≈ 0 atol=1e-10

    @test size(S) == (2883,)
    vars = matread("mat_cases/topflow_unit_tests/DoublePipeBC/S_standard.mat")
    t_S = vars["S"]

    @test norm(t_S - S) / 2883 ≈ 0 atol=1e-10
    @test sum(abs.(t_S - S)) / 2883 ≈ 0 atol=1e-10
end

S = newton(
    nlittot,
    dxv,
    dyv,
    muv, 
    rhov,
    alpha_T,
    S, 
    problem_container.fea, 
    problem_container.bc, 
    problem_container.solver_opts, 
    ND, 
    EN
)

@testset "Newton" begin
    vars = matread("mat_cases/topflow_unit_tests/DoublePipeBC/S_postNewton_standard.mat")
    t_S = vars["S"]

    @test size(S) == size(t_S)
    @test norm(t_S - S) / 2883 ≈ 0 atol=1e-7
    @test sum(abs.(t_S - S)) / 2883 ≈ 0 atol=1e-7
end

obj = sum(SimpleTopOpt.TopFlow.call_PHI(dxv, dyv, muv, alpha_T, S[fea.edofMat']))
change = abs(objOld - obj) / objOld
objOld = obj

L = compute_adjoint_solution(dxv, dyv, muv, rhov, alpha_T, S, fea, bc, ND, EN)

@testset "Adjoint solution" begin
    @test size(L) == (2883,1)

    vars = matread("mat_cases/topflow_unit_tests/DoublePipeBC/L_standard.mat")
    t_L = vars["L"]

    @test size(t_L) == size(L)
    @test norm(t_L - L) / 2883 ≈ 0 atol=1e-8
    @test sum(abs.(t_L - L)) / 2883 ≈ 0 atol=1e-8
end

sens, dV = compute_sensitivities(dxv, dyv, muv, rhov, alpha_T, dalpha_T, S, fea, domain, L)

@testset "Sensitivity computation" begin
    t_dV = ones(domain.nely, domain.nelx) / fea.neltot

    vars = matread("mat_cases/topflow_unit_tests/DoublePipeBC/sens_standard.mat")
    t_sens = vars["sens"]

    @test size(dV) == (30,30)
    @test norm(t_dV - dV) / 900 ≈ 0 atol=1e-10
    @test sum(abs.(t_dV - dV)) / 900 ≈ 0 atol=1e-10

    @test size(sens) == (30,30)
    @test norm(t_sens - sens) / 900 ≈ 0 atol=1e-10
    @test sum(abs.(t_sens - sens)) / 900 ≈ 0 atol=1e-10
end

xPhys = OCUpdate(xPhys, sens, dV, problem_container)

@testset "Optimality Criterion update" begin
    vars = matread("mat_cases/topflow_unit_tests/DoublePipeBC/xPhys_post_OCUpdate_standard.mat")
    t_xPhys = vars["xPhys"]

    @test size(xPhys) == (30,30)
    @test norm(t_xPhys - xPhys) / 900 ≈ 0 atol=1e-10
    @test sum(abs.(t_xPhys - xPhys)) / 900 ≈ 0 atol=1e-10
end

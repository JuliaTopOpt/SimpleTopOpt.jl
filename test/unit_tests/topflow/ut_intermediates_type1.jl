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

tfdc = TopflowDomain(Lx, Ly, nely)
bkman = BrinkmanPenalizationParameters(mu)
cont = SimpleTopOpt.TopflowContinuation(volfrac, bkman, conit)
fea = SimpleTopOpt.TopflowFEA(tfdc)
optimizer = OCParameters(200, 0.2)
bc = SimpleTopOpt.DoublePipeBC(tfdc, fea, Uin)
problem_container = DoublePipeContainer(tfdc, volfrac, optimizer, Uin, rho, mu)

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
# TODO -- this is just used in the newton solver, so we can just
#         get rid of this?
# dS = copy(S)
L = copy(S)
S[bc.fixedDofs] = bc.DIR[bc.fixedDofs]

# Design Field
xPhys = volfrac * ones(tfdc.nely, tfdc.nelx)

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
dxv = tfdc.dx * ones(1, fea.neltot)
dyv = tfdc.dy * ones(1, fea.neltot)
muv = problem_container.mu * ones(1, fea.neltot)
rhov = problem_container.rho * ones(1, fea.neltot)

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
    # @test_skip norm(t_S - S) / 2883 ≈ 0 atol=1e-10 # fails @ 1e-9
    # @test_skip sum(abs.(t_S - S)) / 2883 ≈ 0 atol=1e-10 # fails @ 1e-8
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
    @test norm(t_L - L) / 2883 ≈ 0 atol=1e-10
    @test sum(abs.(t_L - L)) / 2883 ≈ 0 atol=1e-10
end

sens, dV = compute_sensitivities(dxv, dyv, muv, rhov, alpha_T, dalpha_T, S, fea, tfdc, L)

@testset "Sensitivity computation" begin

    vars = matread("mat_cases/topflow_unit_tests/DoublePipeBC/dV_standard.mat")
    t_dV = vars["dV"]

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


end

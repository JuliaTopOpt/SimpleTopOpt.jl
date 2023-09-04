# SimpleTopOpt

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliatopopt.github.io/SimpleTopOpt.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliatopopt.github.io/SimpleTopOpt.jl/dev/)
[![Build Status](https://github.com/juliatopopt/SimpleTopOpt.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/juliatopopt/SimpleTopOpt.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/juliatopopt/SimpleTopOpt.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliatopopt/SimpleTopOpt.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)


A small package implementing basic techniques for classic elastic compliance, thermal compliance,
and fluidic topology optimization problems, following 
[*Efficient topology Optimization in MATLAB using 88 lines of code*](https://link.springer.com/article/10.1007/s00158-010-0594-7), 
the `toph` code listing from [*Topology Optimization*](https://www.amazon.com/Topology-Optimization-Martin-Philip-Bendsoe/dp/3540429921),
and [*A detailed introduction to density-based topology optimization of fluid flow problems with implementation in MATLAB*](https://arxiv.org/abs/2207.13695).

## Features and Usage

`SimpleTopOpt.jl` implements three topology optimization algorithms following
popular educational papers for solutions to problems in structural compliance,
thermal compliance, and fluid flows. The package is primarily accessed through
the `optimize` function and problem containers `Top88Problem`, `TophProblem`,
`DoublePipeProblem`, and `PipeBendProblem`, the latter two of which are both fluid flow problems.


#### Top88

This solves the classical MBB beam problem,
where specifically we seek to minimize the compliance of a structure within a rectangular domain, on which a load is placed
on the top left corner; the left side and bottom right corner are fixed.

We expect a structure composed of triangles eventually connecting the two fixed components.

```julia
domain_1 = Top88Domain(60, 40)

SIMP = ModifiedSIMPParameters(penal=3.0)
optimizer = OptimalityCritera()
filter = SensitivityFilter()

problem = Top88Problem(domain, SIMP, optimizer, filter)

sol = SimpleTopOpt.Top88.optimize(problem)

heatmap(sol.design)
```

#### TopH

This solves a classical thermal compliance
problem, where specifically we seek to minimize the thermal compliance achievable in a unit square domain wherein the top
middle third is removed, allowing heat to
escape, while the domain is constantly and
evenly heated.

We expect a symmetric structure branching out like tree roots as
to maximize the surface area.

```julia
using SimpleTopOpt
using Plots

domain = TophDomain(40, 40)

SIMP = ModifiedSIMPParameters(penal=3.0)
optimizer = OptimalityCriteria()
sensitivity_filter = SensitivityFilter()

problem = TophProblem(domain, SIMP, optimizer, sensitivity_filter)

sol = SimpleTopOpt.TopH.optimize(problem)

heatmap(sol.design)
```

#### Topflow

*Warning:* Topflow currently requires `MATLAB.jl`, which itself requires a
working MATLAB installation before 2022a. See 
[the `MATLAB.jl` page](https://github.com/JuliaInterop/MATLAB.jl) for
more detail.

This implements two fluidic topology optimization problems.

In the first problem, the double pipe problem, we seek to minimize the dissipated energy (with the Brinkman penalty term) over a
rectangular domain with two pairs of cut-outs immediately facing each other. So, 
we expect to see two pipes, each connecting a pair.

In the second problem, the pipe bend problem, we seek to minimize the dissipated
energy in a rectangular domain with a cutout
on the left and bottom sides. So, we expect to see a sort of pipe form connecting immediately these cutouts.

```julia
using SimpleTopOpt
using Plots    

Lx = 1.0; Ly = 1.0; nely = 30
volfrac = 1/3 ; Uin = 1e0; rho = 1e0
mu = 1e0; conit = 50

domain = TopflowDomain(Lx, Ly, nely)
fea = SimpleTopOpt.TopflowFEA(domain)
optimizer = OptimalityCriteria()
pbbc = SimpleTopOpt.PipeBendBC(domain, fea, Uin)
pbp = PipeBendProblem(domain, volfrac, optimizer)

sol = SimpleTopOpt.TopFlow.optimize(pbp)

heatmap(sol.design)
```

## Acknowledgements

Thanks to Mohamed Tarek and Yuqing Zhou for their advice over
summer 2023 as the initial development was underway.


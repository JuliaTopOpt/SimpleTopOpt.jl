# SimpleTopOpt

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliatopopt.github.io/SimpleTopOpt.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliatopopt.github.io/SimpleTopOpt.jl/dev/)
[![Build Status](https://github.com/juliatopopt/SimpleTopOpt.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/juliatopopt/SimpleTopOpt.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/juliatopopt/SimpleTopOpt.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliatopopt/SimpleTopOpt.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)


A small package implementing basic techniques for topology optimization problems.

Current implementations-

1. Classic elastic compliance from [*Efficient topology Optimization in MATLAB using 88 lines of code*](https://link.springer.com/article/10.1007/s00158-010-0594-7).
2. The thermal compliance `toph` code listing from [*Topology Optimization*](https://www.amazon.com/Topology-Optimization-Martin-Philip-Bendsoe/dp/3540429921).
3. Fluidic topology optimization problem from [*A detailed introduction to density-based topology optimization of fluid flow problems with implementation in MATLAB*](https://arxiv.org/abs/2207.13695).
4. 3-D topology optimization from [*An efficient 3D topology optimization code written in Matlab*](https://link.springer.com/article/10.1007/s00158-014-1107-x).

## Features and Usage

`SimpleTopOpt.jl` implements four topology optimization algorithms following
popular educational papers for solutions to problems in structural compliance,
thermal compliance, fluid flows and 3-D topology optimization. The package is primarily accessed through
the `optimize` function and problem containers `Top88Problem`, `TophProblem`,
`DoublePipeProblem`, and `PipeBendProblem`, the latter two of which are both fluid flow problems. 

The 3-D topology optimization implementation `top3d` solves the minimum compliance problem
for the cantilevered beam. It can be modified to solve problems with different boundary conditions,
multiple load cases, active and passive elements, compliant mechanism synthesis, etc. Follow the guidelines
from the corresponding educational paper if you desire to do so. Moreover, the implementation also compares
the time taken to execute different sections of the code with corresponding sections in the original MATLAB code.


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

#### Top3d
This code implements the minimum compliance 3-D topology optimization problem for a cantilevered beam as shown in the following figure-

![Alt text](https://i.postimg.cc/w3N1KKN5/cantilever-beam.png "beam")
<!-- <img src="https://drive.google.com/file/d/1UwWkSvq8HkAYCqM042Jtrhq-Q82qv_94/view?usp=drive_link" alt="Alt text" title="Optional title"> -->

The result is displayed as an interactive visualization in a popup window. 

We have also comapared the code's performance with the original `MATLAB` code and the results are displayed upon execution respectively. If one is more interested in seeing the optimized resultant shape than the time taken by various sections of the code, one should uncomment `250`th line and should comment out `251`st line. Then, upon execution, the popup window will stay on the screen for one to observe the resultant optimized shape closely. The shape should appear as follows-

![Alt text](https://i.postimg.cc/34sWy1sD/shape.png "shape")

The times taken by various sections of the code will still be displayed once you exit from the `Makie` popup window but that will not show accurate results since the time one took to observe the shape on the popup window will also be added to the display_3D section of the code thus increasing the time taken for program execution.

To run the code, simply go to the directory `SimpleTopOpt.jl/src/top_opt_3d/` and run the command `julia top3d.jl` on the terminal.

Also, we have also included `MATLAB` code along with our code for reference. It has been modified to get execution time for various sections of the code.


## Acknowledgements

Thanks to Mohamed Tarek and Yuqing Zhou for their advice over
summer 2023 as the initial development was underway.

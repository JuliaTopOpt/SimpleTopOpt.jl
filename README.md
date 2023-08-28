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

Specifically, we are interested in the numerical minimization of the structural compliance achievable in a rectangular domain.

$$
\begin{align*}
    \min_\mathbf{x} : & c(\mathbf{x}) = \mathbf{U}\mathbf{K}\mathbf{U} = \sum^N_{e=1} E_e(x_e) \mathbf{u}_e^T \mathbf{k}_0 \mathbf{u}_e\\
    \text{subject to } : & \mathbf{K}\mathbf{U} = \mathbf{F}\\
    \mathbf{0} \leq \mathbf{x} \leq \mathbf{1}
\end{align*}
$$

where $\mathbf{x}$ denotes the design field, 

#### TopH

Specifically, we are interested in a numerical minimization 

#### Topflow

## Installation

<> (TODO -- list it on julia package repository and write this?)

## Usage

`SimpleTopOpt.jl` is primarily accessed

<> (TODO -- finish this once the API is finalized into the structures and problem containers as discussed)

## References and Citations

<> (TODO -- fill this out and cite this instead of the above links?)



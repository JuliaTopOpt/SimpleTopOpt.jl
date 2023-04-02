using MAT
using SimpleTopOpt
using Test   
using LinearAlgebra


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
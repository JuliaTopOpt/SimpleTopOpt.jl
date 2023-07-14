using BenchmarkTools
using SimpleTopOpt.Top88
io = IOContext(stdout)

# MATLAB R2022B on 2019 Macbook Pro (Intel i9)
#  4x4   : 0.0009 ≈ 9.0203e-04
#  10x10 : 0.0096...
#  16x16 : 0.0341...
# Julia 1.8.5
#  4x4   : 0.0043
#  10x10 : 0.0288
#  16x16 : 0.1930

# 4x4
println("\n -- Benchmarking for 4x4")
b1 = @benchmark top88(4, 4, 0.4, 3.0, 2.0, true, false);
show(io, MIME("text/plain"), b1)

# 10x10
println("\n -- Benchmarking for 10x10")
b2 = @benchmark top88(10, 10, 0.4, 3.0, 2.0, true, false)
show(io, MIME("text/plain"), b2)

# 16x16
println("\n -- Benchmarking for 16x16")
b3 = @benchmark top88(16, 16, 0.4, 3.0, 2.0, true, false)
show(io, MIME("text/plain"), b3)

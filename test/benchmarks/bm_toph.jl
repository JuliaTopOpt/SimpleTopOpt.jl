using BenchmarkTools
using SimpleTopOpt.TopH
io = IOContext(stdout)

# MATLAB R2022B on 2019 Macbook Pro (Intel i9)
#  20x20 : 0.1474
#  40x40 : 2.54
# Julia 1.8.5
#  20x20 : 0.3187
#  40x40 : 13.962s (!!!) 
#    40x40 benchmark kills itself: over 5801952


# 20x20
println("\n -- Benchmarking for 20x20")
b1 = @benchmark toph(20, 20, 0.4, 3.0, 2.0)
show(io, MIME("text/plain"), b1)

# 40x40
println("\n -- Benchmarking for 40x40")
b2 = @benchmark toph(40, 40, 0.4, 3.0, 2.0)
show(io, MIME("text/plain"), b2)

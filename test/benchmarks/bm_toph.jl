using BenchmarkTools
using SimpleTopOpt

io = IOContext(stdout)

# 20x20
println("\n -- Benchmarking for 20x20")
b1 = @benchmark toph(20, 20, 0.4, 3.0, 2.0)
show(io, MIME("text/plain"), b1)

# 40x40
println("\n -- Benchmarking for 40x40")
b2 = @benchmark toph(40, 40, 0.4, 3.0, 2.0)
show(io, MIME("text/plain"), b2)


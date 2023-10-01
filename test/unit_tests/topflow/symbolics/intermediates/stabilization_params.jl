
#######################################################################################
## NOTE --- here begins the ue tests
@test SymbolicUtils.substitute(
    ue,
    Dict([u1 => 1, u2 => 1, u3 => 1, u4 => 1, u5 => 1, u6 => 1, u7 => 1, u8 => 1]),
) - √2 ≈ 0 atol = 1e-12

#######################################################################################
## NOTE --- here begins the τ tests
@test SymbolicUtils.substitute(
    τ,
    Dict([
        u1 => 1,
        u2 => 1,
        u3 => 1,
        u4 => 1,
        u5 => 1,
        u6 => 1,
        u7 => 1,
        u8 => 1,
        α => 1,
        ρ => 1,
        μ => 1,
        dx => 1,
        dy => 1,
    ]),
) - (1 / √(41)) ≈ 0 atol = 1e-12

@test SymbolicUtils.substitute(
    τ,
    Dict([
        u1 => 1,
        u2 => 2,
        u3 => 3,
        u4 => 4,
        u5 => 5,
        u6 => 6,
        u7 => 7,
        u8 => 8,
        α => 9,
        ρ => 10,
        μ => 11,
        dx => 12,
        dy => 13,
    ]),
) - (3130 / √(13086113)) ≈ 0 atol = 1e-12

@test SymbolicUtils.substitute(
    τ,
    Dict([
        u1 => 0,
        u2 => 1,
        u3 => 0,
        u4 => 1,
        u5 => 0,
        u6 => 1,
        u7 => 0,
        u8 => 1,
        α => 1,
        ρ => 1,
        μ => 1,
        dx => 1,
        dy => 1,
    ]),
) - (sqrt(39)) / 39 ≈ 0 atol = 1e-12

@test SymbolicUtils.substitute(
    τ,
    Dict([
        u1 => 1,
        u2 => 0,
        u3 => 1,
        u4 => 0,
        u5 => 1,
        u6 => 0,
        u7 => 1,
        u8 => 0,
        α => 1,
        ρ => 1,
        μ => 1,
        dx => 1,
        dy => 1,
    ]),
) - (sqrt(39)) / 39 ≈ 0 atol = 1e-12
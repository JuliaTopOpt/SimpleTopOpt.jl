###########################################################################################
# 1
@test norm(
    SymbolicUtils.substitute(
        Ru,
        Dict([
            u1 => 1,
            u2 => 1,
            u3 => 1,
            u4 => 1,
            u5 => 1,
            u6 => 1,
            u7 => 1,
            u8 => 1,
            p1 => 1,
            p2 => 1,
            p3 => 1,
            p4 => 1,
            μ => 1,
            η => 1,
            ξ => 1,
            α => 1,
            ρ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - [
        ((41)^(1 / 2) * (270 * (41)^(1 / 2) - 360)) / 14760,
        ((41)^(1 / 2) * (270 * (41)^(1 / 2) - 360)) / 14760,
        -1 / 4,
        3 / 4,
        -((41)^(1 / 2) * (90 * (41)^(1 / 2) - 360)) / 14760,
        -((41)^(1 / 2) * (90 * (41)^(1 / 2) - 360)) / 14760,
        3 / 4,
        -1 / 4,
    ],
) ≈ 0 atol = 1e-12

###########################################################################################
# 2
@test norm(
    SymbolicUtils.substitute(
        Ru,
        Dict([
            u1 => 0,
            u2 => 0,
            u3 => 1,
            u4 => 1,
            u5 => 1,
            u6 => 1,
            u7 => 0,
            u8 => 0,
            p1 => 1,
            p2 => 1,
            p3 => 1,
            p4 => 1,
            μ => 1,
            η => 1,
            ξ => 1,
            α => 1,
            ρ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - [
        -(38^(1 / 2) * (300 * 38^(1 / 2) + 180)) / 13680,
        (38^(1 / 2) * (60 * 38^(1 / 2) - 180)) / 13680,
        (38^(1 / 2) * (120 * 38^(1 / 2) - 60)) / 13680,
        (38^(1 / 2) * (480 * 38^(1 / 2) - 60)) / 13680,
        (38^(1 / 2) * (480 * 38^(1 / 2) + 300)) / 13680,
        (38^(1 / 2) * (120 * 38^(1 / 2) + 300)) / 13680,
        (38^(1 / 2) * (60 * 38^(1 / 2) - 60)) / 13680,
        -(38^(1 / 2) * (300 * 38^(1 / 2) + 60)) / 13680,
    ],
) ≈ 0 atol = 1e-12

###########################################################################################
# 3
@test norm(
    SymbolicUtils.substitute(
        Ru,
        Dict([
            u1 => 2,
            u2 => 2,
            u3 => 1,
            u4 => 1,
            u5 => 1,
            u6 => 1,
            u7 => 2,
            u8 => 2,
            p1 => 1,
            p2 => 1,
            p3 => 1,
            p4 => 0,
            μ => 0,
            η => 0,
            ξ => 0,
            α => 1,
            ρ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - [
        (10^(1 / 2) * (150 * 10^(1 / 2) - 240)) / 3600,
        (10^(1 / 2) * (120 * 10^(1 / 2) + 360)) / 3600,
        -(10^(1 / 2) * (150 * 10^(1 / 2) + 30)) / 3600,
        (10^(1 / 2) * (150 * 10^(1 / 2) - 60)) / 3600,
        -(10^(1 / 2) * (120 * 10^(1 / 2) - 300)) / 3600,
        -(10^(1 / 2) * (150 * 10^(1 / 2) + 240)) / 3600,
        (10^(1 / 2) * (120 * 10^(1 / 2) - 30)) / 3600,
        -(10^(1 / 2) * (120 * 10^(1 / 2) + 60)) / 3600,
    ],
) ≈ 0 atol = 1e-12

###########################################################################################
# 4
@test norm(
    SymbolicUtils.substitute(
        Ru,
        Dict([
            u1 => 2,
            u2 => 2,
            u3 => 0,
            u4 => 0,
            u5 => 0,
            u6 => 0,
            u7 => 2,
            u8 => 2,
            p1 => 0,
            p2 => 0,
            p3 => 0,
            p4 => 0,
            μ => 0,
            η => 0,
            ξ => 0,
            α => 1,
            ρ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - [
        -(sqrt(5) * (120 * sqrt(5) - 600)) / 1800,
        -(sqrt(5) * (120 * sqrt(5) - 600)) / 1800,
        -(sqrt(5) * (60 * sqrt(5) + 120)) / 1800,
        -(sqrt(5) * (60 * sqrt(5) + 120)) / 1800,
        -(sqrt(5) * (60 * sqrt(5) + 360)) / 1800,
        -(sqrt(5) * (60 * sqrt(5) + 360)) / 1800,
        -(sqrt(5) * (120 * sqrt(5) + 120)) / 1800,
        -(sqrt(5) * (120 * sqrt(5) + 120)) / 1800,
    ],
)/8 ≈ 0 atol = 1e-12

###########################################################################################
# 5
@test norm(
    SymbolicUtils.substitute(
        Ru,
        Dict([
            u1 => 1,
            u2 => 1,
            u3 => 1,
            u4 => 1,
            u5 => 1,
            u6 => 1,
            u7 => 1,
            u8 => 1,
            p1 => 0,
            p2 => 0,
            p3 => 0,
            p4 => 0,
            μ => 0,
            η => 0,
            ξ => 0,
            α => 1,
            ρ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - [
        (5^(1 / 2) * (90 * 5^(1 / 2) - 360)) / 1800,
        (5^(1 / 2) * (90 * 5^(1 / 2) - 360)) / 1800,
        1 / 4,
        1 / 4,
        (5^(1 / 2) * (90 * 5^(1 / 2) + 360)) / 1800,
        (5^(1 / 2) * (90 * 5^(1 / 2) + 360)) / 1800,
        1 / 4,
        1 / 4,
    ],
) / 8 ≈ 0 atol = 1e-12

###########################################################################################
# 6
@test norm(
    SymbolicUtils.substitute(
        Ru,
        Dict([
            u1 => 0,
            u2 => 0,
            u3 => 0,
            u4 => 0,
            u5 => 0,
            u6 => 0,
            u7 => 0,
            u8 => 0,
            p1 => 1,
            p2 => 1,
            p3 => 1,
            p4 => 1,
            μ => 0,
            η => 0,
            ξ => 0,
            α => 1,
            ρ => 1,
            dx => 1,
            dy => 1,
        ]),
    ) - [1 / 2, 1 / 2, -1 / 2, 1 / 2, -1 / 2, -1 / 2, 1 / 2, -1 / 2],
) / 8 ≈ 0 atol = 1e-12


###########################################################################################
# 7
@test norm(
    SymbolicUtils.substitute(
        Ru,
        Dict([
            u1 => 1,
            u2 => 2,
            u3 => 3,
            u4 => 4,
            u5 => 5,
            u6 => 6,
            u7 => 7,
            u8 => 8,
            p1 => 9,
            p2 => 10,
            p3 => 11,
            p4 => 12,
            μ => 13,
            η => 14,
            ξ => 15,
            α => 16,
            ρ => 17,
            dx => 18,
            dy => 19,
        ]),
    ) - [
        (sqrt(7976825)*sqrt(29472388)*((23844516*sqrt(7976825)*sqrt(29472388))/93845 - 7059827304))/3628640410560,
        (sqrt(7976825)*sqrt(29472388)*((74569368*sqrt(7976825)*sqrt(29472388))/319073 - 6485891304))/3628640410560,
        (sqrt(7976825)*sqrt(29472388)*((434057244*sqrt(7976825)*sqrt(29472388))/1595365 + 384725304))/3628640410560,
        (sqrt(7976825)*sqrt(29472388)*((81056904*sqrt(7976825)*sqrt(29472388))/319073 + 339071304))/3628640410560,
        (sqrt(7976825)*sqrt(29472388)*((445944324*sqrt(7976825)*sqrt(29472388))/1595365 + 7131413256))/3628640410560,
        (sqrt(7976825)*sqrt(29472388)*((411211752*sqrt(7976825)*sqrt(29472388))/1595365 + 6582260856))/3628640410560,
        (sqrt(7976825)*sqrt(29472388)*((407820732*sqrt(7976825)*sqrt(29472388))/1595365 - 456311256))/3628640410560,
        (sqrt(7976825)*sqrt(29472388)*((369093432*sqrt(7976825)*sqrt(29472388))/1595365 - 435440856))/3628640410560,
    ],
) / 8 ≈ 0 atol = 1e-12
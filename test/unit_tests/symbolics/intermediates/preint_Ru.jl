###########################################################################################
# 1
@test norm(
    SymbolicUtils.substitute(
        Ru,
        Dict([
            u1 => 2,
            u2 => 2,
            u3 => 2,
            u4 => 2,
            u5 => 2,
            u6 => 2,
            u7 => 2,
            u8 => 2,
            p1 => 2,
            p2 => 2,
            p3 => 2,
            p4 => 2,
            μ => 2,
            η => 2,
            ξ => 2,
            α => 2,
            ρ => 2,
            dx => 2,
            dy => 2,
        ]),
    ) - [
        (8 * sqrt(29)) / 29 + 1 / 2,
        (8 * sqrt(29)) / 29 + 1 / 2,
        -(16 * sqrt(29)) / 29 - 5 / 2,
        -(16 * sqrt(29)) / 29 - 3 / 2,
        (24 * sqrt(29)) / 29 + 15 / 2,
        (24 * sqrt(29)) / 29 + 15 / 2,
        -(16 * sqrt(29)) / 29 - 3 / 2,
        -(16 * sqrt(29)) / 29 - 5 / 2,
    ],
) ≈ 0 atol = 1e-12

###########################################################################################
# 2
# NOTE: Fails. Reaches 1e-8 though
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
        63614347069 / 2736 -
        (15084666695 * sqrt(7976825) * sqrt(8977188)) / 2046798864,
        15946299869 / 684 -
        (45375872405 * sqrt(7976825) * sqrt(8977188)) / 6140396592,
        (16169904317 * sqrt(7976825) * sqrt(8977188)) / 2046798864 -
        72701627557 / 2736,
        (48640353143 * sqrt(7976825) * sqrt(8977188)) / 6140396592 -
        9112228693 / 342,
        27962356445 / 912 - (524906207 * sqrt(7976825) * sqrt(8977188)) / 62024208,
        3504674981 / 114 - (1578959453 * sqrt(7976825) * sqrt(8977188)) / 186072624,
        (5412222403 * sqrt(7976825) * sqrt(8977188)) / 682266288 -
        24467224565 / 912,
        (16280393737 * sqrt(7976825) * sqrt(8977188)) / 2046798864 -
        6133142623 / 228,
    ],
) ≈ 0 atol = 1e-12

###########################################################################################
# 3
# NOTE: Fails. Reaches 1e-7 though
@test norm(
    SymbolicUtils.substitute(
        Ru,
        Dict([
            u1 => 19,
            u2 => 18,
            u3 => 17,
            u4 => 16,
            u5 => 15,
            u6 => 14,
            u7 => 13,
            u8 => 12,
            p1 => 11,
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
        (2210646225 * sqrt(7976825) * sqrt(29472388)) / 407254816 +
        138234168527 / 2736,
        (12143675925 * sqrt(7976825) * sqrt(29472388)) / 2239901488 +
        34515858955 / 684,
        -(2369731999 * sqrt(7976825) * sqrt(29472388)) / 407254816 -
        157981545527 / 2736,
        -(13017576987 * sqrt(7976825) * sqrt(29472388)) / 2239901488 -
        19723290635 / 342,
        (2538508277 * sqrt(7976825) * sqrt(29472388)) / 407254816 +
        60761940895 / 912,
        (13944710601 * sqrt(7976825) * sqrt(29472388)) / 2239901488 +
        7585865179 / 114,
        -(2379422503 * sqrt(7976825) * sqrt(29472388)) / 407254816 -
        53166819895 / 912,
        -(13070809539 * sqrt(7976825) * sqrt(29472388)) / 2239901488 -
        13275302657 / 228,
    ],
) ≈ 0 atol = 1e-12


###########################################################################################
# 4
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
        0,
        0,
        -(sqrt(41)) / 164,
        1 / 4 - sqrt(41) / 164,
        (sqrt(41)) / 82,
        (sqrt(41)) / 82,
        (1 / 4) - (sqrt(41) / 164),
        -sqrt(41) / 164,
    ],
) ≈ 0 atol = 1e-12

###########################################################################################
# 5
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
        0,
        0,
        -38^(1 / 2) / 76 - 1 / 4,
        1 / 4 - 38^(1 / 2) / 76,
        38^(1 / 2) / 38 + 1,
        38^(1 / 2) / 38 + 1 / 2,
        -38^(1 / 2) / 76 - 1 / 4,
        -38^(1 / 2) / 76 - 1 / 4,
    ],
) ≈ 0 atol = 1e-12

###########################################################################################
# 6
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
        3 / 32 - (3 * 10^(1 / 2)) / 160,
        (3 * 10^(1 / 2)) / 160 + 3 / 32,
        -3 / 32,
        3 / 32,
        (3 * 10^(1 / 2)) / 160 - 3 / 32,
        -(3 * 10^(1 / 2)) / 160 - 3 / 32,
        3 / 32,
        -3 / 32,
    ],
) ≈ 0 atol = 1e-12

###########################################################################################
# 7
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
        5^(1 / 2) / 20 - 1 / 16,
        5^(1 / 2) / 20 - 1 / 16,
        -1 / 16,
        -1 / 16,
        -5^(1 / 2) / 20 - 1 / 16,
        -5^(1 / 2) / 20 - 1 / 16,
        -1 / 16,
        -1 / 16,
    ],
) ≈ 0 atol = 1e-12

###########################################################################################
# 8
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
        1 / 16 - 5^(1 / 2) / 20,
        1 / 16 - 5^(1 / 2) / 20,
        1 / 16,
        1 / 16,
        5^(1 / 2) / 20 + 1 / 16,
        5^(1 / 2) / 20 + 1 / 16,
        1 / 16,
        1 / 16,
    ],
) / 8 ≈ 0 atol = 1e-12

###########################################################################################
# 9
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
    ) - [1 / 8, 1 / 8, -1 / 8, 1 / 8, -1 / 8, -1 / 8, 1 / 8, -1 / 8],
) / 8 ≈ 0 atol = 1e-12


###########################################################################################
# 10
@test norm(
    SymbolicUtils.substitute(
        Ru,
        Dict([
            u1 => 1.1,
            u2 => 1.2,
            u3 => 1.3,
            u4 => 1.4,
            u5 => 1.5,
            u6 => 1.6,
            u7 => 1.7,
            u8 => 1.8,
            p1 => 0.7,
            p2 => 0.8,
            p3 => 0.9,
            p4 => 1.0,
            μ => 2,
            η => 3,
            ξ => 4,
            α => 5,
            ρ => 6,
            dx => 7,
            dy => 8,
        ]),
    ) - [
        (246412659 * 9707653^(1 / 2)) / 54362856800 + 188137 / 2800,
        (570735993 * 9707653^(1 / 2)) / 108725713600 + 434999 / 5600,
        -(344355861 * 9707653^(1 / 2)) / 54362856800 - 63379 / 560,
        -(797590047 * 9707653^(1 / 2)) / 108725713600 - 29081 / 224,
        (443853717 * 9707653^(1 / 2)) / 54362856800 + 127423 / 560,
        (1028044959 * 9707653^(1 / 2)) / 108725713600 + 58897 / 224,
        -(69182103 * 9707653^(1 / 2)) / 10872571360 - 378269 / 2800,
        -(160238181 * 9707653^(1 / 2)) / 21745142720 - 881023 / 5600,
    ],
) / 8 ≈ 0 atol = 1e-12
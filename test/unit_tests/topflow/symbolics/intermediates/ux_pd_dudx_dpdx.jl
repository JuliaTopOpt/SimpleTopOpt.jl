@test size(s) == (12,)

#######################################################################################
## NOTE --- here begins the ux tests
@test norm(
    SymbolicUtils.substitute(
        ux,
        Dict([
            u1 => 1,
            u2 => 1,
            u3 => 1,
            u4 => 1,
            u5 => 1,
            u6 => 1,
            u7 => 1,
            u8 => 1,
            ξ => 1,
            η => 1,
        ]),
    ) - [1; 1],
) == 0

@test norm(
    SymbolicUtils.substitute(
        ux,
        Dict([
            u1 => 1,
            u2 => 2,
            u3 => 3,
            u4 => 4,
            u5 => 5,
            u6 => 6,
            u7 => 7,
            u8 => 8,
            ξ => 9,
            η => 10,
        ]),
    ) - [-66; -65],
) == 0

@test norm(
    SymbolicUtils.substitute(
        ux,
        Dict([
            u1 => 2,
            u2 => 4,
            u3 => 6,
            u4 => 8,
            u5 => 10,
            u6 => 12,
            u7 => 14,
            u8 => 16,
            ξ => 18,
            η => 20,
        ]),
    ) - [-632; -630],
) == 0

@test norm(
    SymbolicUtils.substitute(
        ux,
        Dict([
            u1 => 10,
            u2 => 10,
            u3 => 10,
            u4 => 10,
            u5 => 10,
            u6 => 10,
            u7 => 10,
            u8 => 10,
            ξ => 10,
            η => 10,
            dx => 10,
            dy => 10,
        ]),
    ) - [10; 10],
) == 0

#######################################################################################
## NOTE --- here begins the px tests
@test SymbolicUtils.substitute(
    px,
    Dict([p1 => 1, p2 => 1, p3 => 1, p4 => 1, ξ => 1, η => 1]),
) - 1 == 0

@test SymbolicUtils.substitute(
    px,
    Dict([p1 => 1, p2 => 2, p3 => 3, p4 => 4, ξ => 5, η => 6]),
) + 6.5 == 0

@test SymbolicUtils.substitute(
    px,
    Dict([p1 => 1, p2 => 1, p3 => 2, p4 => 3, ξ => 5, η => 8]),
) + 3.5 == 0

#######################################################################################
## NOTE --- here begins the dudx tests
@test norm(
    (SymbolicUtils.substitute(
        dudx,
        Dict([
            u1 => 1,
            u2 => 1,
            u3 => 1,
            u4 => 1,
            u5 => 1,
            u6 => 1,
            u7 => 1,
            u8 => 1,
            ξ => 1,
            η => 1,
            dx => 1,
            dy => 1,
        ]),
    )) - [0 0; 0 0],
) == 0

@test norm(
    (SymbolicUtils.substitute(
        dudx,
        Dict([
            u1 => 1,
            u2 => 2,
            u3 => 3,
            u4 => 4,
            u5 => 5,
            u6 => 6,
            u7 => 7,
            u8 => 8,
            ξ => 9,
            η => 10,
            dx => 11,
            dy => 12,
        ]),
    )) - [-20/11 -7/6; -20/11 -7/6],
) ≈ 0 atol = 1e-12

@test norm(
    (SymbolicUtils.substitute(
        dudx,
        Dict([
            u1 => 1,
            u2 => 1,
            u3 => 2,
            u4 => 3,
            u5 => 5,
            u6 => 8,
            u7 => 13,
            u8 => 21,
            ξ => 34,
            η => 55,
            dx => 89,
            dy => 99,
        ]),
    )) - [-251/89 -97/66; -418/89 -485/198],
) ≈ 0 atol = 1e-12

#######################################################################################
## NOTE --- here begins the dpdx tests
@test norm(
    (SymbolicUtils.substitute(
        dpdx,
        Dict([p1 => 1, p2 => 1, p3 => 1, p4 => 1, ξ => 1, η => 1, dx => 1, dy => 1]),
    )) - [0; 0],
) == 0

@test norm(
    (SymbolicUtils.substitute(
        dpdx,
        Dict([p1 => 1, p2 => 2, p3 => 3, p4 => 4, ξ => 5, η => 6, dx => 7, dy => 8]),
    )) - [-6 / 7; -3 / 8],
) == 0
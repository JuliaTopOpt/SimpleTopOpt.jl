# 3D TOPOLOGY OPITMIZATION CODE BY LIU AND TOVAR (JUL 2013) ORIGINALLY GIVEN IN MATLAB.
# TRANSLATED TO JULIA BY UTKARSH YASHVARDHAN UNDER THE GUIDANCE OF MOHAMED TAREK (2024).
# The authors of the original code gave us permission to translate their code to Julia under MIT license.

using SparseArrays
using Printf
using Colors
using Statistics
using LinearAlgebra
using Makie
using GLMakie
using BenchmarkTools
using TimerOutputs

const to = TimerOutput();

function user_defined_loop_parameters()
    200, 0.01, false;
end

function user_defined_material_properties()
    1, 1e-9, 0.3;
end

function user_defined_load_dofs(nelx, nely, nelz)
    il, jl, kl = meshgrid(nelx, 0, 0:nelz);                 # Coordinates
    loadnid = kl.*(nelx+1).*(nely+1)+il.*(nely+1)+(nely+1 .-jl); # Node IDs
    loaddof = 3*vec(loadnid).-1;                             # DOFs
    loaddof;
end

function user_defined_support_fixed_dofs(nelx, nely, nelz)
    iif, jf, kf = meshgrid(0, 0:nely, 0:nelz);                  # Coordinates
    fixednid = kf.*(nelx+1).*(nely+1).+iif.*(nely+1)+(nely+1 .-jf); # Node IDs (try using broadcast assignment)
    fixeddof = [3*vec(fixednid); 3*vec(fixednid).-1; 3*vec(fixednid).-2]; # DOFs (try using broadcast assignment)
    fixeddof;
end

function prepare_finite_element_analysis(nelx, nely, nelz, loaddof, fixeddof, nu)
    nele = nelx*nely*nelz;
    ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
    F = sparse(loaddof, repeat([1], size(loaddof, 1)), repeat([-1], size(loaddof, 1)), ndof, 1);
    U = zeros(ndof,1);
    freedofs = Int64.(setdiff(1:ndof,fixeddof));
    KE = lk_H8(nu);
    nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
    nodeids = reshape(nodegrd[1:end-1,1:end-1],nely*nelx,1);
    nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
    nodeids = repeat(nodeids, outer=(size(nodeidz, 2), size(nodeidz, 1)))+repeat(reshape(nodeidz, 1, :), outer=size(nodeids));
    edofVec = 3*vec(nodeids).+1;
    edofMat = repeat(edofVec, outer=(1,24)) + repeat([0 1 2 3*nely .+ [3 4 5 0 1 2] -3 -2 -1 3*(nely+1)*(nelx+1).+[0 1 2 3*nely .+ [3 4 5 0 1 2] -3 -2 -1]], outer=(nele, 1));
    iK = Int64.(reshape(kron(edofMat, ones(24, 1))', 24*24*nele));
    jK = Int64.(reshape(kron(edofMat, ones(1, 24))', 24*24*nele));
    nele, F, U, freedofs, KE, edofMat, iK, jK;
end

function prepare_filter(nelx, nely, nelz, nele, rmin)
    iH = ones(Int64(nele*(2*(ceil(rmin)-1)+1)^2));
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;
    for k1 = 1:nelz
        for i1 = 1:nelx
            for j1 = 1:nely
                e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
                for k2 = max(k1-(ceil(rmin)-1), 1):min(k1+(ceil(rmin)-1), nelz)
                    for i2 = max(i1-(ceil(rmin)-1), 1):min(i1+(ceil(rmin)-1), nelx)
                        for j2 = max(j1-(ceil(rmin)-1), 1):min(j1+(ceil(rmin)-1), nely)
                            e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                            k = k+1;
                            if k <= size(iH, 1)
                                iH[k] = e1;
                                jH[k] = e2;
                                sH[k] = max(0, rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                            else
                                append!(iH, e1);
                                append!(jH, e2);
                                append!(sH, max(0, rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2)));
                            end
                        end
                    end
                end
            end
        end
    end
    H = sparse(iH, jH, sH);
    Hs = sum(H, dims=2);
    H, Hs;
end

function initialize_iteration(volfrac, nelx, nely, nelz)
    x = repeat([volfrac], nely, nelx, nelz);
    xPhys = repeat([volfrac], nely, nelx, nelz);
    loop = 0; 
    change = 1;
    x, xPhys, loop, change;
end

function fe_analysis(KE, Emin, xPhys, penal, E0, nele, iK, jK, freedofs, F)
    sK = reshape(vec(KE)*(Emin.+vec(xPhys)'.^penal*(E0-Emin)), 24*24*nele);
    K = sparse(iK, jK, sK);
    K = (K+K')/2;
    KK = cholesky(K[freedofs, freedofs]);
    KK \ F[freedofs];
end

function objective_function_and_sensitivity_analysis(U, KE, edofMat, nelx, nely, nelz, Emin, xPhys, penal, E0)
    ce = reshape(sum((U[edofMat]*KE).*U[edofMat], dims=2), nely, nelx, nelz);
    c = sum((Emin.+xPhys.^penal*(E0-Emin)).*ce);
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely, nelx, nelz);
    c, dc, dv;
end

function filtering_and_modification_of_sensitivities(H, dc, Hs, dv)
    H*(dc[:]./Hs), H*(dv[:]./Hs);
end

function optimality_criteria_update(nelx, nely, nelz, x, dc, dv, xPhys, H, Hs, volfrac, nele)
    l1 = 0; 
    l2 = 1e9; 
    move = 0.2;
    xnew = zeros(nely, nelx, nelz);
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max.(0, max.(x.-move, min.(1, min.(x.+move, x.*(sqrt.(-dc./dv/lmid))))));
        xPhys[:] = (H*xnew[:])./Hs;
        if sum(xPhys) > volfrac*nele
            l1 = lmid; 
        else 
            l2 = lmid; 
        end
    end
    change = findmax(abs.(xnew[:]-x[:]))[1];
    x = xnew;
    change, xPhys, x;
end

function print_results(loop, c, xPhys, change)
    @printf(" It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n", loop, c, mean(vec(xPhys)), change);
end

function optimize(volfrac, nelx, nely, nelz, maxloop, tolx, Emin, penal, E0, nele, KE, iK, jK, freedofs, edofMat, H, Hs, displayflag, F, U)

    # INITIALIZE ITERATION
    @timeit to "initialize_iteration" x, xPhys, loop, change = initialize_iteration(volfrac, nelx, nely, nelz)

    # START ITERATION
    while change > tolx && loop < maxloop
        loop = loop+1;

        # FE-ANALYSIS
        @timeit to "fe_analysis" U[freedofs] = fe_analysis(KE, Emin, xPhys, penal, E0, nele, iK, jK, freedofs, F);

        # OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        @timeit to "objective_function_and_sensitivity_analysis" c, dc, dv = objective_function_and_sensitivity_analysis(U, KE, edofMat, nelx, nely, nelz, Emin, xPhys, penal, E0);

        # FILTERING AND MODIFICATION OF SENSITIVITIES
        @timeit to "filtering_and_modification_of_sensitivities" dc[:], dv[:] = filtering_and_modification_of_sensitivities(H, dc, Hs, dv);

        # OPTIMALITY CRITERIA UPDATE
        @timeit to "optimality_criteria_update" change, xPhys, x = optimality_criteria_update(nelx, nely, nelz, x, dc, dv, xPhys, H, Hs, volfrac, nele);
        
        # PRINT RESULTS
        @timeit to "print_results" print_results(loop, c, xPhys, change);

        # PLOT DENSITIES
            if displayflag
                @timeit to "display_3D" display_3D(xPhys, displayflag);
            end
    end
    xPhys;

end

function top3d(nelx, nely, nelz, volfrac, penal, rmin)

    @timeit to "top3d" begin
        
        # USER-DEFINED LOOP PARAMETERS
        @timeit to "user_defined_loop_parameters" maxloop, tolx, displayflag = user_defined_loop_parameters();

        # USER-DEFINED MATERIAL PROPERTIES
        @timeit to "user_defined_material_properties" E0, Emin, nu = user_defined_material_properties();

        # USER-DEFINED LOAD DOFs
        @timeit to "user_defined_load_dofs" loaddof = user_defined_load_dofs(nelx, nely, nelz);

        # USER-DEFINED SUPPORT FIXED DOFs
        @timeit to "user_defined_support_fixed_dofs" fixeddof = user_defined_support_fixed_dofs(nelx, nely, nelz);

        # PREPARE FINITE ELEMENT ANALYSIS
        @timeit to "prepare_finite_element_analysis" nele, F, U, freedofs, KE, edofMat, iK, jK = prepare_finite_element_analysis(nelx, nely, nelz, loaddof, fixeddof, nu);

        # PREPARE FILTER
        @timeit to "prepare_filter" H, Hs = prepare_filter(nelx, nely, nelz, nele, rmin);

        # OPTIMIZE
        @timeit to "optimize" xPhys = optimize(volfrac, nelx, nely, nelz, maxloop, tolx, Emin, penal, E0, nele, KE, iK, jK, freedofs, edofMat, H, Hs, displayflag, F, U);

        # PLOT
        @timeit to "display_3D" display_3D(xPhys);

    end
end

function meshgrid(x, y, z)
    X = zeros(size(y, 1), size(x, 1), size(z, 1))
    Y = zeros(size(y, 1), size(x, 1), size(z, 1))
    Z = zeros(size(y, 1), size(x, 1), size(z, 1))
    X .= x' .* ones(size(y, 1)) .* ones(1, 1, size(z, 1));
    Y .= ones(size(x, 1))' .* y .* ones(1, 1, size(z, 1));
    Z .= reshape(z, 1, 1, :) .* ones(size(y, 1), size(x, 1));
    X, Y, Z
end

# === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho, displayflag=false)
    nely, nelx, nelz = size(rho);
    hx = 1; hy = 1; hz = 1;            # User-defined unit element size
    scene = Figure();
    flag = false;
    faces = [ 1 2 3; 3 4 1;
    5 6 7; 7 8 5;
    6 2 3; 3 7 6;
    5 1 4; 4 8 5;
    8 7 3; 3 4 8;
    5 6 2; 2 1 5 ]
    for k = 1:nelz
        z = (k-1)*hz;
        for i = 1:nelx
            x = (i-1)*hx;
            for j = 1:nely
                y = nely*hy - (j-1)*hy;
                if (rho[j,i,k] > 0.7)  # User-defined display density threshold
                    vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx; x y-hx z+hx; x+hx y-hx z+hx; x+hx y z+hx];
                    vert[:, 2:3] = vert[:, 3:-1:2]; vert[:, 2, :] = -vert[:, 2, :];
                    if flag == false
                        scene = mesh(vert, faces, color = RGB(0.2+0.8*(1-rho[j,i,k]), 0.2+0.8*(1-rho[j,i,k]), 0.2+0.8*(1-rho[j,i,k])), shading = NoShading; axis=(; show_axis = false))
                        flag = true;
                    else
                        mesh!(vert, faces, color = RGB(0.2+0.8*(1-rho[j,i,k]), 0.2+0.8*(1-rho[j,i,k]), 0.2+0.8*(1-rho[j,i,k])), shading = NoShading)
                    end
                end
            end
        end
    end
    if displayflag
        display(scene);
    else
        # wait(display(scene));
        display(scene);
    end
end

# === GENERATE ELEMENT STIFFNESS MATRIX ===
function lk_H8(nu)
    A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
    k = 1/144*A'*[1; nu];

    K1 = [k[1] k[2] k[2] k[3] k[5] k[5];
    k[2] k[1] k[2] k[4] k[6] k[7];
    k[2] k[2] k[1] k[4] k[7] k[6];
    k[3] k[4] k[4] k[1] k[8] k[8];
    k[5] k[6] k[7] k[8] k[1] k[2];
    k[5] k[7] k[6] k[8] k[2] k[1]];

    K2 = [k[9]  k[8]  k[12] k[6]  k[4]  k[7];
    k[8]  k[9]  k[12] k[5]  k[3]  k[5];
    k[10] k[10] k[13] k[7]  k[4]  k[6];
    k[6]  k[5]  k[11] k[9]  k[2]  k[10];
    k[4]  k[3]  k[5]  k[2]  k[9]  k[12]
    k[11] k[4]  k[6]  k[12] k[10] k[13]];

    K3 = [k[6]  k[7]  k[4]  k[9]  k[12] k[8];
    k[7]  k[6]  k[4]  k[10] k[13] k[10];
    k[5]  k[5]  k[3]  k[8]  k[12] k[9];
    k[9]  k[10] k[2]  k[6]  k[11] k[5];
    k[12] k[13] k[10] k[11] k[6]  k[4];
    k[2]  k[12] k[9]  k[4]  k[5]  k[3]];

    K4 = [k[14] k[11] k[11] k[13] k[10] k[10];
    k[11] k[14] k[11] k[12] k[9]  k[8];
    k[11] k[11] k[14] k[12] k[8]  k[9];
    k[13] k[12] k[12] k[14] k[7]  k[7];
    k[10] k[9]  k[8]  k[7]  k[14] k[11];
    k[10] k[8]  k[9]  k[7]  k[11] k[14]];

    K5 = [k[1] k[2]  k[8]  k[3] k[5]  k[4];
    k[2] k[1]  k[8]  k[4] k[6]  k[11];
    k[8] k[8]  k[1]  k[5] k[11] k[6];
    k[3] k[4]  k[5]  k[1] k[8]  k[2];
    k[5] k[6]  k[11] k[8] k[1]  k[8];
    k[4] k[11] k[6]  k[2] k[8]  k[1]];

    K6 = [k[14] k[11] k[7]  k[13] k[10] k[12];
    k[11] k[14] k[7]  k[12] k[9]  k[2];
    k[7]  k[7]  k[14] k[10] k[2]  k[9];
    k[13] k[12] k[10] k[14] k[7]  k[11];
    k[10] k[9]  k[2]  k[7]  k[14] k[7];
    k[12] k[2]  k[9]  k[11] k[7]  k[14]];

    KE = 1/((nu+1)*(1-2*nu))* [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end

top3d(60, 20, 4, 0.3, 3, 1.5);

# displaying timer outputs of executing above code.
println("───────────────────────────────────────Julia───────────────────────────────────────────");
show(to, allocations = false);
println();

# displaying timer outputs of executing corresponding MATLAB code.
println("───────────────────────────────────────MATLAB───────────────────────────────────────────");
open("MATLAB_results") do f
    while !eof(f)      
       s = readline(f)
       println("$s")
    end
end
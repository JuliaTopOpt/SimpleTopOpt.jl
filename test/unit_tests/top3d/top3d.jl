# 3D TOPOLOGY OPITMIZATION CODE BY LIU AND TOVAR (JUL 2013) ORIGINALLY GIVEN IN MATLAB.
# TRANSLATED TO JULIA BY UTKARSH YASHVARDHAN UNDER THE GUIDANCE OF MOHAMED TAREK (2024).
# The authors of the original code gave us permission to translate their code to Julia under MIT license. Thank you!
using SparseArrays
using Printf
using Colors
using Statistics
using LinearAlgebra
using Makie
using GLMakie
using TimerOutputs
using DelimitedFiles

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
    size_iH = size(iH);
    jH = ones(size_iH);
    sH = zeros(size_iH);
    k = 0;
    for k1 = 1:nelz
        for i1 = 1:nelx
            for j1 = 1:nely
                e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
                for k2 = max(k1-(ceil(rmin)-1), 1):min(k1+(ceil(rmin)-1), nelz)
                    for i2 = max(i1-(ceil(rmin)-1), 1):min(i1+(ceil(rmin)-1), nelx)
                        for j2 = max(j1-(ceil(rmin)-1), 1):min(j1+(ceil(rmin)-1), nely)
                            e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                            val = max(0, rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2))
                            k = k+1;
                            if val != 0
                                if k <= size_iH[1]
                                    iH[k] = e1;
                                    jH[k] = e2;
                                    sH[k] = val;
                                else
                                    append!(iH, e1);
                                    append!(jH, e2);
                                    append!(sH, val);
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    H = sparse(iH, jH, sH);
    Hs = sum(H, dims=2);
    Hs = convert(Matrix{BigFloat},Hs);
    Hs = Float64.(round.(Hs, sigdigits=15));
    H = convert(SparseMatrixCSC{BigFloat, Int64},H);
    H = convert(SparseMatrixCSC{Float64, Int64}, (round.(H, sigdigits=5)));
    H, Hs, iH, jH, sH;
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
    KK \ F[freedofs, :];
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

function optimize(volfrac, nelx, nely, nelz, maxloop, tolx, Emin, penal, E0, nele, KE, iK, jK, freedofs, edofMat, H, Hs, displayflag, F, U, ionm, iom)

    # INITIALIZE ITERATION
    @timeit to "initialize_iteration" x, xPhys, loop, change = initialize_iteration(volfrac, nelx, nely, nelz)
    open("julia_values/initialize_iteration/x.txt", "w") do io
        writedlm(io, x)
    end
    open("julia_values/initialize_iteration/xPhys.txt", "w") do io
        writedlm(io, xPhys)
    end
    open("julia_values/initialize_iteration/loop.txt", "w") do io
        writedlm(io, loop)
    end
    open("julia_values/initialize_iteration/change.txt", "w") do io
        writedlm(io, change)
    end
    m_x = reshape(readdlm("matlab_values/initialize_iteration/x.txt", ',', Float64, '\n'), (nely, nelx, nelz));
    if(m_x != x)
        write(ionm, "initialize_iteration- x\n");
    else
        write(iom, "initialize_iteration- x\n");
    end
    m_xPhys = reshape(readdlm("matlab_values/initialize_iteration/xPhys.txt", ',', Float64, '\n'), (nely, nelx, nelz));
    if(m_xPhys != xPhys)
        write(ionm, "initialize_iteration- xPhys\n");
    else
        write(iom, "initialize_iteration- xPhys\n");
    end
    m_loop = readdlm("matlab_values/initialize_iteration/loop.txt")[1, 1];
    if(m_loop != loop)
        write(ionm, "initialize_iteration- loop\n");
    else
        write(iom, "initialize_iteration- loop\n");
    end
    m_change = readdlm("matlab_values/initialize_iteration/change.txt")[1, 1];
    if(m_change != change)
        write(ionm, "initialize_iteration- change\n");
    else
        write(iom, "initialize_iteration- change\n");
    end

    matlab_iterations = 146;

    # START ITERATION
    while change > tolx && loop < maxloop
        loop = loop+1;

        # FE-ANALYSIS
        @timeit to "fe_analysis" U[freedofs] = fe_analysis(KE, Emin, xPhys, penal, E0, nele, iK, jK, freedofs, F);
        open("julia_values/U/U$loop.txt", "w") do io
            writedlm(io, U)
        end
        if loop <= matlab_iterations
            m_U = readdlm("matlab_values/U/U$loop.txt");
            if(m_U != U)
                write(ionm, "U$loop\n");
            else
                write(iom, "U$loop\n");
            end
        end
        
        # OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        @timeit to "objective_function_and_sensitivity_analysis" c, dc, dv = objective_function_and_sensitivity_analysis(U, KE, edofMat, nelx, nely, nelz, Emin, xPhys, penal, E0);
        open("julia_values/c/c$loop.txt", "w") do io
            writedlm(io, c);
        end
        open("julia_values/dc/dc$loop.txt", "w") do io
            writedlm(io, dc);
        end
        open("julia_values/dv/dv$loop.txt", "w") do io
            writedlm(io, dv, ',');
        end
        if loop <= matlab_iterations
            m_c = readdlm("matlab_values/c/c$loop.txt")[1, 1];
            if m_c != c
                write(ionm, "c$loop\n");
            else
                write(iom, "c$loop\n");
            end
            m_dc = reshape(readdlm("matlab_values/dc/dc$loop.txt", ',', Float64, '\n'), (nely, nelx, nelz));
            if m_dc != dc
                write(ionm, "dc$loop\n");
            else
                write(iom, "dc$loop\n");
            end
            m_dv = reshape(readdlm("matlab_values/dv/dv$loop.txt", ',', Float64, '\n'), (nely, nelx, nelz));
            if m_dv != dv
                write(ionm, "dv$loop\n");
            else
                write(iom, "dv$loop\n");
            end
        end

        # FILTERING AND MODIFICATION OF SENSITIVITIES
        @timeit to "filtering_and_modification_of_sensitivities" dc[:], dv[:] = filtering_and_modification_of_sensitivities(H, dc, Hs, dv);
        open("julia_values/dc(:)/dc$loop.txt", "w") do io
            writedlm(io, dc[:]);
        end        
        open("julia_values/dv(:)/dv$loop.txt", "w") do io
            writedlm(io, dv[:]);
        end
        if loop <= matlab_iterations
            m_dc2 = vec(readdlm("matlab_values/dc(:)/dc$loop.txt"));
            if dc[:] != m_dc2
                write(ionm, "dc[:]$loop\n");
            else
                write(iom, "dc[:]$loop\n");
            end
            m_dv2 = vec(readdlm("matlab_values/dv(:)/dv$loop.txt"));
            if dv[:] != m_dv2
                write(ionm, "dv[:]$loop\n");
            else
                write(iom, "dv[:]$loop\n");
            end
        end
        
        # OPTIMALITY CRITERIA UPDATE
        @timeit to "optimality_criteria_update" change, xPhys, x = optimality_criteria_update(nelx, nely, nelz, x, dc, dv, xPhys, H, Hs, volfrac, nele);
        open("julia_values/change/change$loop.txt", "w") do io
            writedlm(io, change)
        end
        open("julia_values/xPhys/xPhys$loop.txt", "w") do io
            writedlm(io, xPhys)
        end
        open("julia_values/x/x$loop.txt", "w") do io
            writedlm(io, x)
        end
        if loop <= matlab_iterations
            m_change = readdlm("matlab_values/change/change$loop.txt")[1, 1];
            if change != m_change
                write(ionm, "change$loop\n");
            else
                write(iom, "change$loop\n");
            end
            m_xPhys = reshape(readdlm("matlab_values/xPhys/xPhys$loop.txt", ',', Float64, '\n'), (nely, nelx, nelz));
            if xPhys != m_xPhys
                write(ionm, "xPhys$loop\n");
            else
                write(iom, "xPhys$loop\n");
            end
            m_x = reshape(readdlm("matlab_values/x/x$loop.txt", ',', Float64, '\n'), (nely, nelx, nelz));
            if x != m_x
                write(ionm, "x$loop\n");
            else
                write(iom, "x$loop\n");
            end
        end
        
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

    setprecision(BigFloat, 100);
    write("matching.txt", "");
    write("not_matching.txt", "");
    iom = open("matching.txt", "a");
    ionm = open("not_matching.txt", "a");

    @timeit to "top3d" begin
        
        # USER-DEFINED LOOP PARAMETERS
        @timeit to "user_defined_loop_parameters" maxloop, tolx, displayflag = user_defined_loop_parameters();
        open("julia_values/user_defined_loop_parameters.txt", "w") do io
            writedlm(io, [maxloop, tolx, displayflag]);
        end
        m_maxloop, m_tolx, m_displayflag = readdlm("matlab_values/user_defined_loop_parameters.txt");
        if m_maxloop != maxloop
            write(ionm, "maxloop\n");
        else
            write(iom, "maxloop\n");
        end
        if m_tolx != tolx
            write(ionm, "tolx\n");
        else
            write(iom, "tolx\n");
        end
        if m_displayflag != displayflag
            write(ionm, "displayflag\n");
        else
            write(iom, "displayflag\n");
        end

        # USER-DEFINED MATERIAL PROPERTIES
        @timeit to "user_defined_material_properties" E0, Emin, nu = user_defined_material_properties();
        open("julia_values/user_defined_material_properties.txt", "w") do io
            writedlm(io, [E0, Emin, nu]);
        end
        m_E0, m_Emin, m_nu = readdlm("matlab_values/user_defined_material_properties.txt");
        if m_E0 != E0
            write(ionm, "E0\n");
        else
            write(iom, "E0\n");
        end
        if m_Emin != Emin
            write(ionm, "Emin\n");
        else
            write(iom, "Emin\n");
        end
        if m_nu != nu
            write(ionm, "nu\n");
        else
            write(iom, "nu\n");
        end

        # USER-DEFINED LOAD DOFs
        @timeit to "user_defined_load_dofs" loaddof = user_defined_load_dofs(nelx, nely, nelz);
        open("julia_values/loaddof.txt", "w") do io
            writedlm(io, loaddof);
        end
        m_loaddof = vec(readdlm("matlab_values/loaddof.txt", Int));
        if m_loaddof != loaddof
            write(ionm, "loaddof\n");
        else
            write(iom, "loaddof\n");
        end

        # USER-DEFINED SUPPORT FIXED DOFs
        @timeit to "user_defined_support_fixed_dofs" fixeddof = user_defined_support_fixed_dofs(nelx, nely, nelz);
        open("julia_values/fixeddof.txt", "w") do io
            writedlm(io, fixeddof);
        end
        m_fixeddof = vec(readdlm("matlab_values/fixeddof.txt"));
        if fixeddof != m_fixeddof
            write(ionm, "fixeddof\n");
        else
            write(iom, "fixeddof\n");

        end

        # PREPARE FINITE ELEMENT ANALYSIS
        @timeit to "prepare_finite_element_analysis" nele, F, U, freedofs, KE, edofMat, iK, jK = prepare_finite_element_analysis(nelx, nely, nelz, loaddof, fixeddof, nu);
       
        open("julia_values/nele.txt", "w") do io
            writedlm(io, nele);
        end
        (x_F, y_F, v_F) = findnz(F);
        open("julia_values/F.txt", "w") do io
            writedlm(io, [x_F, y_F, v_F]);
        end
        open("julia_values/freedofs.txt", "w") do io
            writedlm(io, freedofs);
        end
        open("julia_values/KE.txt", "w") do io
            writedlm(io, KE);
        end
        open("julia_values/edofMat.txt", "w") do io
            writedlm(io, edofMat);
        end
        open("julia_values/iK.txt", "w") do io
            writedlm(io, iK);
        end
        open("julia_values/jK.txt", "w") do io
            writedlm(io, jK);
        end
       
        m_nele = readdlm("matlab_values/nele.txt");
        if m_nele[1, 1] != nele
            write(ionm, "nele\n");
        else
            write(iom, "nele\n");
        end
        m_FF = readdlm("matlab_values/F.txt");
        if m_FF[1, 1] != size(F, 1) || m_FF[1, 2] != size(F, 2) || m_FF[2:size(m_FF, 1), 1] != x_F || m_FF[2:size(m_FF, 1), 2] != y_F || m_FF[2:size(m_FF, 1), 3] != v_F
            write(ionm, "F\n");
        else
            write(iom, "F\n");
        end
        m_U = readdlm("matlab_values/U.txt");
        if m_U != U
            write(ionm, "U\n");
        else
            write(iom, "U\n");
        end
        m_freedofs = vec(readdlm("matlab_values/freedofs.txt", ',', Int64));
        if m_freedofs != freedofs
            write(ionm, "freedofs\n");
        else
            write(iom, "freedofs\n");
        end
        m_KE = readdlm("matlab_values/KE.txt", ',', Float64, '\n');
        if m_KE != KE
            write(ionm, "KE\n");
        else
            write(iom, "KE\n");
        end
        m_edofMat = readdlm("matlab_values/edofMat.txt", ',', Int64, '\n');
        if m_edofMat != edofMat
            write(ionm, "edofMat\n");
        else
            write(iom, "edofMat\n");
        end
        m_iK = vec(readdlm("matlab_values/iK.txt"));
        if m_iK != iK
            write(ionm, "iK\n");
        else
            write(iom, "iK\n");
        end
        m_jK = vec(readdlm("matlab_values/jK.txt"));
        if m_jK != jK
            write(ionm, "jK\n");
        else
            write(iom, "jK\n");
        end

        # PREPARE FILTER
        @timeit to "prepare_filter" H, Hs, iH, jH, sH = prepare_filter(nelx, nely, nelz, nele, rmin);
        (x_H, y_H, v_H) = findnz(H);
        open("julia_values/H.txt", "w") do io
            writedlm(io, [x_H, y_H, v_H]);
        end
        open("julia_values/Hs.txt", "w") do io
            writedlm(io, Hs);
        end
        m_HH = readdlm("matlab_values/H.txt");
        if m_HH[1, 1] != size(H, 1) || m_HH[1, 2] != size(H, 2) || vec(m_HH[2:size(m_HH, 1), 1]) != x_H || vec(m_HH[2:size(m_HH, 1), 2]) != y_H || vec(m_HH[2:size(m_HH, 1), 3]) != v_H
            write(ionm, "H\n");
        else
            write(iom, "H\n");
        end
        m_Hs = readdlm("matlab_values/Hs.txt");
        if Hs != m_Hs
            write(ionm, "Hs\n");
        else
            write(iom, "Hs\n");
        end

        # OPTIMIZE
        @timeit to "optimize" xPhys = optimize(volfrac, nelx, nely, nelz, maxloop, tolx, Emin, penal, E0, nele, KE, iK, jK, freedofs, edofMat, H, Hs, displayflag, F, U, ionm, iom);
        open("julia_values/xPhys.txt", "w") do io
            writedlm(io, xPhys);
        end
        m_xPhys = readdlm("matlab_values/xPhys.txt", ',', Float64, '\n');
        m_xPhys = reshape(m_xPhys, (nely, nelx, nelz));
        if m_xPhys != xPhys
            write(ionm, "xPhys\n");
        else
            write(iom, "xPhys\n");
        end

        # PLOT
        @timeit to "display_3D" display_3D(xPhys);

        close(ionm);
        close(iom);
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
    A = convert(Matrix{BigFloat},A);
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

    Float64.(round.(KE, sigdigits=15));

end

top3d(60, 20, 4, 0.3, 3, 1.5);

println("The corresponding iterations of MATLAB are as follows-");
open("matlab_iteration_values_displayed.txt") do f
    while !eof(f)      
       s = readline(f)
       println("$s")
    end
end

# displaying timer outputs of executing above code.
println("───────────────────────────────────────Julia───────────────────────────────────────────");
show(to, allocations = false);
println();

# displaying timer outputs of executing corresponding MATLAB code.
println("───────────────────────────────────────MATLAB───────────────────────────────────────────");
open("matlab_results.txt") do f
    while !eof(f)      
       s = readline(f)
       println("$s")
    end
end
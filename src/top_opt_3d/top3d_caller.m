% AN 169 LINE 3D TOPOLOGY OPITMIZATION CODE BY LIU AND TOVAR (JUL 2013).
% MODIFIED BY UTKARSH YASHVARDHAN UNDER THE GUIDANCE OF MOHAMED TAREK (2024).
function top3d_caller()
    profile on;
    top3d(60, 20, 4, 0.3, 3, 1.5);
    profile viewer;
end

function [maxloop, tolx, displayflag] = user_defined_loop_parameters()
    maxloop = 200;    % Maximum number of iterations
    tolx = 0.01;      % Terminarion criterion
    displayflag = 0;  % Display structure flag
end

function [E0, Emin, nu] = user_defined_material_properties()
    E0 = 1;           % Young's modulus of solid material
    Emin = 1e-9;      % Young's modulus of void-like material
    nu = 0.3;         % Poisson's ratio
end

function loaddof = user_defined_load_dofs(nelx, nely, nelz)
    [il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % Coordinates
    loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
    loaddof = 3*loadnid(:) - 1;                             % DOFs
end

function fixeddof = user_defined_support_fixed_dofs(nelx, nely, nelz)
    [iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
    fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
    fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
end

function [nele, F, U, freedofs, KE, edofMat, iK, jK] = prepare_finite_element_analysis(nelx, nely, nelz, loaddof, fixeddof, nu)
    nele = nelx*nely*nelz;
    ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
    F = sparse(loaddof,1,-1,ndof,1);
    U = zeros(ndof,1);
    freedofs = setdiff(1:ndof,fixeddof);
    KE = lk_H8(nu);
    nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
    nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
    nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
    nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
    edofVec = 3*nodeids(:)+1;
    edofMat = repmat(edofVec,1,24)+ ...
        repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
        3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
    iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
    jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
end

function [H, Hs] = prepare_filter(nelx, nely, nelz, nele, rmin)
    iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;
    for k1 = 1:nelz
        for i1 = 1:nelx
            for j1 = 1:nely
                e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
                for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                        for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                            e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                            k = k+1;
                            iH(k) = e1;
                            jH(k) = e2;
                            sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                        end
                    end
                end
            end
        end
    end
    H = sparse(iH,jH,sH);
    Hs = sum(H,2);
end

function [x, xPhys, loop, change] = initialize_iteration(volfrac, nelx, nely, nelz)
    x = repmat(volfrac,[nely,nelx,nelz]);
    xPhys = x;
    loop = 0;
    change = 1;
end

function U = fe_analysis(KE, Emin, xPhys, penal, E0, nele, iK, jK, freedofs, F)
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),24*24*nele,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U = K(freedofs,freedofs)\F(freedofs,:);
end

function [c, dc, dv] = objective_function_and_sensitivity_analysis(U, KE, edofMat, nelx, nely, nelz, Emin, xPhys, penal, E0)
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx,nelz]);
    c = sum(sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce)));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx,nelz);
end

function [dc, dv] = filtering_and_modification_of_sensitivities(H, dc, Hs, dv)
    dc = H*(dc(:)./Hs);  
    dv = H*(dv(:)./Hs);
end

function [change, xPhys, x] = optimality_criteria_update(nelx, nely, nelz, x, dc, dv, xPhys, H, Hs, volfrac, nele)
    l1 = 0; l2 = 1e9; move = 0.2;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
        xPhys(:) = (H*xnew(:))./Hs;
        if sum(xPhys(:)) > volfrac*nele, l1 = lmid; else l2 = lmid; end
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
end

function print_results(loop, c, xPhys, change)
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,mean(xPhys(:)),change);
end

function xPhys = optimize(volfrac, nelx, nely, nelz, maxloop, tolx, Emin, penal, E0, nele, KE, iK, jK, freedofs, edofMat, H, Hs, displayflag, F, U)
    
    % INITIALIZE ITERATION
    tic;
    [x, xPhys, loop, change] = initialize_iteration(volfrac, nelx, nely, nelz);
    toc;
    disp("INITIALIZE ITERATION");

    % START ITERATION
    while change > tolx && loop < maxloop
        loop = loop+1;

        % FE-ANALYSIS
        tic;
        U(freedofs) = fe_analysis(KE, Emin, xPhys, penal, E0, nele, iK, jK, freedofs, F);
        toc;
        disp("FE-ANALYSIS");

        % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        tic;
        [c, dc, dv] = objective_function_and_sensitivity_analysis(U, KE, edofMat, nelx, nely, nelz, Emin, xPhys, penal, E0);
        toc;
        disp("OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS");

        % FILTERING AND MODIFICATION OF SENSITIVITIES
        tic;
        [dc(:), dv(:)] = filtering_and_modification_of_sensitivities(H, dc, Hs, dv);
        toc;
        disp("FILTERING AND MODIFICATION OF SENSITIVITIES");

        % OPTIMALITY CRITERIA UPDATE
        tic;
        [change, xPhys, x] = optimality_criteria_update(nelx, nely, nelz, x, dc, dv, xPhys, H, Hs, volfrac, nele);
        toc;
        disp("OPTIMALITY CRITERIA UPDATE");
        
        % PRINT RESULTS
        tic;
        print_results(loop, c, xPhys, change);
        toc;
        disp("PRINT RESULTS");

        % PLOT DENSITIES
        tic;
        if displayflag, clf; display_3D(xPhys); end %#ok<UNRCH>
        toc;
        disp("PLOT DENSITIES");
    end

end

function top3d(nelx, nely, nelz, volfrac, penal, rmin)

    % USER-DEFINED LOOP PARAMETERS
    tic;
    [maxloop, tolx, displayflag] = user_defined_loop_parameters();
    toc;
    disp("USER-DEFINED LOOP PARAMETERS");

    % USER-DEFINED MATERIAL PROPERTIES
    tic;
    [E0, Emin, nu] = user_defined_material_properties();
    toc;
    disp("USER-DEFINED MATERIAL PROPERTIES");

    % USER-DEFINED LOAD DOFs
    tic;
    loaddof = user_defined_load_dofs(nelx, nely, nelz);
    toc;
    disp("USER-DEFINED LOAD DOFs");

    % USER-DEFINED SUPPORT FIXED DOFs
    tic;
    fixeddof = user_defined_support_fixed_dofs(nelx, nely, nelz);
    toc;
    disp("USER-DEFINED SUPPORT FIXED DOFs");

    % PREPARE FINITE ELEMENT ANALYSIS
    tic;
    [nele, F, U, freedofs, KE, edofMat, iK, jK] = prepare_finite_element_analysis(nelx, nely, nelz, loaddof, fixeddof, nu);
    toc;
    disp("PREPARE FINITE ELEMENT ANALYSIS");

    % PREPARE FILTER
    tic;
    [H, Hs] = prepare_filter(nelx, nely, nelz, nele, rmin);
    toc;
    disp("PREPARE FILTER");

    % OPTIMIZE
    tic;
    xPhys = optimize(volfrac, nelx, nely, nelz, maxloop, tolx, Emin, penal, E0, nele, KE, iK, jK, freedofs, edofMat, H, Hs, displayflag, F, U);
    toc;
    disp("OPTIMIZE");

    % PLOT
    tic;
    clf; display_3D(xPhys);
    toc;
    disp("PLOT");

end

% === GENERATE ELEMENT STIFFNESS MATRIX ===
function [KE] = lk_H8(nu)
    A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
        -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
    k = 1/144*A'*[1; nu];

    K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
        k(2) k(1) k(2) k(4) k(6) k(7);
        k(2) k(2) k(1) k(4) k(7) k(6);
        k(3) k(4) k(4) k(1) k(8) k(8);
        k(5) k(6) k(7) k(8) k(1) k(2);
        k(5) k(7) k(6) k(8) k(2) k(1)];
    K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
        k(8)  k(9)  k(12) k(5)  k(3)  k(5);
        k(10) k(10) k(13) k(7)  k(4)  k(6);
        k(6)  k(5)  k(11) k(9)  k(2)  k(10);
        k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
        k(11) k(4)  k(6)  k(12) k(10) k(13)];
    K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
        k(7)  k(6)  k(4)  k(10) k(13) k(10);
        k(5)  k(5)  k(3)  k(8)  k(12) k(9);
        k(9)  k(10) k(2)  k(6)  k(11) k(5);
        k(12) k(13) k(10) k(11) k(6)  k(4);
        k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
    K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
        k(11) k(14) k(11) k(12) k(9)  k(8);
        k(11) k(11) k(14) k(12) k(8)  k(9);
        k(13) k(12) k(12) k(14) k(7)  k(7);
        k(10) k(9)  k(8)  k(7)  k(14) k(11);
        k(10) k(8)  k(9)  k(7)  k(11) k(14)];
    K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
        k(2) k(1)  k(8)  k(4) k(6)  k(11);
        k(8) k(8)  k(1)  k(5) k(11) k(6);
        k(3) k(4)  k(5)  k(1) k(8)  k(2);
        k(5) k(6)  k(11) k(8) k(1)  k(8);
        k(4) k(11) k(6)  k(2) k(8)  k(1)];
    K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
        k(11) k(14) k(7)  k(12) k(9)  k(2);
        k(7)  k(7)  k(14) k(10) k(2)  k(9);
        k(13) k(12) k(10) k(14) k(7)  k(11);
        k(10) k(9)  k(2)  k(7)  k(14) k(7);
        k(12) k(2)  k(9)  k(11) k(7)  k(14)];
    KE = 1/((nu+1)*(1-2*nu))*...
        [ K1  K2  K3  K4;
        K2'  K5  K6  K3';
        K3' K6  K5' K2';
        K4  K3  K2  K1'];
end

% === DISPLAY 3D TOPOLOGY (ISO-VIEW) ===
function display_3D(rho)
    [nely,nelx,nelz] = size(rho);
    hx = 1; hy = 1; hz = 1;            % User-defined unit element size
    face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
    set(gcf,'Name','ISO display','NumberTitle','off');
    for k = 1:nelz
        z = (k-1)*hz;
        for i = 1:nelx
            x = (i-1)*hx;
            for j = 1:nely
                y = nely*hy - (j-1)*hy;
                if (rho(j,i,k) > 0.7)  % User-defined display density threshold
                    vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                    vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                    patch('Faces',face,'Vertices',vert,'FaceColor',[0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k)),0.2+0.8*(1-rho(j,i,k))]);
                    hold on;
                end
            end
        end
    end
    axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end

% =========================================================================
% === This code was written by K Liu and A Tovar, Dept. of Mechanical   ===
% === Engineering, Indiana University-Purdue University Indianapolis,   ===
% === Indiana, United States of America                                 ===
% === ----------------------------------------------------------------- ===
% === Please send your suggestions and comments to: kailiu@iupui.edu    ===
% === ----------------------------------------------------------------- ===
% === The code is intended for educational purposes, and the details    ===
% === and extensions can be found in the paper:                         ===
% === K. Liu and A. Tovar, "An efficient 3D topology optimization code  ===
% === written in Matlab", Struct Multidisc Optim, 50(6): 1175-1196, 2014, =
% === doi:10.1007/s00158-014-1107-x                                     ===
% === ----------------------------------------------------------------- ===
% === The code as well as an uncorrected version of the paper can be    ===
% === downloaded from the website: http://www.top3dapp.com/             ===
% === ----------------------------------------------------------------- ===
% === Disclaimer:                                                       ===
% === The authors reserves all rights for the program.                  ===
% === The code may be distributed and used for educational purposes.    ===
% === The authors do not guarantee that the code is free from errors, a
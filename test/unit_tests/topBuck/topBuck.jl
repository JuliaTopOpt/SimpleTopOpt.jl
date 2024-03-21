# TOPOLOGY OPTIMIZATION WITH LINEARIZED BUCKLING CRITERIA CODE BY LIU AND TOVAR (APRIL 2021) ORIGINALLY GIVEN IN MATLAB.
# TRANSLATED TO JULIA BY UTKARSH YASHVARDHAN UNDER THE GUIDANCE OF MOHAMED TAREK (2024).
# The authors of the original code gave us permission to translate their code to Julia under MIT license.
using LinearAlgebra
using Match
using Printf
using Images
using Makie
using GLMakie
using SparseArrays
using Statistics
using Roots
import Plots
using Colors
using DelimitedFiles
Plots.default(show = true);

function topBuck250_finalSAMO(nelx,nely,penalK,rmin,ft,ftBC,eta,beta,ocPar,maxit,Lx,penalG,nEig,pAgg,prSel,x0)

    # setprecision(BigFloat, 100);
    write("matching.txt", "");
    write("not_matching.txt", "");
    iom = open("matching.txt", "a");
    ionm = open("not_matching.txt", "a");

    # ---------------------------- PRE. 1) MATERIAL AND CONTINUATION PARAMETERS
    E0,Emin,nu = 1,1e-6,0.3;                                           # Young's moduli & Poisson's ratio
    penalCntK = [25,1,25,0.25];                                                # continuation scheme on K-penal
    penalCntG = [25,1,25,0.25];                                                # " " on G-penal
    betaCnt  = [400,24,25,2];                                                 # " " on beta
    pAggCnt  = [2e5,1,25,2];                                                  # " " on the KS aggregation factor
    cnt(v, vCn, l) = v+(l>=vCn[1]).*(v<vCn[2]).*(mod(l,vCn[3])==0).*vCn[4];    # function applying continuation
    # initialize volume fraction
    if prSel[1][1] == 'V'
        volfrac = 1.0; 
    else
        volfrac = prSel[2][end]; 
    end   
    open("julia_values/E0.txt", "w") do io
        writedlm(io, E0);
    end
    m_E0 = readdlm("matlab_values/E0.txt")[1,1];
    # println("m_E0=$m_E0");
    # println("E0=$E0");
    if m_E0 != E0
        write(ionm, "E0\n");
    else
        write(iom, "E0\n");
    end

    open("julia_values/Emin.txt", "w") do io
        writedlm(io, Emin);
    end
    m_Emin = readdlm("matlab_values/Emin.txt")[1,1];
    # println("m_Emin=$m_Emin");
    # println("Emin=$Emin");
    if m_Emin != Emin
        write(ionm, "Emin\n");
    else
        write(iom, "Emin\n");
    end

    open("julia_values/nu.txt", "w") do io
        writedlm(io, nu);
    end
    m_nu = readdlm("matlab_values/nu.txt")[1,1];
    # println("m_nu=$m_nu");
    # println("nu=$nu");
    if m_nu != nu
        write(ionm, "nu\n");
    else
        write(iom, "nu\n");
    end

    open("julia_values/penalCntK.txt", "w") do io
        writedlm(io, penalCntK);
    end
    m_penalCntK = readdlm("matlab_values/penalCntK.txt", ',', Float64);
    # println("m_penalCntK=$m_penalCntK");
    # println("penalCntK=$penalCntK");
    if m_penalCntK[:] != penalCntK
        write(ionm, "penalCntK\n");
    else
        write(iom, "penalCntK\n");
    end

    open("julia_values/penalCntG.txt", "w") do io
        writedlm(io, penalCntG);
    end
    m_penalCntG = readdlm("matlab_values/penalCntG.txt", ',', Float64);
    if m_penalCntG[:] != penalCntG
        write(ionm, "penalCntG\n");
    else
        write(iom, "penalCntG\n");
    end

    open("julia_values/betaCnt.txt", "w") do io
        writedlm(io, betaCnt);
    end
    m_betaCnt = readdlm("matlab_values/betaCnt.txt", ',', Float64);
    if m_betaCnt[:] != betaCnt
        write(ionm, "betaCnt\n");
    else
        write(iom, "betaCnt\n");
    end

    open("julia_values/pAggCnt.txt", "w") do io
        writedlm(io, pAggCnt);
    end
    m_pAggCnt = readdlm("matlab_values/pAggCnt.txt", ',', Float64);
    if m_pAggCnt[:] != pAggCnt
        write(ionm, "pAggCnt\n");
    else
        write(iom, "pAggCnt\n");
    end

    open("julia_values/volfrac.txt", "w") do io
        writedlm(io, volfrac);
    end
    m_volfrac = readdlm("matlab_values/volfrac.txt")[1, 1];
    if m_volfrac != volfrac
        write(ionm, "volfrac\n");
    else
        write(iom, "volfrac\n");
    end

    # ----------------------------------------- PRE. 2) DISCRETIZATION FEATURES
    Ly = nely/nelx*Lx;                                                         # recover Ly from aspect ratio
    nEl = nelx*nely;                                                           # number of elements
    elNrs = reshape(1:nEl,nely,nelx);                                          # element numbering
    nodeNrs = reshape(1:(1+nely)*(1+nelx),1+nely,1+nelx);               # node numbering (int32)
    # nodeNrs = Int32.(reshape(1:(1+nely)*(1+nelx),1+nely,1+nelx));               # node numbering (int32)
    cMat = reshape(2*nodeNrs[1:end-1,1:end-1].+1,nEl,1).+[0 1 2*nely.+[2 3 0 1] -2 -1];# connectivity matrix
    # cMat = reshape(2*nodeNrs[1:end-1,1:end-1].+1,nEl,1)+Int32.([0,1,2*nely.+[2,3,0,1],-2,-1]);# connectivity matrix
    nDof = (1+nely)*(1+nelx)*2;                                                # total number of DOFs
    # ---------------------------------------------- elemental stiffness matrix
    c1 = [12;3;-6;-3;-6;-3;0;3;12;3;0;-3;-6;-3;-6;12;-3;0;-3;-6;3;12;3;-6;3;-6;12;3;-6;-3;12;3;0;12;-3;12];
    c2 = [-4;3;-2;9;2;-3;4;-9;-4;-9;4;-3;2;9;-2;-4;-3;4;9;2;3;-4;-9;-2;3;2;-4;3;-2;9;-4;-9;4;-4;-3;-4];
    Ke = 1/(1-nu^2)/24*(c1+nu.*c2);                                            # lower symmetric part of Ke
    Ke0 = zeros(1,64);
    Ke0[reshape(tril(ones(8, 8)).==1, 1, 64)] = Ke';
    # Ke0[tril(ones(8, 8)).==1] = Ke';
    Ke0 = reshape(Ke0,8,8);
    Ke0 = Ke0+Ke0'-diagm(diag(Ke0));                                            # recover full elemental matrix
    sI,sII = [0],[0];
    for j = 1:8      # build assembly indices for the lower symmetric part of K
        sI = hcat(sI,reshape(collect(j:8), 1, :));
        sII = hcat(sII, repeat([j],1,8-j+1));
    end
    sI = sI[:, 2:end];
    sII = sII[:, 2:end];
    # println("cmat=$(size(cMat[:,sI]))");
    # println("cmat=$(size(cMat))");
    iK,jK = cMat[:,sI[:]]',cMat[:,sII[:]]';
    Iar = sortslices([iK[:] jK[:]],dims=2,rev=true);                                     # indices for K assembly
    if any(prSel[1]=='B') # >>>>>>>>>>>>>>>> PERFORM ONLY IF BUCKLING IS ACTIVE #B#
        Cmat0 = [1 nu 0; nu 1 0; 0 0 (1-nu)/2]/(1-nu^2);                         # non-dimensional elasticity matrix
        xiG = sqrt(1/3)*[-1,1]; etaG = xiG; wxi = [1,1]; weta = wxi;           # Gauss nodes and weights
        xe = [-1 -1; 1 -1; 1 1; -1 1].*Lx/nelx/2;                                 # dimensions of the elements
        lMat = zeros(3, 4); lMat[1, 1] = 1; lMat[2, 4] = 1; lMat[3, 2:3] = 1;  # placement matrix
        dN(xi, zi) = 0.25*[zi-1 1-zi 1+zi -1-zi; xi-1 -1-xi 1+xi 1-xi];       # shape funct. logical derivatives
        B0(gradN) = lMat * kron(gradN,eye(2));                               # strain-displacement matrix
        indM,t2ind = [1,3,5,7,16,18,20,27,29,34],[ 2,3,4,6,7,9 ];      # auxiliary set of indices (1)
        iG,jG = iK(indM,:),jK(indM,:);                                 # indexing of unique G coefficients
        IkG = sortslices([iG(:), jG(:)],2,rev=true);                                # indexing G entries (lower half)
        a1,a2=reshape(IkG(:,2),10,nEl)', reshape(IkG(:,1),10,nEl)';    # auxiliary set of indices (2)
        dZdu = zeros(10,8);                                                    # build U-derivative of matrix G
        for ii = 1 : 8                    # loop on the displacement components
            tt = 0; Uvec = zeros(8,1); Uvec[ii,1] = 1;                         # set a single displ. component
            se = Cmat0*B0((dN(0,0)*xe)\dN(0,0))*Uvec;                          # stresses at the element center
            for j = 1 : length(xiG)
                for k = 1 : length(etaG)
                    xi = xiG[j]; zi = etaG[k];                                 # current integration points
                    w = wxi[j]*weta[k]*det(dN(xi,zi)*xe);                      # current integration weight
                    gradN = (dN(xi,zi)*xe)\dN(xi,zi);                          # shape funct. physical derivatives
                    B1 = [kron(gradN,[1,0]); kron(gradN,[0,1])];               # deformation gradient
                    tt = tt+(B1'*kron(eye(2)[se[1] se[3]; se[3] se[2]])*B1)*w; # current contribution to dG/du_i
                end
            end
            dZdu[:,ii] = tt([1,3,5,7,19,21,23,37,39,55])';                     # extract independent coefficients
        end
        dZdu[t2ind,:] = 2*dZdu[t2ind,:];                                       # x2 columns for v-m-v product
        fKS(p,v) = max(v)+log(sum(exp(p*(v-max(v)))))/p;                      # KS aggregation function
        dKS(p,v,dv) = sum(exp(p.*(v-max(v)))'.*dv,2)./sum(exp(p.*(v-max(v))));# derivative of the KS function

        open("julia_values/Cmat0.txt", "w") do io
            writedlm(io, Cmat0);
        end
        m_Cmat0 = readdlm("matlab_values/Cmat0.txt");
        if m_Cmat0 != Cmat0
            write(ionm, "Cmat0\n");
        else
            write(iom, "Cmat0\n");
        end

        open("julia_values/xiG.txt", "w") do io
            writedlm(io, xiG);
        end
        m_xiG = readdlm("matlab_values/xiG.txt");
        if m_xiG != xiG
            write(ionm, "xiG\n");
        else
            write(iom, "xiG\n");
        end

        open("julia_values/xe.txt", "w") do io
            writedlm(io, xe);
        end
        m_xe = readdlm("matlab_values/xe.txt");
        if m_xe != xe
            write(ionm, "xe\n");
        else
            write(iom, "xe\n");
        end

        open("julia_values/lMat.txt", "w") do io
            writedlm(io, lMat);
        end
        m_lMat = readdlm("matlab_values/lMat.txt");
        if m_lMat != lMat
            write(ionm, "lMat\n");
        else
            write(iom, "lMat\n");
        end

        open("julia_values/dN.txt", "w") do io
            writedlm(io, dN);
        end
        m_dN = readdlm("matlab_values/dN.txt");
        if m_dN != dN
            write(ionm, "dN\n");
        else
            write(iom, "dN\n");
        end

        open("julia_values/B0.txt", "w") do io
            writedlm(io, B0);
        end
        m_B0 = readdlm("matlab_values/B0.txt");
        if m_B0 != B0
            write(ionm, "B0\n");
        else
            write(iom, "B0\n");
        end

        open("julia_values/indM.txt", "w") do io
            writedlm(io, indM);
        end
        m_indM = readdlm("matlab_values/indM.txt");
        if m_indM != indM
            write(ionm, "indM\n");
        else
            write(iom, "indM\n");
        end

        open("julia_values/t2ind.txt", "w") do io
            writedlm(io, t2ind);
        end
        m_t2ind = readdlm("matlab_values/t2ind.txt");
        if m_t2ind != t2ind
            write(ionm, "t2ind\n");
        else
            write(iom, "t2ind\n");
        end

        open("julia_values/iG.txt", "w") do io
            writedlm(io, iG);
        end
        m_iG = readdlm("matlab_values/iG.txt");
        if m_iG != iG
            write(ionm, "iG\n");
        else
            write(iom, "iG\n");
        end

        open("julia_values/jG.txt", "w") do io
            writedlm(io, jG);
        end
        m_jG = readdlm("matlab_values/jG.txt");
        if m_jG != jG
            write(ionm, "jG\n");
        else
            write(iom, "jG\n");
        end

        open("julia_values/IkG.txt", "w") do io
            writedlm(io, IkG);
        end
        m_IkG = readdlm("matlab_values/IkG.txt");
        if m_IkG != IkG
            write(ionm, "IkG\n");
        else
            write(iom, "IkG\n");
        end

        open("julia_values/a1.txt", "w") do io
            writedlm(io, a1);
        end
        m_a1 = readdlm("matlab_values/a1.txt");
        if m_a1 != a1
            write(ionm, "a1\n");
        else
            write(iom, "a1\n");
        end

        open("julia_values/a2.txt", "w") do io
            writedlm(io, a2);
        end
        m_a2 = readdlm("matlab_values/a2.txt");
        if m_a2 != a2
            write(ionm, "a2\n");
        else
            write(iom, "a2\n");
        end

        open("julia_values/dZdu.txt", "w") do io
            writedlm(io, dZdu);
        end
        m_dZdu = readdlm("matlab_values/dZdu.txt");
        if m_dZdu != dZdu
            write(ionm, "dZdu\n");
        else
            write(iom, "dZdu\n");
        end
    end # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #B#

    open("julia_values/Ly.txt", "w") do io
        writedlm(io, Ly);
    end
    m_Ly = readdlm("matlab_values/Ly.txt")[1,1];
    if m_Ly != Ly
        write(ionm, "Ly\n");
    else
        write(iom, "Ly\n");
    end

    open("julia_values/nEl.txt", "w") do io
        writedlm(io, nEl);
    end
    m_nEl = readdlm("matlab_values/nEl.txt")[1,1];
    if m_nEl != nEl
        write(ionm, "nEl\n");
    else
        write(iom, "nEl\n");
    end

    open("julia_values/elNrs.txt", "w") do io
        writedlm(io, elNrs);
    end
    m_elNrs = readdlm("matlab_values/elNrs.txt", ',', Int64, '\n');
    if m_elNrs != elNrs
        write(ionm, "elNrs\n");
    else
        write(iom, "elNrs\n");
    end

    open("julia_values/nodeNrs.txt", "w") do io
        writedlm(io, nodeNrs);
    end
    m_nodeNrs = readdlm("matlab_values/nodeNrs.txt", ',', Float64, '\n');
    if m_nodeNrs != nodeNrs
        write(ionm, "nodeNrs\n");
    else
        write(iom, "nodeNrs\n");
    end

    open("julia_values/cMat.txt", "w") do io
        writedlm(io, cMat);
    end
    m_cMat = readdlm("matlab_values/cMat.txt", ',', Float64, '\n');
    if m_cMat != cMat
        write(ionm, "cMat\n");
    else
        write(iom, "cMat\n");
    end

    open("julia_values/nDof.txt", "w") do io
        writedlm(io, nDof);
    end
    m_nDof = readdlm("matlab_values/nDof.txt")[1,1];
    if m_nDof != nDof
        write(ionm, "nDof\n");
    else
        write(iom, "nDof\n");
    end

    open("julia_values/c1.txt", "w") do io
        writedlm(io, c1);
    end
    m_c1 = readdlm("matlab_values/c1.txt");
    # println("m_c1.size=$(size(m_c1))");
    # println("c1.size=$(size(c1))");
    if m_c1[:] != c1
        write(ionm, "c1\n");
    else
        write(iom, "c1\n");
    end

    open("julia_values/c2.txt", "w") do io
        writedlm(io, c2);
    end
    m_c2 = readdlm("matlab_values/c2.txt");
    if m_c2[:] != c2
        write(ionm, "c2\n");
    else
        write(iom, "c2\n");
    end

    open("julia_values/Ke.txt", "w") do io # precision issue only
        writedlm(io, Ke);
    end
    m_Ke = readdlm("matlab_values/Ke.txt");
    if m_Ke[:] != Ke
        write(ionm, "Ke\n");
    else
        write(iom, "Ke\n");
    end

    open("julia_values/Ke0.txt", "w") do io
        writedlm(io, Ke0);
    end
    m_Ke0 = readdlm("matlab_values/Ke0.txt", ',', Float64, '\n');
    if m_Ke0 != Ke0
        write(ionm, "Ke0\n");
    else
        write(iom, "Ke0\n");
    end

    open("julia_values/sI.txt", "w") do io
        writedlm(io, sI);
    end
    m_sI = readdlm("matlab_values/sI.txt", ',', Float64, '\n');
    if m_sI != sI
        write(ionm, "sI\n");
    else
        write(iom, "sI\n");
    end
    
    open("julia_values/sII.txt", "w") do io
        writedlm(io, sII);
    end
    m_sII = readdlm("matlab_values/sII.txt", ',', Float64, '\n');
    if m_sII != sII
        write(ionm, "sII\n");
    else
        write(iom, "sII\n");
    end

    # println("size(sI)=$(size(sI))");
    # println("size(m_sI)=$(size(m_sI))");
    # println("size(sII)=$(size(sII))");
    # println("size(m_sII)=$(size(m_sII))");
    
    open("julia_values/iK.txt", "w") do io
        writedlm(io, iK);
    end
    m_iK = readdlm("matlab_values/iK.txt", ',', Float64, '\n');
    if m_iK != iK
        write(ionm, "iK\n");
    else
        write(iom, "iK\n");
    end

    open("julia_values/jK.txt", "w") do io
        writedlm(io, jK);
    end
    m_jK = readdlm("matlab_values/jK.txt", ',', Float64, '\n');
    if m_jK != jK
        write(ionm, "jK\n");
    else
        write(iom, "jK\n");
    end

    open("julia_values/Iar.txt", "w") do io
        writedlm(io, Iar);
    end
    m_Iar = readdlm("matlab_values/Iar.txt", ',', Int64, '\n');
    if m_Iar != Iar
        write(ionm, "Iar\n");
    else
        write(iom, "Iar\n");
    end
    # println("size(Iar)=$(size(Iar))");
    # println("size(m_Iar)=$(size(m_Iar))");
    # ----------------------------- PRE. 3) LOADS, SUPPORTS AND PASSIVE DOMAINS
    fixed = 1:2*(nely+1);                                                      # restrained DOFs (cantilever)
    ind = Int64.(nely/2+1 .+(-8:8));
    lcDof = 2*nodeNrs[ind, end].-1;                                  # loaded DOFs
    modF = 1e-3/Ly/(length(lcDof)-1);                                          # modulus of the force density
    F = sparse(lcDof,repeat([1], length(lcDof)),-modF,nDof,1);                                       # define load vector
    F[lcDof[1]],F[lcDof[end]] = F[lcDof[1]]/2,F[lcDof[end]]/2;         # consistent load on end nodes
    ind = Int64.(nely/2 .+(-9:10));
    # println("ind=$ind");
    pasS,pasV = elNrs[ind,end-9:end],[];                    # define passive domains
    # println("pasS=$pasS");
    free = setdiff(1:nDof, fixed);                                             # set of free DOFs
    act = setdiff((1:nEl)',union(pasS[:],pasV[:]));                            # set of active design variables
    open("julia_values/fixed.txt", "w") do io
        writedlm(io, fixed);
    end
    m_fixed = readdlm("matlab_values/fixed.txt", ',', Float64, '\n');
    if m_fixed[:] != fixed
        write(ionm, "fixed\n");
    else
        write(iom, "fixed\n");
    end
    # println("fixed.size=$(size(fixed))");
    # println("m_fixed.size=$(size(m_fixed))");

    open("julia_values/lcDof.txt", "w") do io
        writedlm(io, lcDof);
    end
    m_lcDof = readdlm("matlab_values/lcDof.txt");
    if m_lcDof[:] != lcDof
        write(ionm, "lcDof\n");
    else
        write(iom, "lcDof\n");
    end
    # println("size(lcDof)=$(size(lcDof))");
    # println("size(m_lcDof)=$(size(m_lcDof))");

    open("julia_values/modF.txt", "w") do io
        writedlm(io, modF);
    end
    m_modF = readdlm("matlab_values/modF.txt");
    if m_modF[1,1] != modF
        write(ionm, "modF\n");
    else
        write(iom, "modF\n");
    end

    (x_F, y_F, v_F) = findnz(F);
    open("julia_values/F.txt", "w") do io
        writedlm(io, [x_F y_F v_F]);
    end
    m_FF = readdlm("matlab_values/F.txt");
    if m_FF[1, 1] != size(F, 1) || m_FF[1, 2] != size(F, 2) || m_FF[2:size(m_FF, 1), 1] != x_F || m_FF[2:size(m_FF, 1), 2] != y_F || m_FF[2:size(m_FF, 1), 3] != v_F
        write(ionm, "F\n");
    else
        write(iom, "F\n");
    end

    open("julia_values/pasS.txt", "w") do io
        writedlm(io, pasS);
    end
    m_pasS = readdlm("matlab_values/pasS.txt", ',', Int64, '\n');
    if m_pasS != pasS
        write(ionm, "pasS\n");
    else
        write(iom, "pasS\n");
    end

    open("julia_values/pasV.txt", "w") do io
        writedlm(io, pasV);
    end
    if filesize("matlab_values/pasV.txt") != 0
        m_pasV = readdlm("matlab_values/pasV.txt");
    else
        m_pasV = [];
    end
    if m_pasV != pasV
        write(ionm, "pasV\n");
    else
        write(iom, "pasV\n");
    end

    open("julia_values/free.txt", "w") do io
        writedlm(io, free);
    end
    m_free = readdlm("matlab_values/free.txt", ',', Int64);
    if m_free[:] != free
        write(ionm, "free\n");
    else
        write(iom, "free\n");
    end
    # println("m_free=$(size(m_free))");
    # println("free=$(size(free))");

    open("julia_values/act.txt", "w") do io
        writedlm(io, act);
    end
    m_act = readdlm("matlab_values/act.txt");
    if m_act[:] != act
        write(ionm, "act\n");
    else
        write(iom, "act\n");
    end
    # println("size(m_act)=$(size(m_act))");
    # println("size(act)=$(size(act))");
    # ------------------------- PRE. 4) PREPARE FILTER AND PROJECTION OPERATORS
    if ftBC == 'N'
        bcF = "symmetric"; 
    else
        bcF = 0; 
    end                      # select filter BC
    dy,dx = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
    h = max.(0,rmin.-sqrt.(dx.^2+dy.^2));                                         # convolution kernel
    Hs = imfilter(ones(nely,nelx),h,bcF);                                      # matrix of weights
    dHs = Hs;
    prj(v,eta,beta) = (tanh(beta*eta).+tanh.(beta*(v[:].-eta)))./(tanh(beta*eta)+tanh(beta*(1-eta)));                                   # relaxed Heaviside projection
    deta(v,eta,beta) = -beta*csch(beta).*sech(beta*(v[:]-eta)).^2 .* sinh(v[:]*beta).*sinh((1-v[:])*beta);                                  # projection eta-derivative
    dprj(v,eta,beta) = beta*(1 .-tanh.(beta*(v.-eta)).^2)./(tanh(beta*eta)+tanh(beta*(1-eta)));                                   # projection x-derivative

    open("julia_values/bcF.txt", "w") do io
        writedlm(io, bcF);
    end
    m_bcF = readdlm("matlab_values/bcF.txt");
    if m_bcF[1,1] != bcF
        write(ionm, "bcF\n");
    else
        write(iom, "bcF\n");
    end

    open("julia_values/dy.txt", "w") do io
        writedlm(io, dy);
    end
    m_dy = readdlm("matlab_values/dy.txt", ',', Int64, '\n');
    if m_dy != dy
        write(ionm, "dy\n");
    else
        write(iom, "dy\n");
    end

    open("julia_values/dx.txt", "w") do io
        writedlm(io, dx);
    end
    m_dx = readdlm("matlab_values/dx.txt", ',', Int64, '\n');
    if m_dx != dx
        write(ionm, "dx\n");
    else
        write(iom, "dx\n");
    end

    open("julia_values/h.txt", "w") do io
        writedlm(io, h);
    end
    m_h = readdlm("matlab_values/h.txt", ',', Float64, '\n');
    if m_h != h
        write(ionm, "h\n");
    else
        write(iom, "h\n");
    end

    open("julia_values/Hs.txt", "w") do io
        writedlm(io, Hs);
    end
    m_Hs = readdlm("matlab_values/Hs.txt", ',', Float64, '\n');
    if m_Hs != Hs
        write(ionm, "Hs\n");
    else
        write(iom, "Hs\n");
    end

    open("julia_values/dHs.txt", "w") do io
        writedlm(io, dHs);
    end
    m_dHs = readdlm("matlab_values/dHs.txt", ',', Float64, '\n');
    if m_dHs != dHs
        write(ionm, "dHs\n");
    else
        write(iom, "dHs\n");
    end

    # ------------------------ PRE. 5) ALLOCATE AND INITIALIZE OTHER PARAMETERS
    x,dsK,dsG,dmKS,dV = zeros(nEl,1),zeros(nEl,1),zeros(nEl,1),zeros(nEl,1),zeros(nEl,1);                                  # initialize vectors of size nElx1
    phiDKphi,phiDGphi,adj = zeros(nEl,nEig),zeros(nEl,nEig),zeros(nEl,nEig);                           # " " of size nElxnEig
    U = zeros(nDof,1); phi = zeros(nDof,nEig); adjL = phi; adjV = phi;         # " " of size nDofx1 & nDofxnEig
    dV[act,1] .= 1/nEl;                                                         # derivative of volume fraction
    xpOld,loop,restartAs,ch,plotL,plotR,muVec = 0,0,0,1,[],[],[];      # misc array & parameters
    # if nargin > 15
    if x0 != 0
        load(x0); x = xInitial;                                                # initialize design from saved data
    else
        # println("x[act]=$((volfrac*(nEl-length(pasV))-length(pasS))/length(act))");
        # println(volfrac);
        # println(nEl);
        # println(length(pasV));
        # println(size(pasS)[1]);
        # println(size(act)[1]);
        x[act] .= (volfrac*(nEl-length(pasV))-size(pasS)[1])/size(act)[1];        # volume fraction on "active" set
        x[pasS] .= 1;                                                           # set x=1 on "passive solid" set
    end
    xPhys = x; 
    # clear iK jK iG jG dx dy;                                        # initialize xPhys and free memory
    open("julia_values/x.txt", "w") do io
        writedlm(io, x);
    end
    m_x = readdlm("matlab_values/x.txt");
    if m_x != x
        write(ionm, "x\n");
    else
        write(iom, "x\n");
    end
    # println("size(x)=$(size(x))");
    # println("size(m_x)=$(size(m_x))");
    open("julia_values/dsK.txt", "w") do io
        writedlm(io, dsK);
    end
    m_dsK = readdlm("matlab_values/dsK.txt");
    if m_dsK != dsK
        write(ionm, "dsK\n");
    else
        write(iom, "dsK\n");
    end

    open("julia_values/dsG.txt", "w") do io
        writedlm(io, dsG);
    end
    m_dsG = readdlm("matlab_values/dsG.txt");
    if m_dsG != dsG
        write(ionm, "dsG\n");
    else
        write(iom, "dsG\n");
    end

    open("julia_values/dmKS.txt", "w") do io
        writedlm(io, dmKS);
    end
    m_dmKS = readdlm("matlab_values/dmKS.txt");
    if m_dmKS != dmKS
        write(ionm, "dmKS\n");
    else
        write(iom, "dmKS\n");
    end

    open("julia_values/dV.txt", "w") do io
        writedlm(io, dV);
    end
    m_dV = readdlm("matlab_values/dV.txt");
    if m_dV != dV
        write(ionm, "dV\n");
    else
        write(iom, "dV\n");
    end
    # println("size(dV)=$(size(dV))");
    # println("size(m_dV)=$(size(m_dV))");

    open("julia_values/phiDKphi.txt", "w") do io
        writedlm(io, phiDKphi);
    end
    # println("size = $(filesize("matlab_values/phiDKphi.txt"))");
    if filesize("matlab_values/phiDKphi.txt") <= 1
        m_phiDKphi = [];
    else
        m_phiDKphi = readdlm("matlab_values/phiDKphi.txt");
    end
    if m_phiDKphi != phiDKphi[:]
        write(ionm, "phiDKphi\n");
    else
        write(iom, "phiDKphi\n");
    end

    open("julia_values/phiDGphi.txt", "w") do io
        writedlm(io, phiDGphi);
    end
    if filesize("matlab_values/phiDGphi.txt") <= 1
        m_phiDGphi = [];
    else
        m_phiDGphi = readdlm("matlab_values/phiDGphi.txt");
    end
    if m_phiDGphi != phiDGphi[:]
        write(ionm, "phiDGphi\n");
    else
        write(iom, "phiDGphi\n");
    end
    # println("size(m_phiDGphi)=$(size(m_phiDGphi))");
    # println("size(phiDGphi)=$(size(phiDGphi[:]))");
    # println("size(phiDKphi)=$(size(phiDKphi[:]))");
    # println("size(m_phiDKphi)=$(size(m_phiDKphi))");
    open("julia_values/adj.txt", "w") do io
        writedlm(io, adj);
    end
    if filesize("matlab_values/adj.txt") <= 1
        m_adj = [];
    else
        m_adj = readdlm("matlab_values/adj.txt");
    end
    if m_adj != adj[:]
        write(ionm, "adj\n");
    else
        write(iom, "adj\n");
    end

    open("julia_values/U.txt", "w") do io
        writedlm(io, U);
    end
    m_U = readdlm("matlab_values/U.txt");
    if m_U != U
        write(ionm, "U\n");
    else
        write(iom, "U\n");
    end

    open("julia_values/phi.txt", "w") do io
        writedlm(io, phi);
    end
    if filesize("matlab_values/phi.txt") <= 1
        m_phi = [];
    else
        m_phi = readdlm("matlab_values/phi.txt");
    end
    if m_phi != phi[:]
        write(ionm, "phi\n");
    else
        write(iom, "phi\n");
    end

    open("julia_values/adjL.txt", "w") do io
        writedlm(io, adjL);
    end
    if filesize("matlab_values/adjL.txt") <= 1
        m_adjL = [];
    else
        m_adjL = readdlm("matlab_values/adjL.txt");
    end
    if m_adjL != adjL[:]
        write(ionm, "adjL\n");
    else
        write(iom, "adjL\n");
    end

    open("julia_values/adjV.txt", "w") do io
        writedlm(io, adjV);
    end
    if filesize("matlab_values/adjV.txt") <= 1
        m_adjV = [];
    else
        m_adjV = readdlm("matlab_values/adjV.txt");
    end
    if m_adjV != adjV[:]
        write(ionm, "adjV\n");
    else
        write(iom, "adjV\n");
    end

    open("julia_values/xpOld.txt", "w") do io
        writedlm(io, xpOld);
    end
    m_xpOld = readdlm("matlab_values/xpOld.txt");
    if m_xpOld[1,1] != xpOld
        write(ionm, "xpOld\n");
    else
        write(iom, "xpOld\n");
    end

    open("julia_values/loop.txt", "w") do io
        writedlm(io, loop);
    end
    m_loop = readdlm("matlab_values/loop.txt");
    if m_loop[1,1] != loop
        write(ionm, "loop\n");
    else
        write(iom, "loop\n");
    end

    open("julia_values/restartAs.txt", "w") do io
        writedlm(io, restartAs);
    end
    m_restartAs = readdlm("matlab_values/restartAs.txt");
    if m_restartAs[1,1] != restartAs
        write(ionm, "restartAs\n");
    else
        write(iom, "restartAs\n");
    end

    open("julia_values/ch.txt", "w") do io
        writedlm(io, ch);
    end
    m_ch = readdlm("matlab_values/ch.txt");
    if m_ch[1,1] != ch
        write(ionm, "ch\n");
    else
        write(iom, "ch\n");
    end

    open("julia_values/plotL.txt", "w") do io
        writedlm(io, plotL);
    end
    if filesize("matlab_values/plotL.txt") <= 1
        m_plotL = [];
    else
        m_plotL = readdlm("matlab_values/plotL.txt");
    end
    if m_plotL != plotL
        write(ionm, "plotL\n");
    else
        write(iom, "plotL\n");
    end

    open("julia_values/plotR.txt", "w") do io
        writedlm(io, plotR);
    end
    if filesize("matlab_values/plotR.txt") <= 1
        m_plotR = [];
    else
        m_plotR = readdlm("matlab_values/plotR.txt");
    end
    if m_plotR != plotR
        write(ionm, "plotR\n");
    else
        write(iom, "plotR\n");
    end

    open("julia_values/muVec.txt", "w") do io
        writedlm(io, muVec);
    end
    if filesize("matlab_values/muVec.txt") <= 1
        m_muVec = [];
    else
        m_muVec = readdlm("matlab_values/muVec.txt");
    end
    if m_muVec != muVec
        write(ionm, "muVec\n");
    else
        write(iom, "muVec\n");
    end

    open("julia_values/xPhys.txt", "w") do io
        writedlm(io, xPhys);
    end
    m_xPhys = readdlm("matlab_values/xPhys.txt");
    if m_xPhys != xPhys
        write(ionm, "xPhys\n");
    else
        write(iom, "xPhys\n");
    end
    # println("size(xPhys)=$(size(xPhys))");
    # println("size(m_xPhys)=$(size(m_xPhys))");

    ## ________________________________________________ START OPTIMIZATION LOOP
    local dg0::Matrix{Float64}, dg1::Matrix{Float64}, g0::Float64, g1::Matrix{Float64}, v0::Float64, c0::Matrix{Float64}, xOld::Vector{Float64}, xOld1::Vector{Float64}, as::Matrix{Float64}, f;
    while loop < maxit && ch > 1e-6
      loop = loop+1;                                                           # update iteration counter
      # --------------------------------- RL. 1) COMPUTE PHYSICAL DENSITY FIELD
      xTilde = imfilter(reshape(x,nely,nelx),h,bcF)./Hs;                       # compute filtered field
      xPhys[act] = xTilde[act];                                                # modify active elements only
      if ft > 1                                                                # apply projection
        #   println("xPhys=$(xPhys)");
        #   prj([1,2,3], eta, beta);
          temp = prj(xPhys,eta,beta);
          f = (mean(temp).-volfrac)*(ft==3);                     # function (volume of x-projected)
          while abs(f) > 1e-6 && prSel[1][1] ~= 'V'                            # Newton loop for finding opt. eta
              eta = eta-f/mean(deta(xPhys(:),eta,beta));
              f = mean(prj(xPhys,eta,beta))-volfrac;
          end
          dHs = Hs./reshape(dprj(xPhys,eta,beta),nely,nelx);                   # modification of the sensitivity
          xPhys = prj(xPhys,eta,beta);                                         # compute projected field
      end
      ch = maximum(abs.(xPhys.-xpOld)); xpOld = xPhys;
      open("julia_values/xTilde/xTilde$loop.txt", "w") do io
        writedlm(io, xTilde);
      end
      m_xTilde = readdlm("matlab_values/xTilde/xTilde$loop.txt", ',', Float64, '\n');
      if m_xTilde != xTilde
          write(ionm, "xTilde$loop\n");
      else
          write(iom, "xTilde$loop\n");
      end
      # -------------------------- RL. 2) SETUP AND SOLVE EQUILIBRIUM EQUATIONS
      sK = (Emin.+xPhys.^penalK*(E0-Emin));                                     # stiffness interpolation
      dsK[act] = penalK*(E0-Emin)*xPhys[act].^(penalK-1);                      # derivative of " "
      sK = reshape(Ke[:]*sK',length(Ke)*nEl,1);
      open("julia_values/sK/sK$loop.txt", "w") do io
        writedlm(io, sK);
      end
      m_sK = readdlm("matlab_values/sK/sK$loop.txt", ',', Float64, '\n');
      if m_sK != sK
          write(ionm, "sK$loop\n");
      else
          write(iom, "sK$loop\n");
      end
    #   m_ = readdlm("matlab_values/U/U$loop.txt", ',', Float64, '\n');
    #   if m_U != U
    #       write(ionm, "U$loop\n");
    #   else
    #       write(iom, "U$loop\n");
    #   end
      K = sparse(Iar[:,1],Iar[:,2],sK[:],nDof,nDof);                           # assemble stiffness matrix
      (x_K, y_K, v_K) = findnz(K);
      open("julia_values/K'/K'$loop.txt", "w") do io
        writedlm(io, [x_K y_K v_K]);
      end
    #   println("Alright");
    #   println("$(size(K))");
      K = K+K'-Diagonal(diag(K));                                                  # symmetrization of K
        #   dK = decomposition(K[free,free],"chol","lower");                         # decompose K and store factor
        (x_K, y_K, v_K) = findnz(K);
        open("julia_values/K/K$loop.txt", "w") do io
            writedlm(io, [x_K y_K v_K]);
        end
        m_KK = readdlm("matlab_values/K/K$loop.txt");
        if m_KK[1, 1] != size(K, 1) || m_KK[1, 2] != size(K, 2) || m_KK[2:size(m_KK, 1), 1] != x_K || m_KK[2:size(m_KK, 1), 2] != y_K || m_KK[2:size(m_KK, 1), 3] != v_K
            write(ionm, "K$loop\n");
        else
            write(iom, "K$loop\n");
        end
    # dK = PositiveFactorizations.cholesky(Positive, K[free, free], true);
      dK = factorize(K[free, free]);
      U[free] = dK \ F[free];                                                  # solve equilibrium system
      open("julia_values/U/U$loop.txt", "w") do io
        writedlm(io, U);
      end
      m_U = readdlm("matlab_values/U/U$loop.txt", ',', Float64, '\n');
      if m_U != U
          write(ionm, "U$loop\n");
      else
          write(iom, "U$loop\n");
      end
      dc = -dsK.*sum((U[cMat]*Ke0).*U[cMat],dims=2);                                # compute compliance sensitivity
      if any(prSel[1]=='B') # >>>>>>>>>>>>>> PERFORM ONLY IF BUCKLING IS ACTIVE #B#
      # ---------------------------------- RL. 3) BUILD STRESS STIFFNESS MATRIX 
      sGP = (Cmat0*B0((dN(0,0)*xe)\dN(0,0))*U(cMat)')';                        # stresses at elements centroids
      Z = zeros(nEl,10);      # allocate array for compact storage of Ge coeff.
      for j = 1:length(xiG)                       # loop over quadrature points
        for k = 1:length(etaG)
            # ---------------------------- current integration point and weight
            xi = xiG(j); zi = etaG(k); w = wxi(j)*weta(k)*det(dN(xi,zi)*xe);
            # - reduced represenation of strain-displacement matrix (see paper)
            gradN = (dN(xi,zi)*xe)\dN(xi,zi);                                  # shape funct. physical derivatives
            a = gradN(1,:); b = gradN(2,:); B = zeros(3,10);
            l = [1 1; 2 1; 3 1; 4 1; 2 2; 3 2; 4 2; 3 3; 4 3; 4 4];
            for jj = 1:10
                B(:,jj) = [a(l(jj,1))*a(l(jj,2));
                           b(l(jj,1))*b(l(jj,2));
                           b(l(jj,2))*a(l(jj,1))+b(l(jj,1))*a(l(jj,2))]; 
            end
            # ----------- current contribution to (unique ~= 0) elements of keG
            Z = Z+sGP*B*w;
        end
      end
      sG0 = E0*xPhys.^penalG;                                                  # stress interpolation
      dsG(act) = penalG*E0*xPhys(act).^(penalG-1);                             # derivative of " "
      sG = reshape((sG0.*Z)',10*nEl,1);
      G = sparse(IkG(:,1)+1,IkG(:,2)+1,sG,[nDof,nDof])+sparse(IkG(:,1),  IkG(:,2),  sG,[nDof,nDof]);                    # assemble global G matrix
      G = G+G'-diag(diag(G));                                                  # symmetrization of G
      # ------------------------------ RL. 4) SOLVE BUCKLING EIGENVALUE PROBLEM
      matFun(x) = dK\(G(free,free)*x);                                       # matrix action function
      eivecs,D = eigs(matFun,length(free),nEig+4,"sa");                      # compute eigenvalues
      mu,ii = sort(diag(-D),"descend");                                      # sorting of eigenvalues (mu=-D(i))
      eivSort = eivecs(:,ii(1:nEig));                                          # sort eigenvectors accordingly
      phi(free,:) = eivSort./sqrt(diag(eivSort'*K(free,free)*eivSort)');       # orthonormalize (phi'*K*phi=1)
      # ----------------------------------- RL. 5) SENSITIVITY ANALYSIS OF BLFs
      dkeG = dsG.*Z;                                                           # x-derivative of Ge
      dkeG(:,t2ind) = 2*dkeG(:,t2ind);                                         # x2 columns for v-m-v product
      for j = 1:nEig     # loop on the eigenvalues included in the optimization
          # 1) ------ Term due to the elastic stiffness matrix (phi'*dK/dx*phi)
          t = phi(:,j);
          phiDKphi(:,j) = dsK.*sum((t(cMat)*Ke0).*t(cMat),2);
          # 2) -------------- Term due to the geometric matrix (phi'*dG/dx*phi)
          p = t(a1).*t(a2)+t(a1+1).*t(a2+1);
          phiDGphi(:,j) = sum(dkeG.*p,2);
          # 3) ----------------------------------------- Setup of adjoint loads
          tmp = zeros(nDof,1);
          for k = 1:8              # contribution of each term dKg/du_i, i=1:nD
              tmp[cMat[:,k]] = tmp(cMat(:,k))+(sG0.*p)*dZdu(:,k);
          end
          adjL[:,j] = tmp;
      end
      # ----------- solve the adjoint problem and compute the term (U'*dK/dx*V)
      adjV[free,:] = dK \ adjL[free,:];               # use the stored K factor
      for j = 1 : nEig
          vv = adjV(:,j);
          adj[:,j] = dsK.*sum((U[cMat]*Ke0).*vv[cMat],2);
      end 
      # --------------- overall sensitivity expression for the "mu" eigenvalues
      dmu = -(phiDGphi+mu(1:nEig )'.*phiDKphi-adj);
      end # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      # ---------------------- RL. 6) SELECT OBJECTIVE FUNCTION AND CONSTRAINTS
      if loop==1
        c0=F'*U; 
        # println("typeof(c0)=$(typeof(c0))");
        v0=mean(xPhys[:]); 
      end # initial compliance & volume fraction
      @match prSel[1] begin                    # select optimization problem to solve
          ['C','V'] => begin          # minimize compliance with volume constraint
              g0 = F'*U/c0;
              dg0 = imfilter(reshape(dc/c0,nely,nelx)./dHs,h,bcF);
              g1 = mean(xPhys(:))/volfrac-1;
              dg1 = imfilter(reshape(dV/volfrac,nely,nelx)./dHs,h,bcF);
          end
          ['V','C'] => begin          # minimize volume with compliance constraint
          g0 = mean(xPhys[:])./v0;
          dg0 = imfilter(reshape(dV/v0,nely,nelx)./dHs,h,bcF);
          g1 = (F'*U)/(prSel[2]*c0).-1;
          # println(typeof(g1));
          dg1 = imfilter(reshape(dc/(prSel[2]*c0),nely,nelx)./dHs,h,bcF);
        end
        ['B','C','V'] => begin     # maximize BLF with compliance & volume constraints (Eq. 13 paper)
        if loop==1
            muKS0=fKS(pAgg,mu(1:nEig)); 
            g0=1; 
            cMax=prSel{2}(1);
        else
            g0=fKS(pAgg,mu(1:nEig))/muKS0; 
        end                         # KS aggregation of mu (=1/lambda)
        dmKS = dKS(pAgg,mu(1:nEig),dmu);                                 # KS aggregation of dmu
        dg0 = imfilter(reshape(dmKS/muKS0,nely,nelx)./dHs,h,bcF);        # back-filter KS sensitivity
        # -- Constraint function: KS aggregation of compliance and volume
        g1Vec = [F'*U;mean(xPhys(:))]./[cMax*c0;volfrac]-1;              # set of constraints ['C','V']
        dg1c = imfilter(reshape(dc/(cMax*c0),nely,nelx)./dHs,h,bcF);     # back-filter compliance derivative
        dg1V = imfilter(reshape(dV/volfrac,nely,nelx)./dHs,h,bcF);       # back-filter volume derivative
        g1 = fKS(pAgg,g1Vec);                                            # aggregate the two constraints
        dg1 = dKS(pAgg,g1Vec,[dg1c(:),dg1V(:)]);                         # sensitivity of the KS constraint
        plotL[loop,:] = [1/g0/muKS0,1/mu(1)]; 
        strL="KS(-),lambda_1(--)";
        plotR[loop,:] = [g1,g1Vec']; 
        strR="g_1(-),gC(--),gV(.-)";
        muVec[loop,:] = mu';
    end
    ['V','C','B'] => begin   # min volume with compliance & BLF constraints (Eq. 14 paper)
    g0 = mean(xPhys(:))./v0;
    dg0 = imfilter(reshape(dV/volfrac,nely,nelx)./dHs,h,bcF);
    # ---- Constraint function: KS aggregation of BLFs and compliance
    muKS = fKS(pAgg,mu(1:nEig));                                     # KS aggregation of mu
    dmKS = dKS(pAgg,mu(1:nEig),dmu);                                 # KS aggregation of dmu
    g1Vec = [prSel{2}(2)*muKS;F'*U]./[1;prSel{2}(1)*c0]-1;           # set of constraints 'B','C'
    dg1l = imfilter(reshape(dmKS*prSel{2}(2),nely,nelx)./dHs,h,bcF); # back-filter dmu
    dg1c = imfilter(reshape(dc/(prSel{2}(1)*c0),nely,nelx)./dHs,h,bcF);# back-filter dc
    g1 = fKS(pAgg,g1Vec);                                            # aggregate the two constraints
    dg1 = dKS(pAgg,g1Vec,[dg1l(:),dg1c(:)]);                         # sensitivity of the KS constraint
    plotL[loop,:] = g0; strL = "g_0";
    plotR[loop,:] = [g1,g1Vec']; strR="g_1(-),gL(--),gC(.-)";
    muVec = cat(1,muVec,mu');
    end
end
# ---------------------------------------- RL. 7) UPDATE DESIGN VARIABLES
if loop==1
    xOld = x[act]; 
    xOld1 = xOld; 
    # println("typeof(xOld)=$(typeof(xOld))");
    as = Array{Float64}(undef, 0, 0); 
end                    # initialize MMA history parameters
x0,as,lmid=ocUpdate(loop,x[act],dg0[act],g1,dg1[act],ocPar,xOld,xOld1,as,beta,restartAs);

    open("julia_values/x0/x0$loop.txt", "w") do io
        writedlm(io, x0);
    end
    m_x0 = readdlm("matlab_values/x0/x0$loop.txt");
    if m_x0 != x0
        write(ionm, "x0$loop\n");
    else
        write(iom, "x0$loop\n");
    end

    open("julia_values/as/as$loop.txt", "w") do io
    writedlm(io, as);
    end
    m_as = readdlm("matlab_values/as/as$loop.txt", ',', Float64, '\n');
    if m_as != as
        write(ionm, "as$loop\n");
    else
        write(iom, "as$loop\n");
    end

    open("julia_values/lmid/lmid$loop.txt", "w") do io
    writedlm(io, lmid);
    end
    m_lmid = readdlm("matlab_values/lmid/lmid$loop.txt");
    if m_lmid[1,1] != lmid
        write(ionm, "lmid$loop\n");
    else
        write(iom, "lmid$loop\n");
    end

    xOld1 = xOld; xOld = x[act]; x[act] = x0;
    # ----------------------------------------- RL. 8) PRINT AND PLOT RESULTS 
    #   println("g1=$g1");
    # f = Figure();
    @printf("It.:%2i g0:%7.4f g1:%0.2e penalK:%7.2f penalG:%7.2f eta:%7.2f beta:%7.1f ch:%0.3e lm:%0.3e\n",loop,g0,g1[1],penalK,penalG,eta,beta,ch,lmid);
    if any(prSel[1] == 'B')  # plot design, g0 & g1 evolution, BLFs evolution
        ax1 = Axis(f[1, 1:2],
        title = "Current design",);
        image!(ax1, 1-reshape(xPhys,nely,nelx));
        ax2 = Axis(f[2, 1],
        title = "Objective and constraint",
        ylabel = strL
        );
        scatter!(ax2, 1:loop, plotL);
        ax3 = Axis(f[2, 2],
        title = "Lowest BLFs",
        ylabel = strR
        );
        scatter!(ax3, 1:loop, plotR);
        display(f);
    else                                       # plot the current design only
        Plots.plot(Gray.(1 .-reshape(xPhys,nely,nelx)), ticks=false, title="plot");
    end
    #  apply continuation on penalization(s), beta & aggregation parameter(s)
    penalKold = penalK; penalGold = penalG; betaOld = beta;
    penalK,penalG,beta,pAgg = cnt(penalK, penalCntK, loop), cnt(penalG,penalCntG,loop), cnt(beta,betaCnt,loop), cnt(pAgg,pAggCnt,loop);
    if (beta-betaOld != 0 || penalK-penalKold != 0 || penalG-penalGold != 0)
        restartAs = 1; 
    else
        restartAs = 0; 
    end                              # restart asymptotes if needed
    open("julia_values/xPhys/xPhys$loop.txt", "w") do io
        writedlm(io, xPhys);
    end
    m_xPhys = readdlm("matlab_values/xPhys/xPhys$loop.txt");
    if m_xPhys != xPhys
        write(ionm, "xPhys$loop\n");
    else
        write(iom, "xPhys$loop\n");
    end
    end
end
    
function ocUpdate( loop, xT, dg0, g1, dg1, ocPar, xOld, xOld1, as, beta, restartAs )
    # -------------------------------- definition of asymptotes and move limits
    xU,xL = min.(xT.+ocPar[1],1), max.(xT.-ocPar[1],0);
    if (loop<2.5 || restartAs==1)
        # println("xT=$(size(xT))");
        # println("xU-xL=$(size(reshape((xU-xL), length(xU-xL), 1)))");
        as = xT.+reshape([-0.5, 0.5], 1, 2).*reshape((xU-xL), length(xU-xL), 1)./(beta+1);
    else
        tmp = (xT-xOld).*(xOld-xOld1);
        gm = ones(length(xT),1);
        gm[tmp.>0] .= ocPar[3]; gm[tmp.<0] .= ocPar[2];
        # println("size(gm)=$(size(gm))");
        # println("size(-(xOld-as[:,1]))=$(size(-(xOld-as[:,1])))");
        as = xT.+gm.*[-(xOld-as[:,1]) (as[:,2]-xOld)];
    end
    # println("typeof(as)=$(typeof(as))");
    xL = max( 0.9*as[:,1]+0.1*xT,xL);                    # adaptive lower bound
    xU = min( 0.9*as[:,2]+0.1*xT,xU);                    # adaptive upper bound
    # ----- split (+) and (-) parts of the objective and constraint derivatives
    p0_0 = (dg0.>0).*dg0; q0_0 = (dg0.<0).*dg0;
    p1_0 = (dg1.>0).*dg1; q1_0 = (dg1.<0).*dg1;
    p0,q0 = p0_0.*(as[:,2]-xT).^2,-q0_0.*(xT-as[:,1]).^2;
    p1,q1 = p1_0.*(as[:,2]-xT).^2,-q1_0.*(xT-as[:,1]).^2;
    # ---------------------------------------- define the primal projection map
    primalProj(lm) = min(xU,max(xL,(sqrt.(p0+lm*p1).*as[:,1]+sqrt.(q0+lm*q1).*as[:,2])./(sqrt.(p0+lm*p1)+sqrt.(q0+lm*q1))));
    psiDual(lm) = g1 .- ( (as[:,2]-xT)'*p1_0 - (xT-as[:,1])'*q1_0 ) .+ sum(p1./(max.(as[:,2]-primalProj(lm),1e-12)) .+ q1./(max.(primalProj(lm)-as[:,1],1e-12)));
    # ------------- compute the Lagrange multiplier and update design variables
    lmUp = 1e6; x = xT; lmid = -1;
    # println("typeof(psiDual( 0 ))= $(typeof(psiDual( 0 )))");
    # println("typeof(psiDual( lmUp ))= $(typeof(psiDual( lmUp )))");
    # println("size(psiDual( 0 ))= $(size(psiDual( 0 )))");
    # println("size(psiDual( lmUp ))= $(size(psiDual( lmUp )))");
    # println("psiDual( 0 )= $(psiDual( 0 ))");
    # println("psiDual( lmUp )= $(psiDual( lmUp ))");
    if (psiDual( 0 ) * psiDual( lmUp ))[1] < 0  # check if LM is within the interval
        lmid = find_zero( psiDual, (0, lmUp));
        x = primalProj( lmid );
    elseif psiDual( 0 )[1] < 0                       # constraint cannot be active
       lmid = 0;
       x = primalProj( lmid );
    elseif psiDual( lmUp )[1] > 0                 # constraint cannot be fulfilled
       lmid = lmUp;
       x = primalProj( lmid );
    end
    x, as, lmid;
end

function meshgrid(x, y)
    X = zeros(size(y, 1), size(x, 1))
    Y = zeros(size(y, 1), size(x, 1))
    # Z = zeros(size(y, 1), size(x, 1), size(z, 1))
    X .= x' .* ones(size(y, 1));
    Y .= ones(size(x, 1))' .* y;
    # Z .= reshape(z, 1, 1, :) .* ones(size(y, 1), size(x, 1));
    X, Y
end

# topBuck250 (320,320,3,3,2, 'N', 0.5,2, [0.1,0.7,1.2], 500, 1, 3, 12, 160, {['V', 'C', 'B'], [2.5, 1.05]});
topBuck250_finalSAMO(480,240,3,4,2, 'N', 0.5,2, [0.1,0.7,1.2], 300, 2, 0, 0, 0, [['V', 'C'], 2.5],0);
# topBuck250(480,240,3,4,2, 'N', 0.5,2, [0.1,0.7,1.2], 750,2,3,12,200, {['B', 'C', 'V'], [2.5,0.25]}, 'IG. mat');
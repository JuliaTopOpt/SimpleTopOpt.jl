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
    Ke0 = reshape(Ke0,8,8);
    Ke0 = Ke0+Ke0'-diagm(diag(Ke0));                                            # recover full elemental matrix
    sI,sII = [0],[0];
    for j = 1:8      # build assembly indices for the lower symmetric part of K
        sI = hcat(sI,reshape(collect(j:8), 1, :));
        sII = hcat(sII, repeat([j],1,8-j+1));
    end
    sI = sI[:, 2:end];
    sII = sII[:, 2:end];
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

    end # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #B#

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
        x[act] .= (volfrac*(nEl-length(pasV))-size(pasS)[1])/size(act)[1];        # volume fraction on "active" set
        x[pasS] .= 1;                                                           # set x=1 on "passive solid" set
    end
    xPhys = x; 
    # clear iK jK iG jG dx dy;                                        # initialize xPhys and free memory
    
    ## ________________________________________________ START OPTIMIZATION LOOP
    local dg0::Matrix{Float64}, dg1::Matrix{Float64}, g0::Float64, g1::Matrix{Float64}, v0::Float64, c0::Matrix{Float64}, xOld::Vector{Float64}, xOld1::Vector{Float64}, as::Matrix{Float64}, f;
    while loop < maxit && ch > 1e-6
      loop = loop+1;                                                           # update iteration counter
      # --------------------------------- RL. 1) COMPUTE PHYSICAL DENSITY FIELD
      xTilde = imfilter(reshape(x,nely,nelx),h,bcF)./Hs;                       # compute filtered field
      xPhys[act] = xTilde[act];                                                # modify active elements only
      if ft > 1                                                                # apply projection
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
      
      # -------------------------- RL. 2) SETUP AND SOLVE EQUILIBRIUM EQUATIONS
      sK = (Emin.+xPhys.^penalK*(E0-Emin));                                     # stiffness interpolation
      dsK[act] = penalK*(E0-Emin)*xPhys[act].^(penalK-1);                      # derivative of " "
      sK = reshape(Ke[:]*sK',length(Ke)*nEl,1);
      K = sparse(Iar[:,1],Iar[:,2],sK[:],nDof,nDof);                           # assemble stiffness matrix
      K = K+K'-Diagonal(diag(K));                                                  # symmetrization of K
    # dK = PositiveFactorizations.cholesky(Positive, K[free, free], true);
      dK = factorize(K[free, free]);
      U[free] = dK \ F[free];                                                  # solve equilibrium system
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
    as = Array{Float64}(undef, 0, 0); 
end                    # initialize MMA history parameters
x0,as,lmid=ocUpdate(loop,x[act],dg0[act],g1,dg1[act],ocPar,xOld,xOld1,as,beta,restartAs);

    xOld1 = xOld; xOld = x[act]; x[act] = x0;
    # ----------------------------------------- RL. 8) PRINT AND PLOT RESULTS 
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
    end
end
    
function ocUpdate( loop, xT, dg0, g1, dg1, ocPar, xOld, xOld1, as, beta, restartAs )
    # -------------------------------- definition of asymptotes and move limits
    xU,xL = min.(xT.+ocPar[1],1), max.(xT.-ocPar[1],0);
    if (loop<2.5 || restartAs==1)
        as = xT.+reshape([-0.5, 0.5], 1, 2).*reshape((xU-xL), length(xU-xL), 1)./(beta+1);
    else
        tmp = (xT-xOld).*(xOld-xOld1);
        gm = ones(length(xT),1);
        gm[tmp.>0] .= ocPar[3]; gm[tmp.<0] .= ocPar[2];
        as = xT.+gm.*[-(xOld-as[:,1]) (as[:,2]-xOld)];
    end
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
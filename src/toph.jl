module TopH

using LinearAlgebra
using SparseArrays
using Statistics

export toph


"""
    toph(nelx, nely, volfrac, penal, rmin, ft)

A direct, naive Julia port of the `toph` code listing from the 
"""
function toph(nelx, nely, volfrac, penal, rmin, write::Bool=false)
    x = volfrac * ones(nely,nelx)
    loop = 0
    change = 1.
    dc = zeros(nely,nelx)
    cValues = []

    while change > 0.01
        loop += 1
        xold = x
        c = 0.
        U = FE(nelx,nely,x,penal)

        KE = lk()
        for ely = 1:nely
            for elx = 1:nelx
                n1 = (nely+1)*(elx-1)+ely
                n2 = (nely+1)* elx   +ely
                Ue = U[[n1; n2; n2+1; n1+1]]

                c += (0.001+0.999*x[ely,elx]^penal)*Ue'*KE*Ue
                dc[ely,elx] = -0.999*penal*(x[ely,elx])^(penal-1)*Ue'*KE*Ue
            end
        end
        push!(cValues,c)
        dc = check(nelx,nely,rmin,x,dc)
        x  = OC(nelx,nely,x,volfrac,dc)
        change = maximum(abs.(x-xold))
        
        write && println("Change = ", change, " c = ", c)
        loop >= 1000 && break
    end

    return x, cValues, loop
end


function OC(nelx,nely,x,volfrac,dc)
    l1 = 0; l2 = 100000; move = 0.2
    xnew = zeros(nely,nelx)

    while (l2-l1) > 1e-4
        lmid = 0.5*(l2+l1)
        RacBe = sqrt.(-dc/lmid)
        XB = x .* RacBe

        for i = 1:nelx
            for j = 1:nely
                xji = x[j,i]
                xnew[j,i]= max(0.001,max(xji-move,min(1,min(xji+move,XB[j,i]))))
            end
        end

        if (sum(sum(xnew)) - volfrac*nelx*nely) > 0
            l1 = lmid
        else
            l2 = lmid
        end
    end


    return xnew
end


function check(nelx,nely,rmin,x,dc)
    dcn=zeros(nely,nelx)

    for i = 1:nelx
      for j = 1:nely
        sum=0.0

        for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
          for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
            l = Int64(l); k = Int64(k)
            fac = rmin-sqrt((i-k)^2+(j-l)^2)
            sum = sum+max(0,fac)
            dcn[j,i] += max(0,fac)*x[l,k]*dc[l,k]
          end
        end

        dcn[j,i] = dcn[j,i]/(x[j,i]*sum)
      end
    end

    return dcn
end

function FE(nelx,nely,x,penal)

    KE = lk()
  
      K = spzeros((nelx+1)*(nely+1), (nelx+1)*(nely+1))
      U = zeros((nely+1)*(nelx+1))
      F = zeros((nely+1)*(nelx+1))
      for elx = 1:nelx
        for ely = 1:nely 
          n1 = (nely+1)*(elx-1)+ely
          n2 = (nely+1)* elx   +ely
          edof = [n1; n2; n2+1; n1+1]
          K[edof,edof] += (0.001+0.999*x[ely,elx]^penal)*KE
        end
      end
  
      F .= 0.01
      fixeddofs = Int64(nely/2+1-(nely/20)):Int64(nely/2+1+(nely/20))
      alldofs = 1:(nely+1)*(nelx+1)
      freedofs = setdiff(alldofs,fixeddofs)
  
      U[freedofs] = K[freedofs, freedofs] \ F[freedofs]
      U[fixeddofs] .= 0
  
      return U
  end

  function check(nelx,nely,rmin,x,dc)
    dcn=zeros(nely,nelx)

    for i = 1:nelx
      for j = 1:nely
        sum=0.0

        for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
          for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
            l = Int64(l); k = Int64(k)

            fac = rmin-sqrt((i-k)^2+(j-l)^2)
            sum = sum+max(0,fac)
            dcn[j,i] += max(0,fac)*x[l,k]*dc[l,k]
          end
        end

        dcn[j,i] = dcn[j,i]/(x[j,i]*sum)
      end
    end
    return dcn
end



function lk()
  return [ 2/3 -1/6 -1/3 -1/6
          -1/6  2/3 -1/6 -1/3
          -1/3 -1/6  2/3 -1/6
          -1/6 -1/3 -1/6  2/3 ]
end


end
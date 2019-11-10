
function cma_es(f, x0 ; tol=1e-5,maxiter=1e3)
      N=length(x0)
      λ=4+floor(3*log(N))
      m=x0
      σ=0.3
      μ=floor(λ/2)
      w=log(μ+0.5).-log.(1:μ)
      w=w/sum(w)
      μ_w=(sum(w)^2)/sum(w.^2)
      c_c=(4+μ_w/N)/(N+4+2*μ_w/N)
      c_σ=(μ_w+2)/(N+μ_w+5)
      c_1=2/((N+1.3)^2+μ_w)
      c_μ=min(1-c_1,
              2*(μ_w-2+1/μ_w)/((N+2)^2+μ_w) )
      damp_σ=1+c_σ+ 2*max(0.0,
                  sqrt( (μ_w-1)/(N+1) )-1    )
      p_c=zeros(N,1)
      p_σ=p_c
      B=@SMatrix(I(N))
      D=SMatrix(ones(N,1))
      C=B*Diagonal(D.^2)*B'
      invsqrtC=B*Diagonal( D.^-1 )*B'
      eigeneval=0

      χ_N=sqrt(N)*(1-1/(4*N) +1/(21*N^2) )

      count=0
      normdiff=Inf
      N=length(x0)

      fx=zeros(Integer(λ))

      C=Diagonal(ones(N))
      X=zeros(Integer(λ),N)

      function loop(X,p_c,p_σ,c_1,c_μ,c_c,σ,C,m,
        B,D,invsqrtC,eigeneval)

          for count in 1:maxiter

          oldxbest=X[1,:]
          for i in 1:Integer(λ)
            X[i,:]=(σ.*B*(D.*randn(N,1))+m)'
            fx[i]=f(X[i,:])
            #count=count+1
          end

          X=X[sortperm(fx),:]
          dm=(X[1:Integer(μ),:]'-repeat(m,1,Integer(μ)))*w
          m=m+dm
          #dm=-dm
          p_σ=(1-c_σ)*p_σ+
                  ((sqrt(c_σ*(2-c_σ)*μ_w)*invsqrtC)*dm)/σ

          hsig=norm(p_σ)/sqrt(1-(1-c_σ)^(2*count/λ))χ_N<1.4+2/(N+1)

          p_c=(1-c_c)*p_c+(hsig*sqrt(c_c*(2-c_c)*μ)/σ) *dm

          aux=(1/σ).*(X[1:Integer(μ),:]'-repeat(m,1,Integer(μ)) )
          C=(1-c_1-c_μ)*C+
                  c_1*(p_c*p_c'+(1-hsig)*c_c*(2-c_c)*C)+
                  c_μ*aux*Diagonal(w)*aux'
          σ=σ*exp( (c_σ/damp_σ)*(norm(p_σ)/χ_N-1) )

          if count-eigeneval>λ/(c_1+c_μ)/N/10
            eigeneval=count
            C=UpperTriangular(C)+UpperTriangular(C)'-Diagonal(C)
            B=eigvecs(C)
            D=eigvals(C).^(0.5)
          end
          invsqrtC=B'*Diagonal( D.^-1 )*B
          count=count+1
          normdiff=norm(X[1,:]-oldxbest)
          #normdiff=norm(dm)
          if normdiff< tol #|| maximum(D)>1e7*minimum(D)
                  break
          end
          end
          return(f(X[1,:]),X[1,:])
      end
      loop(X,p_c,p_σ,c_1,c_μ,c_c,σ,C,m,B,D,invsqrtC,eigeneval)
end

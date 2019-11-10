using Animated_CMA_ES, Plots, Revise
function banana(a,b)
  x->(a-x[1])^2+b*(x[2]-x[1]^2)^2
end
f = banana(1.0,1.0)
#f(x)=x[1]^2+x[2]^2

res=cma_es_animated(f,[0.0,6.0],xrange=[-2., 3.], yrange=[-3.,8.],tol=1e-3)
gif(res[3], "cma_es_banana.gif", fps=2)


res=cma_es_animated(f,[0.0,6.0],xrange=[-2., 3.], yrange=[-3.,8.],tol=1e-3)

g(x)=x[1]^2+x[2]^2
res=cma_es_animated(g,[-3.0,-3.0],xrange=[-5., 3.], yrange=[-5.,3.],tol=1e-3,maxiter=20)

gif(res[3], "cma_es_square.gif", fps=2)

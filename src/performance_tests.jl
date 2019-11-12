using Animated_CMA_ES, BenchmarkTools, Profile

N=50
function h(x)
    h=0.0
    for i in 1:N
        h=h+(x[i]+Float64(i))^2
    end
    return h
end

cma_es(h,zeros(N),tol=1e-5)
new_cma_es(h,zeros(N),tol=1e-5)

@btime cma_es(h,zeros(N),tol=1e-5)
@btime new_cma_es(h,zeros(N),tol=1e-5)


function prof()
    cma_es(h,zeros(100),tol=1e-7)
end

function new_prof()
    new_cma_es(h,zeros(100),tol=1e-7)
end


Profile.clear()
Profile.init(n=10^8,delay=0.000000001)
@profile prof()
Profile.print(noisefloor=2.0)

Profile.clear()
Profile.init(n=10^8,delay=0.000000001)
@profile new_prof()
Profile.print(noisefloor=2.0)

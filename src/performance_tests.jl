using Animated_CMA_ES
using BenchmarkTools
using Profile

#start out with a simple optimization


function h(x)
    h=0.0
    for i in 1:100
        h=h+(x[i]+Float64(i))^2
    end
    return h
end

cma_es(h,zeros(100),tol=1e-7)

@btime cma_es(h,zeros(100),tol=1e-7)

#around 2.1 oseconds, and 326000 allocationsv

function prof()
    cma_es(h,zeros(100),tol=1e-7)
end


Profile.clear()
Profile.init(n=10^8,delay=0.000000001)
@profile prof()
Profile.print(noisefloor=2.0)


@code_warntype cma_es(h,zeros(100),tol=1e-7)

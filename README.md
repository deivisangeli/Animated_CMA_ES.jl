# Animated_CMA_ES

[![Build Status](https://travis-ci.com/deivisangeli/Animated_CMA_ES.jl.svg?branch=master)](https://travis-ci.com/deivisangeli/Animated_CMA_ES.jl)
[![Codecov](https://codecov.io/gh/deivisangeli/Animated_CMA_ES.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/deivisangeli/Animated_CMA_ES.jl)

## This is is a course assignment (ECON622)

For grader of pset 8:

First I created a version of the optimizer without the animation to benchmark against

StaticArrays won't help much as parameter space can be really big (I am testing speed minimizing a function of 50 parameters).

Optimize the way to keep C symmetric (Now C=(C+C')./2 )

@simd when calculating new candidate solutions helped a bit.

I got a big improvement on speed by making the update of the varinace-covariance matrix an in-place operation with an specialized function

@inbounds has improved performance too

I am now decomposing the variance matrix in a single step, which improved speed too.

I still need the eigenvalue decomposition as I need both C and C^-{1/2}

Now most of the time (251/312 profiler tics) are being spent in function calls to the objective function. And that can't be avoided.

There is probably no space for paralelization as each update on the optimum depends on the previous one.

The code is now 3 times faster

Run performance_tests.jl to see benchmarking.

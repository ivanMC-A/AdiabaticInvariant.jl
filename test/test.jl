using AdiabaticInvariant
using Plots

# I'm going to test reconstructing several functions for my fourier reconstruction func.

# First we create an error function

function errdiff(f::Function, fr::Function,x::Number)
    abs(f(x)-fr(x))
end

# Then we create a function that gives me an error in a domain L = xarr[end] - xarr[1]

function contErr(f::Function, fr::Function,xarr::AbstractArray)
    ferr = zeros(length(xarr))
    for i in range(1,length(ferr))
        ferr[i] = errdiff(f,fr,xarr[i])
    end
    return ferr
end

# We are going to define the general domain of our functions and the truncation of our series

xs = range(-π, π; length = 501) # Domain
L = xs[end] - xs[1] # Length of our domain

N = 10 # Truncation order


## Test a for f(x) = sin(x) ##

# We initialize our coefficients

aa = zeros(N)
ba = zeros(N)

# We know that b_1 = 1, thus

ba[2] = 1

# We define our original function

fa(x) = sin(x)

# We define our funcition in fourier form

far(x) = AdiabaticInvariant.evalFS(vcat(aa,ba),x,L)

# We calculate the absolute difference (error) of our functions throughout our domain

faerr = contErr(fa,far,xs)

## Test b for f(x) = x ##
N = 10
# We initialize our coefficients

ab = zeros(N)
bb = zeros(N)

# Our coefficients only exist in b and are of the form 2*(-1)^{n+1}/n
# So we create a function for that

for i in range(2,length(bb))
    bb[i] = 2(-1)^(i)/(i-1)
end




# We define our original function

fb(x) = x

# We define our funcition in fourier form

fbr(x) = AdiabaticInvariant.evalFS(vcat(ab,bb),x,L)

# We calculate the absolute difference (error) of our functions throughout our domain

fberr = contErr(fb,fbr,xs)

plot(xs,fberr)
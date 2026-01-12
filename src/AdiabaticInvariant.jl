module AdiabaticInvariant
export JInvariant, SlabJ, evaluate, deval, fourierFunc

using LinearAlgebra

# -----------------------------------------------
# Definitions
# -----------------------------------------------

"""
    abstract type JInvariant end

        Parent type for all adiabatic invariant representations.

"""

abstract type JInvariant end

"""

    struct SlabJ <: JInvariant

        Spectral representations (coefficients) of an adiabatic invariant.
        a: Array containing Spectral coefficients.

"""

struct SlabJ <: JInvariant

    a::AbstractArray

end

# -----------------------------------------------
# Methods
# -----------------------------------------------

function slabJ(N1,N2,N3)

    a = zeros(N1,N2,N3)
    
    return SlabJ(a)
end
"""
    evalFS(coef::AbstractArray)

Reconstructs, in a domain L, and evaluates, at a point x ∈ L, a function f based on its fourier coefficients (coef). The coef
should be of length N*2 and coef[N + 1] = 0 which corresponds to b0 = 0.
"""

function evalFS(coef::AbstractArray, x::Number, L::Float64)

    n_dims = ndims(coef)

    if n_dims > 1
        return error("coef::AbstractArray dimention is bigger than 1. coef must be a 1D array")
    end

    if iseven(length(coef))
        N = length(coef) ÷ 2
    else
        return error("Length of coef::AbstractArray is not even. coef most be an even array")

    end
    θ = 2*π*x/L

    
    f = 0
    for i in range(start = 0, stop = N-1)
        if i == 0
            f += coef[i + 1]
        else
            f += coef[i + 1]*cos(i*θ) + coef[i + N + 1]*sin(i*θ)
        end
    end
    return f
end


"""
    evalDFS(coef::AbstractArray)

Reconstructs (in a domain L) and evaluates (at a point x ∈ L) the derivative of a function f based on its fourier coefficients (coef). The coef
should be of length N*2 and coef[N + 1] = 0 which corresponds to b0 = 0.
"""

function evalDFS(coef::AbstractArray, x::Float64, L::Float64)

    n_dims = ndims(coef)

    if n_dims > 1
        return error("coef::AbstractArray dimention is bigger than 1. coef must be a 1D array")
    end

    if iseven(length(coef))
        N = length(coef)/2
    else
        return error("Length of coef::AbstractArray is not even. coef most be an even array")

    end
    θ = 2*π*x/L

    
    f = 0
    for i in range(start = 0, stop = N)
        if i == 0
            f += 0
        end
        f += -2*π/L*coef[i + 1]*sin(i*θ) + 2*π/L*coef[i + N + 1]*cos(i*θ)
    end
    return f
end

"""
# Function call

'evaluate(J,fabx,faby,chev,x)'

# Description

Evaluate adiabatic invariant 'J' at a point 'x'.

"""

function evaluate(J::SlabJ,x)
    # TODO: code spcetral reconstruction
    error("evaluate(::$(typeof(J)), ::$(typeof(x))) not implemented")
    
end

"""
    function deval(J::SlabJ, x)

Evaluate adiavatic invariant 'J' derivatives at a point 'x'.

"""

function deval(J::SlabJ,x)
    # TODO: code derivatives w.r.t. X, Y, and r for J.
    # TODO: code spcetral reconstruction for derivatives.
    error("deval(::$(typeof(J)), ::$(typeof(x))) not implemented")

end

###### Fourier Stuff

"""
    function fourA(N::Integer)

Construct the Fourier matrix evaluated at equidistant points.

"""

function fourA(N::Integer)
    A = zeros(ComplexF64,N,N)
    M = Int((N-1)/2)
    for m in -M:M
        for n in 0:(N-1)
            A[m + M + 1,n + 1] = exp(im*(2*π*m*n)/(N))
        end
    end
    return A
end

"""
    function getChebCoef(N::Number, FN::AbstractArray)

Compute Chebyshev polynomial coefficients from function evaluated at Chebyshev points.

"""

function getFourCoef(N::Number, FN::AbstractArray)
    A = fourA(N)
    At = A'
    W = At * A
    b = At * FN
    x = W \ b
    return x
end


###### Chebyshev stuff

"""
    evalChev(coef::AbstractArray, x::Number)
This function evaluates Chebyshev coefficients at a point x∈[-1,1]
"""

function evalChev(coef::AbstractArray, x::Number)
    n_dims = ndims(coef)

    if n_dims > 1
        return error("coef::AbstractArray dimention is bigger than 1. coef must be a 1D array")
    end

    cs = 0
    θ = acos(x)
    for i in range(start = 0, stop = N-1)
        cs += coef[i + 1]*cos(i*θ)
    end
    return cs
end

"""
    evalChev(coef::AbstractArray, x::Number, D::AbstractArray)

Evaluate the Chebyshev series described by `coef` at the physical grid point `x` that lies
in the domain `D`. The point is affinely mapped into [-1, 1] before evaluation. `coef` must
be a 1D array ordered as T0, T1, …, T_{N-1}; returns the scalar value of the series at `x`.
"""
function evalChev(coef::AbstractArray, x::Number, D::AbstractArray)
    n_dims = ndims(coef)

    if n_dims > 1
        return error("coef::AbstractArray dimention is bigger than 1. coef must be a 1D array")
    end

    cs = 0
    xx = grid2chev(D,x)
    θ = acos(xx)
    for i in range(start = 0, stop = N-1)
        cs += coef[i + 1]*cos(i*θ)
    end
    return cs

end

"""
    chp1st(N:Number)
This function generates N Chebyshev points of the first kind
"""

function chept1st(N::Number)
    x = zeros(N)
    
    for k in range(start = 0, stop = N-1)
        x[k+1] = cos((k + 0.5)*π/N)
    end
    return x
end

"""
    chp2nd(N:Number)
This function generates N+1 Chebyshev points of the second kind
"""

function chept2nd(N::Number)
    x = zeros(N+1)

    for i in range(0, stop = N)
        x[i + 1] = cos(i*π/N)
    end
    return x
end

"""
    function chebpt(N::Number, k::String)

Given a number 'N' and a kind 'k', this function calculates the N or N+1 Chebyshev points 
of the first kind k = 'F' or second kind k = 'S'. 
"""

function chept(N::Number, k::String)

    if k == "F"
        x = chept1st(N)
        return x
    elseif k == "S"
         x = chept2nd(N)
        return x
    else
        error("'k' states the kind of chebyshev points you want. k = 'F' is for first kind and k = 'S' is for second kind")
    end

end

"""
    function chebA(N::Number)

Construct the Chebyshev matrix of the first kind evaluated at Chebyshev points.

"""

function chebA(N::Number)
    A = zeros(Float64,N,N)
    xN = chept1st(N)
    for i in 1:N
        θn = acos(xN[i])
        for j in 0:(N-1)
            A[i,j+1] = cos((j) * θn)
        end
    end
    return A
end

"""
    function getChebCoef(N::Number, FN::AbstractArray)

Compute Chebyshev polynomial coefficients from function evaluated at Chebyshev points.

"""

function getChebCoef(N::Number, FN::AbstractArray)
    A = chebA(N)
    At = A'
    W = At * A
    b = At * FN
    x = W \ b
    return x
end

"""
    function chev2grid(D::AbstractArray,y::Number)
Given a chebyshev point y and a domain D, this function maps the point y back to a number x ∈ D
"""
function chev2grid(D::AbstractArray,y::Number)
    a = D[1]
    b = D[end]
    x = 0.5*((a-b)*y + (a+b))
    return x
end

"""
    function grid2chev(D::AbstractArray,y::Number)
Given a point x ∈ D and a Chebyshev domain [-1,1], this function maps the point x back to a y x ∈ [-1,1]
"""
function grid2chev(D::AbstractArray,x::Number)
    a = D[1]
    b = D[end]
    y = 2/(b-a)*(x-a) - 1
    return y
end

end # module

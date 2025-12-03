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

end # moduleS

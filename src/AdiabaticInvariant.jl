module AdiabaticInvariant
export JInvariant, SlabJ, evaluate, deval

using LinearAlgebra, FourierSeries

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

end # module

module AdiabaticInvariant
export AdiabaticInvariant, SlabJ, evaluate, deval

using LinearAlgebra

# -----------------------------------------------
# Definitions
# -----------------------------------------------

"""
    abstract type AdiabaticInvariant end

        Parent type for all adiabatic invariant representations.

"""

abstract type AdiabaticInvariant end

"""

    struct SlabJ <: AdiabaticInvariant

        Spectral representations (coefficients) of an adiabatic invariant.
        a: Array containing Spectral coefficients.

"""

struct SlabJ <: AdiabaticInvariant

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
    evaluate(J::SlabJ,x)

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

module AdiabaticInvariant
export JInvariant, SlabJ, evaluate, deval, fourierFunc, chebA, fourA, get_FFO

using LinearAlgebra
using DifferentialEquations

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

###### Fourier Stuff

"""
    evalFS(coef::AbstractVector)

Reconstructs, in a domain L, and evaluates, at a point x ∈ L, a function f based on its fourier coefficients (coef). The coef
The first coefficients correspond to a_n and the last coefficients correspond to b_n
"""

function evalFS(coef::AbstractVector, x::Number, L::Float64)
    θ = 2*π*x/L
    N = length(coef)
    f = 0
    for i in 0:Int(N/2 - 1/2)
        if i == 0
            f += coef[i + 1]
        else
            f += coef[i + 1]*cos(i*θ) + coef[i + Int(N/2 + 1/2)]*sin(i*θ)
        end
    end
    return f
end

"""
    evalDFS(coef::AbstractVector)

Reconstructs (in a domain L) and evaluates (at a point x ∈ L) the derivative of a function f based on its fourier coefficients (coef). The coef
should be of length N*2 and coef[N + 1] = 0 which corresponds to b0 = 0.
"""

function evalDFS(coef::AbstractVector, x::Float64, L::Float64)
    θ = 2*π*x/L
    N = length(coef)
    f = 0
    for i in 0:Int(N/2 - 1/2)
        if i == 0
            f += 0
        end
        f += -2*π *i/L*coef[i + 1]*sin(i*θ) + 2*π*i/L*coef[i + Int(N/2 + 1/2)]*cos(i*θ)
    end
    return f
end

"""
    function fourA(N::Integer)

Construct the Fourier matrix evaluated at equidistant points.

"""

function fourA(N::Integer)

    M = Int((N-1)/2)

    return [exp(im*(2*π*m*n)/(N)) for m in -M:M, n in 0:N-1]
end


###### Chebyshev stuff

"""
    evalChev(coef::AbstractVector, x::Number)
This function evaluates Chebyshev coefficients at a point x∈[-1,1]
"""

function evalChev(coef::AbstractVector, x::Number)
    N = length(coef)
    cs = 0
    θ = acos(x)
    for i in 0:N-1
        cs += coef[i + 1]*cos(i*θ)
    end
    return cs
end

"""
    evalChev(coef::AbstractVector, x::Number, D::Tuple)

Evaluate the Chebyshev series described by `coef` at the physical grid point `x` that lies
in the domain `D`. The point is affinely mapped into [-1, 1] before evaluation. `coef` must
be a 1D array ordered as T0, T1, …, T_{N-1}; returns the scalar value of the series at `x`.
"""
function evalChev(coef::AbstractVector, x::Number, D::Tuple)
    N = length(coef)
    cs = 0
    xx = grid2chev(D,x)
    θ = acos(xx)
    for i in 0:N-1
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
    
    for k in 0:N-1
        x[k+1] = -cos((k + 0.5)*π/N)
    end
    return x
end

"""
    chp2nd(N:Number)
This function generates N+1 Chebyshev points of the second kind
"""

function chept2nd(N::Number)
    x = zeros(N+1)

    for i in 0:N
        x[i + 1] = cos(i*π/N)
    end
    return x
end

"""
    function chebpt(N::Number, k::String)

Given a number 'N' and a kind 'k', this function calculates the N or N+1 Chebyshev points 
of the first kind k = 'F' or second kind k = 'S'. 
"""

function chept(N::Number; lobatto::Bool = true)

    if lobatto
        return chept2nd(N)
    end

    return chept1st(N)
end

"""
    function chebA(N::Number)

Construct the Chebyshev matrix of the first kind evaluated at Chebyshev points.

"""

function chebA(N::Number)
    [cos(j * θn) for θn in acos.(chept1st(N)), j in 0:N-1]
end

"""
    function chev2grid(D::Tuple,y::Number)
Given a chebyshev point y and a domain D, this function maps the point y back to a number x ∈ D
"""
function chev2grid(D::Tuple,y::Number)
    a = D[1]
    b = D[end]
    x = 0.5*((a-b)*y + (a+b))
    return x
end

"""
    function grid2chev(D::Tuple,y::Number)
Given a point x ∈ D and a Chebyshev domain [-1,1], this function maps the point x back to a y x ∈ [-1,1]
"""
function grid2chev(D::Tuple,x::Number)
    a = D[1]
    b = D[end]
    y = 2/(b-a)*(x-a) - 1
    return y
end

###### coefficients stuff

function getW(A::AbstractMatrix)
    W = Diagonal(A' * A)
    return W
end

function getCoef(A::AbstractMatrix,W::AbstractMatrix,FN::AbstractArray)
    W\(A'*FN)
end

###### Trajectories stuff

function get_FFO(ϵ::Float64,B::Function; tmax = 1e3,  solver=Tsit5(), psec = true)

    # Creating differential equation to solve
    function f!(du,u,p,t)
        # u₁ = x, u₂ = y, # u₃ = z
        # u₄ = vₓ, u₅ = v\_y, u₆ = v\_z
        du[1] = ϵ * u[4]
        du[2] = ϵ * u[5]
        du[3] = ϵ * u[6]
        du[4] = u[5]*B(u)[3] - u[6]*B(u)[2] 
        du[5] = -(u[4]*B(u)[3] - u[6]*B(u)[1])
        du[6] = u[4]*B(u)[2] - u[5]*B(u)[1]
    end

    if psec
            # Write conditioning
        end


    function FFO(u0::AbstractArray; p = nothing)
        prob = ODEProblem(f!,u0, (0.0, tmax), p)
        sol = solve(prob, solver)
        return sol, sol.t[end]

    end

    return FFO
end

end # module

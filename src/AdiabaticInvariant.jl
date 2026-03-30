module AdiabaticInvariant
export JInvariant, slabJ, evaluate, full_J, get_eval_points, deval, numN, numNt, denomD,get_FFO
using OMEinsum
using LinearAlgebra
using DifferentialEquations
using ForwardDiff
using StaticArrays

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

    N::AbstractArray
    bd::Tuple
    ϵ::Float64
    B::Function
    coeff::AbstractArray

end

# -----------------------------------------------
# Methods
# -----------------------------------------------

function slabJ(N1::Int64,N2::Int64,N3::Int64,x1min::Float64,x1max::Float64,x2min::Float64,x2max::Float64,x3min::Float64,x3max::Float64,ϵ::Float64,B::Function)
    N = [N1, N2, N3]
    bd = tuple(x1min,x1max,x2min,x2max,x3min,x3max)
    coeff = coef_J(N,bd,ϵ,B)
    return SlabJ(N,bd,ϵ,B,coeff)
end

function slabJ(N1::Int64,N2::Int64,N3::Int64,x1min::Float64,x1max::Float64,x2min::Float64,x2max::Float64,x3min::Float64,x3max::Float64,F)
    N = [N1, N2, N3]
    bd = tuple(x1min,x1max,x2min,x2max,x3min,x3max)
    coeff = coef_J(N,bd,F)
    return SlabJ(N,bd, 0.1, F, coeff)
end

"""
# Function call

'evaluate(J,fabx,faby,chev,x)'

# Description

Evaluate adiabatic invariant 'J' at a point 'x'.

"""

function evaluate(J::SlabJ,x)
    G1,G2,G3 = get_inv_trans(J.bd)
    vx = fourBase(J.N[1], G1(x[1]))
    vy = conj(fourBase(J.N[2], G2(x[2])))
    vz = chebBase(J.N[3], G3(x[3]))
    return real(ein"i,ijk,j,k->"(vx,J.coeff,vy,vz)[])
end

"""
    function deval(J::SlabJ, x)

Evaluate adiavatic invariant 'J' derivatives at a point 'x'.

"""

function deval(J::SlabJ,x)
    G1,G2,G3 = get_inv_trans(J.bd)
    vx = fourBase(J.N[1], G1(x[1]))
    dvx = devFourBase(J.N[1], G1(x[1]), J.bd[1], J.bd[2])
    vy = conj(fourBase(J.N[2], G2(x[2])))
    dvy = conj(devFourBase(J.N[2], G2(x[2]), J.bd[3], J.bd[4]))
    vr = chebBase(J.N[3], G3(x[3]))
    dvr = DevChebBase(J.N[3], G3(x[3]), J.bd[5], J.bd[6])

    Jx = real(ein"ijk,i,j,k->"(J.coeff,dvx,vy,vr)[])
    Jy = real(ein"ijk,i,j,k->"(J.coeff,vx,dvy,vr)[])
    Jr = real(ein"ijk,i,j,k->"(J.coeff,vx,vy,dvr)[])
    return SVector(Jx,Jy,Jr)
end

function numNt(J::SlabJ,x)
    G1,G2,G3 = get_inv_trans(J.bd)
    vx = fourBase(J.N[1], G1(x[1]))
    dvx = devFourBase(J.N[1], G1(x[1]), J.bd[1], J.bd[2])
    vy = conj(fourBase(J.N[2], G2(x[2])))
    dvy = conj(devFourBase(J.N[2], G2(x[2]), J.bd[3], J.bd[4]))
    vr = chebBase(J.N[3], G3(x[3]))
    dvr = DevChebBase(J.N[3], G3(x[3]), J.bd[5], J.bd[6])

    Jx = real(ein"ijk,i,j,k->"(J.coeff,dvx,vy,vr)[])
    Jy = real(ein"ijk,i,j,k->"(J.coeff,vx,dvy,vr)[])
    Jr = real(ein"ijk,i,j,k->"(J.coeff,vx,vy,dvr)[])

    U = -J.ϵ^2*Jy
    V = J.ϵ^2*Jx
    return [U,V]
end

function numN(J::SlabJ,x)
    G1,G2,G3 = get_inv_trans(J.bd)
    vx = fourBase(J.N[1], G1(x[1]))
    dvx = devFourBase(J.N[1], G1(x[1]), J.bd[1], J.bd[2])
    vy = conj(fourBase(J.N[2], G2(x[2])))
    dvy = conj(devFourBase(J.N[2], G2(x[2]), J.bd[3], J.bd[4]))
    vr = chebBase(J.N[3], G3(x[3]))
    dvr = DevChebBase(J.N[3], G3(x[3]), J.bd[5], J.bd[6])

    

    jx = ein"i,j,k->ijk"(dvx,vy,vr)
    jy = ein"i,j,k->ijk"(vx,dvy,vr)
    U = -J.ϵ^2*jy
    V = J.ϵ^2*jx
    return cat(reshape(U, 1, size(U)...), reshape(V, 1, size(V)...); dims=1)
end

function denomD(J::SlabJ,x)
    G1,G2,G3 = get_inv_trans(J.bd)
    vx = fourBase(J.N[1], G1(x[1]))
    vy = conj(fourBase(J.N[2], G2(x[2])))
    dvy = conj(devFourBase(J.N[2], G2(x[2]), J.bd[3], J.bd[4]))
    vr = chebBase(J.N[3], G3(x[3]))
    dvr = DevChebBase(J.N[3], G3(x[3]), J.bd[5], J.bd[6])

    jy = ein"i,j,k->ijk"(vx,dvy,vr)
    jr = ein"i,j,k->ijk"(vx,vy,dvr)
    return J.B(x)/x[3] * (jr + J.ϵ*J.B(x)^(-1)*jy)
end

function deriv_B(B)
    ∇B  = x -> ForwardDiff.gradient(B,x)
    HB = x -> ForwardDiff.hessian(B,x)
    return ∇B, HB
end

# J for B = B(x,y)e_z

# 0th order J

function perturvedJ(B,x::AbstractVector; perp = true)
    if perp
        vₚ = x[3]
    else
        vₚ = sqrt(x[3]^2 + x[4]^2)
    end
    return vₚ^2/(2*B(x))
end

# 1st order J

function perturvedJ(ϵ::Float64,B,∇B,x::AbstractVector; perp = true)
    if perp
        vₚ = x[3]
    else
        vₚ = sqrt(x[3]^2 + x[4]^2)
    end
    return vₚ^2/(2*B(x)) + ϵ * vₚ^3 * ∇B(x)[2]/(2*B(x)^3)
end

# 2nd order J

function perturvedJ(ϵ::Float64,B,∇B,HB,x::AbstractVector; perp = true)
    if perp
        vₚ = x[3]
    else
        vₚ = sqrt(x[3]^2 + x[4]^2)
    end
    return vₚ^2/(2*B(x)) + ϵ * vₚ^3 * ∇B(x)[2]/(2*B(x)^3) + ϵ^2 * vₚ^4 / (16 * B(x)^5) * (3 * (∇B(x)[1]^2 + 5 * ∇B(x)[2]^2) - B(x) * (HB(x)[1,1] + 5 * HB(x)[2,2]))
end

function get_for_trans(bd::Tuple)
    F1 = x -> AdiabaticInvariant.gTrans.(x, bd[1], bd[2]) # Get transformation function from fourier points to domain
    F2 = y -> AdiabaticInvariant.gTrans.(y, bd[3], bd[4]) # Get transformation function from fourier points to domain
    F3 = r -> AdiabaticInvariant.gTrans.(r, bd[5], bd[6]; F = false) # Get transformation function from chebyshev points to domain
    return F1, F2, F3
end

function get_inv_trans(bd::Tuple)
    InvF1 = x -> AdiabaticInvariant.gInv.(x, bd[1], bd[2]) # Get inverse transformation function from domain to fourier points
    InvF2 = y -> AdiabaticInvariant.gInv.(y, bd[3], bd[4]) # Get inverse transformation function from domain to fourier points
    InvF3 = r -> AdiabaticInvariant.gInv.(r, bd[5], bd[6]; F = false) # Get inverse transformation function from domain to chebyshev points
    return InvF1, InvF2, InvF3
end

function get_trans_func(bd::Tuple)
    F1, F2, F3 = get_for_trans(bd)
    InvF1, InvF2, InvF3 = get_inv_trans(bd)
    return F1, InvF1, F2, InvF2, F3, InvF3
end

function get_col_points(N::AbstractArray;Lx = 2*π,Ly = 2*π)
    x = evenpts(N[1];L = Lx)
    y = evenpts(N[2];L = Ly)
    r = chept1st(N[3])
    return x,y,r
end

function get_dcol_points(x::AbstractArray,y::AbstractArray,r::AbstractArray, bd::Tuple)
    F1, F2, F3 = get_for_trans(bd)
    xt = F1(x)
    yt = F2(y)
    rt = F3(r)
    return xt, yt, rt
end

function get_eval_points(N::AbstractArray,bd::Tuple;Lx = 2*π,Ly = 2*π)
    x, y, r = get_col_points(N;Lx = Lx, Ly = Ly)
    xt, yt, rt = get_dcol_points(x,y,r,bd)
    return x, xt, y, yt, r, rt
end

function full_J(ϵ,B)
    ∇B, HB = deriv_B(B)
    return x -> perturvedJ(ϵ,B,∇B,HB,x)
end

function get_Jn(x::AbstractVector,y::AbstractVector,r::AbstractVector,J)
    Jn = [J(SVector(X,Y,R)) for X in x, Y in y, R in r]
    return Jn
end

function gram_matrices(N::AbstractArray)
    Φ, ΦT, WΦ= fourMatrices(N[1])
    Θ, ΘT, WΘ = fourMatrices(N[2])
    Ψ, ΨT, WΨ = chebMatrices(N[3])
    return ΦT, WΦ, Θ, WΘ, Ψ, WΨ
end

function coef_3d(Fijk,A1T,WA1,A2,WA2,A3,WA3)
    a = ein"ij,jk,kmn->imn"(inv(WA1),A1T,Fijk)
    b = ein"ijk,kl,lm->ijm"(a,A3,inv(WA3))
    c = ein"ijk,jl,lm->imk"(b, A2,inv(WA2))
    return c
end

function coef_J(N::AbstractArray,bd::Tuple,F;Lx = 2*π,Ly = 2*π)
    # Get points and transformations functions
    x, xt, y, yt, r, rt = get_eval_points(N,bd;Lx = Lx, Ly = Ly)
    # Get Gram Matrices for Fourier and Chebyshev
    ΦT, WΦ, Θ, WΘ, Ψ, WΨ = gram_matrices(N)
    # Get J values at transformed points
    Jijk = get_Jn(xt,yt,rt,F)
    Cmnl = coef_3d(Jijk,ΦT,WΦ,Θ,WΘ,Ψ,WΨ)
    return Cmnl
end

function coef_J(N::AbstractArray,bd::Tuple,ϵ::Float64,B::Function;Lx = 2*π,Ly = 2*π)
    # Get points and transformations functions
    x, xt, y, yt, r, rt = get_eval_points(N,bd;Lx = Lx, Ly = Ly)
    # Get Gram Matrices for Fourier and Chebyshev
    ΦT, WΦ, Θ, WΘ, Ψ, WΨ = gram_matrices(N)
    # Get J values at transformed points
    Jijk = get_Jn(xt,yt,rt,full_J(ϵ,B))
    Cmnl = coef_3d(Jijk,ΦT,WΦ,Θ,WΘ,Ψ,WΨ)
    return Cmnl
end
### Transformation functions

"""
    function gTrans(x::Number, a::Number, b::Number; F = true, L = 2*π)
This function transforms a point x in [0, L] to a point in [a, b] if F = true and 
transformsa point x in [-1, 1] to a point in [a, b] if F = false.

"""

function gTrans(x::Number, a::Number, b::Number; F = true, L = 2*π)
    if F
        return (b-a)/L * x + a
    else
        return (b-a)*(x+1)/2 + a
    end
end

"""
    function gTrans(x::AbstractArray, a::Number, b::Number; F = true, L = 2*π)
This function transforms a set of arrays xᵢ ∈ (x, y, r) of domain [0, L], to an array of 
points in [a, b] if F = true and transforms an array of points x in [-1, 1] to an array of
points in [a, b] if F = false. 
"""

function trans_pts(x,y,r,bd;L = 2*π)
    # bd are the boundaries of my domain. It should be an vector (ax,bx,ay,...,br)
    xt = gTrans.(x,bd[1],bd[2], L = L)
    yt = gTrans.(y,bd[3],bd[4], L = L)
    rt = gTrans.(r,bd[5],bd[6]; F = false)
    return xt, yt, rt
end

function gInv(y::Number, a::Number, b::Number; F = true, L = 2*π)
    if F
        return L*(y-a)/(b-a)
    else
        return 2*(y-a)/(b-a) - 1
    end
end

###### Fourier Stuff

"""
    fourBase(N::Number, x::Number)

    This function evaluates the Fourier basis of order N at a point x.
        It returns a vector of length N where the i-th entry corresponds to exp(i*x).

"""

function fourBase(N::Number, x::T) where {T <: Number}
    v = zeros(Complex{T}, N)
    M = Int((N-1)/2)
    for i in -M:M
        v[i + M + 1] = exp(im*(i)*x)
    end
    return v
end

"""
    devFourBase(N::Number, x::Number)

    This function evaluates the derivative of a Fourier basis of order N at a point x.
        It returns a vector of length N where the i-th entry corresponds to im*i*exp(i*x).

"""

function devFourBase(N::Number, x::T,a,b;L = 2π) where {T <: Number}
    v = zeros(Complex{T}, N)
    df =  L/(b - a)
    M = Int((N-1)/2)
    for n in -M:M
        v[n + M + 1] = df*n*exp(im*n*x + im*π/2)
    end
    return v
end


function evenpts(N::Int;L = 2π)
    x = zeros(N)
    for i in 1:N
        x[i] = (i-1)*L/N
    end
    return x
end

"""
    evalFS(coef::AbstractVector)

Reconstructs, in a domain L, and evaluates, at a point x ∈ L, a function f based on its fourier coefficients (coef). The coef
The first coefficients correspond to a_n and the last coefficients correspond to b_n
"""

function evalFS(coef::AbstractVector, x::Number, D::Tuple; Real = true)
    θ = gInv(x, D[1], D[2])
    f = 0
    M = length(coef)
    if Real
        
        for i in 0:Int(M/2 - 1/2)
            if i == 0
                f += coef[i + 1]
            else
                f += coef[i + 1]*cos(i*θ) + coef[i + Int(M/2 + 1/2)]*sin(i*θ)
            end
        end
        return f
    else
        N = Int((M-1)/2)
        for i in -N:N
            f += coef[i + N + 1]*exp(im*i*θ)
        end
        return f
    end
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

function fourA(M::Integer)
    N = (M-1)/2
    return [exp(im*(2*π*m*n)/(M)) for n in 0:M-1, m in -N:N]
end

function fourMatrices(N::Integer)
    A = fourA(N)
    W = getW(A)
    return A, A', W
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

function evalChev(N::Number, x::Number)
    vec = zeros(N)
    θ = acos(x)
    for i in 0:N-1
        vec[i + 1] = cos(i*θ)
    end
    return vec
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

function chebMatrices(N::Number)
    A = chebA(N)
    W = getW(A)
    return A, A', W
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

"""
    chebBase(N::Number, x::Number)

    This function evaluates the Chebyshev basis of the first kind of order N at a point x.
        It returns a vector of length N where the i-th entry corresponds to T_{i-1}(x).

"""

function chebBase(N::Number, x::Number)
    θ = acos(x)
    v = zeros(N)
    for i in 0:N-1
        v[i + 1] = cos(i*θ)
    end
    return v
end

function DevChebBase(N::Integer, x::Number, a,b)
    θ = acos(x)
    v = zeros(N)
    if x == 1 || x == -1
        for i in 0:(N-1)
            v[i + 1] = 2/(b - a)*i^2
        end
        return v
    else
        for i in 0:(N-1)
            v[i + 1] = 2/(b - a)*i * sin(i*θ) / sin(θ)
        end
        return v
    end
end    

###### coefficients stuff

function getW(A::AbstractMatrix)
    W = real(Diagonal(ein"ik,kj->ij"(A', A)))
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

    function affect!(integrator)
        u = integrator.u
        vx = u[4]
        vy = u[5]
        if abs(vy) < 1e-5 && vx > 0
            terminate!(integrator)
        end 
    end

    function affect1!(integrator)
        integrator.u[1] = mod(integrator.u[1], 2π)
    end

    function affect2!(integrator)
        integrator.u[1] = mod(integrator.u[1], 2π)
    end


    if psec
            condition(u, t, integrator) = u[5]

            cb = ContinuousCallback(condition,affect!)

            condition1(u,t,integrator) = u[1] - 2π

            cb1 = ContinuousCallback(condition1, affect1!)

            condition2(u,t,integrator) = u[1] - 2π

            cb2 = ContinuousCallback(condition2, affect2!)

            cbs = CallbackSet(cb, cb1, cb2)
    end


    function FFO(u0::AbstractArray; p = nothing, abstol=1e-10, reltol=1e-10)
        prob = ODEProblem(f!,u0, (0.0, tmax), p)
        if psec
            sol = solve(prob, solver; callback = cbs, abstol = abstol, reltol = reltol)
        else
            sol = solve(prob, solver; abstol = abstol, reltol = reltol)
        end

        return sol, sol.t[end]

    end

    return FFO
end

end # module
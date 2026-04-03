## Packages

using AdiabaticInvariant
using OMEinsum
using CairoMakie
using LinearAlgebra
## Test functions

"""
    function gyroradius(eps,v::Number, B)
        return eps * v / B
    end
    Return the gyroradius given the magnitude of the velocity, the magnetic field, and the small parameter.
"""

function gyroradius(eps,v::Number, B)
    return eps * v / B
end

"""
    function gyroradius(eps,v::AbstractVector, B)
        return eps * norm(v) / B
    end
    Return the gyroradius given the velocity vector, the magnetic field, and the small parameter.
"""

function gyroradius(eps,v::AbstractVector, B)
    return eps * norm(v) / B
end

"""
    function get_solution(B,ϵ,u0)
        ## Magnetic fields
        bint = x -> B(x, scalar = false)

        b = x -> B(x)
        ## Full orbit forward mapping

        FFO = get_FFO(ϵ, Bint)
        
        # Initial conditions

        dx0 = (2* E0/ m)^(0.5)

        u0 = [u0[1], u0[2], u0[3], dx0, 0, 0]

        # Evolving the system

        u, t = FFO(u0)
        return u0, u, t
    end
    It recieves as input the magnetic field, the small parameter, and the initial condition (x,y,z).
    and  time of intersection.
"""

function get_solution(B, ϵ, u0, E0, m)
    ## Full orbit forward mapping

    FFO = get_FFO(ϵ, B)
    
    # Initial conditions

    dx0 = (2* E0/ m)^(0.5)

    u0 = [u0[1], u0[2], u0[3], dx0, 0, 0]

    # Evolving the system

    u, t = FFO(u0)
    return u0, u, t
end

"""
    get_plot(u)
    Plot the solution trajectory in the x-y plane.
"""

function get_plot(u,gr)

    circx = gr.*cos.(0:0.1:2π) .+ u[1,1]

    circy = gr.*sin.(0:0.1:2π) .+ (u[2,1]- gr)

    f = Figure()

    ax = Axis(f[1, 1], xlabel = "x", ylabel = "y", title = "Solution")

    lines!(ax, u[1, :], u[2, : ])

    scatter!(ax, circx, circy, markersize = 5, color = :green)

    display(f)

end

"""
    get_J(dx, B, ϵ, N)
    Construct the adiabatic invariant for the system of interest
"""

function get_J(dx, B, ϵ, N::Tuple)

    ### Setting domain parameters
    x1i, x1f = 0.0, 2π
    x2i, x2f = 0.0, 2π
    x3i, x3f = -5 + dx, 5 + dx

    ### Setting polynomial parameters
    Nx, Ny, Nr = N[1], N[2], N[3]

    J = slabJ(Nx, Ny, Nr, x1i, x1f, x2i, x2f, x3i, x3f, ϵ, B)
    return J
end

"""
    function dotX(x0,x,t)
        function FOtPS(m)
            return [m[1],m[2],m[4]]
        end
        mi = FOtPS(x0)
        mf = FOtPS(x)
        dotx = (mf - mi)./t
        midx = (mf + mi)./2
        return dotx, midx
    end
    Calculate the vector of the time derivative of the poincaré section coordinates using a 
    backward finite difference. The function FOtPS is used to extract the poincaré section
    coordinates from the full orbit coordinates. dotX returns the vector of the time 
    derivative of the poincaré section coordinates and the midpoint of the poincaré section 
    coordinates. The midpoint is used to evaluate the adiabatic invariant at the midpoint of 
    the trajectory between the initial and final points. This is used to get a better estimate
     of the error in the adiabatic invariant.
"""

function dotX(x0,x,t)
    function FOtPS(m)
        return [m[1],m[2],m[4]]
    end
    mi = FOtPS(x0)
    mf = FOtPS(x)
    dotx = (mf - mi)./t
    midx = (mf + mi)./2
    return dotx, midx
end

"""
    Function linV(x,dx,N,D)
        s = size(D(x))
        XD = reshape(dx[1:2],2,1,1,1) .*reshape(D(x),1,s[1],s[2],s[3])
        return XD- N(x)
    end
    Calculate the vector dotX * D - N where D and N are the denominator and numerator of the
    vector field G(x) = N(x)/D(x). G is applied to the J coeffients and this should gave the
    EoM in the poincaré section. This needs to be contracted with the J coeffiecients to get 
    a vector ∈ R^2. If J is a good adiabatic invariant this vector should be close to zero.
"""

function linV(x,dx,N,D)
    s = size(D(x))
    XD = reshape(dx[1:2],2,1,1,1) .*reshape(D(x),1,s[1],s[2],s[3])
    return XD- N(x)
end

"""
    Function Aj(x,dx,N,D,J)
        lin = linV(x,dx,N,D)
        a = real(ein"ijk,ijk->"(J.coeff, lin[1,:,:,:]))[]
        b = real(ein"ijk,ijk->"(J.coeff, lin[2,:,:,:]))[]
        return [a,b]
    end
    Calculate the vector Ai = J * (dotX * D - N) where J is tensor of coefficients of the
    adiabatic invariant, D and N are the denominator and numerator of the vector field G(x) = N(x)/D(x), 
    and dotX is the vector of the time derivative of the poincaré section coordinates. If J is a good
    adiabatic invariant this vector should be close to zero.
"""

function Aj(x,dx,N,D,J)
    lin = linV(x,dx,N,D)
    a = real(ein"ijk,ijk->"(J.coeff, lin[1,:,:,:]))[]
    b = real(ein"ijk,ijk->"(J.coeff, lin[2,:,:,:]))[]
    return [a,b]
end

function linearizationTest(B, ϵ, M, x0; E0 = 3, m = 1)
    ## Magnetic fields
    bint = x -> B(x, scalar = false)

    b = x -> B(x)

    ### getting solution
    
    u0, u, t = get_solution(bint, ϵ, x0, E0, m)

    ## ploting solution

    gr = gyroradius(ϵ, u0[4], b(u0[1:3]))

    get_plot(u,gr)

    ## Constructing J

    J = get_J(u0[4], b, ϵ, M)

    ## getting numerator and denominator

    N = x -> numN(J, x)
    D = x -> denomD(J, x)

    ## Geting derivatives and midpoint in poincaré section

    dxpc, xpc = dotX(u0, u[:, end], t)

    lin = Aj(xpc, dxpc, N, D, J)
    return lin, J, N, D, u, t
end
## State magnetic field and small parameter

# We are going to use the a simple magnetic field model

function mfield(r;scalar = true)
    B0 = 3
    if scalar
        return B0
    else
        return [0,0,B0]
    end
end



Bint = x -> mfield(x, scalar = false)

B = x -> mfield(x)

ϵ = 0.0001

M = tuple(31,31,31)

x0 = [2,3,0]

lin, J, N, D, u, t =linearizationTest(mfield,ϵ,M,x0)

## getting evolution function

FFO = get_FFO(ϵ, Bint)

## Getting first intersection with poincaré section

### Initial velocity in the x direction

### Define mass and initial energy

m, E = 1, 3

vx = (2* E/ m)^(0.5)

u0 = [2,3,0,vx, 0, 0]

### getting solution

u1, t1 = FFO(u0)

## gyrocenter

uc1 = u0 - [0,gyroradius(ϵ, vx, B(u0[1:3])),0,0,0,0]

uc300 = u0 - [0,gyroradius(ϵ, u1[4:6,50], B(u1[1:3,50])),0,0,0,0]

##
circx = gyroradius(ϵ, vx, B(u0[1:3]))*cos.(0:0.1:2π) .+ u0[1]

circy = gyroradius(ϵ, vx, B(u0[1:3]))*sin.(0:0.1:2π) .+ (u0[2]- gyroradius(ϵ, vx, B(u0[1:3])))


## ploting solution

f1 = Figure()

ax1 = Axis(f1[1, 1], xlabel = "x", ylabel = "y", title = "Solution")

lines!(ax1, u1[1, :], u1[2, : ])

scatter!(ax1, circx, circy, markersize = 5, color = :green)

scatter!(ax1, uc1[1], uc1[2], markersize = 20, color = :blue)
scatter!(ax1, uc300[1], uc300[2], markersize = 20, color = :red)

f1


## Constructing J

### Setting domain parameters
x1i, x1f = 0.0, 2π
x2i, x2f = 0.0, 2π
x3i, x3f = -5 + vx, 5 + vx

### Setting polynomial parameters
Nx, Ny, Nr = 31, 31, 31

J = slabJ(Nx, Ny, Nr, x1i, x1f, x2i, x2f, x3i, x3f, ϵ, B)

## Constructing N and D

N = x -> numN(J, x)
D = x -> denomD(J, x)

## Full-orbit to poincaré section

function FOtPS(x)
    return [x[1],x[2],x[4]]
end

X = FOtPS(u1[:, end])
X0 = FOtPS(u1[:, 1])
### Creating backward finite difference

function dotX(x0,x,t)
    dotx = (x - x0)./t
    midx = (x + x0)./2
    return dotx, midx
end
##
dX, midX = dotX(u0, u1[:, end], t1)
## Full-Orbit vs Adiabatic Invariant Dynamics

##
Aj(X,dX,N,D,J)

UU = dotXD(midX,dX,D)

a = real(ein"ijk,ijk->"(J.coeff, UU[1,:,:,:]))[]
b = real(ein"ijk,ijk->"(J.coeff, UU[2,:,:,:]))[]

##

real(ein"ijk,ijk->"(J.coeff,N(X)[2,:,:,:]))

real(ein"ijk,ijk->"(J.coeff,D(X)))


## Test for known B

### Let B equal to 1 + 0.5*cos(x)

function mfield1(r;scalar = true)
    if scalar
        return 1 + 0.5*cos(r[1])
    else
        return [0,0,1 + 0.5*cos(r[1])]
    end
end

ϵ1 = 0.01

M1 = tuple(31,31,31)

lin1, J1, N1, D1, u1, t1 =linearizationTest(mfield1,ϵ,M1,x0)



## Testing if the calculatng D an N first and then calculating Aij gives the same result as calculating Aij directly

function mfield2(r;scalar = true)
    if scalar
        return 1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])
    else
        return [0,0,1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])]
    end
end

ϵ2 = 0.01

M2 = tuple(31,31,31)

lin2, J2, N2, D2, u2, t2 =linearizationTest(mfield2,ϵ2,M2,x0)
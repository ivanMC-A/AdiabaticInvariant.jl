## Packages

using AdiabaticInvariant
using OMEinsum
using CairoMakie
using LinearAlgebra

## State magnetici field and small parameter

# We are going to use the same magnetic field as in the first paper
function mfield(r;scalar = true)
    if scalar
        return 1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])
    else
        return [0,0,1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])]
    end
end

Bint = x -> mfield(x, scalar = false)

B = x -> mfield(x)

ϵ = 0.1

## getting evolution function

FFO = get_FFO(ϵ, Bint)

## Getting first intersection with poincaré section

### Initial velocity in the x direction

### Define mass and initial energy

m, E = 1, 3

vx = 1/ϵ * (2* E/ m)^(0.5)

u0 = [0,0,0,vx, 0, 0]

### getting solution

u1, t1 = FFO(u0)

u2, t2 = FFO(u1[:,end])

u3, t3 = FFO(u2[:,end])

## ploting solution

f1 = Figure()

ax1 = Axis(f1[1, 1], xlabel = "x", ylabel = "y", title = "Solution")

lines!(ax1, u[1, :], u[2, : ])

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
dX, midX = dotX(X0, X, t1[end])
## Full-Orbit vs Adiabatic Invariant Dynamics


function FOvsAI(x,dx,N,D)
    XD = reshape(dx[1:2],2,1,1,1) .*reshape(D(x),1,31,31,31)
    return XD- N(x)
end

function Aij(x,dx,N,D,J)
    lin = FOvsAI(x,dx,N,D)
    a = real(ein"ijk,ijk->"(J.coeff, lin[1,:,:,:]))[]
    b = real(ein"ijk,ijk->"(J.coeff, lin[2,:,:,:]))[]
    return [a,b]
end

Aij(X,X,N,D,J)
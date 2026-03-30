using AdiabaticInvariant
using OMEinsum
using CairoMakie
using LinearAlgebra

# Define domain and parameters

ϵ = 0.1
x1i, x1f = 0.0, 2π
x2i, x2f = 0.0, 2π
x3i, x3f = -2 + sqrt(2)/ϵ, 0.1 + sqrt(2)/ϵ
Nx, Ny, Nr = 31, 31, 31

# Define Magnetic field

# Toy magnetic field

function mfield(r;scalar = true)
    if scalar
        return 1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])
    else
        return [0,0,1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])]
    end
end

Bz = x -> mfield(x)

J = slabJ(Nx, Ny, Nr, x1i, x1f, x2i, x2f, x3i, x3f, ϵ, Bz)

x0, y0, r0 = 0.0, 0.0, 10.0

∇J = deval(J, [x0, y0, r0])

hola = numN(J, [x0, y0, r0])

real(ein"ijk,ijk->"(J.coeff, hola[2]))

numNt(J, [x0, y0, r0])

denomD(J, [x0, y0, r0])
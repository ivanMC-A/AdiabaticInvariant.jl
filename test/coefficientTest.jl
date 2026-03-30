using AdiabaticInvariant
using LinearAlgebra
using OMEinsum
using CairoMakie

# 1D Fourier reconstruction test
# Test function: x^2

# Define domain and parameters

x_min, x_max = 0, 2π
N = 11
bd = tuple(x_min,x_max)

# Define functions to use


fpts = x -> AdiabaticInvariant.evenpts(x)       # Get even spaced points for fourier reconstruction

F = x -> AdiabaticInvariant.gTrans.(x, bd[1], bd[2]) # Get transformation function from fourier points to domain

InvF = x -> AdiabaticInvariant.gInv.(x, bd[1], bd[2]) # Get inverse transformation function from domain to fourier points

y1 = x -> cos(x)                                  # Define function to reconstruct


# Obtain points and transform them

x = fpts(N)                                      # get even spaced points
xt = F(x)                                        # transform points to domain

y1n = y1.(xt)                                    # get function values at transformed points

# Obtain Fourier matrix

ϕ1 = AdiabaticInvariant.fourA(N)
ϕ1T = ϕ1' 
W1 = real(Diagonal(ein"ik,kj->ij"(ϕ1T, ϕ1)))

a1 = W1\ein"ik,k->i"(ϕ1T, y1n)                   # get Fourier coefficients

y1r = x -> AdiabaticInvariant.evalFS(a1,x,bd, Real = false)

xs = range(0, 2π; length = 501)

f1 = Figure()

ax1 = Axis(f1[1, 1], xlabel = "time", ylabel = "error",
    title = "1D Fourier reconstruction of cos(x) error")

lines!(ax1, xs, abs.(y1.(xs) - real(y1r.(xs))), label = "Error")
f1

# Obtain Chebyshev matrix

cpts = x -> AdiabaticInvariant.chept1st(x)

G = x -> AdiabaticInvariant.gTrans.(x, bd[1], bd[2]; F = false)
InvG = x -> AdiabaticInvariant.gInv.(x, bd[1], bd[2]; F = false)

xc = cpts(N)

xct = G(xc)

C = AdiabaticInvariant.chebA(N)
CT = C'

Wc = real(Diagonal(ein"ik,kj->ij"(CT, C)))

y2 = x -> sin(x)

ycn = y2.(xct)

a2 = Wc\ein"ik,k->i"(CT, ycn)

yc = x -> AdiabaticInvariant.evalChev(a2,x,bd)

xs = range(0, 2π; length = 501)

f2 = Figure()

ax2 = Axis(f2[1, 1], xlabel = "time", ylabel = "error",
    title = "1D Chebyshev reconstruction of cos(x) error")

lines!(ax2, xs, abs.(y2.(xs) - yc.(xs)), label = "Error")
f2

# 2D Fourier by Chebyshev reconstruction test
# Test function: cos(x) + sin(y)

# Define domain and parameters

x_min, x_max = 0, 2π
y_min, y_max = 0, 2π
Nx = 31
Ny = 31
bd = tuple(x_min,x_max,y_min,y_max)

# Get even spaced points and transform them for fourier reconstruction

x3 = fpts(Nx)

x3t = F(x3)

# Get Chebyshev points of the first kind and transform them for Chebyshev reconstruction
y3 = cpts(Ny)
y3t = G(y3)

ϕ3 = AdiabaticInvariant.fourA(Nx)
ϕ3T = ϕ3'

ψ3 = AdiabaticInvariant.chebA(Ny)
ψ3T = ψ3'

W3f = real(Diagonal(ein"ik,kj->ij"(ϕ3T, ϕ3)))
W3c = real(Diagonal(ein"ik,kj->ij"(ψ3T, ψ3)))

y3f = (x,y) -> cos(x) + y^2
# y3 = (x,y) -> cos(x)
# y3 = (x,y) -> sin(y)

y3mn = [y3f(x,y) for x in x3t, y in y3t]

ψ3C = W3f\ein"im,mn->in"(ϕ3T,y3mn)

C = ein"ik,kj->ij"(ψ3C,ψ3) / W3c

vc = (N,y) ->AdiabaticInvariant.chebBase(N,InvG(y))

vf = (N,x) -> AdiabaticInvariant.fourBase(N,InvF(x))

x0 = π
y0 = 2

y3r = (x,y) -> real(ein"l,lk,k->"(vf(Nx,x),C,vc(Ny,y))[])
xs = range(0, 2π; length = 501)
ys = range(0, 2π; length = 501)
zs = [abs.(y3f.(x, y) .- y3r.(x, y)) for x in xs, y in ys]

f3 = Figure()
ax3 = Axis(f3[1, 1], xlabel = "x", ylabel = "y",
    title = "2D Fourier by Chebyshev reconstruction of cos(x) + y^2 error")
co = contourf!(ax3, xs, ys, zs, levels = 20)
Colorbar(f3[1, 2], co)
f3

maximum(abs(y3mn[a,b] - y3r(x3t[a], y3t[b]))
        for a in eachindex(x3t), b in eachindex(y3t))

# 3D Fourier by Fourier by Chebyshev reconstruction test
# Test function: cos(x) + sin(y) + cos(z)

x_min, x_max = 0, 2π
y_min, y_max = 0, 2π
z_min, z_max = 0, 2π
Nx = 31
Ny = 33
Nz = 35
N = [Nx, Ny, Nz]
bd = tuple(x_min,x_max,y_min,y_max,z_min,z_max)



# Get even spaced points and transform them for fourier reconstruction

x4 = fpts(Nx)
x4t = F(x4)

# Get even spaced points and transform them for fourier reconstruction

y4 = fpts(Ny)
y4t = F(y4)

# Get Chebyshev points of the first kind and transform them for Chebyshev reconstruction
z4 = cpts(Nz)
z4t = G(z4)

# Get Fourier polinomials evaluated at even spaced points in x
ϕ4 = AdiabaticInvariant.fourA(Nx)
ϕ4T = ϕ4'
W4fx = real(Diagonal(ein"ik,kj->ij"(ϕ4T, ϕ4)))

# Get Fourier polinomials evaluated at even spaced points in y
Θ4 = AdiabaticInvariant.fourA(Ny)
Θ4T = Θ4'
W4fy = real(Diagonal(ein"ik,kj->ij"(Θ4T, Θ4)))

# Get Chebyshev polinomials evaluated at Chebyshev points in z
Ψ4 = AdiabaticInvariant.chebA(Nz)
Ψ4T = Ψ4'
W4c = real(Diagonal(ein"ik,kj->ij"(Ψ4T, Ψ4)))

# Get values of function at transformed points

y4f = (x,y,z) -> cos(x) + sin(y) + cos(z)

y4ft = (r) -> cos(r[1]) + sin(r[2]) + cos(r[3])


y4mnp = [y4f(x,y,z) for x in x4t, y in y4t, z in z4t]

# Get Fourier by Fourier by Chebyshev coefficients

# multiply from the left by the inverse of ψ4

C4Φ4TΘ4T = ein"ji,im,mnp->jnp"(inv(W4fx),ϕ4T,y4mnp)

C4Φ4T = ein"mnp,pi,ij->mnj"(C4Φ4TΘ4T,Ψ4,inv(W4c))

C4 = ein"mnp,ni,ij->mjp"(C4Φ4T,Θ4,inv(W4fy))

C4test = coef_J(N,bd,y4ft;Lx = 2π, Ly = 2π)

# Get Fourier by Fourier by Chebyshev reconstruction
# Create vectors for basis functions

vx4 = (x) -> AdiabaticInvariant.fourBase(Nx, InvF(x))
vy4 = (y) -> conj(AdiabaticInvariant.fourBase(Ny, InvF(y)))
vz4 = (z) -> AdiabaticInvariant.chebBase(Nz, InvG(z))
y4r = (x,y,z) -> real(ein"i,ijk,j,k->"(vx4(x),C4,vy4(y),vz4(z))[])

x0, y0, z0 = 0, π/2, 0

yprint = round(y0, digits = 2)

Err4(x,y,z) = abs.(y4f(x,y,z) - y4r(x,y,z))

# For fixed z
Err4z(x,y) = Err4(x,y,z0)
# For fixed y
Err4y(x,z) = Err4(x,y0,z)
# For fixed x
Err4x(y,z) = Err4(x0,y,z)
xs = range(0, 2π; length = 501)
ys = range(0, 2π; length = 501)
zs = range(0, 2π; length = 501)
Errxy = [Err4z(x,y) for x in xs, y in ys]
Errxz = [Err4y(x,z) for x in xs, z in zs]
Erryz = [Err4x(y,z) for y in ys, z in zs]

fig = Figure()
ax = Axis(fig[1,1][1,1], title = "Error in x-y plane at z = $z0", xlabel = "x", ylabel = "y")
hm = heatmap!(ax, xs, ys, Errxy)
Colorbar(fig[1, 1][1, 2], hm)

ax = Axis(fig[1, 2][1, 1], title = "Error in x-z plane at y = $yprint", xlabel = "x", ylabel = "z")
hm = heatmap!(ax, xs, zs, Errxz)
Colorbar(fig[1, 2][1, 2], hm)

ax = Axis(fig[2, 2][1, 1], title = "Error in y-z plane at x = $x0", xlabel = "y", ylabel = "z")
hm = heatmap!(ax, ys, zs, Erryz)
Colorbar(fig[2, 2][1, 2], hm)

fig

f4 = Figure()

ax = Axis(f4[1, 1], title = "Error in x-z plane at y = yprint")
hm = heatmap!(ax, xs, zs, Errxz)
Colorbar(f4[1, 1][1, 2], hm)
f4

maximum(abs(y4mnp[a,b,c] - y4r(x4t[a], y4t[b], z4t[c]))
        for a in eachindex(x4t), b in eachindex(y4t), c in eachindex(z4t))

maximum(abs(ϕ4[a,i] - vx4(x4t[a])[i]) for a in eachindex(x4t), i in 1:Nx)
maximum(abs(Θ4[a,i] - vy4(y4t[a])[i]) for a in eachindex(y4t), i in 1:Ny)
maximum(abs(Ψ4[a,i] - vz4(z4t[a])[i]) for a in eachindex(z4t), i in 1:Nz)
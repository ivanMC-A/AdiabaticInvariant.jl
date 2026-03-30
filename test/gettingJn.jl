using AdiabaticInvariant
using CairoMakie
using LinearAlgebra
using Dates
using JLD2

function save_slabj(folder::String, base::String, s::AdiabaticInvariant.SlabJ, fig)

    mkpath(folder)

    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    basepath = joinpath(folder, "$(base)_$(timestamp)")

    # save figure
    save(basepath * ".png", fig)

    # save struct data (without the function)
    jldsave(basepath * ".jld2";
        N = s.N,
        bd = s.bd,
        ϵ = s.ϵ,
        coeff = s.coeff
    )

    println("Saved:")
    println(basepath * ".png")
    println(basepath * ".jld2")

end

folder = "test/gettingJnPlots"



# Define domain and parameters

ϵ = 0.1
x1i, x1f = 0.0, 2π
x2i, x2f = 0.0, 2π
x3i, x3f = -1 + sqrt(2)/ϵ, 1 + sqrt(2)/ϵ
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

x0, y0, r0 = 0.0, 0.0, 14.0

JJ = full_J(J.ϵ, J.B)

function Err(x,y,z,f::AdiabaticInvariant.SlabJ,g,pts; abstol = 1e-8)
    fv = evaluate(f,[x,y,z])
    gv = g([x,y,z])
    er = abs.(fv - gv)./ abs.(gv)
    if er > 1
        println("Warning: true value is close to zero, relative error may be large")
        push!(pts, [x,y,z])
        return 1
    else
        return er
    end
end

# Error for fixed x

ErrX = (y,r,pts) ->Err(x0,y,r,J,JJ,pts)

# Error for fixed y

ErrY = (x,r,pts) ->Err(x,y0,r,J,JJ,pts)

# Error for fixed r

ErrR = (x,y,pts) ->Err(x,y,r0,J,JJ,pts)

ptsx = []
ptsy = []
ptsz = []

xs = range(J.bd[1], J.bd[2]; length = 501)
ys = range(J.bd[3], J.bd[4]; length = 501)
rs = range(J.bd[5], J.bd[6]; length = 501)

Erryr = [ErrX(y,r,ptsx) for y in xs, r in rs]
Errxr = [ErrY(x,r,ptsy) for x in xs, r in rs]
Errxy = [ErrR(x,y,ptsz) for x in xs, y in ys]
## plots
fig = Figure()
ax = Axis(fig[1,1][1,1], title = "Error in x-y plane at r = $r0", xlabel = "x", ylabel = "y")
hm = heatmap!(ax, xs, ys, log10.(Errxy))
Colorbar(fig[1, 1][1, 2], hm)

ax = Axis(fig[1, 2][1, 1], title = "Error in x-r plane at y = $y0", xlabel = "x", ylabel = "r")
hm = heatmap!(ax, xs, rs, log10.(Errxr))
Colorbar(fig[1, 2][1, 2], hm)

ax = Axis(fig[2, 2][1, 1], title = "Error in y-r plane at x = $x0", xlabel = "y", ylabel = "r")
hm = heatmap!(ax, ys, rs, log10.(Erryr))
Colorbar(fig[2, 2][1, 2], hm)

fig
##
timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
filename = "image_$timestamp.png"

filepath = joinpath(folder, filename)
save(filepath, fig)
# Trying for a simpler funciton

x1it, x1ft = 0.0, 2π
x2it, x2ft = 0.0, 2π
x3it, x3ft = -20.0, 20.0
Nx, Ny, Nr = 31, 31, 31

y0 = π/2

F = r -> cos(r[1]) + sin(r[2]) + r[3]^2

FR = slabJ(Nx, Ny, Nr, x1it, x1ft, x2it, x2ft, x3it, x3ft, F)

xst = range(FR.bd[1], FR.bd[2]; length = 501)
yst = range(FR.bd[3], FR.bd[4]; length = 501)
rst = range(FR.bd[5], FR.bd[6]; length = 501)

# Error for fixed x

ErrXt(y,r,pts) = Err(x1it,y,r,FR,F,pts)

# Error for fixed y

ErrYt(x,r,pts) = Err(x,y0,r,FR,F,pts)

# Error for fixed r

ErrRt(x,y,pts) = Err(x,y,x3it,FR,F,pts)

ptsxt = []
ptsyt = []
ptszt = []

Erryrt = [ErrXt(y,r,ptsxt) for y in yst, r in rst]
Errxrt = [ErrYt(x,r,ptsyt) for x in xst, r in rst]
Errxyt = [ErrRt(x,y,ptszt) for x in xst, y in yst]

figt = Figure()
ax = Axis(figt[1,1][1,1], title = "Error in x-y plane at r = $x0[3]", xlabel = "x", ylabel = "y")
hm = heatmap!(ax, xst, yst, Errxyt)
Colorbar(figt[1, 1][1, 2], hm)

ax = Axis(figt[1, 2][1, 1], title = "Error in x-r plane at y = $x0[2]", xlabel = "x", ylabel = "r")
hm = heatmap!(ax, xst, rst, Errxrt)
Colorbar(figt[1, 2][1, 2], hm)

ax = Axis(figt[2, 2][1, 1], title = "Error in y-r plane at x = $x0[1]", xlabel = "y", ylabel = "r")
hm = heatmap!(ax, yst, rst, Erryrt)
Colorbar(figt[2, 2][1, 2], hm)

maximum(Errxyt)

figt

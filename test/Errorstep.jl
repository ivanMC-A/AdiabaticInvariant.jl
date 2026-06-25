## Packages

using AdiabaticInvariant
using OMEinsum
using CairoMakie
using LinearAlgebra

## functions

"""
    FOtPS(m)
    Extract the poincaré section coordinates from the full orbit coordinates.
"""

function FOtPS(m)
    return [m[1],m[2],m[4]]
end

"""
    get_J(r, B, ϵ, N)
    Construct the adiabatic invariant for the system of interest. 
"""

function get_J(r, B, ϵ, N::Tuple)

    ### Setting domain parameters
    x1i, x1f = 0.0, 2π
    x2i, x2f = 0.0, 2π
    x3i, x3f = -5 + r, 5 + r

    ### Setting polynomial parameters
    Nx, Ny, Nr = N[1], N[2], N[3]

    J = slabJ(Nx, Ny, Nr, x1i, x1f, x2i, x2f, x3i, x3f, ϵ, B)

    return J
end

""" 
    get_FFO_solution(FFO, x0; E0 = 3, m = 1, tmax = 1e3)
    Calculate the solution trajectory for FO up to N, where N is the number of intersections 
    with the poincaré section.
"""

function get_FFO_solution(FFO,x0,vₓ)
    ## Full orbit forward mapping
    
    # Initial conditions

    u0 = [x0[1], x0[2], x0[3], vₓ, 0, 0]

    # Evolving the system

    u, pspt,pst = FFO(u0)
    return u, pspt, pst
end

"""

    FFOvsFGC(pspt,pst,FGC)
    Compare the solution trajectories obtained from the FFO and FGC integrators.

"""

function FFOvsFGC(u,pspt,tf,FGC,p)
    ## FFO vs FGC

    # Initial conditions

    x0 = FOtPS(u[:,1])

    fof = FOtPS(pspt[1])

    # Evolving the system

    sol = FGC(x0,p; tmax = tf)
    # return norm(sol.u[end]-fof)

    # return norm(fof-sol.u[1])

    return norm(sol.u[end]-fof)/norm(fof-sol.u[1])
end

function FFOvsFGC(u,pspt,tf,FGC)
    ## FFO vs FGC

    # Initial conditions

    x0 = FOtPS(u[:,1])

    fof = FOtPS(pspt[1])

    # Evolving the system

    sol = FGC(x0; tmax = tf)
    # return norm(sol.u[end]-fof)

    # return norm(fof-sol.u[1])

    return norm(sol.u[end]-fof)/norm(fof-sol.u[1])
end

function rel_err(ϵ,vx,x0; pert = true)
    function mfield(r;scalar = true)
        if scalar
            return 1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])
        else
            return [0,0,1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])]
        end
    end

    Bint = x -> mfield(x, scalar = false)

    B = x -> mfield(x)

    if pert
        J = full_J(ϵ,B)
        FGC = get_FGC(J,B)
        p = [ϵ]
    else
        
        M = tuple(31,31,31)

        J = get_J(vx, B, ϵ, M)
        FGC = get_FGC(J)
    end
    # 

    N = 3

    FFO = get_FFO(ϵ, Bint, N)

    u, pspt, pst = get_FFO_solution(FFO,x0,vx)

    if pert
        FFOvsFGC(u, pspt, pst[1], FGC, p)
    else
        FFOvsFGC(u, pspt, pst[1], FGC)
    end
    
end

##


eps = 10 .^LinRange(-6,0,10)

err1 = zeros(length(eps))
err2 = zeros(length(eps))

vx = 0.1

x0 = [3,2,0]

for i in 1:length(eps)
    err1[i] = rel_err(eps[i],vx,x0; pert = false) 
    err2[i] = rel_err(eps[i],vx,x0)
end

## Plotting error

f = Figure()
ax1 = Axis(f[1, 1], xscale = log10, yscale = log10, xlabel = L"\epsilon", ylabel = L"\text{Relative error of FGC}",
    title = L"Screening of $\epsilon$")
lines!(ax1, eps, err1, label = "Interpolation")
lines!(ax1, eps,(1e-6./eps).^2. *10^-2)
lines!(ax1, eps,err2, label = "Perturbative J")
axislegend(ax1)
f

save("relative_error.pdf", f)

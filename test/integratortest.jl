using AdiabaticInvariant
using OMEinsum
using CairoMakie
using LinearAlgebra



"""
    kinetic_energy(u; mass=1.0, vel_indices)

Compute kinetic energy given state vector `u`.
`vel_indices` tells which components are velocities.
"""
function kinetic_energy(u; mass=1.0, vel_indices)
    v = u[vel_indices]
    return 0.5 * mass * sum( ϵ* v.^2)
end


"""
    energy_conservation(sol; mass=1.0, vel_indices)

Compute energy conservation error over a solution `sol`.
Returns:
- energies
- relative error array
"""
function energy_conservation(sol; mass=1.0, vel_indices)
    # Compute energy at each time step
    energies = [kinetic_energy(u; mass=mass, vel_indices=vel_indices) for u in sol.u]

    # Initial energy
    E0 = energies[1]

    # Relative error
    rel_error = abs.((energies .- E0) ./ E0)

    return energies, rel_error
end

""" 
    get_FFO_solution(FFO, x0; E0 = 3, m = 1, tmax = 1e3)
    Calculate the solution trajectory for FO up to N, where N is the number of intersections 
    with the poincaré section.
"""

function get_FFO_solution(FFO,x0, dx0)
    ## Full orbit forward mapping
    
    # Initial conditions

    u0 = [x0[1], x0[2], x0[3], dx0, 0, 0]

    # Evolving the system

    u, pspt,pst = FFO(u0)
    return u, pspt, pst
end

## Initialization of J
ϵ = 0.001
vx = 0.1

function mfield(r;scalar = true)
    if scalar
        return 1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])
    else
        return [0,0,1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])]
    end
end

Bint = x -> mfield(x, scalar = false)

B = x -> mfield(x)

FFO = get_FFO(ϵ, Bint)

x0 = [1,1,0]

u, pspt, pst = get_FFO_solution(FFO,x0,vx)

energies, err = energy_conservation(u; mass=1.0, vel_indices=4:6)

fig = Figure()
ax = Axis(fig[1,1], title="Energy Conservation")

lines!(ax, u.t, err)

fig

save("energyerror.pdf", fig)

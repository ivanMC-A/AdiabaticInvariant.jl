### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 3f2677b0-ca9c-42ae-ac5d-e83e98b020f8
begin
    import Pkg
    Pkg.activate(Base.current_project())

    using AdiabaticInvariant
end

# ╔═╡ 8bbf64f8-4d53-11f1-20d5-77298c3b27b9
begin
	using OMEinsum
	using CairoMakie
	using LinearAlgebra
end

# ╔═╡ 9c52d7d5-d01f-4f9f-a02d-c60828db2b53
md"""
The task at hand is to be able to evaluate `` \mathcal J_\Sigma`` along a streamline of ``\Sigma`` and check if it's invariant, i.e., the relative error of ``\mathcal J_\Sigma`` is 0. To do this, we start by copying a code that I already have to integrate GC trayectories (here, GC refers to the nonperturbative guiding center flow but using the perturbative ``\mathcal J`` for ``\epsilon \ll 1``).
"""

# ╔═╡ e931c995-6934-4654-b809-e5a0c356982f
md"""
Some essential functions are:
- FOtPS: Sends points in ``Z`` to points in ``\Sigma``.
- get\_J: Construct the adiabatic invariant for the system of interest.
- get\_FFO\_solution: Calculate the solution trajectory for FO up to N, where N is the number of intersections with the poincaré section.
"""

# ╔═╡ 32e48647-c26b-4a25-af1e-3510765cc1a7
 begin
	 function FOtPS(m)
    	return [m[1],m[2],m[4]]
	end
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
	 function get_FFO_solution(FFO,x0,vₓ)
	    ## Full orbit forward mapping
	    
	    # Initial conditions
	
	    u0 = [x0[1], x0[2], x0[3], vₓ, 0, 0]
	
	    # Evolving the system
	
	    u, pspt,pst = FFO(u0)
	    return u, pspt, pst
	end
 end

# ╔═╡ 57cbb4a8-e94d-4f6a-9d39-1e94be3ef29c
md"""
## Initialization of ``\mathcal J``
"""

# ╔═╡ 50271395-7f8b-4b4b-858d-da5969e1f7c8
begin
	ϵ = 0.1
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
	
	M = tuple(31,31,31)
	
	J = get_J(vx, B, ϵ, M)
end

# ╔═╡ 85c4ff4b-7ac1-444e-bf8a-fe192906ded6
md"""
## Obtaining first point of intersection w/ ``\Sigma``
Here ``x_0`` is the initial position, and ``N`` is the number of intersection points with ``\Sigma`` we are asking for. **pspt** are the points of intersection and **pst** is the time of intersection of those points. Saving the time and position is essential for integration of FGC and comparison of FGC and FFO.
"""

# ╔═╡ 88cbdbd3-feb3-4b65-949f-cddef4325c69
begin
	x0 = [3,2,0]
	
	N = 3
	
	FFO = get_FFO(ϵ, Bint, N)
	
	u, pspt, pst = get_FFO_solution(FFO,x0,vx)
end

# ╔═╡ fc2d5462-462b-4bec-a384-249e7127e140
pspt

# ╔═╡ 993b175e-4015-4d7d-bad3-0ce6d50c38f7
pst

# ╔═╡ 1ec31a45-c63f-470a-944d-2d08124453c1
md"""
From this information, it feels like the particle is drifting to the right and up every gyroperiod.
"""

# ╔═╡ f0a06083-d787-43ed-a3e4-a52eabfe59bf
md"""
## Obtaining FGC trajectories
"""

# ╔═╡ 95cc26bb-1ed8-470d-b3b2-ce2839c5e566
begin
	FGC = get_FGC(J)

	x00 = [3,2,vx]
	
	sol = FGC(x00; tmax = pst[3])
end

# ╔═╡ 6828e9b9-21df-4eb2-87c7-87034639bae9
md"""
For the FGC trajectories, it feels like the particle is drifting down. Maybe there is a typo in the sign of the FGC DE.
"""

# ╔═╡ 61801143-602d-4e14-b3de-b71f833610e0
md"""
Now that I have the points, I can evaluate ``\mathcal J`` along ``\sigma``. For this, I will create a function that spits back the relative error of ``\mathcal J``.
"""

# ╔═╡ ba5a4e5b-73ad-47fc-87aa-a7a09643ce2d
begin
	function relErrJ(J,x::AbstractArray)
		n = length(x)
		J0 = evaluate(J,x[1])
		err = zeros(n)
		for i in 1:n
			err[i] = abs((evaluate(J,x[i])-J0)/J0)
		end
		return err
	end
end

# ╔═╡ bc8b603c-73cf-462e-8468-74c38136acce
err = relErrJ(J,sol.u)

# ╔═╡ 107517de-a9a9-47a8-b2ba-ecdcd26816f4
begin
	f = Figure()
	ax1 = Axis(f[1, 1], xlabel = L"\text{Time}", ylabel = L"\text{Error}",
    title = L"Relative error of $\mathcal J_\Sigma$ on $\sigma(t)$")
	lines!(ax1, sol.t, err)
	f
end

# ╔═╡ Cell order:
# ╠═3f2677b0-ca9c-42ae-ac5d-e83e98b020f8
# ╠═8bbf64f8-4d53-11f1-20d5-77298c3b27b9
# ╟─9c52d7d5-d01f-4f9f-a02d-c60828db2b53
# ╟─e931c995-6934-4654-b809-e5a0c356982f
# ╠═32e48647-c26b-4a25-af1e-3510765cc1a7
# ╟─57cbb4a8-e94d-4f6a-9d39-1e94be3ef29c
# ╠═50271395-7f8b-4b4b-858d-da5969e1f7c8
# ╟─85c4ff4b-7ac1-444e-bf8a-fe192906ded6
# ╠═88cbdbd3-feb3-4b65-949f-cddef4325c69
# ╠═fc2d5462-462b-4bec-a384-249e7127e140
# ╠═993b175e-4015-4d7d-bad3-0ce6d50c38f7
# ╟─1ec31a45-c63f-470a-944d-2d08124453c1
# ╟─f0a06083-d787-43ed-a3e4-a52eabfe59bf
# ╠═95cc26bb-1ed8-470d-b3b2-ce2839c5e566
# ╟─6828e9b9-21df-4eb2-87c7-87034639bae9
# ╟─61801143-602d-4e14-b3de-b71f833610e0
# ╠═ba5a4e5b-73ad-47fc-87aa-a7a09643ce2d
# ╠═bc8b603c-73cf-462e-8468-74c38136acce
# ╠═107517de-a9a9-47a8-b2ba-ecdcd26816f4

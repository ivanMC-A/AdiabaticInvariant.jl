### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ a6dbc6e7-3849-4e22-be02-a3b9f0eab166
begin
    import Pkg
    Pkg.activate(Base.current_project())

    using AdiabaticInvariant
	using OMEinsum
	using CairoMakie
	using LinearAlgebra
end

# ╔═╡ e7c2f7aa-5551-11f1-03a4-61176487c04e
md"""
# Absolute and Relative error of the Nonperturbative Guiding Center Model

The goal of this notebook is to study the error of the Nonperturbative Guiding Center Model. We start by defining some relevant functions. Let ``\phi_{\text{FO}}:\Omega \times \mathbb R \rightarrow \Omega`` and ``\phi_{\text{GC}}:\Sigma \times \mathbb R \rightarrow \Sigma`` be the flow
operators of full orbit and guiding center respectively, where ``\Omega`` is the domain and ``\Sigma \subset \Omega`` is the Poincaré Section. Then, suppose there exists a minimum time function on the Poincaré section ``T: \Sigma \rightarrow \mathbb R`` so that for ``x = \phi_{\mathrm{FO}}(x,0) \in \Sigma`` , ``\phi_{\mathrm{FO}}(x,T(x)) \in \Sigma``. We define the Poincare maps ``F_j : (x) \mapsto (\phi_j(x,T(x))``, where ``j \in \{\mathrm{FO},\mathrm{GC}\}``. Ideally, we want to find a ``\mathcal J`` so that we have a **relative error**,

```math
E_{\text{1-step}}[\mathcal J] =\frac{ \lVert F_{\mathrm{GC}} - F_{\mathrm{FO}} \rVert_\nu }{ \lVert F_{\mathrm{FO}} - \mathrm{id} \rVert_\nu } << 1,     
```
and an **absolute error**,
```math
E_{\text{abs},\text{1-step}}[\mathcal J] =\lVert F_{\mathrm{GC}} - F_{\mathrm{FO}} \rVert_\nu << 1     
```
"""

# ╔═╡ d5a722fa-5847-4b76-a756-c689d6692dee
md"""
Now we define relevant julia functions.
"""

# ╔═╡ 91f9b88d-735f-4c50-9903-4e27f497ef05
begin
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
	Calculate the solution trajectory for FO up to N, where N is the number of intersections with the poincaré section.
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
	
	    FFOvsFGC(u,pspt,tf,FGC,p)
	Inputs:
		-u: solution from FFO
		-pspt: Point from u at poincaré section.
		-tf: time for integration of FGC. It corresponds to a time where a point of pspt happened.
		-FGC: GC Solver.
		-p: relevant parameters.
	Compare the solution trajectories obtained from the FFO and FGC integrators. Since this function accepts a parameter p, it means that it will take an FGC that works with the analytic J
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

	"""
	    FFOvsFGC(u,pspt,tf,FGC)
	Inputs:
		-u: solution from FFO
		-pspt: Point from u at poincaré section.
		-tf: time for integration of FGC. It corresponds to a time where a point of pspt happened.
		-FGC: GC Solver.
		-p: relevant parameters.
	Compare the solution trajectories obtained from the FFO and FGC integrators. This function works with the approximated J.
	
	"""
	
	function FFOvsFGC(u,pspt,tf,FGC)
	    ## FFO vs FGC
	
	    # Initial conditions
	
	    x0 = FOtPS(u[:,1])
	
	    fof = FOtPS(pspt[1])
	
	    # Evolving the system
	
	    sol = FGC(x0; tmax = tf)
	    # return norm(sol.u[end]-fof)
	
	    # return norm(fof-sol.u[1])
		E_abs = norm(sol.u[end]-fof)
		E_rel = E_abs/norm(fof-sol.u[1])
	
	    return E_abs,E_rel
	end

	"""
		Es(ϵ,vx,x0;pert = true)
	This function condenses the process of calculating the absolute and relative error for the following inputs:
		-ϵ: magnetization of the charged particle.
		-vx: initial velocity of the particle in the x direction.
		-x0: initial position of the particle.
		-pert: tells if you are working with the truncated or approximated J.
		-M: tuple contains the number of modes for the approximated J.

	"""
	
	function Es(ϵ,vx,x0; pert = true, M=(31,31,31))
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
end

# ╔═╡ 590c0030-87da-4081-b0e0-1b96abc57e80
md"""
We proceed to set the values for our inputs:
"""

# ╔═╡ 69704909-1255-4486-b1ab-f193e25407e4
begin
	# Values of ϵ
	
	ϵs = 10 .^LinRange(-6,0,10)


	# Creating arrays to store errors
	
	E_abs_approx = zeros(length(ϵs))
	E_rel_approx = zeros(length(ϵs))

	E_abs_pert = zeros(length(ϵs))
	E_rel_pert = zeros(length(ϵs))

	# Initial velocity
	
	vx = 0.1

	# Initial position
	
	x0 = [3,2,0]
end

# ╔═╡ f9e8dd97-6506-4ade-9d17-f554802aa45e
for i in 1:length(ϵs)
    E_abs_approx[i], E_rel_approx[i] = Es(ϵs[i],vx,x0; pert = false) 
end

# ╔═╡ ad289d64-044b-479c-8deb-53fbf760860f
begin
	f = Figure()
	ax1 = Axis(f[1, 1], xscale = log10, yscale = log10, xlabel = L"\epsilon", ylabel = L"\text{Absolute error of FGC}",
	    title = L"Screening of $\epsilon$")
	lines!(ax1, ϵs, E_abs_approx, label = "Interpolation")
	axislegend(ax1)
	f
end

# ╔═╡ f89c0569-23b4-42e7-a7ac-db2b3d6bb730
begin
	f2 = Figure()
	ax2 = Axis(f2[1, 1], xscale = log10, yscale = log10, xlabel = L"\epsilon", ylabel = L"\text{Relative error of FGC}",
	    title = L"Screening of $\epsilon$")
	lines!(ax2, ϵs, E_rel_approx, label = "Interpolation")
	lines!(ax2, ϵs,(1e-6./ϵs).^2. *10^-2)
	axislegend(ax2)
	f2
end

# ╔═╡ Cell order:
# ╟─e7c2f7aa-5551-11f1-03a4-61176487c04e
# ╠═a6dbc6e7-3849-4e22-be02-a3b9f0eab166
# ╟─d5a722fa-5847-4b76-a756-c689d6692dee
# ╠═91f9b88d-735f-4c50-9903-4e27f497ef05
# ╟─590c0030-87da-4081-b0e0-1b96abc57e80
# ╠═69704909-1255-4486-b1ab-f193e25407e4
# ╠═f9e8dd97-6506-4ade-9d17-f554802aa45e
# ╠═ad289d64-044b-479c-8deb-53fbf760860f
# ╠═f89c0569-23b4-42e7-a7ac-db2b3d6bb730

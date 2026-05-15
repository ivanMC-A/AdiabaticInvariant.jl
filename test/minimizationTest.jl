## Packages

using AdiabaticInvariant
using OMEinsum
using CairoMakie
using LinearAlgebra

## funcitons

"""
    get_J(dxdt, B, ϵ, N)
    Construct the adiabatic invariant for the system of interest. 
"""

function get_J(dxdt, B, ϵ, N::Tuple)

    ### Setting domain parameters
    x1i, x1f = 0.0, 2π
    x2i, x2f = 0.0, 2π
    x3i, x3f = -5 + dxdt, 5 + dxdt

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



function dotX_all(points, times)
    N = length(points) - 1

    dots = Vector{Vector{Float64}}(undef, N)
    mids = Vector{Vector{Float64}}(undef, N)

    for i in 1:N
        if i == 1
            t = times[i]
        else
            t = times[i] - times[i-1]
        end
        dots[i], mids[i] = dotX(points[i], points[i+1], t)
    end

    return dots, mids
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
    function build_A(N,D,X,DX)
        A = []
        for i in getindex(X)
            push!(A, linV(X[i], DX[i], N, D))
        end
        return A
    end
    
    Build the matrix A of vectors a = DX D - N
"""
function build_A(N,D,X,DX)
    A = []
    n = length(X) - 1
    for i in 1:n
        push!(A, linV(X[i], DX[i], N, D))
    end
    return A
end

function build_G(D,X)
    G = []
    n = length(X) - 1
    for i in 1:n
        push!(G, D(X[i]))
    end
    return G
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

function get_A(N,D,xs,ts)
    DX, X = dotX_all(xs, ts)
    A = build_A(N,D,X,DX)
    B = stack(A)
    return B
end

"""
    get_plot(u)
    Plot the solution trajectory in the x-y plane.
"""

function get_plot(u,B)

    xs = 0:0.01:2π
    ys = 0:0.01:2π
    zs = [B([i,j]) for i in xs, j in ys]

    x = [v[1] for v in u]
    y = [v[2] for v in u]

    f = Figure()

    ax = Axis(f[1, 1], xlabel = "x", ylabel = "y", title = "Solution")

    scatter!(ax, x, y)

    contour!(ax, xs, ys, zs)

    display(f)

end

""" 
    get_solution(ϵ, B, N, u0; E0 = 3, m = 1, tmax = 1e3)
    Calculate the solution trajectory for FO up to N, where N is the number of intersections 
    with the poincaré section.
"""

function get_solution(ϵ, B, N, u0; E0 = 3, m = 1, tmax = 1e3)
    ## Full orbit forward mapping

    FFO = get_FFO(ϵ, B, N, tmax = tmax)
    
    # Initial conditions

    dx0 = (2* E0/ m)^(0.5)

    u0 = [u0[1], u0[2], u0[3], dx0, 0, 0]

    # Evolving the system

    u, pspt,pst = FFO(u0)
    return u, pspt, pst, FFO
end

## Testing poincaré points generator

ϵ = 0.1

function mfield(r;scalar = true)
    if scalar
        return 1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])
    else
        return [0,0,1 + 0.03*cos(3*r[1] + 1*r[2]) + 0.03*cos(1*r[1] + 3*r[2])]
    end
end

Bint = x -> mfield(x, scalar = false)

B = x -> mfield(x)

N = 1000

M = tuple(31,31,31)

x0 = [3.6,2,0]

u, pspt, pst, FFO = get_solution(ϵ, Bint, N,x0, tmax = 1e10)

get_plot(pspt,B)

J = get_J(u[4,1],B,ϵ,M)



NN = x -> numN(J,x)
D = x ->denomD(J,x)

A = get_A(NN,D,pspt,pst)

real(ein"ijk,lijkm->lm"(J.coeff, A))
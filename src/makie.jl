import Makie: plot, plot!, scatter, scatter!

@info "Loading some plot recipes for Makie"

# some plotting utilities
export plotif, trimplot, signchart
export arrow!, vectorfieldplot!
export newton_vis

###

#  plot "recipe" for functions


function Makie.plot(f::Function, a::Number, b::Number, args...; kwargs...)
    xs = range(float(a), float(b), length=500)
    Makie.plot(xs, f.(xs), args...; kwargs...)
end

function  Makie.plot!(scene::Makie.AbstractPlotting.Scene, f::Function, args...; kwargs...)
    rect = scene.data_limits[]
    a, b = rect.origin[1],  rect.origin[1] + rect.widths[1]
    xs = range(a, b, length=500)
    Makie.plot!(scene, xs, f.(xs), args...;  kwargs...)
end

Makie.plot!(f::Function, args...; kwargs...) = Makie.plot!(Makie.AbstractPlotting.current_scene(), f, args...; kwargs...)


#  plot "recipe" for parametric functions
function Makie.plot(f::Function, g::Function, a::Number, b::Number, args...; kwargs...)
    xs = range(float(a), float(b), length=500)
    Makie.lines(f.(xs), g.(xs), args...; kwargs...)
end

function  Makie.plot!(scene::Makie.AbstractPlotting.Scene, f::Function, g::Function, a::Number, b::Number, args...; kwargs...)
    xs = range(a, stop=b, length=500)
    Makie.plot!(scene, f.(xs), g.(xs), args...;  kwargs...)
end
Makie.plot!(f::Function, g::Function, a::Number, b::Number, args...; kwargs...) =
    Makie.plot!(Makie.AbstractPlotting.current_scene(), f, g, a, b, args...; kwargs...)



#  plot reciple  for SymPy objects
Makie.plot(ex::Sym, a::Number, b::Number; kwargs...) = Makie.plot(lambdify(ex), a, b; kwargs...)
Makie.plot!(ex::Sym; kwargs...) = Makie.plot!(lambdify(ex); kwargs...)

Makie.plot(ex1::Sym, ex2::Sym,  a::Number, b::Number; kwargs...) =
    Makie.plot(lambdify(ex1), lambdify(ex2), a, b; kwargs...)
Makie.plot!(ex1::Sym, ex2::Sym,  a::Number, b::Number; kwargs...) =
    Makie.plot!(lambdify(ex1), lambdify(ex2), a,b; kwargs...)


# Plot of tuple
# d= 1
# d>1 parametric
function Makie.plot(fs::Tuple, a::Number, b::Number; kwargs...)
    xs = range(a, stop=b, length=500)
    ys = [ [f(x) for x in xs] for f in fs]
    if length(ys) == 1
        Makie.lines(xs, ys[1]; kwargs...)
    else
        Makie.lines(ys...; kwargs...)
    end
end

function Makie.plot!(fs::Tuple, args...; kwargs...)
    xs = range(args[1], stop=args[2], length=500)
    ys = [ [f(x) for x in xs] for f in fs]
    Makie.plot!(ys...; kwargs...)
end

##
## --------------------------------------------------
##

"""
   trimplot(f, a, b, c=20; kwargs...)

Plot f over [a,b] but break graph if it exceeds c in absolute value.
"""
function trimplot(f, a, b, c=20; kwargs...)
    plot(x -> abs(f(x)) <= c ? f(x) : NaN, a, b; kwargs...)
end




"""
    plotif(f, g, a, b)

Plot f colored depending on g >= 0 or not.
"""
function plotif(f, g, a, b, args...; colors=(:red, :blue), linewidth=5, kwargs... )


    xs = a:(b-a)/251:b
    zs = f.(xs)
    p = plot(xs, f.(xs), args...; color=colors[1], linewidth=linewidth, kwargs...)

    ys = g.(xs)
    ys[ys .< 0] .= NaN

    us,vs = Float64[], Float64[]
    for (i,y) in enumerate(ys)
        if isnan(y)
            if length(vs) > 1
                plot!(p, us, vs, color=colors[2], linewidth=5)
            end
            empty!(us)
            empty!(vs)
        else
            push!(us, xs[i])
            push!(vs, zs[i])
        end
    end
    if length(vs) > 1
        plot!(p, us, vs, color=colors[2], linewidth=5)
    end
    p
end

"""
   signchart(f, a, b)

Plot f over a,b with different color when negative.
"""
function signchart(f, a, b)
    p = plotif(f, f, a, b)
    plot!(p, zero)
    p
end



"""
   `arrow!(p, v)`

Add the vector `v` to the plot anchored at `p`.

This would just be a call to `quiver`, but there is no 3-D version of that. As well, the syntax for quiver is a bit awkward for plotting just a single arrow. (Though efficient if plotting many).

```
using Makie
r(t) = [sin(t), cos(t), t]
rp(t) = [cos(t), -sin(t), 1]
lines(unzip(r, 0, 2pi)...)
t0 = 1
arrow!(r(t0), r'(t0))
```
"""
function arrow!(scene::Makie.AbstractPlotting.Scene, p, v; kwargs...)
    P = length(p) == 3 ? Point3f0 : Point2f0
    Makie.arrows!(P.([p]), P.([v]); kwargs...)
end
arrow!(p,v;kwargs...) = arrow!(Makie.AbstractPlotting.current_scene(), p, v; kwargs...)

"""

    vectorfieldplot(V; xlim=(-5,5), ylim=(-5,5), n=10; kwargs...)

V is a function that takes a point and returns a vector (2D dimensions), such as `V(x) = x[1]^2 + x[2]^2`.

The grid `xlim × ylim` is paritioned into (n+1) × (n+1) points. At each point, `pt`, a vector proportional to `V(pt)` is drawn.

This is written to add to an existing plot.

```
plot()  # make a plot
V(x,y) = [x, y-x]
vectorfield_plot!(p, V)
p
```
"""
function vectorfieldplot( V; xlim=(-5,5), ylim=(-5,5), n=10, kwargs...)

    V′(x) = Point2f0(V(x...)...)
    Makie.streamplot(V′, xlim[1]..xlim[2], ylim[1]..ylim[2]; kwargs...)
end
function vectorfieldplot!(scene::Makie.AbstractPlotting.Scene, V; xlim=(-5,5), ylim=(-5,5), n=10, kwargs...)

    V′(x) = Point2f0(V(x...)...)
    Makie.streamplot!(scene, V′, xlim[1]..xlim[2], ylim[1]..ylim[2]; kwargs...)
end
vectorfieldplot!(V; kwargs...) =
    vectorfieldplot!(Makie.AbstractPlotting.current_scene(), V; kwargs...)

## --------------------------------------------------

##
# visualize newtons method
function newton_vis(f, x0, a=Inf,b=-Inf; steps=5, kwargs...)
    xs = Float64[x0]
    for i in 1:steps
        xᵢ = xs[end]
        xᵢ₊₁ = xᵢ - f(xᵢ)/f'(xᵢ)
        push!(xs, xᵢ₊₁)
    end

    m,M = extrema(xs)
    m = min(m, a)
    M = max(M, b)

    p = plot(f, m, M; linewidth=3,  kwargs...)
    plot!([m,M], [0,0])
    plot!(p, zero)
    for i in 1:steps
        plot!(p, [xs[i],xs[i],xs[i+1]], [0,f(xs[i]), 0])
        scatter!(p, xs[i:i],[0], markersize=1)
    end
    scatter!(p, [xs[steps+1]], [0], markersize=1)
    p
end

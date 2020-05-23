#There are a collection of "demos". Run `?visualizations` for a description.
# If can't load into JuliaBox, use this:
#  include(download("https://raw.githubusercontent.com/mth229/MTH229.jl/master/src/MTH229.jl"))
#
# * Run the command `?visualizations` for a description of some interactive features.


"""

`MTH229`: helper functions for using `Julia` with MTH229

This module does two things:

* Install other useful packages with one command (`SimplePlots`, `Roots`,  `ForwardDiff`, `QuadGK`, `SpecialFunctions`, ...)

* Add a number of helper functions.

The helper functions include:

- `secant(f, a, b)`: return a function giving the secant line between ``(a,f(a))`` and ``(b,f(b))``.

- `tangent(f, c)`:  return a function giving the tangent line to ``f(x)`` at the point ``(c,f(c))``.

- `bisection(f, a, b)`: A simple implementation of the bisection
  method. The interval ``[a,b]`` should be a bracketing interval. This
  function makes an illustrative graphic. For real use of the bisection method, the `fzero(f,
  a, b)` function, from the `Roots` package, should be used.

- `adjoint`: This allows the derivative of   a function to be found as with math notation: `f'`.  The notation can be used for higher-order derivatives too: `f''`, `f'''`, ... This uses automatic differentiation from the `ForwardDiff` package.

- `plotif(f, g, a, b)`: Plot the function `f` over the interval `[a,b]` and color differently where ``g(x) > 0`` over ``[a,b]``. By passing in `f` for `g` shows where `f` is positive on `[a,b]`; passing in `f'` shows where `f` is increasing on `[a,b]`; and passing in `f''` shows where `f` is concave up on `[a,b]`.


- `riemann(f, a, b, n; method="right")` An implementation of Riemann sums. The method can be "right" or "left" for Riemann sums, or "trapezoid" or "simpsons" for related approximations.


The package `SimplePlots` provides  a quick to load plottinig  package  with a  syntax  very similar  to the  more  feature   rich  `Plots` package. 

"""
module MTH229


# msg = """
# Loading the `MTH229` package.

# * Run the command `?MTH229` for a short description.


# """

# @info msg



using Reexport
@reexport using SimplePlots
@reexport using Roots
@reexport using SpecialFunctions
@reexport using SymPy
@reexport using QuadGK
@reexport using LinearAlgebra
@reexport using Base.MathConstants

import ForwardDiff


# #using Interact


###
export tangent, secant, D, grad
export lim,  bisection, riemann, fubini
export plotif, trimplot, signchart, newton_vis
export uvec, xs_ys, unzip, arrow!, vectorfieldplot!


" f'(x) will find the derivative of `f` using Automatic Differentation from the `ForwardDiff` package "
Base.adjoint(f::Function) = x -> ForwardDiff.derivative(f, float(x))
D(f, n=1) = n > 1 ? D(D(f), n-1) : x -> ForwardDiff.derivative(f, float(x))
grad(f) = (x, xs...) -> ForwardDiff.gradient(f, vcat(x, xs...))
# could aslo wrap as function ForwardDiff.jacobian, ForwardDiff.hessian


"""
Returns a function describing the tangent line to the graph of f at x=c.

Example. Where does the tangent line intersect the y axis?
```
f(x) = sin(x)
tl(x) = tangent(f, pi/4)(x)  # or tl = tangent(f, pi/3) to use a non-generic function
tl(0)
```

Uses the automatic derivative of `f` to find the slope of the tangent line at `x=c`.

"""
tangent(f,c) = x -> f(c) + f'(c) * (x-c)

"""
Returns a function describing the secant line to the graph of f at x=a and x=b.

Example. Where does the secant line intersect the y axis?
```
f(x) = sin(x)
a, b = pi/4, pi/3
sl(x) = secant(f, a, b)(x)  # or sl = sl(f, a, b) to use a non-generic function
sl(0)
```


"""
secant(f, a, b) = x -> f(a) + (f(b) - f(a)) / (b-a) * (x - a)


"""

`lim(f, c, n, dir="+")`: means to generate numeric table of values of `f` as `h` gets close to `c`.

Example:
```
f(x) = sin(x) / x
lim(f, 0)
```
"""
function lim(f::Function, c::Real; n::Int=6, dir="+")
	 hs = [(1/10)^i for i in 1:n] # close to 0
	 if dir == "+"
	   xs = c .+ hs
	 else
	   xs = c .- hs
	 end
	 ys = map(f, xs)
	 [xs ys]
end


"""

Simple implementation of the bisection method.

Example:

```julia
bisection(sin, 3, 4)
f(x) = x^5 - x^4 - x^3 - x^2 - x - 1
a = bisection(f, 1, 2)
f(a)
```

The display shows a simple graphic illustrating the method's division for the first few steps.

An easier-to-understand alternative to `Roots.find_zero(f, (a,b), Bisection())`.


"""
function bisection(f::Function, a, b)
    a,b = sort([a,b])

    if f(a) * f(b) > 0
        error("[a,b] is not a bracket. A bracket means f(a) and f(b) have different signs!")
    end

    M = a + (b-a) / 2


    i, j = 0, 64
    ss = fill("#", 65)
    ss[i+1]="a"; ss[j+1]="b"
    println("")
    println(join(ss))
    flag = true

    while a < M < b
        if flag && j-i == 1
            ss = fill(" ", 65)
            ss[j:(j+1)] .= "⋮"
            println(join(ss))
            println("")
            flag = false
        end


        if f(M) == 0.0
            println("... exact answer found ...")
	    break
        end
        ## update step
	if f(a) * f(M) < 0
	    a, b = a, M

            if flag
                j = div(i + j, 2)
            end


	else
	    a, b = M, b

            if flag
                i = div(i + j, 2)
            end

	end

        if flag
            ss = fill(".", 65)
            ss[i+1]="a"; ss[j+1]="b"; ss[(i+2):j] .= "#"
            println(join(ss))
        end

        M = a + (b-a) / 2
    end
    M
end

newton(f, fp, x0; kwargs...) = Roots.find_zero((f,fp), x0, Roots.Newton(); kwargs...)
newton(f, x0; kwargs...) = newton(f, D(f), x0; kwargs...)
export newton

#fzero(f, x0; kwargs...) = Roots.find_zero(f, x0; kwargs...)
#fzero(f, a, b; kwargs...) = Roots.find_zero(f, (a, b); kwargs...)
#fzeros(f, a, b; kwargs...) = Roots.find_zeros(f, a, b; kwargs...)


# some plotting utilities

#  plot recipe for functions
function SimplePlots.plot(f::Function, a::Number, b::Number, args...; kwargs...)
    xs = range(float(a), float(b), length=500)
    plot(xs, f.(xs), args...; kwargs...)
end

function  SimplePlots.plot!(plt::SimplePlots.SimplePlot, f::Function, args...; kwargs...)
    xs = plt.data[1]["x"]
    plot!(plt, xs, f.(xs), args...;  kwargs...)
end
SimplePlots.plot!(f::Function, args...; kwargs...) = SimplePlots.plot!(SimplePlots._plot, f, args...; kwargs...)

#  plot reciple  for SymPy objects
SimplePlots.plot(ex::Sym, a::Number, b::Number,  args...; kwargs...) = plot(lambdify(ex), a, b, args...; kwargs...)
SimplePlots.plot!(ex::Sym, args...; kwargs...)    = plot!(lambdify(ex), args..., kwargs...)

##
## --------------------------------------------------
##

"""
   trimplot(f, a, b, c=20; kwargs...)

Plot f over [a,b] but break graph if it exceeds c in absolute value.
"""
function trimplot(f, a, b, c=20; kwargs...)
  xs = range(a, stop=b, length=251)
  ys = f.(xs)

  us, vs = Real[NaN], Real[NaN]
  p = plot(us, vs, xlim=(a, b), legend=false, kwargs...)
  for (x,y) in zip(xs, ys)
    if abs(y) <= c
       push!(us, x); push!(vs, y)
    else
      length(us) > 0 && plot!(p, us, vs, color=1)
      empty!(us); empty!(vs)
    end
 end
 length(us) > 0 && plot!(p, us, vs, color=1)
 p
end


"""
    plotif(f, g, a, b)

Plot f colored depending on g >= 0 or not.
"""
function plotif(f, g, a, b, args...; colors=(1,2), linewidth=5, legend=false,  kwargs... )


    xs = a:(b-a)/251:b
    zs = f.(xs)
    p = plot(xs, f.(xs), args...; color=colors[1], linewidth=linewidth, legend=legend, kwargs...)

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



# visualize newtons method
function newton_vis(f, x0, a=Inf,b=-Inf; steps=5, kwargs...)
    xs = Float64[x0]
    for i in 1:steps
        push!(xs, xs[end] - f(xs[end]) / f'(xs[end]))
    end

    m,M = extrema(xs)
    m = min(m, a)
    M = max(M, b)

    p = plot(f, m, M; linewidth=3, legend=false, kwargs...)
    plot!(p, zero, m, M)
    for i in 1:steps
        plot!(p, [xs[i],xs[i],xs[i+1]], [0,f(xs[i]), 0])
        scatter!(p, xs[i:i],[0])
    end
    scatter!(p, [xs[steps+1]], [0])
    p
end

##
## --------------------------------------------------
##

"""
riemann: compute Riemann sum approximations to a definite integral. As well, implement trapezoid and Simpson's rule.

Example:
```
f(x) = exp(x^2)
riemann(f, 0, 1, 1000)   # default right-Riemann sums
riemann(f, 0, 1, 1000, method="left")       # left sums
riemann(f, 0, 1, 1000, method="trapezoid")  # use trapezoid rule
riemann(f, 0, 1, 1000, method="simpsons")   # use Simpson's rule
```

"""
function riemann(f::Function, a::Real, b::Real, n::Int; method="right")
  if method == "right"
     meth = (f,l,r) -> f(r) * (r-l)
  elseif method == "left"
     meth= (f,l,r) -> f(l) * (r-l)
  elseif method == "trapezoid"
     meth = (f,l,r) -> (1/2) * (f(l) + f(r)) * (r-l)
  elseif method == "simpsons"
     meth = (f,l,r) -> (1/6) * (f(l) + 4*(f((l+r)/2)) + f(r)) * (r-l)
  end

  xs = a .+ (0:n) * (b-a)/n
  as = [meth(f, l, r) for (l,r) in zip(xs[1:end-1], xs[2:end])]
  sum(as)
end

#######
## simplified multivariable integrals

# limits of integration
endpoints(ys,x) = ((f,x) -> isa(f, Function) ? f(x...) : f).(ys, Ref(x))
# avoid specialization in quadgk
struct FWrapper
    f
end
(F::FWrapper)(x) = F.f(x)

"""
fubini(f, dy, dx)
fubini(f, dz, dy, dx)

Computes numeric integral of `f` over region specified by `dz`, `dy`, `dx`. These are a tuple of values of numbers or univariate functions depending on the value of the term on the right (`dy` can depend on `dx`s value).


*Much* slower than `hcubature` from the `HCubature` package, as it refines flat areas too many times, allocates too much, etc. But does allow a more flexible specification of the region to integrate over, as `hcubature` requires box-like regions.

```
f(x,y,z) = x * y^2 * z^3
fubini(f, (0,1), (0,2), (0,3))  # int_0^3 int_0^2 int_0^1 f(x,y,z) dz dy dx
g(v) = f(v...)
hcubature(g, (0,0,0), (3,2,1))  # same. Not order switched

# triangular like region
fubini(f, (0, y->y), (0, x->x), (0,3))
```
"""
fubini(f, dx)     = quadgk(FWrapper(f), dx...)[1]
fubini(f, ys, xs) = fubini(x -> fubini(y -> f(x,y), endpoints(ys, x)), xs)
fubini(f, zs, ys, xs) = fubini(x ->
    fubini(y ->
        fubini(z -> f(x,y,z),
            endpoints(zs, (x,y))),
        endpoints(ys,x)),
    xs)

##################################################

uvec(x) = x / norm(x)

"""
    `unzip(vs)`
    `unzip(v1, v2, ...)`
    `unzip(r::Function, a, b)`

Take a vector of points described by vectors (as returned by, say
`r(t)=[sin(t),cos(t)], r.([1,2,3])`, and return a tuple of collected x
values, y values, and optionally z values.

If the argument is specified as a comma separated collection of vectors, then these are combined and passed along.

If the argument is a function and two end point, then the function is
evaluated at 100 points between `a` and `b`.

This is useful for plotting when the data is more conveniently
represented in terms of vectors, but the plotting interface requires the x and y values collected.

Examples:
```
using SimplePlots
r(t) = [sin(t), cos(t)]
rp(t) = [cos(t), -sin(t)]
plot(unzip(r, 0, 2pi)...)  # calls plot(xs, ys)

t0, t1 = pi/6, pi/4

p, v = r(t0), rp(t0)
plot!(unzip(p, p+v)...)  # connect p to p+v with line

p, v = r(t1), rp(t1)
quiver!(unzip([p])..., quiver=unzip([v]))
```

Based on `unzip` from the `Plots` package.
"""
unzip(vs) = (A=hcat(vs...); Tuple([A[i,:] for i in eachindex(vs[1])]))
unzip(v,vs...) = unzip([v, vs...])
unzip(r::Function, a, b, n=100) = unzip(r.(range(a, stop=b, length=n)))

xs_ys(vs) = (A=hcat(vs...); Tuple([A[i,:] for i in eachindex(vs[1])]))
xs_ys(v,vs...) = xs_ys([v, vs...])
xs_ys(r::Function, a, b, n=100) = xs_ys(r.(range(a, stop=b, length=n)))


"""
   `arrow!(p, v)`

Add the vector `v` to the plot anchored at `p`.

This would just be a call to `quiver`, but there is no 3-D version of that. As well, the syntax for quiver is a bit awkward for plotting just a single arrow. (Though efficient if plotting many).

```
using SimplePlots
r(t) = [sin(t), cos(t), t]
rp(t) = [cos(t), -sin(t), 1]
plot(unzip(r, 0, 2pi)...)
t0 = 1
arrow!(r(t0), r'(t0))
```
"""
function arrow!(plt, p, v; kwargs...)
#function arrow!(plt::SimplePlots.plot, p, v; kwargs...)
  if length(p) == 2
     quiver!(plt, unzip([p])..., quiver=Tuple(unzip([v])); kwargs...)
  elseif length(p) == 3
    # 3d quiver needs support
    # https://github.com/JuliaPlots/Plots.jl/issues/319#issue-159652535
    # headless arrow instead
    plot!(plt, unzip(p, p+v)...; kwargs...)
	end
end
#arrow!(p,v;kwargs...) = arrow!(Plots.current(), p, v; kwargs...)

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
function vectorfieldplot!(plt, V; xlim=(-5,5), ylim=(-5,5), n=10, kwargs...)

    dx, dy = (xlim[2]-xlim[1])/n, (ylim[2]-ylim[1])/n
    xs, ys = xlim[1]:dx:xlim[2], ylim[1]:dy:ylim[2]

    ps = [[x,y] for x in xs for y in ys]
    vs = V.(ps)
    λ = 0.9 * min(dx, dy) /maximum(norm.(vs))

    quiver!(plt, unzip(ps)..., quiver=unzip(λ * vs))

end
vectorfieldplot!(V; kwargs...) = vectorfieldplot!(SimplePlots.current(), V; kwargs...)
###
include("demos.jl")


end

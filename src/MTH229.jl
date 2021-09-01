"""

`MTH229`: helper functions for using `Julia` with MTH229

This module does two things:

* Install other useful packages with one command (`Roots`,  `ForwardDiff`, `QuadGK`, `SpecialFunctions`, ...)

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

- `fisheye(f)` returns the composition `atan ∘ f ∘ tan`, which can be useful to find zeros of a function over the entire range of real numbers.

- `riemann(f, a, b, n; method="right")` An implementation of Riemann sums. The method can be "right" or "left" for Riemann sums, or "trapezoid" or "simpsons" for related approximations.


This package provides some plotting routines for `Plots`, `SimplePlots`, and `Makie`. For the latter two, some "recipes" for plotting functions and symbolic directions; and for all some convenience methods.

The package `SimplePlots` provides  a quick-to-load-plotting  package  with a  syntax  very similar  to the  more  feature   rich  `Plots` package, proving useful with the Binder service.

"""
module MTH229


# msg = """
# Loading the `MTH229` package.

# * Run the command `?MTH229` for a short description.


# """

# @info msg


if isdefined(Base, :Experimental) && isdefined(Base.Experimental, Symbol("@optlevel"))
    @eval Base.Experimental.@optlevel 1
end

using Reexport
@reexport using Roots
@reexport using SpecialFunctions
@reexport using SymPy
@reexport using QuadGK
@reexport using LinearAlgebra
@reexport using ForwardDiff

using Requires

function __init__()
    @require SimplePlots="307c2aad-90be-4152-b348-f51955fac6ce" include("simpleplots.jl")
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plots.jl")
    @require AbstractPlotting="537997a7-5e4e-5d89-9595-2241ea00577e" include("makie.jl")
end

###
export tangent, secant, D, grad, fisheye
export lim
export bisection, newton
export riemann, fubini
export uvec, xs_ys, unzip, parametric_grid
export e

"e is not Base.MathConstants.ℯ, rather `exp(1)` so that it plays nicely with ForwardDiff"
e = exp(1)


 " f'(x) will find the derivative of `f` using Automatic Differentation from the `ForwardDiff` package "
Base.adjoint(f::Function) = x -> ForwardDiff.derivative(f, float(x))
function D(f, n::Int=1)
    n < 0 && throw(ArgumentError("n must be non-negative"))
    n == 0 && return f
    n == 1 && return x -> ForwardDiff.derivative(f, float(x))
    D(D(f), n-1)
end
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
newton(f, x0; kwargs...) = newton(f, f', x0; kwargs...)

#fzero(f, x0; kwargs...) = Roots.find_zero(f, x0; kwargs...)
#fzero(f, a, b; kwargs...) = Roots.find_zero(f, (a, b); kwargs...)
#fzeros(f, a, b; kwargs...) = Roots.find_zeros(f, a, b; kwargs...)

"""
    fisheye(f)

Transform `f` defined on `(-∞, ∞)` to a new function whose domain is in `(-π/2, π/2)` and range is within `(-π/2, π/2)`. Useful for finding all zeros over the real line. For example

```
f(x) = 1 + 100x^2 - x^3
fzeros(f, -100, 100) # empty just misses the zero found with:
fzeros(fisheye(f), -pi/2, pi/2) .|> tan  # finds 100.19469143521222, not perfect but easy to get
```
"""
fisheye(f) = atan ∘ f ∘ tan




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
unzip(vs) = Tuple([vs[j][i] for j in eachindex(vs)] for i in eachindex(vs[1]))
unzip(v,vs...) = unzip([v, vs...])
unzip(r::Function, a, b, n=100) = unzip(r.(range(a, stop=b, length=n)))
# unzip(a) = map(x -> getfield.(a, x), fieldnames(eltype(a)))

# alternate, should deprecate
xs_ys(vs) = (A=hcat(vs...); Tuple([A[i,:] for i in eachindex(vs[1])]))
xs_ys(v,vs...) = xs_ys([v, vs...])
xs_ys(r::Function, a, b, n=100) = xs_ys(r.(range(a, stop=b, length=n)))

"""
    parametric_grid(us, vs, r)

Create matrices for `xs`, `ys`, `zs` from `r(u,v) = [x(u,v), y(u,v), z(u,v)]`

Used to plot parametrically defined surfaces.
"""
function parametric_grid(us, vs, r)
    unzip(r.(us, vs'))
end


# for plotif. This identifies a vector of colors
function identify_colors(g, xs, colors=(:red, :blue, :black))
    F = (a,b) -> begin
        ga,gb=g(a),g(b)
        ga * gb < 0 && return nothing
        ga >= 0 && return true
        return false
    end
    find_colors(F, xs, colors)
end

# F(a,b) returns true, false, or nothing
function find_colors(F, xs, colors=(:red, :blue, :black))
    n = length(xs)
    cols = repeat([colors[1]], n)
    for i in 1:n-1
        a,b = xs[i], xs[i+1]
        val = F(a,b)
        if val == nothing
            cols[i] = colors[3]
        elseif val
            cols[i] = colors[1]
        else
            cols[i] = colors[2]
        end
    end
    cols[end] = cols[end-1]
    cols
end


end

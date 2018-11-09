VERSION < v"0.7.0-" && __precompile__(false)  # Plots issue?

# If can't load into JuliaBox, use this:
#  include(download("https://raw.githubusercontent.com/mth229/MTH229.jl/master/src/MTH229.jl"))
#
# * Run the command `?visualizations` for a description of some interactive features.
msg = """
Loading the `MTH229` package.

* Run the command `?MTH229` for a short description.


"""

using Compat
Compat.@info msg


"""

`MTH229`: helper functions for using `Julia` with MTH229

This module does two things:

* Install other useful packages with one command (`Plots`, `Roots`,  `ForwardDiff`, `QuadGK`, `SpecialFunctions`, ...)

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

"""
module MTH229
#There are a collection of "demos". Run `?visualizations` for a description.

using Reexport
@reexport using Plots
#@reexport using Roots

@reexport using SpecialFunctions

import ForwardDiff
import QuadGK: quadgk
export quadgk

#using Interact



###
export tangent, secant
export lim,  bisection, riemann
export plotif, trimplot, signchart


" f'(x) will find the derivative of `f` using Automatic Differentation from the `ForwardDiff` package "
Base.adjoint(f::Function) = x -> ForwardDiff.derivative(f, float(x))
D(f, n=1) = n > 1 ? D(D(f), n-1) : x -> ForwardDiff.derivative(f, float(x))
export D

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
            ss[j:(j+1)] = "⋮"
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
            ss[i+1]="a"; ss[j+1]="b"; ss[(i+2):j]="#"
            println(join(ss))
        end

        M = a + (b-a) / 2
    end
    M
end

import Roots
import Roots: newton, find_zero, find_zeros
newton(f, fp, x0; kwargs...) = Roots.find_zero((f,fp), x0, Roots.Newton(); kwargs...)
newton(f, x0; kwargs...) = newton(f, D(f), x0; kwargs...)
fzero(f, x0; kwargs...) = Roots.find_zero(f, x0; kwargs...)
fzero(f, a, b; kwargs...) = Roots.find_zero(f, (a, b); kwargs...)
fzeros(f, a, b; kwargs...) = Roots.find_zeros(f, a, b; kwargs...)

export newton, fzero, fzeros

# some plotting utilities

"""
   trimplot(f, a, b, c=20; kwargs...)

Plot f over [a,b] but break graph if it exceeds c in absolute value.
"""
function trimplot(f, a, b, c=20; kwargs...)
  xs = linspace(a, b, 251) #range(a, stop=b, length=251)
  ys = f.(xs)

  us, vs = Real[], Real[]
  p = plot(us, vs, xlim=(a, b), legend=false, kwargs...)
  for (x,y) in zip(xs, ys)
    if abs(y) <= c
       push!(us, x); push!(vs, y)
    else
      length(us) > 0 && plot!(p, us, vs, color=:blue)
      empty!(us); empty!(vs)
    end
 end
 length(us) > 0 && plot!(p, us, vs, color=:blue)
 p
end


"""
    plotif(f, g, a, b)

Plot f colored depending on g < 0 or not.
"""
function plotif(f, g, a, b; colors=(:blue,:red))

    xs = a:(b-a)/251:b
    zs = f.(xs)
    p = plot(xs, f.(xs), color=colors[1], linewidth=5, legend=false)

    ys = g.(xs)
    ys[ys .< 0] .= NaN

    us,vs = Float64[], Float64[]
    for (i,y) in enumerate(ys)
        if isnan(y)
            if !isempty(vs)
                plot!(us, vs, color=colors[2], linewidth=5)
                empty!(us)
                empty!(vs)
            end
        else
            push!(us, xs[i])
            push!(vs, zs[i])
        end
    end
    if !isempty(vs)
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

  xs = a + (0:n) * (b-a)/n
  as = [meth(f, l, r) for (l,r) in zip(xs[1:end-1], xs[2:end])]
  sum(as)
end


#import SymPy: real_roots
#real_roots(f; kwargs...) = PolynomialZeros.poly_roots(f, Over.R, kwargs...)



###
#include("demos.jl")







end

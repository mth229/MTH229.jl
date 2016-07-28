"""

`MTH229`: helper functions for using `Julia` with MTH229

This module does two things:

* Install other useful packages with one command (`Plots`, `Roots`, and `SymPy`)
* add a number of helper functions.

The helper functions include:

- `secant(f, a, b)`: return a function giving the secant line between ``(a,f(a))`` and ``(b,f(b))``.

- `tangent(f, c)`:  return a function giving the tangent line to ``f(x)`` at the point ``(c,f(c))``.

- `bisection(f, a, b)`: A simple implementation of the bisection method. The interval ``[a,b]`` should be a bracketing interval. For real use, the `fzero(f, a, b)` function, from the `Roots` package, should be used.

- `ctranspose`: This allows the derivative of   a function to be found as with math notation: `f'`. It is an alias to `D(f)` from the Roots package. The notation can be used for higher-order derivatives too: `f''`, `f'''`, ... This uses automatic differentiation from the `ForwardDiff` package.

- `plotif(f, g, a, b)`: Plot the function `f` over the interval `[a,b]` and color differently where ``g(x) > 0`` over ``[a,b]``. By passing in `f` for `g` shows where `f` is positive on `[a,b]`; passing in `f'` shows where `f` is increasing on `[a,b]`; and passing in `f''` shows where `f` is concave up on `[a,b]`.


- `riemann(f, a, b, n; method="right")` An implementation of Riemann sums. The method can be "right" or "left" for Riemann sums, or "trapezoid" or "simpsons" for related approximations.
"""
module MTH229

using Reexport
@reexport using Plots
@reexport using Roots
using ForwardDiff
@reexport using SymPy

Plots.plotly()                          # choose as default

### 
export tangent, secant
export lim,  bisection, riemann
export plotif


" f'(x) will find the derivative of `f` using Automatic Differentation from the `ForwardDiff` package "
Base.ctranspose(f::Function) = x -> ForwardDiff.derivative(f, float(x))

"""
tangent
"""
tangent(f,c) = x -> f(c) + f'(c) * (x-c)

"""
secant
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
	   xs = c + hs 
	 else
	   xs = c - hs
	 end
	 ys = map(f, xs)
	 [xs ys]
end


"""

Simple implementation of the bisection method.

Example
```
bisection(sin, 3, 4)
f(x) = x^5 - x^4 - x^3 - x^2 - x - 1
a = bisection(f, 1, 2)
f(a)
```


"""
function bisection(f::Function, a, b)
    a,b = sort([a,b])
    if f(a) * f(b) > 0
        error("[a,b] is not a bracket. A bracket means f(a) and f(b) have different signs/")
    end

    M = a + (b-a) / 2
    
    while a < M < b

        if f(M) == 0.0
	  break
        end
        ## update step
	if f(a) * f(M) < 0 
	   a, b = a, M
	else
	   a, b = M, b
	end

        M = a + (b-a) / 2
    end
    M
end


"""
plotif(f, g, a, b, args...; kwargs...)

Plot the function `f` over the interval `[a,b]`. Replot the same function in a different color when ``g > 0``.

Examples
```
f(x) = x^4 - x^3 - x^2 - x - 1
plotif(f, f,   -1, 2.1)   # where f is positive
plotif(f, f',  -1, 2.1)   # where f is increasing
plotif(f, f'', -1, 2.1)   # where f is concave up
```
"""
function plotif(f, g, a, b, args...; kwargs...)
    p = plot(f, a, b, args...; kwargs..., linewidth=4, legend=false)
    plot!(p,x -> g(x) > 0.0 ? f(x) : NaN, a, b; linewidth=5)
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
     meth(f,l,r) = f(r) * (r-l)
  elseif method == "left"
     meth(f,l,r) = f(l) * (r-l)
  elseif method == "trapezoid"
     meth(f,l,r) = (1/2) * (f(l) + f(r)) * (r-l)
  elseif method == "simpsons"
     meth(f,l,r) = (1/6) * (f(l) + 4*(f((l+r)/2)) + f(r)) * (r-l)
  end

  xs = a + (0:n) * (b-a)/n
  as = [meth(f, l, r) for (l,r) in zip(xs[1:end-1], xs[2:end])]
  sum(as)
end






###
include("demos.jl")







end

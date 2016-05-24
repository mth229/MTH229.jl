"""

`MTH229`: helper functions for using `Julia` with MTH229

This module does two things:

* Install other useful packages with one command (`Plots`, `Roots`, and `SymPy`)
* add a number of helper functions.

The helper functions include:

- `secant(f, a, b)`: return a function giving returning the secant line between $(a,f(a))$ and $(b,f(b))$.

- `tangent(f, c)`:  return a function giving returning the tangent line to $f(x)$ at the point $(c,f(c))$.

- `bisection(f, a, b)`: A simple implementation of the bisection method. The interval $[a,b]$ should be a bracketing interval. For real use, the `fzero(f, a, b)` function, from the `Roots` package, should be used.

- `ctranspose`: This allows the derivative of   a function to be found as with math notation: `f'`. It is an alias to `D(f)` from the Roots package. The notation can be used for higher-order derivatives too: `f''`, `f'''`, ... This uses automatic differentiation from the `ForwardDiff` package.

- `plotif(f, g, a, b)`: Plot the function `f` over the interval `[a,b]` and color differently where $g(x) > 0$ over $[a,b]$. By passing in `f` for `g` shows where `f` is positive on `[a,b]`; passing in `f'` shows where `f` is increasing on `[a,b]`; and passing in `f''` shows where `f` is concave up on `[a,b]`.

- `newton_vis(f, x0, a=Inf,b=-Inf; steps=5, kwargs...)`. A simple function to visualize some `steps` of newton's method. The values of `a` and `b` are optional, but if not set the x-viewing window will be determined by the points the sequence.

- `riemann(f, a, b, n; method="right")` An implementation of Riemann sums. The method can be "right" or "left" for Riemann sums, or "trapezoid" or "simpsons" for related approximations.
"""
module MTH229

using Plots
using Roots
using SymPy



export tangent, secant
export plotif, newton_viz, riemann

" f'(x) will find the derivative of `f` using Automatic Differentation from the `ForwardDiff` package "
Base.ctranspose(f::Function) = D(f)

"""
tangent
"""
tangent(f,c) = x -> f(c) + f'(c) * (x-c)

"""
secant
"""
secant(f, a, b) = x -> f(a) + (f(b) - f(a)) / b-a * (x - a)

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

Plot the function `f` over the interval `[a,b]`. Replot the same function in a different color when $g > 0$.

Examples
```
f(x) = x^4 - x^3 - x^2 - x - 1
plotif(f, f,   -1, 2.1)   # where f is positive
plotif(f, f',  -1, 2.1)   # where f is increasing
plotif(f, f'', -1, 2.1)   # where f is concave up
```
"""
function plotif(f, g, a, b, args...; kwawrgs...)
  plot([f, x -> g(x) > 0.0 ? f(x) : NaN], a, b, args...; linewidth=5, kwargs...)
end


"""
newton_viz

Simple visualization of Newton's method

```
f(x) = log(x) - 1/10
newton_vis(f, 1)
"""
function newton_vis(f, x0, a=Inf,b=-Inf; steps=5, kwargs...)
    xs = Float64[x0]
    for i in 1:5
        push!(xs, xs[end] - f(xs[end]) / f'(xs[end]))
    end
    
    m,M = extrema(xs)
    m = min(m, a)
    M = max(M, b)
    
    plot(f, m, M; linewidth=3, legend=false, kwargs...)
    plot!(zero)
    for i in 1:4
        plot!([xs[i],xs[i],xs[i+1]], [0,f(xs[i]), 0])
        scatter!(xs[i:i],[0])
    end
    scatter!(xs[5:5], [0])
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














end

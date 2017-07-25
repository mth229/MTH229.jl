using Reactive
using Interact

export trim_viz, function_viz, fzeros_viz,
       bisection_viz, limit_viz, 
       derivative_viz, critical_pts_viz,
       newton_viz, riemann_viz

"""

Some visualizations in the `MTH229` package:

- `trim_viz(f, a, b). Shows a graph of `f` with controls to trim large negative and positive `y` values.

- `functions_viz(f, a, b)`. Shows where a function is positive, increasing, or concave up over `[a,b]`.

- `fzeros_viz(f, a, b)`. Shows zeros of a function over an interval using `fzeros`

- `bisection_viz(f, a , b)`. For a bracketing interval, shows first few bisection steps.

- `limit_viz(f, c, dir="+")`. Shows limiting value as `x` approaches `c`.

- `derivative_viz(f, c)`. Shows the secant line and its approach to the tangent line for different values of "`h`."

- `critical_pts_viz(f, a, b)`. Shows zeros of the derivative found with `fzeros`. (So *may* show some places where derivative is undefined too.)

- `newton_viz(f, x0). A simple function to visualize some `steps` of Newton's method. The values of `a` and `b` are optional, but, if not set, the x-viewing window will be determined by the points of the sequence of approximations

- `riemann_viz(f, a, b)`. Show the Riemann sum partition in action.


"""
visualizations = nothing
export visualizations


"""

`trim_viz(f::Function, a, b)`: plot `f` over `[a,b]` with trimming of large and small values

"""
function trim_viz(f::Function, a=-5, b=5)
    trim(f, lo, hi) = x -> f(x) < lo ? NaN : (f(x) > hi ? NaN : f(x))
    @manipulate for lo in -[2^i for i in 1:10], hi in [2^i for i in 1:10]
        plot(trim(f, lo, hi), a, b)
    end
end
        


"""

Simple illustration of how to find positive, increasing or concavity graphically. Uses `plotif` helper.

```
function_viz(sin, a, b)
```
"""
function function_viz(f::Function, a=-5, b=5)
    @manipulate for u in ["positive", "increasing", "concave up"]
        if u == "positive"
            p = plotif(f, f, a, b, title="Highlighting where f > 0")
        elseif u == "increasing"
            p = plotif(f, D(f), a, b, title="Showng where f' > 0")
        else
            p = plotif(f, D(f,2), a, b, title="Showing where f'' > 0")
        end
        p
    end
end


"""

fzeros_viz: show zeros of a function over `[a,b]`. Uses `fzeros`, so may miss some.

```
fzeros_viz(f, a, b)
```
"""
function fzeros_viz(f::Function, a, b)
    zs = fzeros(f, a, b)
    plot(f, a, b, legend=false)
    scatter!(zs, 0*zs)
end

"""

`bisection_viz(f, a, b)`: for a bracketing interval `[a,b]`, show first few steps of the bisection method

```
bisection_viz(sin, 3, 4)
```
"""
function bisection_viz(f, a, b)
    lo,hi = a, b
    rt = fzero(f, lo, hi)

    @manipulate for n=slider(1:8, value=1, label="n")
        cvged = false
        a,b = lo, hi
        for i in 1:n
            c = (a + b)/2
            for j in 2:i
                if f(c) == 0
                    cvged = true
                    break 
                elseif f(a)*f(c) < 0
                    b = c
                else
                    a = c
                end
            end
        end
        

        
        p = plot(f, lo, hi, linewidth=3, legend=false)
        scatter!(p, [rt], [0], markersize=3)
        if cvged
            scatter!(p, [rt], [0], markersize=6)
        else
            plot!(p, [a, (a + b)/2], [0,0], linewidth=3)
            plot!(p, [(a + b)/2, b], [0,0], linewidth=3)
        end
        p
    end
end


"""
`limit_viz(f, c)`. Show graph of `f` and how `x` approaches `c`

```
f(x) = sin(x) / x
limit_viz(f, 0)
```
"""
function limit_viz(f::Function, c, dir="+")
    L = N(limit(f, c))
    @manipulate for i in slider(0:5, value=0, label="n")
        h = (dir == "+" ? (1/2)^i : -(1/2)^i)
        s = (dir == "+" ? "+" : "-")
        val = f(c + h)
        plot(f, c-1, c + 1, linewidth=3, title="At c $s (1/2)^$i f is $(round(val,3))", legend=false)
        plot!([c-1, c+1], L * [1,1])
        plot!([c-1, c+1], val * [1,1])
        scatter!([c+h], [f(c+h)], markersize=3)
    end
end
    





"""
Simple illustration of how secant line approaches tangent line as `h` goes to `0`. This plots `f` from `[c-1/2, c + 1/2]` and draws a secant line for different values of `h`.

Example:
```
derivative_viz(sin, pi/4)
```
"""
function derivative_viz(f::Function, c, a=c-1/2, b=c+1/2)
    
    @manipulate for j in slider(1:8, value=1, label="h=(1/2)^i")
	h = (1/2)^j
        m = (f(c + h) - f(c)) / h
        err = abs(f'(c) - m)
        e = Int(round(log2(1 / err), 0))

	p = plot(f, a, b, title="Plotting secant line with h=(1/2)^$j, error is about (1/2)^(-$e)", linewidth=3, legend=false)
	m = (f(c + h) - f(c)) / h
	plot!(p, x -> f(c) + m * (x-c), a, b, linewidth=3)
        scatter!(p, [c, c+h], [f(c), f(c+h)], markersize=3)
#	plot!(p, x -> f(c) + f'(c)*(x-c))
        p
    end
    
end


"""

critical_pts_viz: show critical points of a function over `[a,b]`. Uses `fzeros`, so may miss some.

```
critical_pts_viz(f, a, b)
```
"""
function critical_pts_viz(f::Function, a, b)
    zs = fzeros(f', a, b)
    plot(f, a, b, legend=false)
    scatter!(zs, f.(zs))
end


"""
newton_viz

Simple visualization of Newton's method

```
f(x) = log(x) - 1/10
newton_viz(f, 1) # newton_viz(f, 1, steps=20)
"""
function newton_viz(f::Function, x0, a=Inf, b = -Inf; steps=10)
    @manipulate for nsteps = slider(1:steps, value=1, label="no. steps")
        newton_vis(f, x0, a, b; steps=nsteps)
    end
end


"""

Riemann visualization.

Show right-riemann approximation for different values of `n`.

```
riemann_viz(sin, 0, pi)
```

"""
function riemann_viz(f::Function, a, b)
    @manipulate for j in slider(1:8, value=1, label="no. rectangles")
        n = 2^j
        val = riemann(f, a, b, n)
        act_val = quadgk(f, a, b)[1]
        err = abs(val - act_val)
        e = Int(round(log2(1 / abs(val - act_val)), 0))
        p = plot(f, a, b, linewidth=3, legend=false, title="Riemann sum with n=2^$j: $(round(val,3)). Error is about 2^(-$e)")
        delta = (b-a)/n
        for i in 1:n
            plot!(p, a + delta * [i-1, i, i, i-1, i-1], f(a + i * delta) * [0,0,1,1,0], color=:blue)
        end
        p
    end
end

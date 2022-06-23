# A single script to contain most of this package.
# For use within a Pluto notebook
using Plots, SymPy, Roots, ForwardDiff, LinearAlgebra, SpecialFunctions, QuadGK
using Printf, PlutoUI

begin
    Base.adjoint(f::Function) = x -> ForwardDiff.derivative(f, float(x))
    const e = exp(1)
    fisheye(f)=x->atan(f(tan(x)))
    tangent(f,c) = x -> f(c) + f'(c) * (x-c)
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

    "`newton(f,x; kwargs)` run Newton's method. Use `with_terminal` to display trace"
    function newton(f, x; kwargs...)
    	p = Roots.ZeroProblem((f,f'),x)
	i = init(p, Roots.Newton(); kwargs...)
    	println(
            @sprintf(
	    	"%s = % 18.16f,\t %s = % 18.16f",
	    	"x0", float(x), "f(x0)", float(f(x)))
	)
	α = x
	for (i,xᵢ) ∈ enumerate(i)
	    α = xᵢ
    	    println(
                @sprintf(
                    "%s = % 18.16f,\t %s = % 18.16f",
                    "x$i", float(xᵢ), "f(x$i)", float(f(xᵢ)))
            )
	end
	α
    end

    """
   sign_chart(f, a, b; atol=1e-4)

Create a sign chart for `f` over `(a,b)`. Returns a tuple with an identified zero or vertical asymptote and the corresponding sign change. The tolerance is used to disambiguate numerically found values.

# Example

```
julia> sign_chart(x -> x/(x-1)^2, -5, 5)
2-element Vector{NamedTuple{(:∞0, :sign_change), Tuple{Float64, String}}}:
 (∞0 = 0.0, sign_change = "- → +")
 (∞0 = 1.0000000000000002, sign_change = "+ → +")
```

"""
function sign_chart(f, a, b; atol=1e-6)
    pm(x) = x < 0 ? "-" : x > 0 ? "+" : "0"
    summarize(f,cp,d) = (DNE_0_∞=cp, sign_change=pm(f(cp-d)) * " → " * pm(f(cp+d)))

    # check endpoint
    if min(abs(f(a)), abs(f(b))) <= max(max(a,b)*eps(), atol)
        return "Sorry, the endpoints must not be zeros for the function"
    end

    zs = find_zeros(f, a, b)
    pts = vcat(a, zs, b)
    for (u,v) ∈ zip(pts[1:end-1], pts[2:end])
        zs′ = find_zeros(x -> 1/f(x), u, v)
        for z′ ∈ zs′
            flag = false
            for z ∈ zs
                if isapprox(z′, z, atol=atol)
                    flag = true
                    break
                end
            end
            !flag && push!(zs, z′)
        end
    end


    if isempty(zs)
	fc = f(a + (b-a)/2)
	return "No sign change, always " * (fc > 0 ? "positive" : iszero(fc) ? "zero" : "negative")
    end

    sort!(zs)
    m,M = extrema(zs)
    d = min((m-a)/2, (b-M)/2)
    if length(zs) > 1
        d′ = minimum(diff(zs))/2
        d = min(d, d′ )
    end
    summarize.(f, zs, d)
end


    "`riemann(f,a,b,n,[method])` compute Riemann summs"
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

        xs = range(a, b, length=n+1)
        as = [meth(f, l, r) for (l,r) in zip(xs[1:end-1], xs[2:end])]
        sum(as)
    end

    # Visualization
    "`bisection(f,a,b; n=0)` visualize `n` steps of the bisection method"
    function bisection(f, a, b; n=0)
	p=plot(f, a , b, legend=false)
        a′,b′ = a,b
	show_segment=true
	for i in 1:n
	    c = (a + b)/2
	    if f(c) == 0
		scatter!([c], [f(c)])
		show_segment=false
		break
	    end
	    if f(c) * f(a) < 0
		a,b = a,c
	    else
		a,b = c,b
	    end
	end
	if show_segment
            plot!([a′,b′], [0,0], linewidth=1, color=:gray)
	    plot!([a,b], [0,0], linewidth=5)
	end
	p
    end


    "`trimplot(f,a,b, c= 20;...)` trim large values of `f(x)`, then plot"
    function trimplot(f, a, b, c=20; color=:black, legend=false, kwargs...)
        F = (a,b) -> begin
            fa, fb = f(a), f(b)
            M = max(fa, fb)
            m = min(fa, fb)
            m < -c && return false
            M > c && return false
            true
        end
        xs = range(a, b, length=251)
        cols = find_colors(F, xs, (color, :transparent, :red))
        Plots.plot(xs, f.(xs), colors=cols, legend=legend, kwargs...)
    end



    "`newton_vis` visualize Newton's method"
    function newton_vis(f, x0, a=Inf,b=-Inf; steps=5, kwargs...)
        xs = Float64[x0]
        for i in 1:steps
            push!(xs, xs[end] - f(xs[end]) / f'(xs[end]))
        end

        m,M = extrema(xs)
        m = min(m, a)
        M = max(M, b)

        p = Plots.plot(f, m, M; linewidth=3, legend=false, kwargs...)
        Plots.plot!(p, zero)
        for i in 1:steps
            Plots.plot!(p, [xs[i],xs[i],xs[i+1]], [0,f(xs[i]), 0])
            Plots.scatter!(p, xs[i:i],[0])
        end
        Plots.scatter!(p, [xs[steps+1]], [0])
        p
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

    "`plotif(f, g, a, b; kwargs...)` Plot `f`, using red to indicate where `g(x)>0`."
    function plotif(f, g, a, b; kwargs...)
    	xs = range(a, b, length=251)
	cols = identify_colors(g, xs)
    	Plots.plot(xs, f; color=cols, legend=false, kwargs...)
    end

    nothing
end

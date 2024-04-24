"""

`MTH229`: helper functions for using `Julia` with MTH229.

This package reexports functions from the `CalculusWithJulia` package.

Use `mth229()` to download course notebooks.

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
@reexport using CalculusWithJulia
@reexport using QuadGK
@reexport using SimpleExpressions
@reexport using SymPy


## simpleexpressions
import CalculusWithJulia.Roots.CommonSolve: solve
function solve(ex::SimpleExpressions.SymbolicEquation, x₀, args...; kwargs...)
    find_zero(ex, x₀, args...; kwargs...)
end
function solve(ex::SimpleExpressions.SymbolicEquation, I::Interval; kwargs...)
    find_zeros(ex, I; kwargs...)
end
Base.adjoint(f::SimpleExpressions.AbstractSymbolic) = SimpleExpressions.D(f)

###
export fisheye
export bisection, newton
export fubini


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


"""
    fisheye(f)

Transform `f` defined on `(-∞, ∞)` to a new function whose domain is in `(-π/2, π/2)` and range is within `(-π/2, π/2)`. Useful for finding all zeros over the real line. For example

```
f(x) = 1 + 100x^2 - x^3
fzeros(f, -100, 100) # empty just misses the zero found with:
fzeros(fisheye(f), -pi/2, pi/2) .|> tan  # finds 100.19469143521222, not perfect but easy to get
```

By Gunter Fuchs.
"""
fisheye(f) = atan ∘ f ∘ tan


# ---

"""
    fubini(f, [zs], [ys], xs; rtol=missing, kws...)

Integrate `f` of 1, 2, or 3 input variables.

The zs may depend (x,y), the ys may depend on x

## Examples
```
# integrate over the unit square
fubini((x,y) -> sin(x-y), (0,1), (0,1))

# integrate over a triangle
fubini((x,y) -> 1, (0,identity), (0,1 ))

#
f(x,y,z) = x*y^2*z^3
fubini(f, (0,(x,y) ->  x+ y), (0, x -> x), (0,1))
```


!!! Note
    This uses nested calls to `quadgk`. The use of `hcubature` is recommended, typically after a change of variables to make a rectangular domain. The relative tolerance increases at each nested level.
"""
fubini(@nospecialize(f), dx; rtol=missing, kws...) =
    quadgk(f, dx...; rtol=first(skipmissing((rtol, nothing))), kws...)[1]

fubini(@nospecialize(f), ys, xs; rtol=missing, kws...) =
    fubini(x -> fubini(y -> f(x,y), endpoints(ys, x); rtol=rtol), xs;
           rtol = 10*rtol, kws...)

fubini(@nospecialize(f), zs, ys, xs; rtol=missing, kws...) =
    fubini(x ->
           fubini(y ->
                  fubini(z -> f(x,y,z), endpoints(zs, (x,y));
                         rtol=10*10*rtol, kws...),
                  endpoints(ys,x);
                  rtol = 10*rtol),

           xs;
           rtol=rtol)

endpoints(ys,x) = ((f,x) -> isa(f, Function) ? f(x...) : f).(ys, Ref(x))


#=
import ZipFile

"""
    mth229(dirnm=homedir())

Entry point to install projects for MTH229.
"""
function mth229(dirnm=homedir())
    if !isfile(joinpath(dirnm, "01-calculator.ipynb"))
        @warn "installing projects in $dirnm"

        zf = "https://www.github.com/mth229/229-projects/archive/master.zip"
        zarchive = ZipFile.Reader(download(zf))
        !isdir(dirnm) && mkdir(dirnm)
        cd(dirnm)

        @show zarchive.files
        for f in zarchive.files
            nm = basename(f.name)
            occursin("ipynb", nm) || continue
            @info "installing $nm"
            open(nm, "w") do io
                write(io, read(f, String))
            end
        end
    end
end
export mth229
=#

end

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
@reexport using SymPy


###
export bisection, newton

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

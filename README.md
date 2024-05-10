# MTH229

[![Run on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mth229/229-projects/lite?labpath=blank-notebook.ipynb)


Helper files for using `Julia` with MTH229.

This package and the accompanying projects can be accessed online through [binder](https://mybinder.org/v2/gh/mth229/229-projects/master).


Documentation is available at [mth229.github.io](https://mth229.github.io/).

----

To use this package (and a plotting package) with `Julia` issue the command:

```noeval
using MTH229
using Plots
plotly()
```

----

To install `Julia` you can use [Juliaup](https://github.com/JuliaLang/juliaup). For windows, try:

> `winget install julia -s msstore`

For other operating systems this command can be run at the shell:

> `curl -fsSL https://install.julialang.org | sh`


The `MTH229` package can be installed like other `Julia` packages. For example:

```noeval
import Pkg
Pkg.add("MTH229")
```

This package also installs and re-exports several other packages we make use of (`Roots`, `SymPy`, etc/) in  MTH229 at the College of Staten Island.

----

This package does not install a plotting package. The `Plots` package is suggested. Here are some commands to ensure the interactive `plotly` backend is available:

```noeval
Pkg.add("Plots")
Pkg.add("PlotlyBase")
Pkg.add("PlotlyKaleido")
```

In class, this package is used within a notebook environment, started with:

```noeval
using IJulia
notebook()
```

The notebook opens in the home directory. This behavior can be modified by passing a value to the `dir` argument.

If `IJulia` is not installed, it can be done with these commands:

```noeval
Pkg.add("PyCall")
Pkg.add("IJulia")
Pkg.add("Conda")

using Conda
Conda.add("notebook=6.5.6") ## this is a workaround with Windows
Pkg.build("IJulia")
```

## Projects

MTH229 at CSI has several "projects." There are accompanying `ipynb` notebooks to be used from within `IJulia`.

These notebooks can be installed locally by copying and pasting then executing the following commands

```
Pkg.add("ZipFile")
using ZipFile
zf = "https://www.github.com/mth229/229-projects/archive/main.zip"
zarchive = ZipFile.Reader(download(zf))
dirnm = "./229-projects-main"
isdir(dirnm) && error("Directory $dirnm already exists")

mkdir(dirnm)

for f in zarchive.files
    nm = f.name
    occursin("ipynb", nm) || continue
    @show nm
    open(nm, "w") do io
        write(io, read(f, String))
    end
end
```

## Running remotely

Here are three options to run the `MTH229` package on a remote `Julia` installation:

* If you are registered in a class at CSI, we have been providing a service we call `JuliaBox`.

* You can run freely through the resource-limited binder: [![Run on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mth229/229-projects/lite?labpath=blank-notebook.ipynb)

* [Google colab](https://colab.research.google.com/) offers a free service. To use thi package with `colab`, you need to execute a block of commands that downloads `Julia` and installs this package and a plotting package:

```
# Installation cell
%%capture
%%shell
if ! command -v julia 3>&1 > /dev/null
then
    wget -q 'https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.2-linux-x86_64.tar.gz' \
        -O /tmp/julia.tar.gz
    tar -x -f /tmp/julia.tar.gz -C /usr/local --strip-components 1
    rm /tmp/julia.tar.gz
fi
julia -e 'using Pkg; pkg"add IJulia MTH229; precompile;"'
julia -e 'using Pkg; Pkg.add(url="https://github.com/mth229/BinderPlots.jl")' # lighter weight alternative to `Plots+plotly()`
echo 'Now change the runtime type'
```

After this executes (which can take quite some time, as in a few minutes) under the `Runtime` menu select `Change runtime type` and then select `Julia`.

After that, in a cell execute the commands to set up this package and a plotting package:

```
using MTH229
using BinderPlots
```

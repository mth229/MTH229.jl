# MTH229

[![Run on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mth229/229-projects/lite?labpath=blank-notebook.ipynb) Run on binder (with SymPy)

[![Run on Binder w/o SymPy](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mth229/MTH229Lite.jl/main?labpath=blank-notebook.ipynb) Run on binder (without SymPy)




These can be accessed online through [binder](https://mybinder.org/v2/gh/mth229/229-projects/master).


Helper files for using `Julia` with MTH229.

Documentation is available at [mth229.github.io](https://mth229.github.io/).

To use this package (and a plotting package) issue the command:

```noeval
using MTH229
using Plots
```


This package can be installed like other `Julia` packages. For example:

```noeval
import Pkg
Pkg.add("MTH229")
```

This package also installs and re-exports several other packages we make use of (`Roots`, `SymPy`, etc/) in  MTH229 at the College of Staten Island.

This package does not install a plotting package. The `Plots` package is suggested. Here are some commands to ensure the interactive `plotly` backend is available:

```noeval
Pkg.add("Plots")
Pkg.add("PlotlyBase")
Pkg.add("PlotlyKaleido")
```

In class, this package is used within a notebook environment, as installed and opened with:

```noeval
Pkg.add("PyCall")
Pkg.add("IJulia")
Pkg.add("Conda")

using Conda
Conda.add("notebook=6.5.6") ## this is a workaround with Windows
Pkg.build("IJulia")
using IJulia
notebook()
```

## Projects

MTH229 at CSI has several "projects." There are `ipynb` notebooks to be used from within `IJulia`.

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

To run the `MTH229` package remotely, we have three options:

* If registered in a class, we have been provided a service we call `JuliaBox`

* You can run freely through the resource-limited binder: [![Run on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mth229/229-projects/lite?labpath=blank-notebook.ipynb)

* [Google colab](https://colab.research.google.com/) offers a free service. To run this, you need to execute a command that downloads `Julia` and installs this package and a plotting package:

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
julia -e 'using Pkg; Pkg.add(url="https://github.com/mth229/BinderPlots.jl")'
echo 'Now change the runtime type'
```

After this executes (which can take quite some time, as in a few minutes) under the `Runtime` menu select `Change runtime type` and then select `Julia`.

After that, in a cell execute the commands to set up this package and a plotting package:

```
using MTH229
using BinderPlots
```

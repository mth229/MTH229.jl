# MTH229

[![Run on Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mth229/229-projects/lite?labpath=blank-notebook.ipynb)


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

This package does not install a plotting package. The `Plots` package is suggested.

## Projects

MTH229 at CSI has several "projects." There are `ipynb` notebooks to be used from within `IJulia`.

These notebooks can be installed locally by copying and pasting then executing the following commands

```
using ZipFile
zf = "https://www.github.com/mth229/229-projects/archive/main.zip"
zarchive = ZipFile.Reader(download(zf))
dirnm = "./229-projects-master"
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

# MTH229

Helper files for using `Julia` with MTH229.

To use package (and a plotting package) issue the command:

```noeval
using MTH229
using Plots
```


If using `Julia` version `1.7` or higher, or using `Pluto` as an interface, this package will be installed on demand. The installation may take a bit of time to download the necessary files, but this is only done the first time.

For other uses, the package can be installed with commands like:

```noeval
import Pkg
Pkg.add("MTH229")
Pkg.add("Plots")
```

(This command can be issued at a command line *or* just by itself within an IJulia cell.)

In addition to `MTH229` and `Plots`, this should also install and re-export several other packages we make use of (`Roots`, `SymPy`, etc/) in  MTH 229 at the College of Staten Island.


----

To use this package we have to load it into a session with the command:

```
using MTH229
using Plots
```

That also loads a plotting package.


## Projects

MTH229 at CSI has several "projects."

These can be accessed online through [binder](https://mybinder.org/v2/gh/mth229/229-projects/master).

These can be installed locally by copying and pasting then executing the following commands

```
using MTH229.ZipFile
zf = "https://www.github.com/mth229/229-projects/archive/master.zip"
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

There are also `Pluto` notebooks available.

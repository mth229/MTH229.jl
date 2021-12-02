# MTH229

Helper files for using `Julia` with MTH229.

This is installed with commands like:

```noeval
import Pkg
Pkg.add("MTH229")
```

(This command can be issued at a command line *or* just by itself within an IJulia cell.)

This should also install several other packages we make use of (`Roots`, `SymPy`) in  MTH 229 at the College of Staten Island.

A plotting package must be installed, among `Plots`, `SimplePlots`, or `Makie`. For example, this command will install `Plots`:

```noeval
Pkg.add("Plots")
```

----

To use this package we have to load it into a session with the command:

```
using MTH229, Plots
```

That also loads a plotting package.

To find out what is in the package read the help page for the package:

```
?MTH229
```


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

There are also `Pluto` noteboooks available.

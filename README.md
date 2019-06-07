# MTH229

Helper files for using `Julia` with MTH229.

This is installed with

```noeval
] add https://github.com/mth229/MTH229.jl
```

This should also install several other packages we make use of (`Roots`, `Plots`, `SymPy`)


To use this package we have:

```
using MTH229
```

To find out what is in the package read the help page for the package:

```
?MTH229
```

If installation fails, then this command will pull in the necessary functions:

```noeval
include(download("https://raw.githubusercontent.com/mth229/MTH229.jl/master/src/229.jl"))
```

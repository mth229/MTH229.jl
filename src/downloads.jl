using ZipFile


"""
    download_projects()

Download MTH229 projects into directory `./229-projects-master`.
If the directory already exists, this will exit.
"""
function download_projects()
    
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
end

export download_projects

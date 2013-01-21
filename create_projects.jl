## create the projects

using TeachersHelpers, Gadfly, Compose, KnitJ

projects = ["projects.jmd", "calculator.jmd", "functions.jmd", "graphing.jmd", "zeroes.jmd",
            "limits.jmd", "derivatives.jmd", "newton.jmd", "extrema.jmd", "integration.jmd"]
wd = expanduser("~/pmg/JuliaNotes/Tutorial.Rmd")

function run_project(p)
    q = joinpath(wd, p)
    println("---- Working on $q ----")
    knitj_html(q)
    println("---- all done with $q ----")
end
    

function run_projects()
    map(run_project, projects)
end

function scp_gauss()
    html_projects = [replace(x, r"\.jmd$", ".html") for x in projects]
    files = [joinpath(wd, x) for x in html_projects]
    run(`scp $files verzani@wiener:public_html/tmp/julia`)
end

using Pkg

curr_prj = Pkg.API.Context().env.project_file

try
    if startswith(curr_prj, pwd())
        mv(curr_prj, curr_prj * "_orig")
    end

    Pkg.add(PackageSpec(url="https://github.com/JuliaHEP/RadiationDetectorSignals.jl.git", rev="master"))

finally
    if startswith(curr_prj, pwd()) && isfile(curr_prj * "_orig")
        mv(curr_prj * "_orig", curr_prj, force = true)
    end
    Pkg.resolve()
end

# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [fixdoctests]
#
# for local builds.

using Documenter
using RadiationDetectorDSP

makedocs(
    sitename = "RadiationDetectorDSP",
    modules = [RadiationDetectorDSP],
    format = :html,
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    html_prettyurls = !("local" in ARGS),
    html_canonical = "https://JuliaHEP.github.io/RadiationDetectorDSP.jl/stable/",
)

deploydocs(
    repo = "github.com/JuliaHEP/RadiationDetectorDSP.jl.git"
)

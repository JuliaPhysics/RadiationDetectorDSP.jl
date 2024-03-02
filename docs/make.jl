# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using RadiationDetectorDSP

# Doctest setup
DocMeta.setdocmeta!(
    RadiationDetectorDSP,
    :DocTestSetup,
    :(using RadiationDetectorDSP);
    recursive=true,
)

makedocs(
    sitename = "RadiationDetectorDSP",
    modules = [RadiationDetectorDSP],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://JuliaPhysics.github.io/RadiationDetectorDSP.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/JuliaPhysics/RadiationDetectorDSP.jl.git",
    forcepush = true,
    push_preview = true,
)

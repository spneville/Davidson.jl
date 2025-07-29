using Davidson
using Documenter
using Pkg: Pkg

# Where to get files from and where to build them
SRCPATH = joinpath(@__DIR__, "src")
BUILDPATH = joinpath(@__DIR__, "build")
ROOTPATH = joinpath(@__DIR__, "..")
CONTINUOUS_INTEGRATION = get(ENV, "CI", nothing) == "true"

# Setup julia dependencies for docs generation if not yet done
Pkg.activate(@__DIR__)
Pkg.develop(Pkg.PackageSpec(; path = ROOTPATH))
Pkg.instantiate()

# Recursively set the documentation metadata for the Davidson module
DocMeta.setdocmeta!(Davidson, :DocTestSetup,
                    :(using Davidson); recursive=true)

# Generate the docs in the BUILDPATH
makedocs(;
    modules=[Davidson],
    authors="Simon Neville",
    sitename="Davidson.jl",
    format=Documenter.HTML(;
        canonical="https://spneville.github.io/Davidson.jl",
        edit_link="main",
        assets=String[],
        # Use clean URLs, unless built as a "local" build
        prettyurls = CONTINUOUS_INTEGRATION,
    ),
    pages=[
        "Home" => "index.md",
        "solver" => "solver.md",
        "Matrix-Vector Multiplication" => "matvec.md"
    ],
)

# Deploy docs to the gh-pages branch
deploydocs(;
    repo="github.com/spneville/Davidson.jl",
    devbranch="main",
)

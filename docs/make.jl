using Davidson
using Documenter
using Pkg: Pkg

# Where to get files from
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
        "Solvers" => "solver.md",
        "Matrix-vector multiplication function" => "matvec.md",
        "Allowed types" => "allowed_types.md",
        "Creating work arrays" => "work_arrays.md"
    ],
)

# Deploy docs to the gh-pages branch
deploydocs(;
    repo="github.com/spneville/Davidson.jl",
    devbranch="main",
)

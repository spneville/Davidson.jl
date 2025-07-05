using Davidson
using Documenter

DocMeta.setdocmeta!(Davidson, :DocTestSetup, :(using Davidson); recursive=true)

makedocs(;
    modules=[Davidson],
    authors="Simon Neville",
    sitename="Davidson.jl",
    format=Documenter.HTML(;
        canonical="https://spneville.github.io/Davidson.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/spneville/Davidson.jl",
    devbranch="main",
)

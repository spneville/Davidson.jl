using GenDav
using Documenter

DocMeta.setdocmeta!(GenDav, :DocTestSetup, :(using GenDav); recursive=true)

makedocs(;
    modules=[GenDav],
    authors="Simon Neville",
    sitename="GenDav.jl",
    format=Documenter.HTML(;
        canonical="https://spneville.github.io/GenDav.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/spneville/GenDav.jl",
    devbranch="main",
)

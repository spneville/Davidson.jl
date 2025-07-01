module Davidson

using LinearAlgebra
using UnPack

AllowedFloat = Union{Float32, Float64}
AllowedComplex = Union{ComplexF32, ComplexF64}

abstract type Cache end

export solver
export DavidsonCache

include("cache.jl")
include("wrapper.jl")
include("solver.jl")

end

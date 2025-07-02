module Davidson

using LinearAlgebra
using UnPack

AllowedFloat = Union{Float32, Float64}
AllowedComplex = Union{ComplexF32, ComplexF64}
AllowedTypes = Union{AllowedFloat, AllowedComplex}

abstract type Cache end

export solver
export DavidsonCache
export worksize

include("cache.jl")
include("wrapper.jl")
include("solver.jl")

end

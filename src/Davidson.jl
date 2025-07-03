module Davidson

using LinearAlgebra
using UnPack

AllowedFloat = Union{Float32, Float64}
AllowedComplex = Union{ComplexF32, ComplexF64}
Allowed64 = Union{Float64, ComplexF64}
AllowedTypes = Union{AllowedFloat, AllowedComplex}


abstract type Cache end

export solver, solver!
export DavidsonCache
export workarrays

include("cache.jl")
include("wrapper.jl")
include("solver.jl")

end

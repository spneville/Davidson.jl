module Davidson

using LinearAlgebra
using UnPack
using Printf
using InteractiveUtils

AllowedFloat = Union{Float32, Float64}
AllowedComplex = Union{ComplexF32, ComplexF64}
Allowed64 = Union{Float64, ComplexF64}

"""
    AllowedTypes = Union{Float32, Float64, ComplexF32, ComplexF64}

The allowed types `T` of the matrix whose eigenpairs are to be computed.
"""
AllowedTypes = Union{AllowedFloat, AllowedComplex}

abstract type Cache end

export solver, solver!
export DavidsonCache
export workarrays

export AllowedTypes

include("cache.jl")
include("wrapper.jl")
include("solver.jl")

end

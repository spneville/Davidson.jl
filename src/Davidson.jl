module Davidson

using LinearAlgebra
using UnPack
using Printf
using InteractiveUtils

"""
    AllowedFloat = Union{Float32, Float64}

The allowed types `R` of the eigenvalues and on-diagonal matrix elements.
"""
AllowedFloat = Union{Float32, Float64}

AllowedComplex = Union{ComplexF32, ComplexF64}

"""
    AllowedTypes = Union{Float32, Float64, ComplexF32, ComplexF64}

The allowed types `T` of the matrix whose eigenpairs are to be computed.
"""
AllowedTypes = Union{AllowedFloat, AllowedComplex}

Allowed64 = Union{Float64, ComplexF64}

abstract type Cache end

export AllowedTypes
export AllowedFloat

export solver, solver!
export DavidsonCache
export workarrays

include("cache.jl")
include("wrapper.jl")
include("solver.jl")

end

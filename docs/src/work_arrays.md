# [Creating work arrays](@id WorkArrays)

## Required work arrays

The in-place `solver!` function can be made _almost_ allocation free
via the passing of optional pre-allocated work arrays. This can be
advantageous if a large number of different diagonalisations are to be
performed and the maximum dimensions across all calculations are known
ahead of time. The work arrays have a _minimum_ dimension, but not a
maximum one. Thus, they can be allocated once, using the maximum dimensions
across all calculations to be performed, and then used in every calculation.

Two work vectors have to be supplied:

* `Twork::Vector{T}`, where [`T<:AllowedTypes`](@ref AllowedTypes)
* `Rwork::Vector{R}`, where [`R<AllowedFloat`](@ref AllowedFloat)

where the types `T` and `R` are compatible (see [`Allowed types`](@ref AllowedTypes))

## Generating the work arrays

The required dimensions of the `Twork` and `Rwork` vectors are dependent on
the size of the matrix being diagonalised, the block size, and the maximum
subspace dimension. There exists a function, `workarrays`, that will return
the correctly dimensioned arrays given these values.

```@docs
workarrays
```

# [Matrix-Vector Multiplication](@id matvec)
This section details the required form of the in-place matrix-vector
multiplication function `f` that the `solver` and `solver!` functions
take as an argument.

## Arguments of the function `f`
* `v::AbstractMatrix{T}`: Vectors to be multiplied by the matrix `A`, where [`T<:AllowedTypes`](@ref AllowedTypes)
* `Av::AbstractMatrix{T}`: Matrix-vector products, where [`T<:AllowedTypes`](@ref AllowedTypes)

Here, the columns of the matrices `v` and `Av` contain the vectors in question.

The function must take the following form, where [`T<:AllowedTypes`](@ref AllowedTypes):

```julia
function f!(v::AbstractMatrix{T}, Av::AbstractMatrix{T})

# The contents of the function will be application-dependent,
# but the matrix A will be applied to the vectors contained in
# v to yield the vectors in Av

end
```
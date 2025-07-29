# [solver](@id solver)
This section details the use of the `solver` function, which computes
a given number of the lowest-lying eigenpairs of a given matrix `A` in a
'matrix-free' manner. That is, only a function returning matrix-vector
products is required as input, as opposed to the matrix itself.

```@docs
solver
```
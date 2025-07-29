# [Davidson.jl](@id man-davidson)
A Julia package implementing the Davidson algorithm for the
calculation of the lowest-lying eigenpairs of large, sparse,
diagonally dominant Hermitian matrices.

The implementation is 'matrix-free'. That is, for a given matrix `A`,
the matrix itself does not have to be provided, only a function `f`
that returns the action of `A` on a set of input vectors.

Multiple roots are computed simultaneously using a block variant of the
algorithm.

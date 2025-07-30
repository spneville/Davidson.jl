# Davidson

A Julia package implementing the Davidson algorithm for the
calculation of the lowest-lying eigenpairs of large, sparse,
diagonally dominant Hermitian matrices.

The implementation is 'matrix-free'. That is, for a given matrix `A`,
the matrix itself does not have to be provided, only a function `f`
that returns the action of `A` on a set of input vectors.

Multiple roots are computed simultaneously using a block variant of the
algorithm.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://spneville.github.io/Davidson.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://spneville.github.io/Davidson.jl/dev/)
[![Build Status](https://github.com/spneville/Davidson.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/spneville/Davidson.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/spneville/Davidson.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/spneville/Davidson.jl)

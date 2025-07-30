# Davidson

A Julia package implementing the Davidson algorithm for the
calculation of the lowest-lying eigenpairs of large, sparse,
diagonally dominant Hermitian matrices.

The implementation is 'matrix-free'. That is, for a given matrix `A`,
the matrix itself does not have to be provided, only a function `f`
that returns the action of `A` on a set of input vectors.

Multiple roots are computed simultaneously using a block variant of the
algorithm.

| **Documentation** | **Build Status** | **License** |
|:-----------------:|:----------------:|:---------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![CI][github-img]][github-url] [![][codecov-img]][codecov-url] | [![license][license-img]][license-url] |

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://spneville.github.io/Davidson.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://spneville.github.io/Davidson.jl/stable/

[github-img]: https://github.com/spneville/Davidson.jl/actions/workflows/CI.yml/badge.svg?branch=main
[github-url]: https://github.com/spneville/Davidson.jl/actions/workflows/CI.yml?query=branch%3Amain

[codecov-img]: https://codecov.io/gh/spneville/Davidson.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/spneville/Davidson.jl

[license-img]: http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat
[license-url]: LICENSE.md

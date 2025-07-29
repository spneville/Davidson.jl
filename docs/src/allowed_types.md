# Allowed types

## [Allowed matrix types](@id AllowedTypes)

```@docs
AllowedTypes
```

## [Allowed real types](@id AllowedFloat)

```@docs
AllowedFloat
```

## Compatible types
The arrays of type `R<:AllowedFloat` supplied to the `solver` and `solver!` functions
must be compatible with the arrays of type `R<:AllowedTypes`. For example, if `T==ComplexF64`, then
it is required that `R==Float64`. The following table lists the compatible type pairs:

| `T`           | `R`         |
| ------------  | ----------- |
| `Float32`     | `Float32`   |
| `Float64`     | `Float64`   |
| `ComplexF32`  | `Float32`   |
| `ComplexF64`  | `Float64`   |

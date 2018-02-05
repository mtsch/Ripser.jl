# Ripser.jl

[![Build Status](https://travis-ci.org/mtsch/Ripser.jl.svg?branch=master)](https://travis-ci.org/mtsch/Ripser.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/le4fbrk5hsgnf3ax?svg=true)](https://ci.appveyor.com/project/mtsch/ripser-jl)
[![Coverage Status](https://coveralls.io/repos/github/mtsch/Ripser.jl/badge.svg?branch=master)](https://coveralls.io/github/mtsch/Ripser.jl?branch=master)

Simple wrapper to [Ripser](https://github.com/Ripser/ripser) in Julia.

## Installation

This package requires
[PersistenceBarcodes.jl](https://github.com/mtsch/PersistenceBarcodes.jl).

To install, run:

```
Pkg.clone("https://github.com/mtsch/PersistenceBarcodes.jl")
Pkg.clone("https://github.com/mtsch/Ripser.jl")
Pkg.build("Ripser")
```

## Usage

```
ripser(mat; dim_max = 1, thresh = Inf, modulus = 2)
```

Run Ripser on a distance matrix, returning a `PersistenceDiagram`. Arguments:

* `mat`: the distance matrix. Only the lower triangle of the matrix is used,
  unless the type of matrix is `UpperTriangular`.
* `dim_max`: compute persistent homology up to this dimension.
* `threshold`: compute Rips complexes up to this diameter.
* `modulus`: compute homology with coefficients in the primer field Z/*p*Z,
  where *p* is the value given.

```
read_lowertridist(filename)
```

Read lower diagonal matrix in comma-separated format. See
[`examples`](examples) for example files.

## Notes

The current implementation relies on capturing Ripser's `STDOUT` and parsing it.

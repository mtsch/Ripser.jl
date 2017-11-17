# Ripser.jl

[![Build Status](https://travis-ci.org/mtsch/Ripser.jl.svg?branch=master)](https://travis-ci.org/mtsch/Ripser.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/le4fbrk5hsgnf3ax?svg=true)](https://ci.appveyor.com/project/mtsch/ripser-jl)
[![Coverage Status](https://coveralls.io/repos/github/mtsch/Ripser.jl/badge.svg?branch=master)](https://coveralls.io/github/mtsch/Ripser.jl?branch=master)

Simple wrapper to [Ripser](https://github.com/Ripser/ripser) in Julia.

## Usage

```
ripser(mat; dim_max = 1, thresh = Inf, modulus = 2)
```

Run Ripser on a distance matrix, returning a `PersistenceDiagram`. Arguments:

* `mat`: the distance matrix. Only the lower triangle of the matrix is used,
  unless the type of matrix is `UpperTriangular`.
* `dim_max`: compute persistent homology up to this dimension.
* `threshold`: compute Rips complexes up to this diameter.
* `modulus`: compute homology with coefficients in the primer field Z/_p_Z,
  where _p_ is the value given.

```
PersistenceDiagram
```

The type returned by `ripser`. Contains a vector of vectors of persistence pairs
(`Tuple{Float64, Float64}`). Operations:

* `pd[i]`: indexing (like an array). Note: the indexing is 1-based.
* `dim(pd)`: get the maximum dimension of `pd`.
* `print(pd)`: show all persistence pairs.
* `plot(pd)`: plot a persistence diagram using
  [Plots.jl](https://github.com/JuliaPlots/Plots.jl).

```
read_lowertridist(filename)
```

Read lower diagonal matrix in comma-separated format. See
[`examples`](examples) for example files.

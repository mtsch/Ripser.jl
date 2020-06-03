# Ripser.jl

[![Build Status](https://travis-ci.org/mtsch/Ripser.jl.svg?branch=master)](https://travis-ci.org/mtsch/Ripser.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/le4fbrk5hsgnf3ax?svg=true)](https://ci.appveyor.com/project/mtsch/ripser-jl)
[![codecov](https://codecov.io/gh/mtsch/Ripser.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mtsch/Ripser.jl)
<!---
[![Coverage -->
<!--Status](https://coveralls.io/repos/github/mtsch/Ripser.jl/badge.svg?branch=master)](https://coveralls.io/github/mtsch/Ripser.jl?branch=master)
-->

Simple wrapper to [Ripser](https://github.com/Ripser/ripser) in Julia with almost no
features and minimal dependencies.

Ripser.jl _should_ work, but is not maintained. If you are interested in this project, it
might be worth taking a look at [Ripserer.jl](https://github.com/mtsch/Ripserer.jl).

## Installation

This package is unregistered. To install, run it with:

```
(v1.0) pkg> add https://github.com/mtsch/Ripser.jl#master
(v1.0) pkg> test Ripser
```

## Usage

```
ripser(dists; modulus = 2, dim_max = 1, threshold = Inf, cocycles = false)
```

Run Ripser on a dists, and a square (sparse) matrix of `T<:AbstractFloat`.

* `dists`: the distance matrix. Matrix can be sparse.
* `dim_max`: compute persistent homology up to this dimension.
* `threshold`: compute Rips complexes up to this diameter.
* `modulus`: compute homology with coefficients in the prime field Z/*p*Z,
  where *p* is the value given.
* `cocycles`: return representative cocycles.

Returns a `Vector{Vector{Tuple{T, T}}}` and if `cocycles` is set to `true`, a
`Vector{Vector{Tuple{Vector{Int64},Int64}}}`.

The first part of the result contains a persistence diagram per each dimension (and is of
length `dim_max + 1`). Each persistence diagram is a vector of tuples, where a tuple is the
persistence pair `(birth, death)`. A single cocycle is represented as a vector of `(simplex,
value)` pairs. They are only included if the `value` is nonzero. Each persistence pair gets
a corresponding cocycle.

Note that representative cocycles for the zeroth dimension are not computed.

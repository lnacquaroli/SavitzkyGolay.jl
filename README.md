# SavitzkyGolay.jl
Implementation of the 1D Savitzky-Golay filter in [JuliaLang](https://julialang.org/).

Simple-plain implementation of the Savitzky-Golay filter in Julia.

## Installation

This package is not yet registered. It can be installed in Julia with the following ([see further](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html#Adding-unregistered-packages-1)):
```julia
julia> ]
pkg> add https://github.com/lnacquaroli/SavitzkyGolay.jl
```

## Usage

After installation, we load the package:
```julia
using SavitzkyGolay
using Plots # This is for visualization purposes, not required in the SG package itself
```

Suppose we have a signal with noise that want to smooth out. The function for this is `savitzky_golay`, that accepts the following arguments:

```julia
savitzky_golay(y::Vector, window_size::Int, order::Int; deriv::Int=0, rate::Real=1.0)    
```

where `y` is the data vector with noise, and `window_size`, `order`, `deriv` and `rate` are options for the Savitzky-Golay filter. Notice that `y`, `window_size` and `order` are mandatory arguments to input, while `deriv=0` and `rate=1.0` by default.

## Examples

```julia
t = LinRange(-4, 4, 500)
y = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
sg = savitzky_golay(y, 21, 4)
plot(t, [y sg.y], label=["Original signal" "Filtered signal"])
```
![Example 2: SG2](https://github.com/lnacquaroli/SavitzkyGolay.jl/blob/main/examples/Figure_2.png "Example 2: SG2")

Another simpler example:
```julia
t = 0:20
y = collect(0:20)
sg = savitzky_golay(y, 11, 2)
plot(t, [y sg.y], label=["Original signal" "Filtered signal"])
```
![Example 1: SG1](https://github.com/lnacquaroli/SavitzkyGolay.jl/blob/main/examples/Figure_1.png "Example 1: SG1")

The solution `sg` is a `SGolayResults` type that contains four fields: 
- `y` with the filtered signal,
- `params` type `SGolay` with the initial parameters
- `coeff` with the computed coefficients
- `Vdm` with the Vandermonde matrix

## Constructor

There is an option to call the constructor `SGolay` to build the filter and then use it in different places. To call the constructor you need to specify at least two parameters, the full window size and the polynomial order. The constructor accepts the following arguments:

```julia
SGolay(window_size, polynomial_order, derivative, rate)
```

For instance:
```julia
sgfilter = SGolay(11, 2)

sgfilter = SGolay(11, 2, 1)

sgfilter = SGolay(11, 2, 1, 0.1)
```

By default, if not specified, `deriv=0` and `rate=1.0`.

The same examples above with constructors are as follow:

```julia
t = 0:20
y = collect(0:20)
sgfilter1 = SGolay(11, 2)
y1 = sgfilter1(y)
```

```julia
t = LinRange(-4, 4, 500)
y = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
sgfilter2 = SGolay(21, 4)
y2 = sgfilter2(y)
```

The solutions `y1` and `y2` are the same type as the `SGolayResults`.

## References

- [SciPyCookbook](https://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html)
- [StackOverflow - 48421852](https://stackoverflow.com/a/48421852)
- [Gist - jiahao](https://gist.github.com/jiahao/b8b5ac328c18b7ae8a17)

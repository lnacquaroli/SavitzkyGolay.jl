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

Suppose we have a signal with noise that want to smooth out.

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

The solution `sg` is a NamedTuple that contains three fields: 
- `y` with the filtered signal,
- `params` with the initial parameters
- `coeff` with the computed coefficients


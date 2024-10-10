# SavitzkyGolay.jl
Implementation of the 1D Savitzky-Golay filter in [JuliaLang](https://julialang.org/).

Simple-plain implementation of the Savitzky-Golay filter in Julia.

## Installation

This package is registered and can be installed in Julia with the following:
```julia
julia> ]
pkg> add SavitzkyGolay
```

## Usage

After installation, we load the package:
```julia
using SavitzkyGolay
using Plots # This is for visualization purposes, not required in the SG package itself
```

Suppose we have a signal with noise that want to smooth out. The function for this is `savitzky_golay`,
which accepts the following arguments:

```julia
sg = savitzky_golay(y::AbstractVector, window_size::Int, order::Int; deriv::Int=0, rate::Real=1.0, haswts=false)    

sg = savitzky_golay(y::AbstractVector, wts::AbstractVector, window_size::Int, order::Int; deriv::Int=0, rate::Real=1.0, haswts=true)    
```

- `y`: The data vector with noise to be filtered.
- `wts` a non-negative weights vector of length `window_size` (optional)
- `window_size`: The length of the filter window (i.e., the number of coefficients). Must be an odd number.
- `order`: The order of the polynomial used to fit the samples. Must be less than `window_size`.
- `deriv`: The order of the derivative to compute. This must be a non-negative integer. The default is 0, which means to filter the data without differentiating. If `deriv > 0` it may need scaling which can be achieved using the `rate` optional argument. (optional) 
- `rate`: Scaling real number when using the derivative. (optional)
- `haswts` (Bool whether a weight vector is to be used, defaults to false if no `wts` argument given)

The solution `sg` is a `SGolayResults` type that contains four fields: 

- `y` with the filtered signal,
- `params` type `SGolay` with the initial parameters
- `coeff` with the computed coefficients
- `Vdm` with the Vandermonde matrix
- `haswts` (Bool whether a weight vector is to be used)

## Examples

```julia
t = LinRange(-4, 4, 500)
y2 = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
y2_sg = savitzky_golay(y2, 21, 4)
plot(t, [y2 y2_sg.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)
```
![Example 1: SG1](https://github.com/lnacquaroli/SavitzkyGolay.jl/blob/main/examples/output-sg-exp.png "Example 1: SG1")

Another simpler example:
```julia
t = 0:20
y1 = collect(0:20)
y1_sg = savitzky_golay(y1, 11, 2)
plot(t, [y1 y1_sg.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)
```
![Example 2: SG2](https://github.com/lnacquaroli/SavitzkyGolay.jl/blob/main/examples/output-sg-line.png "Example 2: SG2")

Example with derivatives:

```julia
x = LinRange(-5, 15, 200)
data = 0.15*x.^3 - 2*x.^2 + x  .+ randn(length(x))
data_derivative = 0.45*x.^2 - 4*x .+ 1
sg = savitzky_golay(data, 21, 3, deriv=1)
sg_rate = savitzky_golay(data, 21, 3, deriv=1, rate=200/(15-(-5)))
plot(x, [data data_derivative sg.y sg_rate.y ], label=["Data" "Exact Derivative" "SG" "SG with rate"])
```

![Example 3: SG3](https://github.com/lnacquaroli/SavitzkyGolay.jl/blob/main/examples/output-sg-deriv.png "Example 3: SG3 with derivative")

This is filtering with a constant weights vector which is the same as the un-weighted Savitzky-Golay
filtering above in example 2:

```julia
y1 = collect(0:20)
wts_11 = ones(11)
y1_sg_w1 = savitzky_golay(y1,wts_11,11,2)
plot(y1, [y1 y1_sg_w1.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)
```
This demonstrates filtering with a triangle weights vector and the figure shows the difference between the
un-weighted SG and the weighted SG:

```julia
t = LinRange(-4, 4, 500)
y2 = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
y2_sg = savitzky_golay(y2, 21, 4)
tri_21 = Float64.(vcat( 1:11, 10:-1:1 ))
y2_sg_tri = savitzky_golay(y2, tri_21, 21, 4)
plot(t, [y2 y2_sg.y y2_sg_tri.y], label=["Original signal" "Filtered (no weights)" "Filtered (triangle weights)"],
     lc=[RGBA(0.3,0.5,0.7,0.3) 2 3], ylabel="", xlabel="t", legend=:topleft)
```

![Example 4: SG4](https://github.com/devel-chm/SavitzkyGolay.jl/blob/main/examples/output-sg-tri-wts.png "Example 4: SG4 with derivative")

## Constructor

There is an option to call the constructor `SGolay` to build the filter and then use it in different places. To call the constructor you need to specify at least two parameters, the full window size, and the polynomial order. The constructor accepts the following arguments:

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

The same examples above with constructors are as follows:

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

## With window weights


## References

- [SciPyCookbook](https://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html), https://stackoverflow.com/a/49330505
- [StackOverflow - 48421852](https://stackoverflow.com/a/48421852)
- [Gist - jiahao](https://gist.github.com/jiahao/b8b5ac328c18b7ae8a17)

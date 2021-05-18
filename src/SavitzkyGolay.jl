module SavitzkyGolay

#
# References:
# https://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html
# https://stackoverflow.com/a/48421852
#

using LinearAlgebra

export savitzky_golay

function savitzky_golay(y::AbstractVector, window_size::T, order::T; deriv::T=0, rate::T=1) where T <: Real

    p = _check_inputs_sg(y, window_size, order, deriv, rate)

    # Determine order range and half window size
    order_range = 0 : p.order
    hw = Int((p.w - 1) / 2)

    # Compute coefficients
    c = [k^i for k in -hw:hw, i in order_range]
    m = pinv(c)[p.deriv+1,:] * (p.rate)^p.deriv * factorial(p.deriv)

    # Pad the signal at the extremes with values taken from the signal itself
    initvals = p.y[1] .- abs.(reverse(p.y[2:hw+1]) .- p.y[1] )
    endvals = p.y[end] .+ abs.(reverse(p.y[end-hw:end-1] .- p.y[end]) )
    y_ = vcat(initvals, p.y, endvals)

    return (y=_convolve_1d(y_, m), params=p, coeff=m)
end

function _check_inputs_sg(y, w, o, d, r)
    isodd(w) || throw(ArgumentError("w must be an even number."))
    w ≥ 1 || throw(ArgumentError("w must greater than or equal to 1."))
    w ≥ o + 2 || throw(ArgumentError("w too small for the polynomial order chosen (w ≥ order + 2)."))
    length(y) > 1 || throw(ArgumentError("vector x must have more than one element."))
    return (y=Float64.(y), w=Int64(w), order=Int64(o), deriv=Int64(d), rate=Float64(r))
end

function _convolve_1d(u::Vector, v::Vector)
    m = length(u)
    n = length(v)
    m > n || throw(ArgumentError("length of signal u must be greater than length of kernel v."))
    w = zeros(m + n - 1)
    @inbounds for j in 1:m, k in 1:n
        w[j+k-1] += u[j]*v[k]
    end
    return w[n:end-n+1]
end

end  # module SavitzkyGolay

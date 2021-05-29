module SavitzkyGolay

export savitzky_golay, SGolay, SGFilter

using LinearAlgebra

Base.@kwdef struct SGolay
    w::Int64          # Window size
    order::Int64      # Polynomial order
    deriv::Int64 = 0  # Derivative order
    rate::Real = 1.0  # Rate
end

struct SGFilter{T <: Float64}
    y::Array{T}
    params::NamedTuple
    coeff::Array{T}
    Vdm::Matrix{Float64}
end

function (p::SGolay)(y::AbstractVector)
    return savitzky_golay(y, p.w, p.order; deriv=p.deriv, rate=p.rate)
end

function savitzky_golay(
    y::AbstractVector, window_size::T0, order::T0;
    deriv::T0=0, rate::T1=1.0,
    ) where {T0 <: Int64, T1 <: Real}

    p = _check_input_sg(y, window_size, order, deriv, rate)

    # Determine order range and half window size
    order_range = 0 : p.order
    hw = Int((p.w - 1) / 2)

    # Build Vandermonde matrix
    V = _vandermonde(hw, order_range)

    # Compute coefficients
    c = _coefficients(V, order_range, p)

    # Pad the signal at the extremes with values taken from the signal itself
    y_ = _padding_signal(p, hw)

    # Convolve signal and kernel
    y_conv = _convolve_1d(y_, c)

    # return (y=y_conv, params=p, coeff=c, Vdm=V)
    return SGFilter(y_conv, p, c, V)
end

function _check_input_sg(y::Vector, w, o, d, r)
    isodd(w) || throw(ArgumentError("w must be an even number."))
    w ≥ 1 || throw(ArgumentError("w must greater than or equal to 1."))
    w ≥ o + 2 || throw(ArgumentError("w too small for the polynomial order chosen (w ≥ order + 2)."))
    length(y) > 1 || throw(ArgumentError("vector x must have more than one element."))
    return (y=Float64.(y), w=w, order=o, deriv=d, rate=Float64(r))
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

function _vandermonde(hw, order_range)
    V = zeros(2*hw + 1, length(order_range))
    @inbounds for i in -hw:hw, j in order_range
        V[i+hw+1, j+1] = i^j
    end
    return V
end

function _coefficients(V, order_range, p)
    Vqr = qr(V')
    c = Vqr.R \ (Vqr.Q' * _onehot(p.deriv + 1, length(order_range)))
    c .*= (p.rate)^p.deriv * factorial(p.deriv)
    return c
end

function _onehot(i, m)
    m > i || throw(ArgumentError("length of vector must be greater than the position"))
    oh = zeros(m)
    oh[i] = 1.0
    return oh
end

function _padding_signal(p, hw)
    initvals = p.y[1] .- abs.(reverse(p.y[2:hw+1]) .- p.y[1])
    endvals = p.y[end] .+ abs.(reverse(p.y[end-hw:end-1] .- p.y[end]))
    return vcat(initvals, p.y, endvals)
end

end  # module SavitzkyGolay

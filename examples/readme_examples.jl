
using SavitzkyGolay
using Plots

# Plain call
t = LinRange(-4, 4, 500)
y2 = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
y2_sg = savitzky_golay(y2, 21, 4)
plot(t, [y2 y2_sg.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

t = 0:20
y1 = collect(0:20)
y1_sg = savitzky_golay(y1, 11, 2)
plot(t, [y1 y1_sg.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

# With derivative
x = LinRange(-5, 15, 200)
data = 0.15*x.^3 - 2*x.^2 + x  .+ randn(length(x))
data_derivative = 0.45*x.^2 - 4*x .+ 1
sg = savitzky_golay(data, 21, 3, deriv=1)
sg_rate = savitzky_golay(data, 21, 3, deriv=1, rate=200/(15-(-5)))
plot(x, [data data_derivative sg.y sg_rate.y ], label=["Data" "Exact Derivative" "SG" "SG with rate"])

# with uniform window weights (same is un-weighted Savitzky-Golay filtering)
y1 = collect(0:20)
wts_11 = ones(11)
y1_sg_w1 = savitzky_golay(y1,wts_11,11,2)
plot(y1, [y1 y1_sg_w1.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

# with triangle window weights (a good default weighting---near optimum)
t = LinRange(-4, 4, 500)
y2 = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
y2_sg = savitzky_golay(y2, 21, 4)
tri_21 = Float64.(vcat( 1:11, 10:-1:1 ))
y2_sg_tri = savitzky_golay(y2, tri_21, 21, 4)
plot(t, [y2 y2_sg.y y2_sg_tri.y], label=["Original signal" "Filtered (no weights)" "Filtered (triangle weights)"],
     lc=[RGBA(0.3,0.5,0.7,0.3) 2 3], ylabel="", xlabel="t", legend=:topleft)

# With constructor
t = 0:20
y1 = collect(0:20)
f1 = SGolay(11, 2)
y1_sg = f1(y1)
plot(t, [y1 y1_sg.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

t = LinRange(-4, 4, 500)
y2 = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
f2 = SGolay(21, 4)
y2_sg = f2(y2)
plot(t, [y2 y2_sg.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

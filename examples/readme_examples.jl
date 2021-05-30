
using SavitzkyGolay
using Plots

# Plain call
t = 0:20
y = collect(0:20)
y_ = savitzky_golay(y, 11, 2)
plot(t, [y y_.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

t = LinRange(-4, 4, 500)
y = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
y2_ = savitzky_golay(y, 21, 4)
plot(t, [y y2_.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)


# With constructor
t = 0:20
y = collect(0:20)
f1 = SGolay(11, 2)
y1 = f1(y)
plot(t, [y y1.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

t = LinRange(-4, 4, 500)
y = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
f2 = SGolay(21, 4)
y2 = f2(y)
plot(t, [y y2.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

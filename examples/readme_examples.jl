using SavitzkyGolay
using Plots

t = 0:20
y = collect(0:20)
sg = savitzky_golay(y, 11, 2)
plot(t, [y y_.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

t = LinRange(-4, 4, 500)
y = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
sg = savitzky_golay(y, 21, 4)
plot(t, [y sg.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

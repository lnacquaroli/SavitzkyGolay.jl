
cd(raw"C:\Users\flvisa\Desktop\lna\numerical_analysis\github\savitzky_golay")

include("SavitzkyGolay.jl")
import .SavitzkyGolay as sgolay
import Plots as plt
# plt.pyplot()

t = 0:20
y = collect(0:20)
y_ = sgolay.savitzky_golay(y, 11, 2)
plt.plot(t, [y y_.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

t = LinRange(-4, 4, 500)
y = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
y_ = sgolay.savitzky_golay(y, 21, 4)
plt.plot(t, [y y_.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

# With constructor
t = 0:20
y = collect(0:20)
f1 = sgolay.SGolay(w=11, order=2)
y1 = f1(y)
plt.plot(t, [y y1.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)

t = LinRange(-4, 4, 500)
y = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
f2 = sgolay.SGolay(w=21, order=4)
y2 = f2(y)
plt.plot(t, [y y2.y], label=["Original signal" "Filtered signal"], ylabel="", xlabel="t", legend=:topleft)


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

using SavitzkyGolay
using Test

@testset begin
t = LinRange(-4, 4, 500)
y1 = collect(0:20)
y2 = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))

# Plain call
y1_ = savitzky_golay(y1, 11, 2)
y2_ = savitzky_golay(y2, 21, 4)

# With constructor
f1 = SGolay(11, 2)
y1c_ = f1(y1)

f2 = SGolay(21, 4)
y2c_ = f2(y2)
end

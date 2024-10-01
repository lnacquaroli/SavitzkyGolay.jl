using SavitzkyGolay
using Test

@testset "SGolay input argument checks" begin
    @test_throws ArgumentError SGolay(0,2)                    # window size must be non-zero
    @test_throws ArgumentError SGolay(2,2)                    # window size must be odd
    @test_throws ArgumentError SGolay(3,2)                    # windows size must be > 2+order
    @test_throws ArgumentError SGolay(5,2,-2)                 # derivative must be non-negative
end

    t = LinRange(-4, 4, 500)

@testset "savitzky_golay input argument checks" begin
    @test_throws ArgumentError savitzky_golay(t,0,2)          # window size must be non-zero
    @test_throws ArgumentError savitzky_golay(t,2,2)          # window size must be odd
    @test_throws ArgumentError savitzky_golay(t,3,2)          # windows size must be > 2+order
    @test_throws ArgumentError savitzky_golay(t,5,2;deriv=-1) # derivative must be non-negative
end

@testset "SGolay and savitzky_golay output" begin
    # Plain call
    y1 = collect(0:20)
    y1_ = savitzky_golay(y1, 11, 2)

    y2 = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
    y2_ = savitzky_golay(y2, 21, 4)

    # With constructor
    f1 = SGolay(11, 2)
    y1c_ = f1(y1)

    f2 = SGolay(21, 4)
    y2c_ = f2(y2)

    @test sum(abs.(y1_.y .- y1)) ≈ 0     atol = 1e-12
    @test sum(abs.(y1c_.y .- y1_.y)) ≈ 0 atol = 1e-12

end

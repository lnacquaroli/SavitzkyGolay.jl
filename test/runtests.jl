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

# Plain call
y1 = collect(0:20)
y1_sg = savitzky_golay(y1, 11, 2)
wts_11 = ones(11)

y2 = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
y2_sg = savitzky_golay(y2, 21, 4)
wts_21 = ones(21)

# With constructor
f1 = SGolay(11, 2)
y1c_ = f1(y1)

f2 = SGolay(21, 4)
y2c_ = f2(y2)

@testset "SGolay and savitzky_golay output" begin
    @test sum(abs.(y1_sg.y .- y1)) ≈ 0     atol = 1e-12
    @test sum(abs.(y1c_.y .- y1_sg.y)) ≈ 0 atol = 1e-12
end

@testset "wts handling checks" begin
    # The following test throws the wrong king of error in
    # running the tests but seems to work from the REPL
    # directly
    #
    @test_throws "wts vector argument unusable" f1(y1, wts_11)
    @test_throws "Missing wts vector as second argument" (SGolay(11,1,0,1.0,haswts=true))(y1)
end


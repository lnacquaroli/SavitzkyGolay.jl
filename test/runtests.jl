using SavitzkyGolay
using Test

@testset "SGolay input argument checks" begin
    @test_throws ArgumentError SGolay(0,2)                    # window size must be non-zero
    @test_throws ArgumentError SGolay(2,2)                    # window size must be odd
    @test_throws ArgumentError SGolay(3,2)                    # windows size must be > 2+order
    @test_throws ArgumentError SGolay(5,2,-2)                 # derivative must be non-negative
    @test_throws ArgumentError SGolay(5,2,3)                  # derivative must be ≤ order
end

t = LinRange(-4, 4, 500)

@testset "savitzky_golay input argument checks" begin
    @test_throws ArgumentError savitzky_golay(t,0,2)          # window size must be non-zero
    @test_throws ArgumentError savitzky_golay(t,2,2)          # window size must be odd
    @test_throws ArgumentError savitzky_golay(t,3,2)          # windows size must be > 2+order
    @test_throws ArgumentError savitzky_golay(t,5,2;deriv=-1) # derivative must be non-negative
    @test_throws ArgumentError savitzky_golay(t,5,2;deriv=3)  # derivative must be ≤ order
end

# Plain call
y1 = collect(0:20)
y1_sg = savitzky_golay(y1, 11, 2)

wts_10 = ones(10)
wts_11 = ones(11)
wts_11x = [1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

y2 = exp.(-t.^2) .+ 0.05 .* (1.0 .+ randn(length(t)))
y2_sg = savitzky_golay(y2, 21, 4)

wts_21 = ones(21)

tri_11 = Float64.(vcat( 1:6, 5:-1:1 ))
tri_21 = Float64.(vcat( 1:11, 10:-1:1 ))

# With constructor
f1 = SGolay(11, 2)
y1c_ = f1(y1)

f2 = SGolay(21, 4)
y2c_ = f2(y2)

@testset "SGolay and savitzky_golay output" begin
    @test sum(abs.(y1_sg.y .- y1)) ≈ 0     atol = 1e-12
    @test sum(abs.(y1c_.y .- y1_sg.y)) ≈ 0 atol = 1e-12
end

@testset "wts argument checking" begin
    @test_throws "wts vector argument unusable" f1(y1, wts_11)
    @test_throws "Missing wts vector as second argument" (SGolay(11,1,0,1.0,haswts=true))(y1)
    @test_throws "wts vector length must equal window size"  savitzky_golay(y1, wts_10, 11, 2)
    @test_throws "wts vector must be positive" savitzky_golay(y1, wts_11x, 11, 2)
end

@testset "weighted SG check" begin
    y1_sg_w1 = savitzky_golay(y1,wts_11,11,2)
    @test sum(abs.(y1_sg.y .- y1_sg_w1.y)) ≈ 0  atol = 1e-12
    y1_sg_tri = savitzky_golay(y1, tri_11, 11, 2)
    y2_sg_tri = savitzky_golay(y2, tri_21, 21, 4)
end


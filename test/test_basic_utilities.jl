using Test

@testset "central_diff" begin
    # Simple evenly spaced values should give the same spacing back.
    @test GSG.central_diff([1.0, 3.0, 5.0, 7.0]) == [2.0, 2.0, 2.0, 2.0]

    # Uneven spacing should use forward/backward differences at the edges and
    # the mean of adjacent differences in the interior.
    @test GSG.central_diff([1.0, 2.0, 4.0, 7.0]) == [1.0, 1.5, 2.5, 3.0]
end

@testset "non_nan_indices" begin
    values = [NaN, 1.0, NaN, 2.5, 3.5]

    @test GSG.non_nan_indices(values) == [2, 4, 5]
end

@testset "non_nan_index_pairs" begin
    values = [NaN, 1.0, 2.0, NaN, 5.0]
    sigma = [0.1, NaN, 0.3, 0.4, 0.5]

    @test GSG.non_nan_index_pairs(values, sigma) == [3, 5]
end

@testset "search_z_func" begin
    @test GSG.search_z_func(0) == 50
    @test GSG.search_z_func(499.9) == 50
    @test GSG.search_z_func(500) == 100
    @test GSG.search_z_func(1500) == 250
    @test GSG.search_z_func(2500) == 1000
end

@testset "grid_horz_dist" begin
    # Duplicate station locations should be collapsed before the distance is computed.
    latitudes = [0.0, 0.0, 1.0, 1.0, 2.0]
    longitudes = [0.0, 0.0, 0.0, 0.0, 0.0]
    ll_grid = [0.0, 1.0, 2.0, 3.0]

    expected = fill(83.4, length(ll_grid))

    @test isapprox(GSG.grid_horz_dist(latitudes, longitudes, ll_grid), expected, atol=1e-6)
end

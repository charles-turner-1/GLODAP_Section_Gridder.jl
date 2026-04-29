using Test
using Interpolations

@testset "corrlen_scale_factor" begin
    @test GSG.corrlen_scale_factor(horzCoordinate="latitude", meanLatitude=45.0) == 111.2

    expected = 111.2 * cos(pi / 3)
    @test GSG.corrlen_scale_factor(horzCoordinate="longitude", meanLatitude=60.0) ≈ expected

    # latitudes intentionally overwrite meanLatitude, matching legacy behaviour.
    @test GSG.corrlen_scale_factor(
        horzCoordinate="longitude",
        latitudes=[0.0, 60.0],
        meanLatitude=0.0,
    ) ≈ 111.2 * cos(pi / 6)

    @test_throws ErrorException GSG.corrlen_scale_factor(horzCoordinate="depth", meanLatitude=0.0)
end

@testset "horzCorrDistanceKilometres" begin
    @test GSG.horzCorrDistanceKilometres(
        [1.0, 2.0];
        meanLatitude=60.0,
        horzCoordinate="longitude",
    ) ≈ [111.2 * cos(pi / 3), 2 * 111.2 * cos(pi / 3)]

    @test GSG.horzCorrDistanceKilometres(
        [1.0, 2.0];
        meanLatitude=30.0,
        horzCoordinate="latitude",
    ) == [111.2, 222.4]
end

@testset "corrlen_mean_latitude" begin
    @test GSG.corrlen_mean_latitude([10.0, NaN, 20.0]) == 15.0
end

@testset "corrlen_finalize_lenz" begin
    @test GSG.corrlen_finalize_lenz(10_000, [1.0, 2.0, 3.0]) == fill(1000, 3)
    @test GSG.corrlen_finalize_lenz([500.0, 1500.0], [1.0, 2.0]) == [500.0, 1000.0]
end

@testset "corrlen_interpolate_to_grid" begin
    @test GSG.corrlen_interpolate_to_grid([0.0, 2.0], [10.0, 20.0], [0.0, 1.0, 2.0]; extrapolation_bc=NaN) == [10.0, 15.0, 20.0]
    @test GSG.corrlen_interpolate_to_grid([0.0, 2.0], [10.0, 20.0], [0.0, 1.0, 2.0, 3.0]; extrapolation_bc=Flat()) == [10.0, 15.0, 20.0, 20.0]
end

@testset "corrlen_retry_search" begin
    searchf = GSG.corrlen_retry_search(x -> x + 1, 3)
    @test searchf(2) == 9
end

using Test
using DataFrames

@testset "legacy_partial_cruise_deltas" begin
    deltas = GSG.legacy_partial_cruise_deltas([0.0, 5.0, 10.0], [1.0, 4.0, 9.0])

    @test deltas.ΔminVal == -1.0
    @test deltas.ΔmaxVal == 1.0
    @test deltas.minObsVal == 1.0
    @test deltas.maxObsVal == 9.0
end

@testset "legacy_truncate_mask" begin
    mask = trues(4, 3)
    horz_grid = [0.0, 5.0, 10.0]

    truncated = GSG.legacy_truncate_mask(mask, horz_grid, 2.0, 8.0)

    @test truncated == Bool[
        0 1 0
        0 1 0
        0 1 0
        0 1 0
    ]
    @test mask == trues(4, 3)
end

@testset "legacy_variable_exception_good_indices" begin
    station = [1.0, 1.0, 2.0, 2.0, 3.0]
    pressure = [10.0, 60.0, 20.0, 80.0, 40.0]
    station_list = [1, 2]
    min_pressure_list = [0, 10]
    max_pressure_list = [50, 70]
    col_idx = [1, 2]

    good_idx = GSG.legacy_variable_exception_good_indices(
        5,
        station,
        pressure,
        station_list,
        min_pressure_list,
        max_pressure_list,
        col_idx,
    )

    @test good_idx == [2, 4, 5]
end

@testset "legacy_horzlen_factor" begin
    exception_df = DataFrame(
        Expocode=["EXPO1", "EXPO2"],
        Variable=["G2salinity", "G2oxygen"],
        Gridding=["isobaric", "isopycnic"],
        Factor=[1.5, 2.5],
    )

    @test GSG.legacy_horzlen_factor(exception_df, "EXPO1", "G2salinity", "isobaric") == 1.5
    @test GSG.legacy_horzlen_factor(exception_df, "missing", "G2salinity", "isobaric") == 1
    @test GSG.legacy_horzlen_factor(exception_df, "missing", "G2salinity", "isopycnic") == 1.0
end

@testset "legacy_obs_xvals validation" begin
    @test GSG.legacy_obs_xvals("latitude", [0.0, 1.0], [10.0, 11.0], [120.0, 121.0]) == [10.0, 11.0]
    @test GSG.legacy_obs_xvals("longitude", [0.0, 1.0], [10.0, 11.0], [120.0, 121.0]) == [120.0, 121.0]
    @test_throws ErrorException GSG.legacy_obs_xvals("depth", [0.0, 1.0], [10.0, 11.0], [120.0, 121.0])
end

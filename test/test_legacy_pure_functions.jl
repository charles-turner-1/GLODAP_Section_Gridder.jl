using Test
using DataFrames
using CSV

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

@testset "partial cruise wrappers" begin
    @test GSG.checkPartialCruise(
        [0.0, 5.0, 10.0];
        horzCoordinate="latitude",
        obsLat=[1.0, 4.0, 9.0],
        obsLon=[120.0, 121.0, 122.0],
    ) == false

    @test GSG.checkPartialCruise(
        [0.0, 5.0, 10.0];
        horzCoordinate="latitude",
        obsLat=[4.0, 5.0, 6.0],
        obsLon=[120.0, 121.0, 122.0],
    ) == true

    mask = Matrix{Bool}(trues(4, 3))
    truncated = GSG.maskPartialCruise(
        mask;
        obsLat=[10.0, 11.0],
        obsLon=[2.0, 8.0],
        horzGrid=[0.0, 5.0, 10.0],
        horzCoordinate="longitude",
    )

    @test truncated == Bool[
        0 1 0
        0 1 0
        0 1 0
        0 1 0
    ]
end

@testset "exception-handling wrappers" begin
    mktempdir() do tmpdir
        horzlen_path = joinpath(tmpdir, "horzlen.csv")
        CSV.write(
            horzlen_path,
            DataFrame(
                Expocode=["EXPO1"],
                Variable=["G2salinity"],
                Gridding=["isobaric"],
                Factor=[1.5],
            ),
        )

        @test GSG.checkHorzLenFactor(
            expocode="EXPO1",
            variableName="G2salinity",
            griddingType="isobaric",
            HORZLEN_EXCEPTIONS=horzlen_path,
        ) == 1.5

        variable_exceptions_path = joinpath(tmpdir, "variable_exceptions.csv")
        write(
            variable_exceptions_path,
            "expocode,variable,station,minPressure,maxPressure\n" *
            "EXPO1,G2oxygen,1,0,50\n" *
            "EXPO1,G2oxygen,2,10,70\n",
        )

        @test GSG.checkVariableExceptions(
            expocode="EXPO1",
            variableName="G2oxygen",
            variable=[1.0, 2.0, 3.0, 4.0, 5.0],
            station=[1.0, 1.0, 2.0, 2.0, 3.0],
            pressure=[10.0, 60.0, 20.0, 80.0, 40.0],
            EXCEPTIONS_DIR=tmpdir,
            EXCEPTIONS_FILENAME="variable_exceptions.csv",
        ) == [2, 4, 5]
    end
end

@testset "legacy_obs_xvals validation" begin
    @test GSG.legacy_obs_xvals("latitude", [0.0, 1.0], [10.0, 11.0], [120.0, 121.0]) == [10.0, 11.0]
    @test GSG.legacy_obs_xvals("longitude", [0.0, 1.0], [10.0, 11.0], [120.0, 121.0]) == [120.0, 121.0]
    @test_throws ErrorException GSG.legacy_obs_xvals("depth", [0.0, 1.0], [10.0, 11.0], [120.0, 121.0])
end

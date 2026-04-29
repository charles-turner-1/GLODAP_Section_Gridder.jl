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

function GSG.fithorzlen(_x::Tuple{Vector{Int},Vector{Int},Vector{Int}}, data::Vector{Int}, pr_grid::Vector{Int}; searchz, kwargs...)
    value = searchz(2)
    value < 4 && error("retry")
    return ([value], nothing)
end

function GSG.fitvertlen(_x::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}, data::Vector{Float64}, pr_grid::Vector{Float64}; kwargs...)
    if haskey(kwargs, :smoothz)
        return ([0.2, 0.4], nothing)
    end
    return ([200.0, 1200.0], nothing)
end

function GSG.fithorzlen(_x::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}, data::Vector{Float64}, pr_grid::Vector{Float64}; kwargs...)
    if haskey(kwargs, :smoothz)
        error("force prescribed density-space lenx")
    end
    return ([1.0, 3.0], nothing)
end

@testset "fithorzlen_increasing_search" begin
    lenx, _ = GSG.fithorzlen_increasing_search(([1], [2], [3]), [4], [5]; searchz=x -> x, multiplier=1)
    @test lenx == [4]

    @test_throws ErrorException GSG.fithorzlen_increasing_search(nothing, nothing, nothing; searchz=x -> x, multiplier=16)
end

@testset "calcCorrLengths wrapper" begin
    variable = [1.0, 2.0, 3.0]
    lenz, lenxkm = GSG.calcCorrLengths(
        variable=variable,
        obsLat=[10.0, 20.0, 30.0],
        obsLon=variable,
        obsPres=[0.0, 10.0, 20.0],
        presGrid=[0.0, 10.0, 20.0],
        pressureStepNumber=2,
        lenxFactor=0,
    )

    @test length(lenz) == 3
    @test length(lenxkm) == 3
    @test all(isfinite, lenz)
    @test all(==(0.0), lenxkm)
end

@testset "calcDensityCorrLengths wrapper" begin
    lenz, lenxkm = GSG.calcDensityCorrLengths(
        [1.0, 2.0, 3.0];
        obsLat=[10.0, 20.0, 30.0],
        obsLon=[100.0, 101.0, 102.0],
        obsSigma=[0.0, 1.0, 2.0],
        sigGrid=[0.0, 1.0, 2.0],
        lenxPrescribed=2.0,
        sigmaStepNumber=2,
    )

    @test lenz ≈ [0.2, 0.3, 0.4]
    @test length(lenxkm) == 3
    @test lenxkm[1] ≈ lenxkm[2]
    @test lenxkm[2] ≈ lenxkm[3]
end

using Test
using ZipFile

@testset "rm_flagged_data" begin
    data = [1.0, 2.0, 3.0, 4.0, 5.0]
    flags = [0, 1, 1, 0, 1]
    expected_result = [1.0, NaN, NaN, 4.0, NaN]

    result = GSG.rm_flagged_data(data, flags)

    @test isapprox(result, expected_result, nans=true)

    # Test with all flagged data.
    data = [1.0, 2.0, 3.0, 4.0, 5.0]
    flags = [2, 1, 1, 0, 1]
    expected_result = [1.0, NaN, NaN, NaN, NaN]

    result = GSG.rm_flagged_data(data, flags)

    @test isapprox(result, expected_result, nans=true)

    # Test with no bad data.
    data = [1.0, 2.0, 3.0]
    flags = [0, 0, 0]
    expected_result = [1.0, 2.0, 3.0]

    result = GSG.rm_flagged_data(data, flags)

    @test result == expected_result
end

@testset "adjust_tco2" begin
    # We only test with the default adjustment table bundled in the repo.
    @test GSG.adjust_tco2(expocode="dud value") == 0.0
    @test GSG.adjust_tco2(expocode="06AQ20021124") == 0.0
    @test GSG.adjust_tco2(expocode="06AQ19860627") == 12.0
    @test GSG.adjust_tco2(expocode="06AQ19840719") == 0.0
end

@testset "load_section_info" begin
    # Smoke test with one real bundled section mask so CI exercises the current loader path.
    section_info = GSG.load_section_info("A05")

    # Check the results.
    @test length(section_info.ll_grid) == 670
    @test length(section_info.pr_grid) == 651
    @test size(section_info.mask) == (670, 651)
    @test sum(section_info.mask) == 322623
    @test section_info.horz_coord == "longitude"
end

@testset "resolve_glodap_db_path uses mocked env var" begin
    mktempdir() do tmpdir
        csv_path = joinpath(tmpdir, "fixture_glodap.csv")

        open(csv_path, "w") do io
            write(io, "G2expocode,G2cruise,G2theta,G2pressure,G2longitude\n")
            write(io, "33RO19980123,1,1.5,10,-30\n")
            write(io, "33RO19980123,1,1.6,20,-31\n")
            write(io, "OTHER,2,9.9,999,999\n")
        end

        withenv("GLODAP_DB" => csv_path) do
            resolved = GSG.resolve_glodap_db_path(nothing; allow_download=false)
            @test resolved == csv_path

            vars = GSG.load_glodap_vars(["G2theta", "G2pressure", "G2longitude"], "33RO19980123", resolved)
            @test names(vars) == ["G2theta", "G2pressure", "G2longitude"]
            @test size(vars, 1) == 2
            @test vars[!, "G2theta"] == [1.5, 1.6]
            @test vars[!, "G2pressure"] == [10, 20]
        end
    end
end

@testset "resolve_glodap_db_path uses stubbed download source" begin
    mktempdir() do tmpdir
        csv_name = "fixture_glodap.csv"
        csv_path = joinpath(tmpdir, csv_name)
        zip_path = joinpath(tmpdir, "fixture_glodap.zip")

        open(csv_path, "w") do io
            write(io, "G2expocode,G2cruise,G2theta,G2pressure,G2longitude\n")
            write(io, "33RO19980123,1,1.5,10,-30\n")
            write(io, "33RO19980123,1,1.6,20,-31\n")
            write(io, "OTHER,2,9.9,999,999\n")
        end

        writer = ZipFile.Writer(zip_path)
        try
            file = ZipFile.addfile(writer, csv_name)
            write(file, read(csv_path))
        finally
            close(writer)
        end

        cache_dir = joinpath(tmpdir, "cache")
        resolved = withenv("GLODAP_DB" => nothing, "GLODAP_DB_URL" => nothing, "GLODAP_DB_FILENAME" => nothing) do
            GSG.resolve_glodap_db_path(
                nothing;
                cache_dir=cache_dir,
                filename=csv_name,
                download_url="file://$(zip_path)",
            )
        end

        @test isfile(resolved)
        @test resolved == joinpath(cache_dir, csv_name)
    end
end

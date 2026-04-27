using Test

@testset "rm_flagged_data" begin
    data = [1.0, 2.0, 3.0, 4.0, 5.0]
    flags = [0, 1, 1, 0, 1]
    expected_result = [1.0, NaN, NaN, 4.0, NaN]

    result = GSG.rm_flagged_data(data, flags)

    @test isapprox(result, expected_result, nans=true)

    data = [1.0, 2.0, 3.0, 4.0, 5.0]
    flags = [2, 1, 1, 0, 1]
    expected_result = [1.0, NaN, NaN, NaN, NaN]

    result = GSG.rm_flagged_data(data, flags)

    @test isapprox(result, expected_result, nans=true)

    data = [1.0, 2.0, 3.0]
    flags = [0, 0, 0]
    expected_result = [1.0, 2.0, 3.0]

    result = GSG.rm_flagged_data(data, flags)

    @test result == expected_result
end

@testset "adjust_tco2" begin
    @test GSG.adjust_tco2(expocode="dud value") == 0.0
    @test GSG.adjust_tco2(expocode="06AQ20021124") == 0.0
    @test GSG.adjust_tco2(expocode="06AQ19860627") == 12.0
    @test GSG.adjust_tco2(expocode="06AQ19840719") == 0.0
end

@testset "load_section_info" begin
    section_info = GSG.load_section_info("A05")

    @test length(section_info.ll_grid) == 670
    @test length(section_info.pr_grid) == 651
    @test size(section_info.mask) == (670, 651)
    @test sum(section_info.mask) == 322623
    @test section_info.horz_coord == "longitude"
end

@testset "load_glodap_vars placeholder" begin
    @test true
end

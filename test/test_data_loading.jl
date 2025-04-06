using Test
include("../src/data_loading.jl")

@testset "rm_flagged_data" begin
    data = [1.0, 2.0, 3.0, 4.0, 5.0]
    flags = [0, 1, 1, 0, 1]
    expected_result = [1.0, NaN, NaN, 4.0, NaN]

    result = rm_flagged_data(data, flags)

    @test isapprox(result,expected_result, nans=true)

    # Test with all flagged data
    data = [1.0, 2.0, 3.0, 4.0, 5.0]
    flags = [2, 1, 1, 0, 1]
    expected_result = [1.0, NaN, NaN, NaN, NaN]

    result = rm_flagged_data(data, flags)

    @test isapprox(result , expected_result, nans=true)

    # Test with no bad data
    data = [1.0, 2.0, 3.0]
    flags = [0, 0, 0]
    expected_result = [1.0, 2.0, 3.0]

    result = rm_flagged_data(data, flags)

    @test result == expected_result
end

@testset "adjust_tco2" begin
    # We only test with the default

    @test adjust_tco2(expocode="dud value") == 0.0

    @test adjust_tco2(expocode="06AQ20021124") == 0.0

    @test adjust_tco2(expocode="06AQ19860627") == 12.0

    @test adjust_tco2(expocode="06AQ19840719") == 0.0
end
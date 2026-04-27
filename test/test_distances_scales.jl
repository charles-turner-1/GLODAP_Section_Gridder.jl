using Test

@testset "grid_vert_dist Tests" begin
    pr_grid = [0.0, 90.0, 180.0, 270.0]

    @test GSG.grid_vert_dist(pr_grid) == [10, 10, 10, 10]
end

@testset "create_sigma_grid Tests" begin
    x = collect(0:6000)
    gridded_sig = GSG.create_sigma_grid(x)
    expected = collect(0:10:6000)

    @test gridded_sig == expected

    gridded_sig = GSG.create_sigma_grid(x, 1200)
    expected = collect(0:5:6000)

    @test gridded_sig == expected

    x = collect(0:6001)
    gridded_sig = GSG.create_sigma_grid(x)
    expected = collect(0:10:6000)

    @test gridded_sig == expected

    x = collect(0:0.1:59.9)
    expected = collect(0:0.1:59.9)
    gridded_sig = GSG.create_sigma_grid(x)

    @test gridded_sig == expected
end

@testset "grid_sigma_distance Tests" begin
    sigma_grid = [0.0, 10.0, 20.0, 30.0]

    @test GSG.grid_sigma_distance(sigma_grid) == fill(7.5, 4)
end

@testset "calc_scale_factors Tests" begin
    vert_dist = [10.0, 20.0]
    horz_dist = [100.0, 200.0, 400.0]

    scale_factors = GSG.calc_scale_factors(vert_dist, horz_dist; print_scales=false)

    @test size(scale_factors.vert) == (2, 3)
    @test size(scale_factors.horz) == (2, 3)
    @test scale_factors.vert[:, 1] == [0.1, 0.05]
    @test scale_factors.horz[1, :] == [0.01, 0.005, 0.0025]
end

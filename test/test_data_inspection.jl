using Test

@testset "list_section_expocodes" begin
    expocodes = GSG.list_section_expocodes("A05")

    @test names(expocodes) == ["Year", "GLODAP Expocode", "GO-SHIP Expocode"]
    @test size(expocodes, 1) > 0
    @test expocodes[1, "GLODAP Expocode"] == "29HE19920714"

    # Cover the manual expocode_dir branch as well as the default packaged-data path.
    expocodes_from_dir = GSG.list_section_expocodes("A05", joinpath(GSG.root, "data", "SectionExpocodes"))
    @test expocodes_from_dir == expocodes
end

using Test

@testset "list_section_expocodes" begin
    expocodes = GSG.list_section_expocodes("A05")

    @test names(expocodes) == ["Year", "GLODAP Expocode", "GO-SHIP Expocode"]
    @test nrow(expocodes) > 0
    @test expocodes[1, "GLODAP Expocode"] == "29HE19920714"
end

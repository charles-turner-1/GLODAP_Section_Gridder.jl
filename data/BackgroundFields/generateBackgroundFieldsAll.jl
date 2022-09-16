include("../../src/DIVA_WrapperFunctions.jl")
using Plots

const GLODAP_DIR = "/home/ct/MATLAB/GLODAP"
const GOSHIP_DIR = "/home/ct/MATLAB/GO_SHIP"
const MASK_MATFILE = "/home/ct/Julia/GLODAP_Easy_Ocean/GOSHIP_MaskStruct.mat"

const gridDir, repDir, convDir = loadGoShipDirectories(GOSHIP_DIR)
##

variablesToGrid = ["G2theta","G2salinity","G2tco2"]

horzGridFile = "/home/ct/Julia/GLODAP_Easy_Ocean/SectionDataFiles/horzCoords.csv"

horzCoords = CSV.read(horzGridFile,DataFrame)
dropmissing!(horzCoords)

elements = collect(eachcol(horzCoords))
startVal=30
sectionNames = elements[1][startVal:end]
horzCoordinates = elements[2][startVal:end]


##
for sectionName in enumerate(reverse(sectionNames)), variable in variablesToGrid
    secName = String(sectionName[2])
    horzCoordinate = String(horzCoordinates[sectionName[1]])
    llGrid, prGrid, _ = loadSectionInfo(MASK_MATFILE,secName,gridDir)

    titleStr = "Background field: " * variable * ", " * secName

    # Test if the background field file exists and skip it if so.
    fileStr =  "/home/ct/Julia/GLODAP_Easy_Ocean/BackgroundFields/" *secName *
     "_" * variable * ".nc"
     isfile(fileStr) ? continue : nothing

    numSectionExpocodes = size(listSectionExpocodes(secName,convDir))[1]

    numSectionExpocodes > 0 ? nothing : continue


    bgField = createBackgroundField(variable=variable,sectionName=secName
                                   ,horzCoordinate=horzCoordinate)

    display(heatmap(llGrid,prGrid,bgField,yflip=true,c=:jet1,title=titleStr
    ,ylabel="Pressure", xlabel=horzCoordinate))
    println("Are you happy with this field? Y/n")
    fieldIsGood = readline(stdin)
    while fieldIsGood == "n"
        println("Specify new horizontal correlation length factor (default = 100): ")
        horzCorrLen = readline(stdin); horzCorrLen = parse(Float64,horzCorrLen)
        println("Specify new vertical correlation length factor (default = 1): ")
        vertCorrLen = readline(stdin); vertCorrLen = parse(Float64,vertCorrLen)

        bgField = createBackgroundField(variable=variable,sectionName=secName
                                       ,horzCoordinate=horzCoordinate
                                       ,horzLenFactor=horzCorrLen
                                       ,vertLenFactor=vertCorrLen)

        display(heatmap(llGrid,prGrid,bgField,yflip=true,c=:jet1,title=titleStr
        ,ylabel="Pressure", xlabel=horzCoordinate))
        println("Are you happy with this field? Y/n")
        fieldIsGood = readline(stdin)

    end
    willWrite = "n"
    println("Confirm you wish to write background field: y/N")
    willWrite = readline(stdin)
    if willWrite == "y"
        writeBackgroundField(variable=bgField,sectionName=secName,variableName=variable)
    end
end

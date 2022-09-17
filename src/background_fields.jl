function createBackgroundField(;variable::String,sectionName::String
                              ,GOSHIP_DIR::String="/home/ct/MATLAB/GO_SHIP"
                              ,GLODAP_DIR::String="/home/ct/MATLAB/GLODAP"
                              ,horzCoordinate::String
                              ,horzLenFactor::Number=100
                              ,vertLenFactor::Number=1
                              ,MASK_MATFILE::String="/home/ct/Julia/GLODAP_Easy_Ocean/GOSHIP_MaskStruct.mat")
    # Create a background field using GLODAP climatology data
    if variable in ["G2tco2","DIC","TCO2"]
        climatology = Dataset("/home/ct/Julia/GLODAP_Easy_Ocean/Climatologies/GLODAPv2.2016b.TCO2.nc")
        clim_var = climatology["TCO2"][:]
    elseif variable in ["G2theta","Temperature","TMP"]
        climatology = Dataset("/home/ct/Julia/GLODAP_Easy_Ocean/Climatologies/GLODAPv2.2016b.temperature.nc")
        clim_var = climatology["temperature"][:]
    elseif variable in ["G2salinity","Salinity","SAL"]
        climatology = Dataset("/home/ct/Julia/GLODAP_Easy_Ocean/Climatologies/GLODAPv2.2016b.salinity.nc")
        clim_var = climatology["salinity"][:]
    end

    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    replace!(clim_var,missing=>NaN)

    clim_lon = climatology["lon"][:]
    clim_lat = climatology["lat"][:]
    clim_dep = climatology["Depth"][:]

    llGrid, prGrid, mask = loadSectionInfo(MASK_MATFILE,sectionName,gridDir)
    _, _, convDir = loadGoShipDirectories(GOSHIP_DIR)
    expocodes = listSectionExpocodes(sectionName,convDir)[:,"GLODAP Expocode"]

    vars = loadGLODAPVariables(["G2theta","G2tco2","G2salinity"
    ,"G2gamma","G2latitude","G2longitude","G2pressure"]
    ,GLODAP_DIR,GLODAP_expocodes=expocodes)



    lon = vars["G2longitude"]
    lat = vars["G2latitude"]

    if sectionName[1] == 'P'
         lon = matchLonConventions(clim_lon,lon)
         lonOffset = 0
    else
         lonOffset = 200
    end

    ll = unique!(DataFrame(lon=lon,lat=lat)) # Find unique locations, ie. stations
    obsLonVals = ll[!,"lon"]
    obsLatVals = ll[!,"lat"]

    cLonGrid = fill(0.0,length(clim_lon),length(obsLonVals))
    [cLonGrid[:,i] = clim_lon .- obsLonVals[i] .- lonOffset for i in eachindex(obsLonVals)]

    k = abs.(cLonGrid)
    lonIdx = fill(0,length(obsLonVals))
    [lonIdx[i] = argmin(k[:,i]) for i in 1:length(obsLonVals)]

    cLatGrid = fill(0.0,length(clim_lat),length(obsLatVals))
    [cLatGrid[:,i] = clim_lat .- obsLatVals[i] for i in eachindex(obsLatVals)]

    k = abs.(cLatGrid)
    latIdx = fill(0,length(obsLatVals))
    [latIdx[i] = argmin(k[:,i]) for i in eachindex(obsLatVals)]


    llG = unique!(DataFrame(lon=lonIdx,lat=latIdx)) # Find unique locations, ie. grid cells
    lons = llG[!,"lon"]
    lats = llG[!,"lat"]
    # These are the same length so we can continue to just use lons as a counter

    # Grid is on 33 levels
    bgField = convert(Matrix{Union{Float64,Missing}},fill(NaN,length(lons),33))

    [bgField[i,:] = clim_var[lons[i],lats[i],:] for i in eachindex(lons)]

    #return (clim_lon[lons],clim_lat[lats],bgField)

    oLonGrid = fill(NaN,size(bgField))
    oLatGrid = fill(NaN,size(bgField))
    oPrsGrid = fill(NaN,size(bgField))

    [oLonGrid[:,i] = clim_lon[lons] for i in 1:length(clim_dep)]
    [oLatGrid[:,i] = clim_lat[lats] for i in 1:length(clim_dep)]
    [oPrsGrid[i,:] = clim_dep for i in 1:length(lons)]

    zVar = bgField[:];  goodIdx = findNonNaNIndices(zVar)
    zVar = convert(Vector{Float64},zVar)[goodIdx]

    w = oLatGrid[:][goodIdx]
    x = oLonGrid[:][goodIdx]
    y = oPrsGrid[:][goodIdx]

    P_grid,L_grid = ndgrid(clim_dep,lons)

    horzDist = gridHorzDistance(w,x,llGrid)
    vertDist = gridVertDistance(prGrid)
    scaleVert, scaleHorz = calcScaleFactors(vertDist,horzDist)

    (clVert, clHorz) = calcCorrLengths(zVar,obsLat=w,obsLon=x
    ,obsPres=y,presGrid=prGrid)

    # There should be some sort of call to matchLonConventions here I think
    # instead of just subtracting 200 from all longitudes carelessly to undo our fuckery.
    if horzCoordinate == "longitude"
        divaGridClim = easyDIVAGrid(variable=zVar,vertVar=y,latLon=x.-lonOffset
        ,horzCoordinate=horzCoordinate,meanValue="scalar",vertGrid=prGrid,horzGrid=llGrid
        ,horzScale=scaleHorz,vertScale=scaleVert,mask=mask,Epsilon=0.1,horzCorrLength=horzLenFactor*clHorz
        ,vertCorrLength = clVert*vertLenFactor)
    elseif horzCoordinate == "latitude"
        divaGridClim = easyDIVAGrid(variable=zVar,vertVar=y,latLon=w
        ,horzCoordinate=horzCoordinate,meanValue="scalar",vertGrid=prGrid,horzGrid=llGrid
        ,horzScale=scaleHorz,vertScale=scaleVert,mask=mask,Epsilon=0.1,horzCorrLength=horzLenFactor*clHorz
        ,vertCorrLength = clVert*vertLenFactor)
    end

    return divaGridClim

end

function writeBackgroundField(;variable::Matrix{Float64},sectionName::String
                             ,variableName::String
                             ,outputPathStr::String="/home/ct/Julia/GLODAP_Easy_Ocean/BackgroundFields/")
    # This function takes a background field on the same grid as the mask we
    # plan on using and writes it to disk.

    # I've only implemented functionality to use a climatological background field
    # for DIC, salinity and temperature. As a result, we can only write if the
    # variable name is one of these. For convenience, I'm going to force them to
    # be written in the GLODAP nomenculature ie. G2tco2.
    if variableName ∉ ["G2tco2","G2theta","G2salinity"]
        error("variableName must be either \"G2tco2\", \"G2theta\" or \"G2salinity\"")
    end

    fileName = sectionName * "_" * variableName * ".nc"

    fileNameFull = joinpath(outputPathStr,fileName)

    ds = Dataset(fileNameFull,"c")

    defDim(ds,"Vertical_Cell_No",size(variable,1))
    defDim(ds,"Horizontal_Cell_No",size(variable,2))

    titleStr= "Background field for " * variableName * " for " * sectionName
    ds.attrib["Title"] = titleStr

    v1 = defVar(ds,variableName,Float64,("Horizontal_Cell_No","Vertical_Cell_No"))

    v1[:,:] = variable'

    ds.attrib["Comments"] = "Generated using the GLODAP Climatology"
    ds.attrib["Creator"] = "Charles Turner"

    close(ds)

    return nothing

end

function readBackgroundField(;sectionName::String,variableName::String)
    # This function will go looking for the background field and load it
    if variableName != "G2tco2"  && variableName != "G2theta" && variableName != "G2salinity"
        error("variableName must be either \"G2tco2\", \"G2theta\" or \"G2salinity\"")
    end

    fileName = sectionName * "_" * variableName * ".nc"
    pathStr = "/home/ct/Julia/GLODAP_Easy_Ocean/BackgroundFields/"

    fileNameFull = joinpath(pathStr,fileName)

    ds = Dataset(fileNameFull)

    field = ds[variableName][:]
    field = convert(Matrix{Float64},field')

    return field

end


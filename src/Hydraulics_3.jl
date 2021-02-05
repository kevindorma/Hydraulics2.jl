# November, 2020
# Kevin Dorma
# module for common hydraulic calculations
# rev 0, December 2020
# rev 1, February 2021, added network hydraulics


# I need to make some decisions about how to interace with these functions
# I assume that I use dataframes to store information about sizing lines and general hydraulics
# then I will execute a function and it will return a result, but not change the raw information
# size lines
    # Input dataframe for describing the lines
    # Output line descriptor, pipe schedule, size criteria DP/100 and C, IDmm for each criteria and the chosen IDmm
# pick NPS will return the next size larger for NPS.This will be the chosen size.
# we will have utility functions for increasing or decreasing the NPS by 1. 


module Hydraulics

using DataFrames
using CSV
using ExcelFiles # more convenient than CSV for user input

export calcMoodyF
export addPipeProperties, addFluidProperties, getReynolds, getSegmentDP
export getLineSize, getListLargerNPS, getLargerNPS
# and the network hydraulic files
export indexLookup, doNetworkHydraulics, extractResults, combineResults

# these are the reference data.
# and this needs to be loaded correctly when we use the code in a module
include("npsList.jl")
include("schedList.jl")
include("pipeRoughnessList.jl")
include("fitting3K.jl")
include("pipeIDlist.jl")

#schedList = CSV.File("schedList.csv") |> DataFrame
#pipeRoughness = CSV.File("pipeRoughness.csv") |> DataFrame
#fitting3K = CSV.File("fitting3K.csv") |> DataFrame
#IDmm = CSV.File("pipeIDlist.csv") |> DataFrame
## not used
#pipeTable = CSV.read("pipeIDtable.csv");


function packageInfo()
    # return information about the package as a string
    return("Hydraulics package, Kevin Dorma. Written in Julia, rev 1 February 2021.")
end

# need function to append piping data
# later


function calcMoodyF(Reynolds,eD)
    # Moody friction factor, this is Halaand
    # eD is e/D relative roughness
    # Reynolds is Reynolds number
    invSqrtF = -1.8 .* log10.((eD ./ 3.70) .^ 1.11 + 6.9 ./ Reynolds)
    moodyF = (1 ./ (invSqrtF .^2) )
    return (moodyF)
end

function checkFittingList(fittingList)
    # fittingList is our list of fittings, fittingReference is the generic data
    # go through the list of fittings (fittingList) and compare with the reference list (possibly our fitting3K global list)
    # return items that are not found
    errorList = DataFrame(message = String[], entry=Int64[], Segment=String[], fittingType=String[])
    for i = 1:(size(fittingList)[1])
        thisFitting = fittingList[i,:fittingType]
        match = fitting3K[fitting3K.fittingType .== thisFitting,:]
        if ((size(match)[1]) == 0)
            push!(errorList, ("Fitting not found in row", i, fittingList[i,:Segment], fittingList[i,:fittingType] ))
        end
    end
    return (errorList)
end



function checkLineList(lines,fluidList)
    # go through our list of lines and compare with our reference lists
    # return the line items that are not found
    errorList = DataFrame(message = String[], entry=Int64[], Segment=String[], NPS=Float64[], Schedule=String[], material=String[], fluidName=String[])
    for i = 1:(size(lines)[1])
        # NPS, must check if this field exists
	if (indexLookup("NPS",names(lines)) > 0)
        thisItem = lines[i,:NPS]
        match = npsList[npsList.NPS .== thisItem,:]
        if ((size(match)[1]) == 0)
            push!(errorList, ("NPS not found in row", i, lines[i,:Segment], lines[i,:NPS], lines[i,:Schedule], lines[i,:Material], lines[i,:fluidName] ));
        end;
	end;
        # schedule
        thisItem = lines[i,:Schedule]
        match = schedList[schedList.Schedule .== thisItem,:]
        if ((size(match)[1]) == 0)
            push!(errorList, ("Schedule not found in row", i, lines[i,:Segment], lines[i,:NPS], lines[i,:Schedule], lines[i,:Material], lines[i,:fluidName] ))
        end
        # roughness material
        thisItem = lines[i,:Material]
        match = pipeRoughness[pipeRoughness.Material .== thisItem,:]
        if ((size(match)[1]) == 0)
            push!(errorList, ("Material not found in row", i, lines[i,:Segment], lines[i,:NPS], lines[i,:Schedule], lines[i,:Material], lines[i,:fluidName] ))
        end
        # fluid
        thisItem = lines[i,:fluidName]
        match = fluidList[fluidList.fluidName .== thisItem,:]
        if ((size(match)[1]) == 0)
            push!(errorList, ("FluidName not found in row", i, lines[i,:Segment], lines[i,:NPS], lines[i,:Schedule], lines[i,:Material], lines[i,:fluidName] ))
        end
    end
    return (errorList)
end


function addPipeProperties(df)
    # roughness, IDmm, eD ratio
    # this appends columns to df
    # used either line sizing or hydraulic calculations

    df[:,:roughnessMM] .= 0.0
    df[:,:IDmm] .= 0.0

    for i = 1:(size(df)[1])
        df[i,:roughnessMM]=pipeRoughness[df[i,:Material] .== pipeRoughness[:,:Material],:roughnessMM][1]
        df[i,:IDmm]= IDmm[(df[i,:NPS] .== IDmm[:,:NPS]) .& (df[i,:Schedule] .== IDmm[:,:Schedule]),:Idmm][1]
    end

    return (0.0)
end

function addFluidProperties(df,fluidList)
    # extract the density, viscosity and whatever else is needed from the fluidList
    # Add this to the df
    # used either line sizing or hydraulic calculations

    

    df[:,:rho_kgm3] .= 0.0
    df[:,:mu_mPas] .= 0.0
    df[:,:roughnessMM] .= 0.0


    for i = 1:(size(df)[1])
        df[i,:rho_kgm3] = fluidList[df[i,:fluidName] .== fluidList[:,:fluidName], :rho_kgm3][1]
        df[i,:mu_mPas]  = fluidList[df[i,:fluidName] .== fluidList[:,:fluidName], :mu_mPas][1]
        df[i,:roughnessMM]  = pipeRoughness[df[i,:Material] .== pipeRoughness[:,:Material], :roughnessMM][1]
    end

    return (0.0)
end



function lineSizeDP(df)
    # this is for sizing lines based on pressure drop kPa per 100 m, return the ID in mm
    # very convenient form straight out of Perrys handbook
    # uses friction factor correlation from Chen
    # Do not modify the original data in df

    g = 9.80665
    Sf = df.kPaPer100m*1000 ./ (df.rho_kgm3*g*100)
    q = df.massFlow .* df.margin ./(df.rho_kgm3*3600)
    eps = (df.roughnessMM/1000)
    kinVisc = (df.mu_mPas/1000) ./ df.rho_kgm3
    termA = ((eps.^5)*g .* Sf ./ (q.^2)).^0.25
    termB = ((kinVisc .^ 5)./((q.^3) .* Sf * g )).^0.2
    termC = 0.125*(termA .+ termB).^0.2
    returnValue = (1000*(termC .* (q.^2) ./ (g*Sf)).^0.2)
    return (returnValue)
end

function lineSizeErosion(df)
    # df is the dataframe with hydraulic information
    # this is the erosion C factor method
    # C = v/sqrt(rho), with v in m/s and rho in kg/m3
    # for most service, C = 120, this is similar to C = 100 for imperial units
    # return the ID in mm

    g = 9.80665
    maxVeloc = df.frictionCsi ./ sqrt.(df.rho_kgm3)
    q = df.massFlow .* df.margin ./(df.rho_kgm3*3600)
    area = q ./ maxVeloc
    returnValue = (sqrt.(4*area ./ pi )*1000.0)
    return (returnValue)
end

function getLineSize(df,fluidList)
    # from the calculated line size for the two different methods, determine the required line size
    
    addFluidProperties(df,fluidList)
    
    theSchedule = df[:,:Schedule]

    theSegment = df[:,:Segment]

    tempDP100 = lineSizeDP(df)
    tempErosion = lineSizeErosion(df)
    tempNeeded = tempErosion[:] + tempDP100[:]
    for j in 1:size(tempDP100)[1]
        tempNeeded[j] =max(tempDP100[j], tempErosion[j])
    end
    returnValue = DataFrame(Segment = theSegment, mmDP100 = tempDP100, mmErosion = tempErosion, mmNeeded = tempNeeded, Schedule = theSchedule)
    return (returnValue)
end

# I need a function to pick the next larger line size given the pipe schedule


function getLargerNPS(ourIDmm, ourSchedule)
    # from the common list of line sizes, pick the pipe size that matches the required schedule
    # and is the next larger that the required ID
    # this works for a single line size
    ourScheduleList = IDmm[IDmm[:,:Schedule] .== ourSchedule, :];
    refinedList = ourScheduleList[(ourScheduleList[:,:Idmm] .- ourIDmm) .> 0.0,:];
    theMinID = minimum(refinedList[:,:Idmm])
    ourNPS = refinedList[(refinedList[:,:Idmm] .== theMinID), :NPS]

    return (ourNPS)
end


function getListLargerNPS(lineSizeDF)
    # we will use the required line size and schedule and return the NPS and actual ID
    diamList = DataFrame(Segment=String[], NPS=Float64[], Schedule=String[], IDmm=Float64[])
    for i = 1:(size(lineSizeDF)[1])
        ourSchedule = lineSizeDF[i,:Schedule];
        ourIDmm = lineSizeDF[i,:mmNeeded]
        ourScheduleList = IDmm[IDmm[:,:Schedule] .== ourSchedule, :];
        refinedList = ourScheduleList[(ourScheduleList[:,:Idmm] .- ourIDmm) .> 0.0,:];
        theMinID = minimum(refinedList[:,:Idmm])
        ourNPS = refinedList[(refinedList[:,:Idmm] .== theMinID), :NPS][1]

        push!(diamList, (lineSizeDF[i,:Segment], ourNPS, ourSchedule, theMinID))
    end
    return (diamList)
end

function sizing2hydraulics(sizingDF, chosenLineSizeDF)
    # given the line sizing DF, copy the info into hydraulics dataframe
    # chosenSizeDF contains the NPS and Schedule
    # return the hydraulics DF
    # we also need a simple fittingDF, but this is done with a separate function
    hydraulicsDF = DataFrame(Segment=String[], Description=String[], LineTag=String[], PnID=String[], NPS=Float64[], Schedule=String[],
		Material=String[],	fluidName=[],	inletP_kPaa=[],	massFlow=Float64[])
    for i = 1:(size(sizingDF)[1])
        push!(hydraulicsDF, (sizingDF[i,:Segment], sizingDF[i,:Description],
                sizingDF[i,:LineTag], sizingDF[i,:PnID], chosenLineSizeDF[i,:NPS],
                chosenLineSizeDF[i,:Schedule], sizingDF[i,:Material], sizingDF[i,:fluidName],
                100.0, (sizingDF[i,:massFlow]*sizingDF[i,:margin])))
    end
    return(hydraulicsDF)
end

function sizing2fittingList(sizingDF, chosenLineSizeDF)
    # given the line sizing DF, copy the info into the fitting list
    # chosenSizeDF contains the NPS and Schedule
    # return the pre-populsted fittign list DF
    fittingListDF = DataFrame(Segment=String[], fittingType=String[], num_length_m=Float64[], comment=String[], revision=String[])
    
    for i = 1:(size(sizingDF)[1])
        push!(fittingListDF, (sizingDF[i,:Segment], "PIPE", 100.0, "from line sizing", "A"))
    end
    return(fittingListDF)
end




function incrementLineSize(rowNum, lineSizeDF, increment)
    # increment the specified row number needed line size by +1 or -1
    # not the cleanest function, but it works
    
    neededDF = DataFrame(Segment=String[], mmNeeded=Float64[], Schedule=String[])
    
    i = rowNum
    push!(neededDF, (lineSizeDF[i,:Segment], lineSizeDF[i,:IDmm] + increment, lineSizeDF[i,:Schedule]))

    incrementedLines = getListLargerNPS(neededDF)

    lineSizeDF[i,:NPS] = incrementedLines[1,:NPS]
    lineSizeDF[i,:IDmm] = incrementedLines[1,:IDmm]

    return (3.14)
end

function getReynolds(df)
    # preliminary calculations for pipe hydraulics
    # return the velocity, Reynolds, friction factor, density, diameterMM, diameterInch

    theSegment = df.Segment
    diamMM = df.IDmm
    diamInch = diamMM / 25.4
    rho = df.rho_kgm3
    volFlow = df.massFlow ./ df.rho_kgm3
    area = 0.25 * pi * (df.IDmm / 1000.0).^2
    velocity = volFlow ./ area / 3600
    Reynolds = df.rho_kgm3 .* velocity .* (df.IDmm / 1000.0) ./ (df.mu_mPas/1000.0)
    eD = df.roughnessMM ./ df.IDmm
    moodyF = calcMoodyF(Reynolds, eD)
    returnValue = DataFrame(Segment = theSegment, velocity_ms = velocity, rho_kgm3 = rho, Re = Reynolds, frictF = moodyF, IDmm = diamMM, IDinch = diamInch)
    return (returnValue)
end


function elementDP(lines, fittingList)
    # given the hydraulic elements in df (long list of pipe and fittings)
    # and the preliminary values in prelim
    # calculate the DP in each element
    # dp = rho.f. (L/D) v2/2
    # dp = rho K v2/2
    # where I need to calculate K = K1/Re + Kinf*(1 + Kd/Dinch^0.3)
    # can I get all of the pipe segments and use K = Kp * f L/D, where Kp = 1 for pipe
    # return the value for Kp (1 for pipe, 0 for fitting) and the pressure drop in the item
    returnValue = DataFrame(Segment = String[], Kp = Float64[], fittingK=Float64[], pipeK=Float64[], elementK = Float64[], DPkpa = Float64[])
    
    # start with getting Reynolds and other important things
    prelim = getReynolds(lines)

    for i = 1:(size(fittingList)[1])
        theSegment = fittingList[i,:Segment]

        # first we get all of the K values
        valK1 = fitting3K[fitting3K[:,:fittingType] .== fittingList[i,:fittingType],:K1][1]
        valKinf = fitting3K[fitting3K[:,:fittingType] .== fittingList[i,:fittingType],:Kinf][1]
        valKd = fitting3K[fitting3K[:,:fittingType] .== fittingList[i,:fittingType],:Kd][1]
        valKp = fitting3K[fitting3K[:,:fittingType] .== fittingList[i,:fittingType],:Kp][1]

        # now we get the hydraulic properties
        ff = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:frictF][1]
        rho = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:rho_kgm3][1]
        idInch = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:IDinch][1]
        idMM = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:IDmm][1]
        veloc = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:velocity_ms][1]
        Re = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:Re][1]
        
        Kfitting = valK1/Re + valKinf*(1.0 + valKd/(idInch^0.3))

        # pressure drop in kPa
        thisKfitting = fittingList[i,:num_length_m] * Kfitting
        thisKpipe = valKp * ff * fittingList[i,:num_length_m] /(idMM/1000.0) 
        thisK = thisKfitting + thisKpipe
        thisDP = thisK * (0.001 * rho * 0.5 * (veloc)^ 2)

        push!(returnValue, (theSegment, valKp, thisKfitting, thisKpipe, thisK, thisDP))
    end
    return (returnValue)
end

function compileDP(lines, fittingDP)
    # for each entry in lines, sum all of the fittingDP
    returnValue = DataFrame(Segment = String[], segmentK = Float64[], DPkpa = Float64[], inletP = Float64[], outletP = Float64[])

    for i = 1:(size(lines)[1])
        theSegment = lines[i,:Segment]
        theSumK = sum(fittingDP[fittingDP[:,:Segment] .== theSegment,:elementK])
        theSumDP = sum(fittingDP[fittingDP[:,:Segment] .== theSegment,:DPkpa])
        theInletP = lines[i,:inletP_kPaa]
        theOutletP = theInletP - theSumDP
        push!(returnValue, (theSegment, theSumK, theSumDP, theInletP, theOutletP))
    end
    return (returnValue)
end

function getSegmentDP(lines, fittingList, fluidList)
    # lines describes the hydraulics, fittingList has the piping details, fluidList is the fluid

    addFluidProperties(lines,fluidList)
    addPipeProperties(lines)
    fittingDP = elementDP(lines, fittingList)
    returnValue = compileDP(lines, fittingDP)
    return (returnValue)
end


# functions added for network hydraulics
# february, 2021
# exported
function indexLookup(findMe, inThis)
    # find the index where findMe can be found inThis, which is a vector
    maxI = size(inThis)[1]
    foundIdx = 0
    for i = 1:maxI
        if (findMe == inThis[i]) 
            foundIdx = i
        end
    end
    return(foundIdx)
end

function getConnectivity(connectivity, nodeList)
    # connectivity, lineList and nodeList are dataframe
#    maxLine = size(lineList)[1]

    maxLine = size(connectivity)[1]

    maxNode = size(nodeList)[1]
    C = zeros(Float64,maxLine,maxNode)
    
    maxConnect = size(connectivity)[1]
    
    for i = 1:maxConnect
        theSegment = connectivity[i,:Segment]
        inNode = connectivity[i,:inNode]
        outNode = connectivity[i,:outNode]
#        theRow = indexLookup(theSegment,lineList[:,:Segment])
        theRow = indexLookup(theSegment,connectivity[:,:Segment])
        inCol = indexLookup(inNode,nodeList[:,:Node])
        outCol = indexLookup(outNode,nodeList[:,:Node])
        C[theRow,inCol] = 1.0
        C[theRow,outCol] = -1.0
    end
    return(C)
end

function removeTrivialRows(A)
    # remove rows from matrix A that contain only one entry
    numRows = size(A)[1]
    numColumns = size(A)[2]
    indexList = zeros(Int64, numRows,1)
    numNonTrivialRows = 0
    for i = 1:numRows
        countNonZero = 0
        for j = 1:numColumns
            if (A[i,j] != 0.0)
                countNonZero = countNonZero + 1
            end
        end
        if (countNonZero >= 2)
            numNonTrivialRows += 1
            indexList[numNonTrivialRows] = i
        end
                
    end
    returnMatrix = zeros(numNonTrivialRows,numColumns)
    for i = 1:numNonTrivialRows
        for j = 1:numColumns
            returnMatrix[i,j] = A[indexList[i],j]
        end
    end
    
    return(returnMatrix)
end

function initializeHydraulicParameters(setRe, df)
    # set all of the Reynolds numbers in lines df
    # return the hydraulicParameters (Reynolds, friction factor, density, diameterMM, diameterInch)

    theSegment = df.Segment
    diamMM = df.IDmm
    rho_kgm3 = df.rho_kgm3
    diamInch = diamMM / 25.4

    volFlow = df.massFlow ./ df.rho_kgm3
    area = 0.25 * pi * (df.IDmm / 1000.0).^2
    velocity = volFlow ./ area / 3600
    Reynolds = setRe .* (df.rho_kgm3 ./ df.rho_kgm3)
    eD = df.roughnessMM ./ df.IDmm
    moodyF = calcMoodyF(Reynolds, eD)

    returnValue = DataFrame(Segment = theSegment, rho_kgm3 = rho_kgm3, Re = Reynolds, frictF = moodyF, IDmm = diamMM, IDinch = diamInch)
    return (returnValue)
end


function initializeMassFlow(df, connectivity)
    # set mass flow in df to 1000 kg/h
    # but only if the flow rate is undefined or zero
    # I want to create a df with fields Segment and massFlow
    # we have viscosity, density and IDmm in the df
    # we will use a frictional pressure gradient of 20 kPa per 100 m for initial guess on mass flow
    # and friction factor 0.02
    # but then overwrite if we have something filled in for massFlow in the df

    returnValue = DataFrame(Segment = String[], massFlow = Float64[])

    # first, set all mass flows equal to 1000 kg/h
    numRows = size(connectivity)[1]
    for i = 1:numRows
        push!(returnValue, (connectivity[i,:Segment], 1000.0) )
    end


    # now, plow thorugh the hydraulic information
    # if there is no information here, then use DP per 100m to calculate a mass flow
    # if there is information here, use this
    numHydraulicRows = size(df)[1]

    dp100 = 20.0
    ff = 0.02 # estimate for friction factor
    
    for i = 1:numHydraulicRows
        ii = indexLookup(df[i,:Segment],connectivity[:,:Segment])
        # calculate this just in case
        dd = df[i,:IDmm]
        rho = df[i,:rho_kgm3]
        visc = df[i,:mu_mPas]
        mm = 3600 * sqrt(20*1000*rho*(pi^2 / (8*ff))*((dd^5)/100)*(1/1000)^5)
        if (ismissing(df[i,:massFlow]) || (df[i,:massFlow] == 0.0) ) 
            df[i,:massFlow] = mm
            returnValue[ii,:massFlow] = mm
        else
            returnValue[ii,:massFlow] = df[i,:massFlow]
        end
    end
    return (returnValue)
end

function getResistanceMatrix(hydraulicParameters, elementK, connectivity)
    # given the hydraulic parameters, a list of all segments
    # and the long list of K values for every element in the network
    # calculate the sum of the overall K for each segment
    # return a diagonal matrix with the K values on the diagonal
    # wait, the resistance K value is in terms of velocity in m/s
    # I need mass flow in kg/h
    #
    # I think I need to add the connectivity matrix to get the correct value for the variable number
    
    numRows = size(hydraulicParameters)[1]
    numColumns = size(connectivity)[1]

    numElements = size(elementK)[1]
#    matrixK = zeros(Float64,numRows,numRows)
#     this needs to be the number of variables in the connectivity matrix
    matrixK = zeros(Float64,numRows,numColumns)

    
    for i = 1:numElements
#        theRow = indexLookup(elementK[i,:Segment], hydraulicParameters[:,:Segment])
#        theRow = indexLookup(elementK[i,:Segment], hydraulicParameters[:,:Segment])

        theRow = indexLookup(elementK[i,:Segment], hydraulicParameters[:,:Segment])
        theColumn = indexLookup(elementK[i,:Segment], connectivity[:,:Segment])

        matrixK[theRow,theColumn] += elementK[i,:elementK]
    end
    # now we multiply each row element by 8 x 1000 / (rho pi^2 D^4) (1000/3600)^2 
    # where D is in mm, mass flow in kg/h and DP in kPa
    for i = 1:numRows
        multK = 8.0*1000.0 * ((1000/3600)^2) / ((pi^2) * hydraulicParameters[i,:rho_kgm3] * (hydraulicParameters[i,:IDmm]^4)   )
        for j = 1:numColumns
            matrixK[i,j] *= multK
        end
    end
    return (matrixK)
end

function getSegLength(connectivity, fittingList)
    
    returnValue = DataFrame(Segment = String[], pipeLength=Float64[])
    
# initialize the segLength data
    for i=1:(size(connectivity)[1])
	push!(returnValue, (connectivity[i,:Segment], 0.0))
    end
    for i = 1:(size(fittingList)[1])
        valKp = fitting3K[fitting3K[:,:fittingType] .== fittingList[i,:fittingType],:Kp][1]
	ii = indexLookup(fittingList[i,:Segment], connectivity[:,:Segment])
	if (valKp == 1) 
		returnValue[ii,:pipeLength] += fittingList[i,:num_length_m]
	end
    end
    return(returnValue)
end




function getConstraintMatrix(constraintDF, constraintType, segmentDF, nodeDF)
    # given the constraint data, the type (massFlow or pressure or value), the list of segments and nodes
    # create a matrix with values of 1 in the row, or the actual value
    
    numConstraints = size(constraintDF)[1]
    numSegments = size(segmentDF)[1]
    numNodes = size(nodeDF)[1]
    
    # we need the dimensions for the return matrix
    # assume a value for now
    numCols = 2
    if (constraintType == "flowRate")
        numCols = numSegments
    end
    if (constraintType == "pressure")
        numCols = numNodes
    end
    if (constraintType == "value")
        numCols = 1
    end
    returnMatrix = zeros(numConstraints, numCols)
    
    # cycle through all of the constraints
    # of the constraintType matches the :variable in the constraintDF, then place a 1 in the matrix entry
    for i = 1:numConstraints
        if (constraintType == "value")
            returnMatrix[i] = constraintDF[i,:value]
        elseif (constraintType == "flowRate")
            varIndexLocal = indexLookup(constraintDF[i,:object], segmentDF[:,:Segment])
            if (varIndexLocal > 0)
                returnMatrix[i,varIndexLocal] = 1.0
            end
        elseif (constraintType == "pressure")
            varIndexLocal = indexLookup(constraintDF[i,:object], nodeDF[:,:Node])
            if (varIndexLocal > 0)
                returnMatrix[i,varIndexLocal] = 1.0
            end
        end
    end

    return(returnMatrix)
end

function getQuadraticMatrix(hydraulicDF, which)
    # given the estimates of mass flow rate in the hydraulic matrix, fill out the
    # linearized estimate to the quadratic mass flow^2
    # m^2 ~= 2 M m - M*M,
    # or
    # m2 - 2 M m = - M*M, where we are looking for m2 (m^2) and m, and M is the
    # previous estiamte of m
    # the value of which is either quadratic, linear, value

    # given the constraint data, the type (massFlow or pressure or value), the list of segments and nodes
    # create a matrix with values of 1 in the row, or the actual value
    
    numEquations = size(hydraulicDF)[1]
    numCols = numEquations
    
    if (which == "rhs")
        numCols = 1
    end

    returnMatrix = zeros(numEquations, numCols)
    
    # cycle through all of the segments
    # fill in the value for the matrix equation
    # m2 - 2 M m = - M*M, where we are looking for m2 (m^2) and m, and M is the
    # if quadratic rM[i,i] = 1.0
    # if linear rm[i,i] = -2*hydraulicDF[i,:massFlow]
    # if value rm[i] = -hydraulicDF[i,:massFlow]*hydraulicDF[i,:massFlow]

    for i = 1:numEquations
        if (which == "rhs")
            returnMatrix[i] = -hydraulicDF[i,:massFlow]*hydraulicDF[i,:massFlow]
        elseif (which == "quadratic")
            returnMatrix[i,i] = 1.0
        elseif (which == "linear")
            returnMatrix[i,i] = -2.0 * hydraulicDF[i,:massFlow]
        end
    end

    return(returnMatrix)
end

function getSegmentIMatrix(hydraulicParameters, connectivity)
    # given the Identity matrix for the segments
    # return Identity matrix
    
    numRows = size(hydraulicParameters)[1]
    numColumns = size(connectivity)[1]

    matrixI = zeros(Float64,numRows,numColumns)
    
    for i = 1:numRows
        findMe = hydraulicParameters[i,:Segment]
        inThis = connectivity[:,:Segment]
        j = indexLookup(findMe, inThis)
        matrixI[i,j] = 1.0
    end

    return (matrixI)
end

function getElementK(prelim, fittingList)
    # given the hydraulic elements in df (long list of pipe and fittings)
    # and the preliminary friction properties (Re, eD, ff)
    # calculate the loss coefficient K in each element
    # where I need to calculate K = K1/Re + Kinf*(1 + Kd/Dinch^0.3)
    # can I get all of the pipe segments and use K = Kp * f L/D, where Kp = 1 for pipe
    # return the value for Kp (1 for pipe, 0 for fitting) and the pressure drop in the item
    returnValue = DataFrame(Segment = String[], Kp = Float64[], fittingK=Float64[], pipeK=Float64[], elementK = Float64[], pipeLength=Float64[])
    

    for i = 1:(size(fittingList)[1])
        theSegment = fittingList[i,:Segment]

        # first we get all of the K values
        valK1 = Hydraulics.fitting3K[Hydraulics.fitting3K[:,:fittingType] .== fittingList[i,:fittingType],:K1][1]
        valKinf = Hydraulics.fitting3K[Hydraulics.fitting3K[:,:fittingType] .== fittingList[i,:fittingType],:Kinf][1]
        valKd = Hydraulics.fitting3K[Hydraulics.fitting3K[:,:fittingType] .== fittingList[i,:fittingType],:Kd][1]
        valKp = Hydraulics.fitting3K[Hydraulics.fitting3K[:,:fittingType] .== fittingList[i,:fittingType],:Kp][1]

        # now we get the hydraulic properties
        ff = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:frictF][1]
        rho = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:rho_kgm3][1]
        idInch = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:IDinch][1]
        idMM = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:IDmm][1]
#        veloc = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:velocity_ms][1]
        Re = prelim[prelim[:,:Segment] .== fittingList[i,:Segment],:Re][1]
        
        Kfitting = valK1/Re + valKinf*(1.0 + valKd/(idInch^0.3))

        # pressure drop in kPa
        thisKfitting = fittingList[i,:num_length_m] * Kfitting
        thisKpipe = valKp * ff * fittingList[i,:num_length_m] /(idMM/1000.0) 
	if (valKp == 1) 
		thisPipeLength = fittingList[i,:num_length_m]
	else
		thisPipeLength = 0.0
	end

        thisK = thisKfitting + thisKpipe
#        thisDP = thisK * (0.001 * rho * 0.5 * (veloc)^ 2)

        push!(returnValue, (theSegment, valKp, thisKfitting, thisKpipe, thisK, thisPipeLength))
    end
    return (returnValue)
end

function addToBlockMatrix(gBM, localMatrix, gR, gC, eqnType, varType, numLocalRows, numLocalColumns)
    # we will traverse the localMatrix and place the entries in the global block matrix
    # gR and gC are the table of indices that map gR[eqnType,i] to the row in global matrix
    # and gC[varType,j] to the column in the global matrix
    # eqnType and varType refer to which block
    # numLocalRows and Columns tells us how many local entries are avalable

    for i = 1:numLocalRows
        for j = 1:numLocalColumns
            gi = gR[eqnType,i]
            gj = gC[varType,j]
            gBM[gi,gj] = localMatrix[i,j]
        end
    end
    return(3.14)
end



function addToBlockRHS(gBv, localVec, gR, eqnType, numLocalRows)
    # we will traverse the localMatrix and place the entries in the global block matrix
    # gR and gC are the table of indices that map gR[eqnType,i] to the row in global matrix
    # and gC[varType,j] to the column in the global matrix
    # eqnType and varType refer to which block
    # numLocalRows and Columns tells us how many local entries are avalable

    for i = 1:numLocalRows
        gi = gR[eqnType,i]
        gBv[gi] = localVec[i]
    end
    return(3.14)
end

# create routines to extract results
# extractFromBlockSolution(testSolution, globalColumn, j (which type is it),columnsInSegment[j])
function extractFromBlockSoln(gBx, gR, varType, numLocalRows)
    # we will traverse the appropriate entries in the global block matrix
    # gR is the vector of indices that map gR[varType,i] to the row in global solution vector
    # varType refer to which block
    # numLocalRows tells us how many local entries are avalable

    localX=zeros(numLocalRows,1)
    for i = 1:numLocalRows
        gi = gR[varType,i]
        localX[i] = gBx[gi]
    end
    return(localX)
end

# create routines to extract results
# extractFromBlockSolution(testSolution, globalColumn, j (which type is it),columnsInSegment[j])
function catIt(dfHydraulic1, dfHydraulic2, dfAll)
    # dfAll has all of the rows, we need to fill in blank rows for the other dfs

    localX=zeros(numLocalRows,1)
    for i = 1:numLocalRows
        gi = gR[varType,i]
        localX[i] = gBx[gi]
    end
    return(localX)
end

# create routines to extract results
function extractResults(gBx, lines, connectivity, nodeList, fittingList)
    # we will traverse the solutions in the block x vector
    # lines is the df with the segments defined
    # connectivity is the connectivity df
    # we will find all of the from entries in the connectivity matrix
    numLines = getSize(connectivity)
    numHydraulicLines = getSize(connectivity)    
    numNodes = getSize(nodeList)

    columnsInSegment = getColumnsInSegment(numLines,numNodes)
    globalColumn = variableMapping(numLines, numNodes, columnsInSegment)

    ourDetails = (Hydraulics.getReynolds(lines))
	# note that fittingList is a really long vector with every fitting
    ourElementK = getElementK(ourDetails, fittingList)
    Km = getResistanceMatrix(ourDetails, ourElementK, connectivity)  # this is an actual matrix
    segLength = getSegLength(connectivity,fittingList)  # dataframe for length and other stuff

    segList = copy(connectivity[:,:Segment])
    nList = copy(nodeList[:,:Node]) 

    # connectivity, used for pressure continuity
    localA = getConnectivity(connectivity,nodeList)
    
    # extract all of the entries with +1 and -1
    fromA = max.(localA, 0)
    toA = -min.(localA, 0)

    varType = 1        # mass flow
    xMassFlow = extractFromBlockSoln(gBx, globalColumn, varType, columnsInSegment[varType])

    varType = 2     # mass flow^2
    xMassFlow2 = extractFromBlockSoln(gBx,globalColumn,varType,columnsInSegment[varType])


    varType = 3     # DP
    xDP = extractFromBlockSoln(gBx,globalColumn,varType,columnsInSegment[varType])

    varType = 4     # node pressure
    xPnode = extractFromBlockSoln(gBx,globalColumn,varType,columnsInSegment[varType])
    fromP = fromA*xPnode
    toP = toA*xPnode


    # messy, but it works
#    allResults = [segList xMassFlow xDP fromP toP Konly]
    allResults = [segList xMassFlow xDP fromP toP]


#    returnVal = DataFrame(Segment=allResults[:,1], massFlow=allResults[:,2], dp=allResults[:,3], fromP=allResults[:,4], toP=allResults[:,5], elementK=allResults[:,6])

    returnVal = DataFrame(Segment=allResults[:,1], massFlow=allResults[:,2], dp=allResults[:,3], fromP=allResults[:,4], toP=allResults[:,5], fromNode=connectivity[:,:inNode], toNode=connectivity[:,:outNode], segmentLength=segLength[:,:pipeLength])


    return(returnVal)
end

function getSize(df)
    # get the number of rows
    return(size(df)[1])
end

function getColumnsInSegment(numLines, numNodes)
    # define the number of variables in each of the different types of variables
    numEqnTypes = 5 # not used
    numVarTypes = 4 # m, m2, dp, p

    columnsInSegment = zeros(Int64,numVarTypes,1)
    columnsInSegment[1] = numLines
    columnsInSegment[2] = numLines
    columnsInSegment[3] = numLines
    columnsInSegment[4] = numNodes
    
    return(columnsInSegment)
end

function variableMapping(numLines, numNodes, columnsInSegment)
    # create the variable mapping as a 2 D array
    # map[row,column] = global index
    # row is the type of variable (m, m2, dp, p)
    # column contains the index for each of the available variables in the global vector
    # function to count number of equations and variables
    numEqnTypes = 5 # not used
    numVarTypes = 4 # m, m2, dp, p

    maxOtherIndex = numLines + numNodes
    globalColumn = zeros(Int64,numVarTypes,maxOtherIndex)

    index = 1
    for i = 1:size(columnsInSegment)[1]
        for j = 1:columnsInSegment[i]
            globalColumn[i,j] = index
            index += 1
        end
    end
    
    return(globalColumn)
end

function doNetworkHydraulics(lineHydraulics, fluidList, connectivity, dofList, nodeList, fittingList)
    # there is a lot of code here. I do not want the user typing this. It should only be function calls
    #
    # check input data
        # numLines = checkSegments(lineHydraulics)
        # numNodes = checkNodes
        # I have a DOF check in the code
        # or numLines = getSize(lineHydraulics)
        # or numNodes = getSize(nodeList)
        # variableMapping = getVariableMapping(numLines, numNodes)
    #
    # resultVector = doHydraulics(lineHydraulics, fluidList,dofList,nodeList)
    # Then we need to create a new table with the segment results
    # and we need another table with the pressures
    #
    # then do this
    # segmentResults = extractResults("segment", lineHydraulics, resultVector)
    # nodeResults = extractResults("node", lineHydraulics, resultVector)
    # I could supply a list of dataframe names to include in the results
    #

    # driver script
    addFluidProperties(lineHydraulics,fluidList)
    addPipeProperties(lineHydraulics)

    # and this is where we need to initialize
    # return the matrix of mass flow and m2
    # I will need to provide the connectivity df as well to get all of the fields
    massFlow = initializeMassFlow(lineHydraulics,connectivity) # if we do not have a mass flow, make one up

    # pass the mass flow matrix
    #and set the Reynolds number and friction factor to a reasonable value
    hydraulicParameters = initializeHydraulicParameters(40000.0, lineHydraulics)


    numLines = size(connectivity)[1]
    numResistance = size(lineHydraulics)[1]
    numNodes = size(nodeList)[1];
    #C = zeros(Float64,numLines,numNodes)

    segList = copy(connectivity[:,:Segment])
#    segList = copy(lineHydraulics[:,:Segment])
    nList = copy(nodeList[:,:Node]) 

    # connectivity, used for pressure continuity
    # each row defines a segment: we have the inlet and outlet node
    # I don't think I need the lineHydraulics df
    Amatrix = getConnectivity(connectivity,nodeList)
    # I think I only need to assemble the block matrix.
    Ia = zeros(Float64,numLines,numLines)
    for i = 1:numLines
        Ia[i,i] = 1.0; 
    end;


    # create the mass balance matrix
    Bprelim = Amatrix'
    Bmatrix = removeTrivialRows(Bprelim);

    # I need to create the constraint equations
    # use the connectivity df instead of lineHydraulics
    Cm = getConstraintMatrix(dofList, "flowRate", connectivity, nodeList)
    Cp = getConstraintMatrix(dofList, "pressure", connectivity, nodeList)
    xc = getConstraintMatrix(dofList, "value", connectivity, nodeList);
    
    # need to define where each matrix will go.
    # sure I could move this outside the loop, but it should not slow things down much
    # need vectors for the number of rows and columns in each segment
    # I could put most of this batch of code in a function
    # let 
    # l = num lines
    # n = num nodes
    # k = resistance calcs
    # d = degrees of freedom

    # num unknonws = 3*l + n
    # count equations
    # num of mass balance relations = num rows in B
    # quadratic equations = l
    # resistance eqns = k
    # pressure connectivity = l
    # degree of freedom = d
    
        # the block matrix is
        #  1.  2.  3.  4.  
        #    m, m2, dp,  p | rhs
        # 1  B,  0,  0,  0 | 0
        # 2 Qm,Qm2,  0,  0 | xq
        # 3  0,-Km,Iseg, 0 | 0
        # 4  0,  0, -I, Am | 0
        # 5 Cm,  0,  0, Cp | xc



    # function to count number of equations and variables
    # edit this to account fo rnumber of resistances
    numEqnTypes = 5
    numVarTypes = 4
    rowsInSegment = zeros(Int64,numEqnTypes,1)
    rowsInSegment[1] = size(Bmatrix)[1]
    rowsInSegment[2] = numLines
    rowsInSegment[3] = numResistance
    rowsInSegment[4] = numLines
    rowsInSegment[5] = size(dofList)[1]

    columnsInSegment = zeros(Int64,numVarTypes,1)
    columnsInSegment[1] = numLines
    columnsInSegment[2] = numLines
    columnsInSegment[3] = numLines
    columnsInSegment[4] = numNodes

    numEquations = sum(rowsInSegment)
    numVariables = sum(columnsInSegment)

    # I need this screening as a preliminary check
    (numEquations != numVariables) ? error("error: ", numEquations, " equations and ", numVariables, " variables") : println("OK: ", numEquations, " equations and variables")

    # this dozen lines could go into a function
    # I could create a matrix with global indices
    # globalRow[segment][localRow]
    # globalColumn[segment][localColumn]
    maxOtherIndex = numLines + numNodes
    globalRow = zeros(Int64,numEqnTypes,maxOtherIndex)
    globalColumn = zeros(Int64,numVarTypes,maxOtherIndex)

    index = 1
    for i = 1:size(rowsInSegment)[1]
        for j = 1:rowsInSegment[i]
            globalRow[i,j] = index
            index += 1
        end
    end

    index = 1
    for i = 1:size(columnsInSegment)[1]
        for j = 1:columnsInSegment[i]
            globalColumn[i,j] = index
            index += 1
        end
    end

    # define the global block matrix and vector

    nGlobalRow = maximum(globalRow)
    nGlobalColumn = maximum(globalColumn)
    gBM = zeros(Float64,nGlobalRow,nGlobalColumn) # global block matrix
    gBv = zeros(Float64,nGlobalRow,1) # global block vector
    gBx = zeros(Float64,nGlobalRow,1) # global block solution

    epsM = 0.01 # mass flow kg/h
    maxIter = 15

    errM = 100.0
    iter = 0

    # and this is where we need to iterate

    while ((errM > epsM) && (iter < maxIter) )
        # the resistance matrix is not applied to segments that are control valves
        ourElementK = getElementK(hydraulicParameters, fittingList)
        Km = getResistanceMatrix(hydraulicParameters, ourElementK, connectivity)
        Iseg = getSegmentIMatrix(hydraulicParameters, connectivity)

        # we need the quadratic matrix.
        # must realize that our lineHydraulic table is not the definative source of infomration
        # for this
        # we should use the connectivity table
        # Or do we create a simple table for only mass flow and m2?
        # m2 = 2 m M - M M
        # m2 - 2 m M = -M M
        # change this from lineHydraulics to massFlow
        # why does massFlow not work?
#        Qm2 = getQuadraticMatrix(lineHydraulics, "quadratic")
#        Qm  = getQuadraticMatrix(lineHydraulics, "linear")
#        xq  = getQuadraticMatrix(lineHydraulics, "rhs");

        Qm2 = getQuadraticMatrix(massFlow, "quadratic")
        Qm  = getQuadraticMatrix(massFlow, "linear")
        xq  = getQuadraticMatrix(massFlow, "rhs");



        # the block matrix is
        #  1.  2.  3.  4.  
        #    m, m2, dp,  p | rhs
        # 1  B,  0,  0,  0 | 0
        # 2 Qm,Qm2,  0,  0 | xq
        # 3  0,-Km,Iseg, 0 | 0
        # 4  0,  0, -I, Am | 0
        # 5 Cm,  0,  0, Cp | xc
        #
        # P1 - P2 = K.m2, where P1 is in, P2 is out
        # DP = K.m2
        # P1 - P2 = DP
        # or
        # -K.m2 + DP = 0
        # -DP + P1 - P2 = 0

        # once I have the answer, how do I get inlet and outlet pressure for each segment.
            # The A matrix contains the information
            # The rows in A are for each segment.
            # columns with +1 are inlet, -1 are outlet.
            # Transform the A matrix by extracting the entries with +1. Then matrix multiply.
            # then get the entries with -1 and matrix multiply.



        # assemble the global block matrix
        i = 1
        j = 1
        addToBlockMatrix(gBM,Bmatrix,globalRow,globalColumn,i,j,rowsInSegment[i],columnsInSegment[j])

        i = 2
        j = 1
        addToBlockMatrix(gBM,Qm,globalRow,globalColumn,i,j,rowsInSegment[i],columnsInSegment[j])

        i = 2
        j = 2
        addToBlockMatrix(gBM,Qm2,globalRow,globalColumn,i,j,rowsInSegment[i],columnsInSegment[j])
        addToBlockRHS(gBv, xq, globalRow, i, rowsInSegment[i])


        i = 3
        j = 2
        addToBlockMatrix(gBM,-Km,globalRow,globalColumn,i,j,rowsInSegment[i],columnsInSegment[j])

        i = 3
        j = 3
        addToBlockMatrix(gBM,Iseg,globalRow,globalColumn,i,j,rowsInSegment[i],columnsInSegment[j])

        i = 4
        j = 3
        addToBlockMatrix(gBM,-Ia,globalRow,globalColumn,i,j,rowsInSegment[i],columnsInSegment[j])

        i = 4
        j = 4
        addToBlockMatrix(gBM,Amatrix,globalRow,globalColumn,i,j,rowsInSegment[i],columnsInSegment[j])

        i = 5
        j = 1
        addToBlockMatrix(gBM,Cm,globalRow,globalColumn,i,j,rowsInSegment[i],columnsInSegment[j])

        i = 5
        j = 4
        addToBlockMatrix(gBM,Cp,globalRow,globalColumn,i,j,rowsInSegment[i],columnsInSegment[j])
        addToBlockRHS(gBv, xc, globalRow, i, rowsInSegment[i])

        # and solve
        testSoln = gBM \ gBv

        # convergnece
#        errM = maximum(abs.(lineHydraulics[:,:massFlow] - testSoln[1:numLines]))

        errM = maximum(abs.(massFlow[:,:massFlow] - testSoln[1:numLines]))

        iter += 1


        # update the mass flow rate solutions
        massFlow[:,:massFlow] = testSoln[1:numLines]
        for i = 1:numResistance
            ii = indexLookup(lineHydraulics[i,:Segment], connectivity[:,:Segment])
            lineHydraulics[i,:massFlow] = testSoln[ii]
        end

        # then I need to update the hydraulic parameters
        hydraulicParameters = Hydraulics.getReynolds(lineHydraulics)
        # and I need to update the pressures and DP estimates

        gBx = testSoln
    end
    return (gBx)
end

# create routines to extract results
# extractFromBlockSolution(testSolution, globalColumn, j (which type is it),columnsInSegment[j])
function combineResults(dfAllResults, dfLines)
    # dfAll has all of the rows
    # dfLines has the pressur drop hydraulic info
    # do an innerjoin to produce a single table

    hydraulicDetails = (Hydraulics.getReynolds(dfLines));
    combined1=outerjoin(dfAllResults, hydraulicDetails, on = intersect(names(dfAllResults), names(hydraulicDetails)), makeunique=false,indicator=nothing,validate=(false,false));
    combined2=outerjoin(combined1, dfLines, on = intersect(names(combined1), names(dfLines)), makeunique=false,indicator=nothing,validate=(false,false));

    return(combined2)
end


end # module


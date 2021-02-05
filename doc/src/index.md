<!-- Index -->

## Hydraulics.jl documentation

Units of measure used in the package are:

* Flow rate, kg/h
* Pressure, kPaa
* Diameter, mm
* Piping length, m

Three sets of example files are provided to demonstrate functionality.

* Line Sizing
	* demoHydraulicsSizingRev1.ipynb
	* sizingHydraulicsSampleRev1.xls
* Simple hydraulics (pressure drop)
	* demoHydraulicsSimpleRev2.ipynb
	* simpleHydraulicsSampleRev2.xls
* Network hydraulics
	* demoHydraulicsNetworkRev1.ipynb
	* networkHydraulicsSampleRev1.xls

## Limitations

* Gravity not yet implemented.



## Line Sizing

Line sizing is based on two parameters:

1. Pressure drop per 100 m.
2. Erosional velocity determined by C factor (SI). v_max = C / sqrt(rho), where C = 120 is typical (equivalent to C = 100 for Imperial units in API 14E).

The calculation method comes from Perry's Handbook, 6th edition, page xxx. This presents the correlation for friction factor (by Chen) where the pressure gradient is used as the explicit variable.

## Simple Hydraulics

Simple hydraulic calculations are standard pressure drop for piping and fittings.

* Halaand correlation for friction factor (explicit in pipe Reynolds number).
* 3K method used for calculating the loss coefficient for fittings.
* A very large list of piping components can be included for each line segment.

## Network Hydraulics

Network hydraulic calculations are used where there is a large network of piping that supplies and removes fluid. These calculations are more complex than simple hydraulic calculations.

The following calculations are performed:
1. Mass balance at every pipe junction.
2. Pressure drop for each pipe segment.
3. Pressure continuity for each segment (inletPressure - outletPressure = Pressure drop).

The full set of nonlinear hydraulic equations (nominally quadrtic) are solved using a Newton iteration. The resulting matrix equation (with nominally 4 unknowns for each pipe segment) are solved with a full matrix solver at each iteration. This means that it is not necessary to specify the mass flow for each piping segment. Specifying the pressure at each end of a piping segment will result in the correct mass flow rate to balance the pressure drop.

## How to use the tools

Performing line sizing or hydraulic (pressure drop) calculations requires standard piping data (tabulated with the Hydraulic module source code):

* pipe material (roughness)
* pipe schedule (for estimating the nominal pipe size)

Each of the three hydraulic calculations has a sample input file and example source code.

## Line sizing

The line sizing application is described with the following examples for calculation methodology and problem definitions (inputs):

* demoHydraulicsSizingRev1.ipynb
* sizingHydraulicsSampleRev1.xls

The spreadsheet contains the following tabs.

* fluidList
	* Defines physical properties (density, viscosity) for the fluids that are used for line sizing
* sizing
	* Defines the piping properties Schedule, material type (for roughness), and the fluid contained.
	* Defines the normal mass flowrate and applied margin.
	* Defines the erosional velocity sizing criteria (C value, in SI units)
	* Defines the pressure gradient (kPa per 100 m) as the sizing criteria
* Other tabs are provided for reference and to provide the list of appropriate options for lookup tables for NPS, Schedule, piping material

You could use a CSV file for the input data, but then you would need to ensure that all of the piping data matches the available choices (dropdown menus in a spreadsheet is more convenient).

A line sizing calculation requires a spreadsheet table (or CSV file) with the following data:

* Segment: Descriptor of the segment. Must be a unique in this table. This could be alphanumeric (A, B, C, or 1, 2, 3). A segment must have the same fluid flow rate, fluid properties, line size.
* Description: optional, text description of the line.
* LineTag: optional, detailed text description of the line. Typically the line number from the PnID.
* PnID: optional, text description of the drawing reference
* fluidName: discrete field (taken from fluid property table) to identify the fluid
* Schedule: discrete field (taken from the pipe schedule table) to identify the pipe wall thickness
* Material: discrete field (taken from the pipe material table) to identify the pipe wall rooughness
* massFlow: nominal flowrate in kg/h. Consider this to be the normal flow rate.
* margin: margin applied to the flowrate for the line sizing calculation.
* fricionCsi: friction C factor, in SI units. C = 120 is considered typical, and is equivalent to C = 100 for Imperial units from API-12. v_max = C_si / sqrt(rho)
* kPaPer100m: design pressure gradient in kPa per 100 m.


### Line sizing functions

* getLineSize(lines,fluidList)
    * given the lines (line sizing) dataframe, return line sizing values in dataframe
    * Segment: segment in df
    * mmDP100: ID based on kPa per 100 m
    * mmErosion: ID based on erosion C value
    * mmNeeded: larger of the two values
    * Schedule: the line schedule to be used
* getLargerNPS(ourIDmm, ourSchedule)
    * return the NPS of the next larger pipe (given the schedule)
* incrementLineSize(rowNum, lineSizeDF, increment)
    * increases (or decreases) the line NPS (increment is 1 or -1)
* sizing2hydraulics(sizingDF, chosenLineSizeDF)
    * copy the line sizing dataframe to a hydraulics dataframes.
    * all line lengths pre-populated at 100 m.
* sizing2fittingList(sizingDF, chosenLineSizeDF)
    * populates a simple fitting list with each line specified in the sizing DF, and makes each line 100 m long

### Sample use

Check for consistancy errors.  
```checkLineList(lineSizing,fluidList)```

Find the minimum required pipe ID.
```myLineSize = getLineSize(lineSizing, fluidList)```

Get the next larger NPS from the line list.
```chosenLineSizes = getListLargerNPS(myLineSize)```

Increment/decrement a line NPS to the next larger/smaller size.
```incrementLineSize(3,chosenLineSizes,1)
chosenLineSizes```

Combine information to a common table.
```
 combined1=outerjoin(lineSizing, chosenLineSizes, on = intersect(names(lineSizing), names(chosenLineSizes)), makeunique=false,indicator=nothing,validate=(false,false))
```
```
 combined2=outerjoin(combined1, myLineSize, on = intersect(names(combined1), names(myLineSize)), makeunique=false,indicator=nothing,validate=(false,false))
```
Get list of names for result headings.  
``` names(combined2)```

Create table of results.
```
 resultTable1 = select(combined2, [:Segment, :fluidName, :massFlow, :margin, :mmDP100, :mmErosion, :NPS, :Schedule])
```


## Simple hydraulic calculations

The simply hydraulic calculation application is described with the following examples for calculation methodology and problem definitions (inputs):

* demoHydraulicsSimpleRev2.ipynb
* simpleHydraulicsSampleRev2.xls

The spreadsheet is more complex than the line sizing application, and contains the following tabs.

* fluidList
	* Defines physical properties (density, viscosity) for the fluids that are used for line sizing.
* hydraulics
	* Defines the piping properties Nominal Pipe Size, Schedule, material type (for roughness) and the fluid contained.
	* Defines the mass flowrate for each segment.
* fittingList
	* A large list to define the number of fittings and pipe length for each pipe segment.
* Other tabs are provided for reference and to provide the list of appropriate options for lookup tables for NPS, Schedule, piping material, fitting types.

CSV files could also be used to describe the inputs.

The hydraulics input sheet requires the following information:

* Segment: Descriptor of the segment, typically a PFD stream number or a PnID line number. Must be unique in this table. This could also be alphanumeric (A, B, C, or 1, 2, 3). A segment must have the same fluid flow rate, fluid properties, line size.
* Description: text descrption of the line.
* LineTag: detailed text description of the line. Typically the line number from the PnID.
* NPS: nominal pipe diameter (taken from pipe NPS table) to identify the pipe diameter.
* Schedule: discrete field (taken from the pipe schedule table) to identify the pipe wall thickness.
* Material: discrete field (taken from the pipe material table) to identify the pipe wall rooughness
* PnID: optional, text description of the drawing reference.
* fluidName: discrete field (taken from fluid property table) to identify the fluid.
* inletP_kPaa: pressure at pipe inlet
* massFlow: actual flowrate in kg/h. 

Pipe lengths and fitting counts are tabulated separately for each line segment:

* Segment: Descriptor of the segment, typically a PFD stream number or a PnID line number. Must be a unique in this table. This could also be alphanumeric (A, B, C, or 1, 2, 3). A segment must have the same fluid flow rate, fluid properties, line size.
* fittingType: discrete field (taken from the fitting3K list) to identify the type of fitting. Example, PIPE or EL45-THD-STD
* num_length_m: either the total number of fittings of the given type, or the total length of pipe specified.
* comment: optional comment
* revision: optional revision comment

Defining the pipe lengths and fitting counts in a large table permits easy checking of the input data. The numbers can come directly from piping isometrics. For example, we may specify a number of PIPE entries for a line segment, where each PIPE length matches a value on an isometric.


### Hydraulic functions 

* getSegmentDP(lines, fittingList)
	* This is the main function for hydraulic calculations
    * lines describes the hydraulics, fittingList has the piping details, fluidList is the fluid
    * this calculates the pressure drop through each item listed in the fittingList
    * then the pressure drop for each segment is compiled
    * function returns a vector of pressure drops
* getReynolds(df)
    * given the hydraulic dataframe df, return dataframe with
        * velocity, Reynolds, friction factor, density, diameterMM, diameterInch

### Sample use
Check input sheet for consistancy errors.
```Hydraulics.checkLineList(lineHydraulics,fluidList)```

Bulk of the simple hydraulic calculations.
```myDP = getSegmentDP(lineHydraulics,fittingList,fluidList)```

## Hydraulic network calculations

Hydraulic network calculations require more preparation than calculating the pressure drop in several lines. A network diagram is recommended before setting up any calculations:

* Label piping segments.
* Label nodes (source, mix point, sink).
* Label degrees of freedom (constraints: flowrate through a segment, pressure at a point).
* Note segments that are "wild", such as a control valve that maintains a specified pressure downstream. The pressure drop in this node is unknown.

The hydraulic network calculation application is described with the following examples for calculation methodology and problem definitions (inputs):

* demoHydraulicsNetworkRev1.ipynb
* networkHydraulicsSampleRev1.xls

The spreadsheet is more complex than the simple hydraulic application, and contains the following tabs.

* fluidList
	* Defines physical properties (density, viscosity) for the fluids that are used for line sizing.
* nodeList
	* Defines each of the nodes (sources, junction points, sinks) in the piping network.
* connectivity
	* Defines each piping segment, and the connecting nodes (inlet and outlet)
* DoF (Degrees of Freedom)
	* Defines the known values for flowrate through a segment, or the pressure at a node.
* hydraulics
	* Defines the piping properties Nominal Pipe Size, Schedule, material type (for roughness) and the fluid contained.
	* Provides an initial guess for the mass flow rate for each segment (if known).
* fittingList
	* A large list to define the number of fittings and pipe length for each pipe segment.
* Other tabs are provided for reference and to provide the list of appropriate options for lookup tables for NPS, Schedule, piping material, fitting types.

CSV files could also be used to describe the inputs.

The connectivity input sheet requires the following information for describing the inlet and outlet for each piping segment, even if frictional pressure drop is not calculated for a segment.

* Segment: Descriptor of the segment, typically a PFD stream number or a PnID line number. Must be unique in this table. This could also be alphanumeric (A, B, C, or 1, 2, 3). A segment must have the same fluid flow rate, fluid properties, line size.
* inNode: Discrete field to define which node is associated with the upstream pressure.
* outNode: Discrete field to define which node is associated with the downstream pressure.
* segmentType (I don't think this is used, but it is convenient for our reference.

The spreadsheet provides a number of checks to display information about the nodes (from the nodeList sheet) and the piping segments (from the hydraulics sheet).

The network hydraulics input sheet requires the following information for calculating frictional pressure drop in each segment:

* Segment: Descriptor of the segment, typically a PFD stream number or a PnID line number. Must be unique in this table. This could also be alphanumeric (A, B, C, or 1, 2, 3). A segment must have the same fluid flow rate, fluid properties, line size.
* Description: text descrption of the line.
* LineTag: detailed text description of the line. Typically the line number from the PnID.
* NPS: nominal pipe diameter (taken from pipe NPS table) to identify the pipe diameter.
* Schedule: discrete field (taken from the pipe schedule table) to identify the pipe wall thickness.
* Material: discrete field (taken from the pipe material table) to identify the pipe wall rooughness.
* PnID: optional, text description of the drawing reference.
* fluidName: discrete field (taken from fluid property table) to identify the fluid.
* massFlow: estimate of the actual flowrate in kg/h. If the mass flow rate is not provided, then an estimate is calculated based on the pressure gradient 20 kPa per 100 m.

Pipe lengths and fitting counts are tabulated separately for each line segment that requires a frictional pressure drop calculation:

* Segment: Descriptor of the segment, typically a PFD stream number or a PnID line number. Must be a unique in this table. This could also be alphanumeric (A, B, C, or 1, 2, 3). A segment must have the same fluid flow rate, fluid properties, line size.
* fittingType: discrete field (taken from the fitting3K list) to identify the type of fitting. Example, PIPE or EL45-THD-STD
* num_length_m: either the total number of fittings of the given type, or the total length of pipe specified.
* comment: optional comment
* revision: optional revision comment

Defining the pipe lengths and fitting counts in a large table permits easy checking of the input data. The numbers can come directly from piping isometrics. For example, we may specify a number of PIPE entries for a line segment, where each PIPE length matches a value on an isometric.


### Network hydraulic functions 

Several functions are used for network hydraulic calculations. Most of the calculations are used for defining the block matrix for each of the linear equations, assembling/solving the full matrix, and extracting results.


### Sample use

Check for consistancy errors.
```
 checkLineList(lineHydraulics,fluidList) # check for consistancy errors
```

Bulk of network hydraulic calculations.
``` 
 gBx = doNetworkHydraulics(lineHydraulics, fluidList, connectivity, dofList, nodeList, fittingList)
```

Extract results from global solution vector.
```
 ourResults = extractResults(gBx, lineHydraulics, connectivity, nodeList, fittingList)
```

Combine hydraulic inputs with extracted results.  
```
 allResults = combineResults(ourResults,lineHydraulics)
```

Get a list of all available information.  
```
 names(allResults)
```

Create tables for dumping to csv files.

```
# overall pressure profile
resultTable1 = select(allResults, [:Segment, :massFlow, :dp, :fromP, :toP])
# piping dimensions
resultTable2 = select(allResults, [:Segment, :NPS, :Schedule, :IDmm, :roughnessMM, :segmentLength])
# friction information
resultTable3 = select(allResults, [:Segment, :fluidName, :velocity_ms, :rho_kgm3, :mu_mPas, :Re, :frictF])
```

## Utility DataFrames

These tables are provided for your convenience. Both CSV files and Julia code are provided for defining these tables.

* npsList: list of available NPS. Used for dropdownlist in spreadsheet
* schedList: list of different pipe schedules. Used for dropdown list in spreadsheet
* pipeRoughness: tabulated roughness (mm) for different piping materials. Used for dropdownlist in spreadsheet.
* fitting3K: 3K factors for fittings
* IDmm: list of pipe NPS, pipe Schedule, ID in mm and pipe OD.

## Utility Functions

The following functions are available for use:

* calcMoodyF(Reynolds,eD)
    * Reynolds number and e/D, where e is roughness (in mm) and D is ID in mm
    * return floating (or array of float) for Moody friction factor
    * Halaand correlation is used for the Moody friction factor.
    * sorry, there is no check that the pipe Reynolds number is in the turbulent regime. (Re > 2500)


# GLODAP_Section_Gridder.jl

This is a package which will enable you to compute gridded sections of bottle 
data from GLODAP. It glues together 3 data resourses into an easy to use package:

1. [GLODAP](https://www.glodap.info/): The Global Ocean Data Analysis Project (GLODAP) is a synthesis activity for ocean surface to bottom biogeochemical data collected through chemical analysis of water samples. Data for 13 core variables (salinity, oxygen, phosphate, nitrate, silicate, dissolved inorganic carbon, total alkalinity, pH, CFC-11, CFC-12, CFC-113, CCl4 and SF6) are subjected to primary and secondary quality control to identify outliers and correct for measurement biases.

    GLODAP is publicly available, discoverable, and citable. GLODAP enables quantification of the ocean carbon sink, ocean acidification and evaluation of ocean biogeochemical models.

    GLODAP was first published in 2004. The second version, GLODAPv2, was released in 2016. This forms the basis for regular updates, containing new data sets as well as updates of older ones. The first such regular update was GLODAPv2.2019.

    GLODAPv2.2022 adds data from 96 cruises to the previous release and extends coverage in time until 2021. GLODAPv2.2022 contains data from almost 1.4 million water samples collected on 1085 cruises. An accompanying manuscript has been submitted for consideration for publication in Earth System Science Data.
2. [GO-SHIP Easy Ocean](https://github.com/kkats/GO-SHIP-Easy-Ocean):
    Despite technological advances over the last several decades, ship-based hydrography remains the only method for obtaining high-quality, high spatial and vertical resolution measurements of physical, chemical, and biological parameters over the full water column essential for physical, chemical, and biological oceanography and climate science. The Global Ocean Ship-based Hydrographic Investigations Program (GO-SHIP) coordinates a network of globally sustained hydrographic sections. These data provide a unique data set that spans four decades, comprised of more than 40 cross-ocean transects. The section data are, however, difficult to use owing to inhomogeneous format. The purpose of this new temperature, salinity, and dissolved oxygen data product is to combine, reformat and grid these data measured by Conductivity-Temperature-Depth-Oxygen (CTDO) profilers in order to facilitate their use by a wider audience. The product is machine readable and readily accessible by many existing visualisation and analysis software packages. The data processing can be repeated with modifications to suit various applications such as analysis of deep ocean, validation of numerical simulation, and calibration of autonomous platforms.
    
 3. [DIVAnd](https://github.com/gher-uliege/DIVAnd.jl): DIVAnd (Data-Interpolating Variational Analysis in n dimensions) performs an n-dimensional variational analysis/gridding of arbitrarily located observations. Observations will be interpolated/analyzed on a curvilinear grid in 1, 2, 3 or more dimensions. In this sense it is a generalization of the original two-dimensional DIVA version (still available here https://github.com/gher-ulg/DIVA but not further developed anymore).

    The method bears some similarities and equivalences with Optimal Interpolation or Krigging in that it allows to create a smooth and continous field from a collection of observations, observations which can be affected by errors. The analysis method is however different in practise, allowing to take into account topological features, physical constraints etc in a natural way. The method was initially developped with ocean data in mind, but it can be applied to any field where localized observations have to be used to produce gridded fields which are "smooth".

***This package acts as a wrapper, enabling easy gridding of GLODAP data onto sections provided by GO-SHIP Easy Ocean, using DIVAnd.***

## Getting started

- Presently, the easiest way to install this toolbox is to simply clone this 
repository (`git clone https://github.com/charles-turner-1/GLODAP_Section_Gridder.jl`).
    It will soon be available through the native Julia package manager.
- This package requires that you have a MATLAB licence, and a copy of MATLAB 
installed, as it interfaces with the MATLAB versions of the GLODAP and GO-SHIP 
Easy Ocean datasets. It will also write exception data to .mat files. (As MATLAB 
is not open source, it is intended that this MATLAB dependency will be removed 
eventually. However, this is not currently a priority.)
- Default settings (data location, filenames etc.) are saved in `defaults.toml`.
These are explained in detail below.
- Top level functionality is called using one of two functions: `gridCruisePipeline` 
or `gridSectionPipeline`. 
    - `gridCruisePipeline` will grid a variable from a single cruise.
    - `gridSectionPipeline` will grid a single variable, for all repeat occupations 
    of a hydrographic section.
- Gridding functionality can also be called in more granular detail, exceptions 
or correlation lengths manually specified etc. These will be detailed at a later time.

### Reading defaults and data inspection
`src/data_inspection.jl` contains a number of functions which allows a user to 
inspect the datasets they are using. 
- `readDefaults()` parses the `defaults.toml` file. This allows the user to 
manually inspect the defaults settings.
    - `GOSHIP_DIR` specifies the directory to which the GO-SHIP Easy Ocean toolbox 
    has been cloned to.
    - `GLODAP_DIR` specifies the directory which the GLODAP dataset has been saved
    in.
    - `MASK_MATFILE` specifies the `.mat` file in which all masks being used have 
    been saved. If a user wishes to specify their own mask, they may do so 
    either by editing this file or by changing `MASK_MATFILE`.
    - `GLODAP_FILENAME` specifies the name of the GLODAP dataset. This allows 
    the user to update the GLODAP version easily.
    - `EXCEPTIONS_DIR` contains all exception data: this is where data has been 
    manually adjusted and the adjustment saved.
    - `EXCEPTIONS_FILENAME` contains a list of expocodes (both used in GLODAP and 
    GO-SHIP, see https://github.com/charles-turner-1/WOCE_GLODAP_Conversions 
    for more details) for which exceptions have been specified.
    - `VARIABLE_EXCEPTIONS` specifies bottle data which the user wishes to 
    manually exclude from gridding.
    - `HORZLEN_EXCEPTIONS` specifies data for which the user wishes to manually
    adjust the horizontal correlation length. This is particularly useful for 
    cruises where a number of stations lack profiles for a given variable, for 
    example.
    - `USERNAME`: If you wish to write a gridded to field to a netCDF file using 
    the built in functionality, changing this will change the Author field to 
    your name.
    - `OUTPUTS_DIR`: If using the built in functionality to write output fields, 
    this will set the directory to which these fields are written.
- In nearly all use cases, defaults can be manually changed to whatever the user
wishes in function calls. However, it is more straightforward to change them in 
the `defaults.toml` file. A `writeDefaults()` functionality will be implemented 
soon, however, manual adjustment of `defaults.toml` is still necessary. 
- `listAvailableMasks()` will list all masks contained in `MASK_MATFILE`. If 
called with no arguments, it will read the default mask file. If called with an 
argument, a mask file can be manually specified, eg. `listAvailableMasks(myMaskFile)`.
- `listSectionExpocodes(sectionName)` allows you to view all the expocodes found
for a given section. Alternatively, it may be called as `listSectionExpocodes(mySectionName,myExpocodeDir)`,
which allows the user to specify a different group of expocodes for a given section 
by specifying the section name and the location of the `.csv` file storing the 
expocodes.
- `listAvailableGLODAPVariables()` will list all the variables contained in the 
GLODAP dataset loaded. This saves you opening MATLAB to do so. It can alternatively 
be called on non defaults using `listAvailableGLODAPVariables(glodap_dir,glodap_file)`
to specify a different directory and GLODAP version.

## Gridding cruises and repeat hydrographic sections

### gridCruisePipeline

`gridCruisePipeline` allows the user to grid a variable from a single hydrographic
section. All arguments are positional, and all but four are optional.

##### Mandatory Arguments
- `sectionName`: A string describing the hydrographic section which a cruise 
occupied, for example "A05", "P06".
- `horzCoordinate`: A string describing the major direction of the cruise. Currently,
this must be either `longitude` or `latitude`.
- `expocode`: A string telling the function which cruise you wish to grid. For 
example, to grid the 1992 occupation of A05, set `expocode="29HE19920714"`. 
Expocodes for a section can be found using `listSectionExpocodes(sectionName)`.
- `variableName`: A string telling the function which expocode we wish to grid.
For example, set `variableName = "G2theta"` to grid temperature.

##### Optional Arguments
- `gridding`: A string specifying whether to use isobaric (default) or isopycnic 
gridding. Beware, isopycnic gridding is still in testing.
- `meanValue`: A string, specifying which mean value to use as a background. Options 
are `scalar` (default), `horzMean` (horizontal mean value), or `climatology` 
(climatological mean value). `climatology` is derived from the GLODAP climatology, 
and only available for `G2theta`, `G2salinity`, and `G2tco2` (at present).
- `epsilonVal`: A float, specifying the value Epsilon_0 (default 0.1). 
Smaller values force the gridded fields to fit more closely to the observations
- `plotResults`: A boolean (default false). Settings `plotResults=true` will 
plot out the returned field as a heatmap.
- `autoTruncateMask`: A boolean (default false). Setting `autoTruncateMask=true`
will automatically truncate the mask such that regions with no observations are 
masked out. This can be useful for partial cruises.
- `crossValidate`: A boolean (default false). Setting `crossValidate=true` will 
enable a cross validation of the gridding to determine whether a better `epsilonVal`
can be obtained, and auto adjust to return the field gridded with this value. This
is computationally expensive and does not necessarily return a better field.
- `crossValidationNum`: An integer (default 5). If `crossValidate` is set to true,
`crossValidateNum` controls the number of cross validations performed. A larger 
number should return a better final field but is more computationally expensive.

### gridSectionPipeline

`gridSectionPipeline` allows the user to grid a variable from all hydrographic 
occupations of a given section. All arguments are positional, and all but three
 are optional.

## Installation

Presently, this package can be installed by cloning this repository (`git clone https://github.com/charles-turner-1/GLODAP_Section_Gridder.jl`). It will soon
also be available through the native Julia package manager.

## Bugs
# `grid_cruise(...)`

`grid_cruise(...)` is the cleaner single-cruise entrypoint in the newer snake_case path.

## Signature

```julia
grid_cruise(
    expocode::String,
    section_name::String,
    varname::String,
    gridding::String="isobaric";
    glodap_db::Union{Nothing,String}=nothing,
)
```

## What it does

At a high level it:

1. loads section metadata from the bundled NetCDF masks
2. resolves a usable GLODAP CSV path
3. loads the requested cruise data and core gridding coordinates
4. removes a scalar mean from the target variable
5. does an initial DIVAnd fit
6. estimates correlation lengths from residual structure
7. performs a second fit using those learned correlation lengths
8. returns a `GriddedCruise` struct containing the result and intermediate metadata

## Data resolution

If `glodap_db` is omitted, `grid_cruise(...)` tries to find or create a local CSV cache automatically. See [`../data-sources.md`](../data-sources.md) for the current resolution order.

## Current limitations

This path is promising, but it is still early.

Known caveats include:

- the second-fit `pmn` settings are still hardcoded
- the broader package config story is not unified yet
- the newer path is cleaner than the old one, but not fully rolled through the whole package

## Example

```julia
gridded = grid_cruise(
    "33RO20230101",
    "A05",
    "G2theta";
    glodap_db="/data/GLODAPv2.2023_Merged_Master_File.csv",
)
```

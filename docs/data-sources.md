# Data sources

## Section masks

The newer section-loading path uses NetCDF files bundled in the repository under:

- `data/masks/`

`load_section_info(section_name)` reads one of these mask files and returns:

- horizontal grid
- pressure grid
- section mask
- the horizontal coordinate name stored in the file metadata

This avoids relying on the older MATLAB mask-loading path for the new workflow.

## GLODAP data

There are currently two GLODAP stories in the repo.

### Legacy path

A lot of the older code still expects:

- `GLODAP_DIR`
- `GLODAP_FILENAME`
- a MAT-format GLODAP dataset

That path is still present because a fair bit of the older pipeline code depends on it.

### New CSV path

The newer `grid_cruise(...)` workflow uses a CSV-format merged master file.

You can pass that explicitly with:

```julia
grid_cruise(expocode, section_name, varname; glodap_db="/path/to/GLODAPv2.2023_Merged_Master_File.csv")
```

If `glodap_db` is not supplied, the current resolution order is:

1. `ENV["GLODAP_DB"]`
2. a CSV path implied by `defaults.toml` if one is configured there
3. `~/.glodap/GLODAPv2.2023_Merged_Master_File.csv`
4. if still missing, download the default zipped archive and unpack it into `~/.glodap/`

## Why the cache exists

The earlier version of the new path still had a hardcoded machine-local CSV path inside `grid_cruise(...)`, which made it useless as a portable public API.

Caching into `~/.glodap/` is the first pass at fixing that:

- it gives a predictable local location
- it keeps the user from having to thread the path manually every time
- it still allows explicit override when needed

## Current default download

The current default archive target is the zipped GLODAPv2.2023 merged-master CSV hosted by GLODAP.

That default is convenient, but it should still be treated as an implementation detail rather than the final long-term config story.

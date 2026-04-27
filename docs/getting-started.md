# Getting started

## Current state of the package

`GLODAP_Section_Gridder.jl` is in the middle of a transition.

There are currently two broad usage styles in the repo:

1. **Older pipeline functions** like `gridCruisePipeline(...)` and `gridSectionPipeline(...)`
   - these still lean on older defaults/config patterns
   - some of that path still expects MATLAB-era data layout and MAT files
2. **Newer snake_case functions** like `load_section_info(...)`, `load_glodap_vars(...)`, and `grid_cruise(...)`
   - this is the cleaner direction for future work
   - it uses bundled NetCDF section masks and a CSV-format GLODAP merged master file

## Install

Right now the simplest setup is still to clone the repo:

```bash
git clone https://github.com/charles-turner-1/GLODAP_Section_Gridder.jl
cd GLODAP_Section_Gridder.jl
```

Then activate it in Julia and instantiate dependencies in the normal way.

## Data you need

At minimum, useful gridding runs need:

- the package checkout itself
- the bundled section mask files under `data/masks/`
- access to a GLODAP merged-master dataset

Depending on which code path you use, that GLODAP dataset may be:

- a legacy MAT file
- or a CSV merged master file

## Recommended direction

If you are starting new work, prefer the newer snake_case path where possible.

In particular:

- use `load_section_info(...)` instead of the old MATLAB-backed section loader
- use `load_glodap_vars(...)` for the CSV path
- use `grid_cruise(...)` as the cleaner single-cruise entrypoint

The rest of these docs focus mostly on that newer direction.

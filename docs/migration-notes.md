# Migration notes

This repo is gradually moving from older MATLAB/MAT-oriented internals toward a cleaner Julia-native path.

## Broad direction

The rough direction now looks like:

- **older path**: defaults-heavy pipeline functions + MATLAB/MAT-backed loaders
- **newer path**: bundled NetCDF masks + CSV-based GLODAP access + smaller snake_case functions

## Recent changes captured here

Recent work has started to make the newer path something that can actually be used outside one developer machine.

### `grid_cruise(...)` portability

The new single-cruise path now:

- accepts `glodap_db=` as an explicit keyword argument
- no longer hardcodes a private local CSV path
- can resolve a local CSV automatically
- can download and cache the default CSV archive into `~/.glodap/`

### Docstrings

The main `src/` files now have much broader docstring coverage so the code itself is less opaque while the higher-level docs are still being built out.

## What still needs work

This is not the end state.

Still to sort out properly:

- one coherent config story across legacy and new entrypoints
- clearer user-facing docs for which APIs are preferred
- examples that exercise the newer path directly
- deciding how much of the old camelCase API should be retained versus deprecated

# demo-land-model

**Summary:** a minimal Fortran land model with the same structure as JULES, CLM and ICON-land, but toy physics.

- **Structure:** a grid of cells. Each cell has one soil type, one shared 1D soil column, and several surface tiles (e.g. broadleaf, grass, bare soil). Each tile has its own canopy.
- **Physics:**
  - Soil water: 1D Richards equation (van Genuchten, implicit, modified Picard).
  - Soil heat: 1D diffusion (implicit).
  - Canopy: prescribed water and heat interception.
- **Coupling per step:** tile canopies → fraction-weighted aggregation → soil water → soil heat.
- **Inputs are files.** The namelist gives the run length, time step and file locations. Soil types, tile types, grid cells and hourly forcing are CSV files. Types are referenced by integer ID.
- **Style:** pure functions with explicit arguments and no module-level state. The driver is a library module (`src/driver_mod.f90`), and `src/main.f90` is a thin entry point.
- **Build and run** (from the project root, because the example namelist uses relative paths):

```bash
cmake -S . -B build                 # Debug by default (-O0 -g -fcheck=all, FP traps on)
cmake --build build
ctest --test-dir build              # 8 test programs
./build/land_driver config/land.nml # writes output/timeseries.csv and output/profiles.csv
# production: cmake -S . -B build_release -DCMAKE_BUILD_TYPE=Release   (-O3)
```

## Input files

All CSV files have a header line, are read by column name (column order does not matter), and may contain `#` comment lines and blank lines.

**Namelist** (`config/land.nml`): group `&run_config` is required, and `&soil_config` sets the soil column grid.

| Variable | Meaning |
|---|---|
| `soil_parameters_file`, `tile_parameters_file`, `cells_file` | locations of the three parameter CSVs |
| `precipitation_file`, `net_radiation_file`, `air_temperature_file`, `potential_evaporation_file` | locations of the forcing CSVs |
| `run_days`, `time_step_seconds` | run length and time step (required; the time step must divide 3600 s) |
| `output_directory`, `output_interval_steps` | output location and frequency (defaults `output`, 1) |
| `n_layers`, `top_layer_thickness`, `layer_growth_factor`, `initial_temperature` (in `&soil_config`) | soil column layers (default 10, 0.05 m, 1.3) and initial temperature (default 285 K) |

**Soil types** (`soil_types.csv`): `soil_id, name, theta_sat, theta_res, van_genuchten_alpha, van_genuchten_n, k_sat, specific_storage, heat_capacity_dry, conductivity_dry, conductivity_sat`.

**Tile types** (`tile_types.csv`): `tile_id, name, leaf_area_index, radiation_extinction, max_water_storage, interception_fraction, heat_capacity, canopy_air_exchange, canopy_ground_exchange`.

**Grid cells** (`cells.csv`): `cell_id, soil_id, water_table_depth, tile_fraction_<tile_id>...`, with one fraction column per tile ID. Each row's fractions must sum to 1.

**Forcing** (one CSV per variable): `hour, <cell_id>, <cell_id>, ...`, with one row per hour starting at hour 0. A record applies for the whole hour that starts at its `hour`. Column headers must match the `cell_id` values in `cells.csv`. Units are SI: precipitation and potential evaporation in kg m-2 s-1, net radiation in W m-2, air temperature in K. The files must cover the whole run.

**Types are referenced by integer ID, never by name.** The `name` column in the type files is only a record of what each ID means. The example config in `config/` uses these IDs:

| Soil ID | Name | | Tile ID | Name |
|---|---|---|---|---|
| 1 | sand | | 1 | broadleaf |
| 2 | loam | | 2 | grass |
| 3 | clay_loam | | 3 | bare_soil |

IDs must be unique but need not be contiguous. Invalid input stops the run with a message naming the file and the problem. Examples: unknown ID, duplicate ID, fractions not summing to 1, unphysical parameter values, missing column, or forcing shorter than the run.

The example forcing in `config/forcing/` is synthetic and identical for every cell: a diurnal cycle with 10 mm/h rain for 3 h on day 2 and 4 mm/h for 2 h on day 4, stored as hourly means.

## Layout

| Path | Purpose |
|---|---|
| `src/main.f90` | program entry point (`land_driver [namelist]`) |
| `src/driver_mod.f90` | `read_run_settings`, `run_land_model`: reads inputs, builds the grid, steps the model, writes output |
| `src/csv_mod.f90`, `src/parameters_mod.f90`, `src/forcing_mod.f90` | CSV reading, soil/tile/cell input with validation, hourly forcing |
| `src/types_mod.f90` | `GridCell`, `SoilColumn`, `SoilType`, `Tile`, `CanopyParameters`, `CanopyState`, `TileFluxes`, `AtmosphericForcing`, `CellDefinition`, `CellDiagnostics` |
| `src/soil_properties_mod.f90` | van Genuchten θ(ψ), K(ψ), C(ψ); thermal conductivity and heat capacity |
| `src/soil_water_mod.f90`, `src/soil_heat_mod.f90`, `src/canopy_mod.f90` | the three physics components |
| `src/land_mod.f90` | `build_grid_cell`, `aggregate_tile_fluxes`, `cell_step`, `land_step` |
| `src/tridiagonal_mod.f90`, `io_mod.f90`, `errors_mod.f90`, `constants_mod.f90`, `kinds_mod.f90` | support |
| `config/` | example namelist, parameter CSVs and forcing |
| `tests/` | one program per module, plus an end-to-end driver test |

## Model details

**Grid cell structure.**
- A `GridCell` holds an index into the soil library, one `SoilColumn`, and an array of `Tile`s whose fractions sum to 1.
- Each `Tile` holds an index into the tile library (a `CanopyParameters` set) and its own `CanopyState`.
- This follows the JULES approach: tiles act at the surface and canopy, and the soil column is shared.
- IDs from the input files are resolved to library positions when the files are read.

**Canopy** (`canopy_step`).
- Rain interception: a fraction of the rain over the canopy is caught in a store capped at `max_water_storage`. The overflow becomes throughfall.
- Store evaporation: the store evaporates in proportion to its wetness.
- Radiation: the canopy absorbs a Beer's-law fraction of net radiation. The rest reaches the ground.
- Canopy temperature: relaxes implicitly towards air temperature.
- Ground heat flux: transmitted radiation plus a canopy-ground exchange term.
- Bare soil: LAI and storage are zero, so all rain and radiation pass through.

**Soil water.**
- Equation: mixed-form Richards equation, with matric head as the unknown. Cell-centred finite volumes, implicit Euler, modified Picard iteration.
- Top boundary: throughfall minus soil evaporation. Infiltration is capped by the surface capacity, and the excess is reported as surface runoff.
- Bottom boundary: free drainage.
- Saturated soil: water content grows linearly with head through a specific-storage term, so mass is conserved even when saturated.
- Dry soil: each Picard update is limited to `max(1 m, 0.5·|ψ|)`. Without this, the iteration diverges in very dry soil, where the capacity is tiny.
- Non-convergence: if Picard fails, the step is halved recursively.

**Soil heat.**
- Equation: implicit diffusion with θ-dependent conductivity and heat capacity.
- Top boundary: aggregated ground heat flux, minus latent heat of soil evaporation.
- Bottom boundary: zero flux.
- No phase change.
- The top boundary coupling is explicit in the soil temperature. It is stable while `exchange · dt` is below the heat capacity of the top layer.

**Forcing in time.** Forcing is piecewise constant per hour: a model step uses the record of the hour containing its start time. The time step must divide 3600 s so that no step straddles two records.

## Tests

| Test | Checks |
|---|---|
| `tridiagonal` | known solution |
| `csv` | field splitting, comments, CRLF, column lookup by name; failure on a missing column, non-numeric value, ragged row or missing file |
| `parameters` | reading soil, tile and cell files, ID-to-index mapping with non-contiguous IDs; failure on duplicate or unknown IDs, unphysical values, missing columns, or fractions not summing to 1 |
| `soil_water` | retention inversion, per-step mass conservation, steady infiltration (drainage → supply) |
| `soil_heat` | energy conservation, and the constant-flux half-space analytic solution (3% tolerance) |
| `canopy` | per-step water conservation, store bounds, bare-soil pass-through |
| `land` | fraction-weighted aggregation, 3-day multi-cell run (water budget < 1e-6 kg/m² per step, bounded T and θ, solver converged) |
| `driver` | namelist reading, a full run from CSV files with forcing columns matched by cell ID, water budget; failure on a bad time step, short forcing or missing namelist values |

Failure cases are tested by re-running the test executable with a case name and checking that it exits with an error. Water budget of a cell step: `P = Δcanopy + E_canopy + E_soil + runoff + drainage + Δsoil`. It is reported as `water_balance_error` in the output.

## Limitations (deliberate)

- No snow, soil freezing or phase change.
- No vegetation transpiration or root uptake.
- Energy is not fully closed. Canopy evaporation does not cool the canopy, and rain does not carry heat.
- Forcing is hourly only, with no interpolation.
- All cells share the same layer structure and initial temperature.
- Serial code only.

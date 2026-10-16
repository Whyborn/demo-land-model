
**AI Disclaimer**: This toy model was made using Claude Sonnet 5.5, as the intention for this is not to do science, but to be a target for technical demonstrations. The precise contents are not important, just that it *looks* somewhat like a classical land model. It was built with the prompt:

"I want to make a minimal land model in Fortran that structurally looks similar to production land models like JULES, CLM or ICON-land, in that it follows the same fundamental model structure- 1D soil columns, with each grid cell having a soil type and a series of tiles, representing different surface types. Implement minimal physics- a basic 1D richards model for soil water, a simple 1D diffusion soil heat model and a basic canopy model that prescribes some amount of water and heat interception."

# demo-land-model

**Summary:** a minimal Fortran land model with the same structure as JULES, CLM and ICON-land, but toy physics.

- **Structure:** a grid of cells. Each cell has one soil type, one shared 1D soil column, and several surface tiles (broadleaf, grass, bare soil). Each tile has its own canopy.
- **Physics:**
  - Soil water: 1D Richards equation (van Genuchten, implicit, modified Picard).
  - Soil heat: 1D diffusion (implicit).
  - Canopy: prescribed water and heat interception.
- **Coupling per step:** tile canopies → fraction-weighted aggregation → soil water → soil heat.
- **Style:** pure functions with explicit arguments and no module-level state. All state lives in derived types that are passed in and returned.
- **Build and run:**

```bash
cmake -S . -B build                 # Debug by default (-O0 -g -fcheck=all, FP traps on)
cmake --build build
ctest --test-dir build              # 5 unit-test programs
./build/land_driver config/land.nml # writes output/timeseries.csv and output/profiles.csv
# production: cmake -S . -B build_release -DCMAKE_BUILD_TYPE=Release   (-O3)
```

## Layout

| Path | Purpose |
|---|---|
| `src/types_mod.f90` | `GridCell`, `SoilColumn`, `SoilType`, `Tile`, `CanopyParameters`, `CanopyState`, `TileFluxes`, `AtmosphericForcing`, `CellDiagnostics` |
| `src/soil_properties_mod.f90` | van Genuchten θ(ψ), K(ψ), C(ψ); thermal conductivity and heat capacity; soil library (sand, loam, clay loam) |
| `src/soil_water_mod.f90` | Richards solver (see below) |
| `src/soil_heat_mod.f90` | heat diffusion solver |
| `src/canopy_mod.f90` | `canopy_step`, plus tile library (broadleaf, grass, bare soil) |
| `src/land_mod.f90` | `build_grid_cell`, `aggregate_tile_fluxes`, `cell_step`, `land_step` |
| `src/forcing_mod.f90` | synthetic diurnal forcing and two rain events |
| `src/tridiagonal_mod.f90`, `io_mod.f90`, `constants_mod.f90`, `kinds_mod.f90` | support |
| `app/main.f90` | namelist-driven driver |
| `config/land.nml` | grid, soil, tile fractions and run settings |
| `tests/` | one program per physics module, plus whole-model tests |

## Model details

**Grid cell structure.**
- A `GridCell` holds a `soil_type_index` (into a soil library), one `SoilColumn`, and an array of `Tile`s whose fractions sum to 1.
- Each `Tile` holds a `tile_type_index` (into a tile library of `CanopyParameters`) and its own `CanopyState`.
- This follows the JULES approach: tiles act at the surface and canopy, and the soil column is shared.

**Canopy** (`canopy_step`).
- Rain interception: a fraction of the rain over the canopy is caught in a store capped at `max_water_storage`. The overflow becomes throughfall.
- Store evaporation: the store evaporates in proportion to its wetness.
- Radiation: the canopy absorbs a Beer's-law fraction of net radiation. The rest reaches the ground.
- Canopy temperature: relaxes implicitly towards air temperature.
- Ground heat flux: transmitted radiation plus a canopy-ground exchange term.
- Bare soil: LAI and storage are zero, so all rain and radiation pass through and the canopy temperature equals the air temperature.

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
- The top boundary coupling is explicit in the soil temperature. It is stable while `exchange · dt` is below the heat capacity of the top layer (about 0.2 of it for the defaults).

**Forcing.** Synthetic: a diurnal temperature, radiation and PET cycle, plus rain at 10 mm/h for 3 h on day 2 and 4 mm/h for 2 h on day 4.

## Tests

| Test | Checks |
|---|---|
| `tridiagonal` | known solution |
| `soil_water` | retention inversion, per-step mass conservation, steady infiltration (drainage → supply) |
| `soil_heat` | energy conservation, and the constant-flux half-space analytic solution (3% tolerance) |
| `canopy` | per-step water conservation, store bounds, bare-soil pass-through |
| `land` | fraction-weighted aggregation, 3-day multi-cell run (water budget < 1e-6 kg/m² per step, bounded T and θ, solver converged) |

Water budget of a cell step: `P = Δcanopy + E_canopy + E_soil + runoff + drainage + Δsoil`. It is reported as `water_balance_error` in the output. The reference run in `config/land.nml` closes it to 5e-8 kg/m² or better per step.

## Limitations (deliberate)

- No snow, soil freezing or phase change.
- No vegetation transpiration or root uptake.
- Energy is not fully closed. Canopy evaporation does not cool the canopy, and rain does not carry heat.
- Atmospheric forcing is the same for every cell.
- Serial code only.
- The only ground heat flux is the one passed to the soil heat solver.

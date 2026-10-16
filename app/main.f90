! Driver: reads a namelist, builds a small grid, and steps the land model with synthetic forcing.
! Usage: land_driver [path/to/land.nml]
program land_driver
  use kinds_mod, only: wp
  use constants_mod, only: seconds_per_day
  use types_mod, only: SoilType, CanopyParameters, GridCell, AtmosphericForcing
  use soil_properties_mod, only: default_soil_library
  use canopy_mod, only: default_tile_library
  use forcing_mod, only: synthetic_forcing
  use land_mod, only: geometric_layer_thicknesses, build_grid_cell, land_step
  use io_mod, only: open_timeseries, write_timeseries_rows, write_profiles
  implicit none

  integer, parameter :: max_cells = 16
  integer, parameter :: n_tile_types = 3

  ! Namelist contents
  integer :: n_cells = 4
  integer :: soil_type_index(max_cells) = 2
  real(wp) :: water_table_depth(max_cells) = 2.0_wp
  real(wp) :: tile_fractions(n_tile_types, max_cells) = 0.0_wp
  integer :: n_layers = 10
  real(wp) :: top_layer_thickness = 0.05_wp
  real(wp) :: layer_growth_factor = 1.3_wp
  real(wp) :: initial_temperature = 285.0_wp
  real(wp) :: time_step_seconds = 1800.0_wp
  real(wp) :: run_days = 5.0_wp
  integer :: output_interval_steps = 1
  character(len=256) :: output_directory = 'output'
  namelist /grid_config/ n_cells, soil_type_index, water_table_depth, tile_fractions
  namelist /soil_config/ n_layers, top_layer_thickness, layer_growth_factor, initial_temperature
  namelist /run_config/ time_step_seconds, run_days, output_interval_steps, output_directory

  type(SoilType) :: soil_library(3)
  type(CanopyParameters) :: tile_library(n_tile_types)
  type(GridCell), allocatable :: cells(:)
  real(wp), allocatable :: layer_thickness(:)
  character(len=256) :: namelist_path
  integer :: namelist_unit, io_status, timeseries_unit, cell_index, step_index, n_steps
  type(AtmosphericForcing) :: forcing

  namelist_path = 'config/land.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, namelist_path)

  open(newunit=namelist_unit, file=trim(namelist_path), status='old', action='read', iostat=io_status)
  if (io_status /= 0) error stop 'land_driver: cannot open namelist file'
  read(namelist_unit, nml=grid_config, iostat=io_status)
  if (io_status /= 0) error stop 'land_driver: error reading &grid_config'
  read(namelist_unit, nml=soil_config, iostat=io_status)
  if (io_status /= 0) error stop 'land_driver: error reading &soil_config'
  read(namelist_unit, nml=run_config, iostat=io_status)
  if (io_status /= 0) error stop 'land_driver: error reading &run_config'
  close(namelist_unit)

  if (n_cells < 1 .or. n_cells > max_cells) error stop 'land_driver: n_cells out of range'

  soil_library = default_soil_library()
  tile_library = default_tile_library()
  layer_thickness = geometric_layer_thicknesses(n_layers, top_layer_thickness, layer_growth_factor)

  allocate(cells(n_cells))
  do cell_index = 1, n_cells
    cells(cell_index) = build_grid_cell(soil_type_index(cell_index), soil_library(soil_type_index(cell_index)), &
                                        tile_fractions(:, cell_index), layer_thickness, &
                                        water_table_depth(cell_index), initial_temperature)
  end do

  call execute_command_line('mkdir -p ' // trim(output_directory))
  call open_timeseries(trim(output_directory) // '/timeseries.csv', timeseries_unit)

  n_steps = nint(run_days * seconds_per_day / time_step_seconds)
  do step_index = 1, n_steps
    forcing = synthetic_forcing((step_index - 1) * time_step_seconds)
    cells = land_step(cells, soil_library, tile_library, forcing, time_step_seconds)
    if (mod(step_index, output_interval_steps) == 0) then
      call write_timeseries_rows(timeseries_unit, step_index * time_step_seconds, cells)
    end if
  end do

  close(timeseries_unit)
  call write_profiles(trim(output_directory) // '/profiles.csv', cells)
  write(*, '(a,i0,a,i0,a)') 'land_driver: ran ', n_steps, ' steps over ', n_cells, ' cells; output in ' &
      // trim(output_directory)
end program land_driver

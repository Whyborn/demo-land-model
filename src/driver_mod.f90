! Model driver: reads the run namelist and the input CSV files, builds the grid and steps the model.
module driver_mod
  use kinds_mod, only: wp
  use constants_mod, only: seconds_per_day
  use errors_mod, only: fail_with_message
  use types_mod, only: SoilType, CanopyParameters, CellDefinition, GridCell, AtmosphericForcing
  use parameters_mod, only: read_soil_library, read_tile_library, read_cell_definitions
  use forcing_mod, only: ForcingSeries, read_forcing_series, forcing_at_time, forcing_hours_available, &
                         seconds_per_forcing_record
  use land_mod, only: geometric_layer_thicknesses, build_grid_cell, land_step
  use io_mod, only: open_timeseries, write_timeseries_rows, write_profiles
  implicit none
  private

  public :: RunSettings, read_run_settings, run_land_model

  integer, parameter :: path_length = 256

  type :: RunSettings
    character(len=path_length) :: soil_parameters_file = ''
    character(len=path_length) :: tile_parameters_file = ''
    character(len=path_length) :: cells_file = ''
    character(len=path_length) :: precipitation_file = ''
    character(len=path_length) :: net_radiation_file = ''
    character(len=path_length) :: air_temperature_file = ''
    character(len=path_length) :: potential_evaporation_file = ''
    character(len=path_length) :: output_directory = 'output'
    real(wp) :: run_days = -1.0_wp
    real(wp) :: time_step_seconds = -1.0_wp
    integer :: output_interval_steps = 1
    integer :: n_layers = 10
    real(wp) :: top_layer_thickness = 0.05_wp
    real(wp) :: layer_growth_factor = 1.3_wp
    real(wp) :: initial_temperature = 285.0_wp
  end type RunSettings

contains

  ! Read the &run_config and &soil_config namelist groups. File locations, run length and time step
  ! have no defaults and must be given.
  function read_run_settings(namelist_path) result(settings)
    character(len=*), intent(in) :: namelist_path
    type(RunSettings) :: settings

    character(len=path_length) :: soil_parameters_file, tile_parameters_file, cells_file
    character(len=path_length) :: precipitation_file, net_radiation_file, air_temperature_file
    character(len=path_length) :: potential_evaporation_file, output_directory
    real(wp) :: run_days, time_step_seconds
    integer :: output_interval_steps
    integer :: n_layers
    real(wp) :: top_layer_thickness, layer_growth_factor, initial_temperature
    namelist /run_config/ soil_parameters_file, tile_parameters_file, cells_file, precipitation_file, &
        net_radiation_file, air_temperature_file, potential_evaporation_file, run_days, time_step_seconds, &
        output_directory, output_interval_steps
    namelist /soil_config/ n_layers, top_layer_thickness, layer_growth_factor, initial_temperature

    integer :: namelist_unit, io_status

    soil_parameters_file = settings%soil_parameters_file
    tile_parameters_file = settings%tile_parameters_file
    cells_file = settings%cells_file
    precipitation_file = settings%precipitation_file
    net_radiation_file = settings%net_radiation_file
    air_temperature_file = settings%air_temperature_file
    potential_evaporation_file = settings%potential_evaporation_file
    output_directory = settings%output_directory
    run_days = settings%run_days
    time_step_seconds = settings%time_step_seconds
    output_interval_steps = settings%output_interval_steps
    n_layers = settings%n_layers
    top_layer_thickness = settings%top_layer_thickness
    layer_growth_factor = settings%layer_growth_factor
    initial_temperature = settings%initial_temperature

    open(newunit=namelist_unit, file=namelist_path, status='old', action='read', iostat=io_status)
    if (io_status /= 0) call fail_with_message('cannot open namelist file ' // namelist_path)
    read(namelist_unit, nml=run_config, iostat=io_status)
    if (io_status /= 0) call fail_with_message('cannot read &run_config from ' // namelist_path)
    rewind(namelist_unit)
    read(namelist_unit, nml=soil_config, iostat=io_status)
    if (io_status /= 0) call fail_with_message('cannot read &soil_config from ' // namelist_path)
    close(namelist_unit)

    settings = RunSettings(soil_parameters_file=soil_parameters_file, tile_parameters_file=tile_parameters_file, &
                           cells_file=cells_file, precipitation_file=precipitation_file, &
                           net_radiation_file=net_radiation_file, air_temperature_file=air_temperature_file, &
                           potential_evaporation_file=potential_evaporation_file, output_directory=output_directory, &
                           run_days=run_days, time_step_seconds=time_step_seconds, &
                           output_interval_steps=output_interval_steps, n_layers=n_layers, &
                           top_layer_thickness=top_layer_thickness, layer_growth_factor=layer_growth_factor, &
                           initial_temperature=initial_temperature)
    call require_valid_settings(settings)
  end function read_run_settings

  subroutine require_valid_settings(settings)
    type(RunSettings), intent(in) :: settings
    real(wp) :: steps_per_record

    if (len_trim(settings%soil_parameters_file) == 0) call fail_with_message('soil_parameters_file is required')
    if (len_trim(settings%tile_parameters_file) == 0) call fail_with_message('tile_parameters_file is required')
    if (len_trim(settings%cells_file) == 0) call fail_with_message('cells_file is required')
    if (len_trim(settings%precipitation_file) == 0) call fail_with_message('precipitation_file is required')
    if (len_trim(settings%net_radiation_file) == 0) call fail_with_message('net_radiation_file is required')
    if (len_trim(settings%air_temperature_file) == 0) call fail_with_message('air_temperature_file is required')
    if (len_trim(settings%potential_evaporation_file) == 0) then
      call fail_with_message('potential_evaporation_file is required')
    end if
    if (settings%run_days <= 0.0_wp) call fail_with_message('run_days must be given and positive')
    if (settings%time_step_seconds <= 0.0_wp) call fail_with_message('time_step_seconds must be given and positive')
    steps_per_record = seconds_per_forcing_record / settings%time_step_seconds
    if (abs(steps_per_record - nint(steps_per_record)) > 1.0e-9_wp) then
      call fail_with_message('time_step_seconds must divide the hourly forcing interval (3600 s)')
    end if
    if (settings%output_interval_steps < 1) call fail_with_message('output_interval_steps must be at least 1')
    if (settings%n_layers < 1) call fail_with_message('n_layers must be at least 1')
    if (settings%top_layer_thickness <= 0.0_wp .or. settings%layer_growth_factor < 1.0_wp) then
      call fail_with_message('need top_layer_thickness > 0 and layer_growth_factor >= 1')
    end if
  end subroutine require_valid_settings

  function build_cells(definitions, soil_library, layer_thickness, initial_temperature) result(cells)
    type(CellDefinition), intent(in) :: definitions(:)
    type(SoilType), intent(in) :: soil_library(:)
    real(wp), intent(in) :: layer_thickness(:), initial_temperature
    type(GridCell), allocatable :: cells(:)

    integer :: cell_index
    allocate(cells(size(definitions)))
    do cell_index = 1, size(definitions)
      cells(cell_index) = build_grid_cell(definitions(cell_index)%soil_index, &
                                          soil_library(definitions(cell_index)%soil_index), &
                                          definitions(cell_index)%tile_fractions, layer_thickness, &
                                          definitions(cell_index)%water_table_depth, initial_temperature)
      cells(cell_index)%cell_id = definitions(cell_index)%cell_id
    end do
  end function build_cells

  subroutine run_land_model(namelist_path)
    character(len=*), intent(in) :: namelist_path

    type(RunSettings) :: settings
    type(SoilType), allocatable :: soil_library(:)
    type(CanopyParameters), allocatable :: tile_library(:)
    type(CellDefinition), allocatable :: definitions(:)
    type(ForcingSeries) :: forcing_series
    type(GridCell), allocatable :: cells(:)
    real(wp), allocatable :: layer_thickness(:)
    integer, allocatable :: cell_ids(:)
    integer :: n_steps, step_index, timeseries_unit

    settings = read_run_settings(namelist_path)
    soil_library = read_soil_library(trim(settings%soil_parameters_file))
    tile_library = read_tile_library(trim(settings%tile_parameters_file))
    definitions = read_cell_definitions(trim(settings%cells_file), soil_library, tile_library)
    cell_ids = definitions%cell_id
    forcing_series = read_forcing_series(trim(settings%precipitation_file), trim(settings%net_radiation_file), &
                                         trim(settings%air_temperature_file), &
                                         trim(settings%potential_evaporation_file), cell_ids)

    n_steps = nint(settings%run_days * seconds_per_day / settings%time_step_seconds)
    if (n_steps * settings%time_step_seconds > forcing_hours_available(forcing_series) &
                                               * seconds_per_forcing_record + 1.0e-6_wp) then
      call fail_with_message('forcing files cover fewer hours than the requested run length')
    end if

    layer_thickness = geometric_layer_thicknesses(settings%n_layers, settings%top_layer_thickness, &
                                                  settings%layer_growth_factor)
    cells = build_cells(definitions, soil_library, layer_thickness, settings%initial_temperature)

    call execute_command_line('mkdir -p ' // trim(settings%output_directory))
    call open_timeseries(trim(settings%output_directory) // '/timeseries.csv', timeseries_unit)
    do step_index = 1, n_steps
      cells = land_step(cells, soil_library, tile_library, &
                        forcing_at_time(forcing_series, (step_index - 1) * settings%time_step_seconds), &
                        settings%time_step_seconds)
      if (mod(step_index, settings%output_interval_steps) == 0) then
        call write_timeseries_rows(timeseries_unit, step_index * settings%time_step_seconds, cells)
      end if
    end do
    close(timeseries_unit)
    call write_profiles(trim(settings%output_directory) // '/profiles.csv', cells)

    write(*, '(a,i0,a,i0,a,a)') 'land_driver: ran ', n_steps, ' steps over ', size(cells), &
        ' cells; output in ', trim(settings%output_directory)
  end subroutine run_land_model

end module driver_mod

! End-to-end test: write a complete configuration, run the driver module, check the output.
program test_driver
  use kinds_mod, only: wp
  use driver_mod, only: RunSettings, read_run_settings, run_land_model
  use csv_mod, only: CsvTable, read_csv_table, csv_integer_column, csv_real_column
  use test_utilities_mod, only: assert_true, assert_close, write_text_file, expect_failure_of_subprocess
  implicit none

  character(len=*), parameter :: work_directory = 'test_driver_files'
  character(len=1), parameter :: newline = achar(10)
  character(len=*), parameter :: namelist_file = work_directory // '/land.nml'
  character(len=32) :: case_name

  if (command_argument_count() >= 1) then
    call get_command_argument(1, case_name)
    call run_failing_case(trim(case_name))
    stop
  end if

  call execute_command_line('rm -rf ' // work_directory // ' && mkdir -p ' // work_directory)
  call write_configuration(1800.0_wp, 12)
  call test_read_settings()
  call test_run()
  call expect_failure_of_subprocess('time_step_not_dividing_hour', 'a time step that does not divide 3600 s stops the program')
  call expect_failure_of_subprocess('forcing_too_short', 'forcing shorter than the run stops the program')
  call expect_failure_of_subprocess('missing_namelist_value', 'a namelist without file locations stops the program')

contains

  ! Two cells (IDs 5 and 6) with different soils; hourly forcing CSVs with rain in hours 10 and 11.
  subroutine write_configuration(time_step_seconds, forcing_hours)
    real(wp), intent(in) :: time_step_seconds
    integer, intent(in) :: forcing_hours
    character(len=:), allocatable :: precipitation_text, radiation_text, temperature_text, evaporation_text
    character(len=64) :: step_text
    integer :: hour

    write(step_text, '(f10.1)') time_step_seconds
    call write_text_file(namelist_file, &
        '&run_config' // newline &
        // "  soil_parameters_file = '" // work_directory // "/soil_types.csv'" // newline &
        // "  tile_parameters_file = '" // work_directory // "/tile_types.csv'" // newline &
        // "  cells_file = '" // work_directory // "/cells.csv'" // newline &
        // "  precipitation_file = '" // work_directory // "/precipitation.csv'" // newline &
        // "  net_radiation_file = '" // work_directory // "/net_radiation.csv'" // newline &
        // "  air_temperature_file = '" // work_directory // "/air_temperature.csv'" // newline &
        // "  potential_evaporation_file = '" // work_directory // "/potential_evaporation.csv'" // newline &
        // '  run_days = 0.5' // newline &
        // '  time_step_seconds = ' // trim(adjustl(step_text)) // newline &
        // "  output_directory = '" // work_directory // "/output'" // newline &
        // '/' // newline // '&soil_config' // newline // '  n_layers = 8' // newline // '/' // newline)
    call write_text_file(work_directory // '/soil_types.csv', &
        'soil_id,name,theta_sat,theta_res,van_genuchten_alpha,van_genuchten_n,k_sat,specific_storage,' &
        // 'heat_capacity_dry,conductivity_dry,conductivity_sat' // newline &
        // '1,sand,0.43,0.045,14.5,2.68,8.25e-5,1e-4,1.3e6,0.3,2.2' // newline &
        // '2,clay_loam,0.41,0.095,1.9,1.31,6.4e-7,1e-4,1.1e6,0.2,1.4' // newline)
    call write_text_file(work_directory // '/tile_types.csv', &
        'tile_id,name,leaf_area_index,radiation_extinction,max_water_storage,interception_fraction,' &
        // 'heat_capacity,canopy_air_exchange,canopy_ground_exchange' // newline &
        // '1,forest,5,0.5,1.0,0.7,2e4,20,5' // newline // '2,bare,0,0.5,0,0,0,10,10' // newline)
    call write_text_file(work_directory // '/cells.csv', &
        'cell_id,soil_id,water_table_depth,tile_fraction_1,tile_fraction_2' // newline &
        // '5,1,2.0,0.5,0.5' // newline // '6,2,1.0,0.0,1.0' // newline)

    precipitation_text = 'hour,6,5' // newline
    radiation_text = 'hour,5,6' // newline
    temperature_text = 'hour,5,6' // newline
    evaporation_text = 'hour,5,6' // newline
    do hour = 0, forcing_hours - 1
      ! Column order differs between files on purpose: columns are matched by cell_id.
      if (hour >= 10 .and. hour < 12) then
        ! cell 6 first, then cell 5 (header order is 6,5)
        precipitation_text = precipitation_text // integer_text(hour) // ',1.0e-3,2.0e-3' // newline
      else
        precipitation_text = precipitation_text // integer_text(hour) // ',0,0' // newline
      end if
      radiation_text = radiation_text // integer_text(hour) // ',150,150' // newline
      temperature_text = temperature_text // integer_text(hour) // ',288,288' // newline
      evaporation_text = evaporation_text // integer_text(hour) // ',2e-5,2e-5' // newline
    end do
    call write_text_file(work_directory // '/precipitation.csv', precipitation_text)
    call write_text_file(work_directory // '/net_radiation.csv', radiation_text)
    call write_text_file(work_directory // '/air_temperature.csv', temperature_text)
    call write_text_file(work_directory // '/potential_evaporation.csv', evaporation_text)
  end subroutine write_configuration

  function integer_text(value) result(text)
    integer, intent(in) :: value
    character(len=:), allocatable :: text
    character(len=16) :: buffer
    write(buffer, '(i0)') value
    text = trim(buffer)
  end function integer_text

  subroutine test_read_settings()
    type(RunSettings) :: settings
    settings = read_run_settings(namelist_file)
    call assert_close(settings%time_step_seconds, 1800.0_wp, 1.0e-12_wp, 'time step is read from the namelist')
    call assert_close(settings%run_days, 0.5_wp, 1.0e-12_wp, 'run length is read from the namelist')
    call assert_true(trim(settings%soil_parameters_file) == work_directory // '/soil_types.csv', &
                     'soil parameter file location is read from the namelist')
    call assert_true(settings%n_layers == 8 .and. settings%output_interval_steps == 1, &
                     'unset optional values keep their defaults')
  end subroutine test_read_settings

  subroutine test_run()
    type(CsvTable) :: table
    integer, allocatable :: cell_ids(:)
    real(wp), allocatable :: precipitation(:), balance_error(:), runoff(:)
    logical, allocatable :: is_cell_5(:), is_cell_6(:)
    real(wp), parameter :: time_step = 1800.0_wp

    call run_land_model(namelist_file)

    table = read_csv_table(work_directory // '/output/timeseries.csv')
    cell_ids = csv_integer_column(table, 'cell_id')
    precipitation = csv_real_column(table, 'precip_kg_m2_s')
    balance_error = csv_real_column(table, 'water_balance_error_kg_m2')
    runoff = csv_real_column(table, 'runoff_kg_m2_s')
    is_cell_5 = cell_ids == 5
    is_cell_6 = cell_ids == 6

    call assert_true(count(is_cell_5) == 24 .and. count(is_cell_6) == 24, '12 h at 1800 s gives 24 rows per cell')
    call assert_close(sum(precipitation, is_cell_5) * time_step, 14.4_wp, 1.0e-9_wp, &
                      'cell 5 receives its own forcing column (2 h at 2e-3 kg m-2 s-1)')
    call assert_close(sum(precipitation, is_cell_6) * time_step, 7.2_wp, 1.0e-9_wp, &
                      'cell 6 receives its own forcing column despite different column order')
    call assert_true(maxval(abs(balance_error)) < 1.0e-6_wp, 'water balance closes in every cell and step')
    call assert_true(all(runoff >= 0.0_wp), 'runoff is never negative')
  end subroutine test_run

  subroutine run_failing_case(name_of_case)
    character(len=*), intent(in) :: name_of_case
    select case (name_of_case)
    case ('time_step_not_dividing_hour')
      call write_configuration(1000.0_wp, 12)
    case ('forcing_too_short')
      call write_configuration(1800.0_wp, 3)
    case ('missing_namelist_value')
      call write_text_file(namelist_file, '&run_config' // newline // '  run_days = 1.0' // newline // '/' // newline &
                           // '&soil_config' // newline // '/' // newline)
    end select
    call run_land_model(namelist_file)
  end subroutine run_failing_case

end program test_driver

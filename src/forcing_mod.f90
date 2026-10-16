! Hourly atmospheric forcing read from one CSV file per variable.
! Each file has the columns: hour, <cell_id>, <cell_id>, ... with one row per hour starting at hour 0.
! Values are SI: precipitation and potential evaporation in kg m-2 s-1, net radiation in W m-2,
! air temperature in K. A record applies for the whole hour that starts at its hour value.
module forcing_mod
  use kinds_mod, only: wp
  use errors_mod, only: fail_with_message
  use types_mod, only: AtmosphericForcing
  use csv_mod, only: CsvTable, read_csv_table, csv_real_column, number_of_rows, integer_to_string
  implicit none
  private

  public :: ForcingSeries, read_forcing_series, forcing_at_time, forcing_hours_available

  real(wp), parameter, public :: seconds_per_forcing_record = 3600.0_wp

  type :: ForcingSeries
    real(wp), allocatable :: precipitation(:, :)          ! (hour, cell)
    real(wp), allocatable :: net_radiation(:, :)
    real(wp), allocatable :: air_temperature(:, :)
    real(wp), allocatable :: potential_evaporation(:, :)
  end type ForcingSeries

contains

  ! Values of one forcing file for the requested cells, as an array (hour, cell).
  function read_forcing_variable(file_path, cell_ids) result(values)
    character(len=*), intent(in) :: file_path
    integer, intent(in) :: cell_ids(:)
    real(wp), allocatable :: values(:, :)

    type(CsvTable) :: table
    real(wp), allocatable :: hours(:)
    integer :: cell_index, hour_index

    table = read_csv_table(file_path)
    hours = csv_real_column(table, 'hour')
    do hour_index = 1, size(hours)
      if (abs(hours(hour_index) - real(hour_index - 1, wp)) > 1.0e-9_wp) then
        call fail_with_message('hour column in ' // trim(file_path) // ' must count 0, 1, 2, ... (row ' &
                               // integer_to_string(hour_index) // ')')
      end if
    end do

    allocate(values(number_of_rows(table), size(cell_ids)))
    do cell_index = 1, size(cell_ids)
      values(:, cell_index) = csv_real_column(table, integer_to_string(cell_ids(cell_index)))
    end do
  end function read_forcing_variable

  function read_forcing_series(precipitation_file, net_radiation_file, air_temperature_file, &
                               potential_evaporation_file, cell_ids) result(series)
    character(len=*), intent(in) :: precipitation_file, net_radiation_file, air_temperature_file
    character(len=*), intent(in) :: potential_evaporation_file
    integer, intent(in) :: cell_ids(:)
    type(ForcingSeries) :: series

    series%precipitation = read_forcing_variable(precipitation_file, cell_ids)
    series%net_radiation = read_forcing_variable(net_radiation_file, cell_ids)
    series%air_temperature = read_forcing_variable(air_temperature_file, cell_ids)
    series%potential_evaporation = read_forcing_variable(potential_evaporation_file, cell_ids)

    if (any(series%precipitation < 0.0_wp) .or. any(series%potential_evaporation < 0.0_wp)) then
      call fail_with_message('precipitation and potential evaporation must be non-negative')
    end if
    if (any(series%air_temperature < 150.0_wp) .or. any(series%air_temperature > 350.0_wp)) then
      call fail_with_message('air temperature must be given in K (found a value outside 150-350 K)')
    end if
  end function read_forcing_series

  ! Number of complete hours for which every variable has a record.
  pure function forcing_hours_available(series) result(n_hours)
    type(ForcingSeries), intent(in) :: series
    integer :: n_hours
    n_hours = min(size(series%precipitation, 1), size(series%net_radiation, 1), &
                  size(series%air_temperature, 1), size(series%potential_evaporation, 1))
  end function forcing_hours_available

  ! Forcing for every cell at a time (seconds since the start of the series), piecewise constant per hour.
  ! The caller must ensure the time lies within forcing_hours_available.
  pure function forcing_at_time(series, time_seconds) result(forcing)
    type(ForcingSeries), intent(in) :: series
    real(wp), intent(in) :: time_seconds
    type(AtmosphericForcing) :: forcing(size(series%precipitation, 2))

    integer :: hour_index
    hour_index = floor(time_seconds / seconds_per_forcing_record) + 1
    forcing%precipitation = series%precipitation(hour_index, :)
    forcing%net_radiation = series%net_radiation(hour_index, :)
    forcing%air_temperature = series%air_temperature(hour_index, :)
    forcing%potential_evaporation = series%potential_evaporation(hour_index, :)
  end function forcing_at_time

end module forcing_mod

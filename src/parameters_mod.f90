! Reading of soil types, tile types and grid cell definitions from CSV files.
! Types are identified by integer IDs (columns soil_id and tile_id). The name column is a
! human-readable record of what each ID means and is never used for lookup.
module parameters_mod
  use kinds_mod, only: wp
  use errors_mod, only: fail_with_message
  use types_mod, only: SoilType, CanopyParameters, CellDefinition
  use csv_mod, only: CsvTable, read_csv_table, column_index, csv_real_column, csv_integer_column, &
                     csv_string_column, number_of_rows, integer_to_string, field_length
  implicit none
  private

  public :: read_soil_library, read_tile_library, read_cell_definitions, find_index_by_id

  character(len=*), parameter :: tile_fraction_prefix = 'tile_fraction_'
  real(wp), parameter :: fraction_sum_tolerance = 1.0e-6_wp

contains

  ! Position of an ID in a list of IDs, or 0 when absent.
  pure function find_index_by_id(ids, id) result(position)
    integer, intent(in) :: ids(:), id
    integer :: position
    integer :: candidate
    position = 0
    do candidate = 1, size(ids)
      if (ids(candidate) == id) then
        position = candidate
        return
      end if
    end do
  end function find_index_by_id

  subroutine require_unique_ids(ids, id_column_name, file_path)
    integer, intent(in) :: ids(:)
    character(len=*), intent(in) :: id_column_name, file_path
    integer :: first, second
    do first = 1, size(ids)
      do second = first + 1, size(ids)
        if (ids(first) == ids(second)) then
          call fail_with_message('duplicate ' // id_column_name // ' ' // integer_to_string(ids(first)) &
                                 // ' in ' // trim(file_path))
        end if
      end do
    end do
  end subroutine require_unique_ids

  function read_soil_library(file_path) result(library)
    character(len=*), intent(in) :: file_path
    type(SoilType), allocatable :: library(:)

    type(CsvTable) :: table
    integer, allocatable :: ids(:)
    character(len=field_length), allocatable :: names(:)
    real(wp), allocatable :: theta_sat(:), theta_res(:), alpha(:), n_exponent(:), k_sat(:)
    real(wp), allocatable :: specific_storage(:), heat_capacity_dry(:), conductivity_dry(:), conductivity_sat(:)
    integer :: row

    table = read_csv_table(file_path)
    ids = csv_integer_column(table, 'soil_id')
    call require_unique_ids(ids, 'soil_id', file_path)
    names = csv_string_column(table, 'name')
    theta_sat = csv_real_column(table, 'theta_sat')
    theta_res = csv_real_column(table, 'theta_res')
    alpha = csv_real_column(table, 'van_genuchten_alpha')
    n_exponent = csv_real_column(table, 'van_genuchten_n')
    k_sat = csv_real_column(table, 'k_sat')
    specific_storage = csv_real_column(table, 'specific_storage')
    heat_capacity_dry = csv_real_column(table, 'heat_capacity_dry')
    conductivity_dry = csv_real_column(table, 'conductivity_dry')
    conductivity_sat = csv_real_column(table, 'conductivity_sat')

    allocate(library(number_of_rows(table)))
    do row = 1, size(library)
      library(row) = SoilType(id=ids(row), name=names(row), theta_sat=theta_sat(row), theta_res=theta_res(row), &
                              van_genuchten_alpha=alpha(row), van_genuchten_n=n_exponent(row), k_sat=k_sat(row), &
                              specific_storage=specific_storage(row), heat_capacity_dry=heat_capacity_dry(row), &
                              conductivity_dry=conductivity_dry(row), conductivity_sat=conductivity_sat(row))
      call require_valid_soil(library(row), file_path)
    end do
  end function read_soil_library

  subroutine require_valid_soil(soil, file_path)
    type(SoilType), intent(in) :: soil
    character(len=*), intent(in) :: file_path
    character(len=:), allocatable :: context
    context = ' for soil_id ' // integer_to_string(soil%id) // ' in ' // trim(file_path)
    if (soil%theta_res < 0.0_wp .or. soil%theta_res >= soil%theta_sat .or. soil%theta_sat > 1.0_wp) then
      call fail_with_message('need 0 <= theta_res < theta_sat <= 1' // context)
    end if
    if (soil%van_genuchten_alpha <= 0.0_wp) call fail_with_message('van_genuchten_alpha must be positive' // context)
    if (soil%van_genuchten_n <= 1.0_wp) call fail_with_message('van_genuchten_n must exceed 1' // context)
    if (soil%k_sat <= 0.0_wp) call fail_with_message('k_sat must be positive' // context)
    if (soil%specific_storage <= 0.0_wp) call fail_with_message('specific_storage must be positive' // context)
    if (soil%heat_capacity_dry <= 0.0_wp .or. soil%conductivity_dry <= 0.0_wp &
        .or. soil%conductivity_sat < soil%conductivity_dry) then
      call fail_with_message('need heat_capacity_dry > 0 and 0 < conductivity_dry <= conductivity_sat' // context)
    end if
  end subroutine require_valid_soil

  function read_tile_library(file_path) result(library)
    character(len=*), intent(in) :: file_path
    type(CanopyParameters), allocatable :: library(:)

    type(CsvTable) :: table
    integer, allocatable :: ids(:)
    character(len=field_length), allocatable :: names(:)
    real(wp), allocatable :: leaf_area_index(:), radiation_extinction(:), max_water_storage(:)
    real(wp), allocatable :: interception_fraction(:), heat_capacity(:), canopy_air_exchange(:)
    real(wp), allocatable :: canopy_ground_exchange(:)
    integer :: row

    table = read_csv_table(file_path)
    ids = csv_integer_column(table, 'tile_id')
    call require_unique_ids(ids, 'tile_id', file_path)
    names = csv_string_column(table, 'name')
    leaf_area_index = csv_real_column(table, 'leaf_area_index')
    radiation_extinction = csv_real_column(table, 'radiation_extinction')
    max_water_storage = csv_real_column(table, 'max_water_storage')
    interception_fraction = csv_real_column(table, 'interception_fraction')
    heat_capacity = csv_real_column(table, 'heat_capacity')
    canopy_air_exchange = csv_real_column(table, 'canopy_air_exchange')
    canopy_ground_exchange = csv_real_column(table, 'canopy_ground_exchange')

    allocate(library(number_of_rows(table)))
    do row = 1, size(library)
      library(row) = CanopyParameters(id=ids(row), name=names(row), leaf_area_index=leaf_area_index(row), &
                                      radiation_extinction=radiation_extinction(row), &
                                      max_water_storage=max_water_storage(row), &
                                      interception_fraction=interception_fraction(row), &
                                      heat_capacity=heat_capacity(row), canopy_air_exchange=canopy_air_exchange(row), &
                                      canopy_ground_exchange=canopy_ground_exchange(row))
      call require_valid_tile(library(row), file_path)
    end do
  end function read_tile_library

  subroutine require_valid_tile(tile_parameters, file_path)
    type(CanopyParameters), intent(in) :: tile_parameters
    character(len=*), intent(in) :: file_path
    character(len=:), allocatable :: context
    context = ' for tile_id ' // integer_to_string(tile_parameters%id) // ' in ' // trim(file_path)
    if (tile_parameters%leaf_area_index < 0.0_wp .or. tile_parameters%radiation_extinction < 0.0_wp &
        .or. tile_parameters%max_water_storage < 0.0_wp .or. tile_parameters%heat_capacity < 0.0_wp) then
      call fail_with_message('leaf_area_index, radiation_extinction, max_water_storage and heat_capacity ' &
                             // 'must be non-negative' // context)
    end if
    if (tile_parameters%interception_fraction < 0.0_wp .or. tile_parameters%interception_fraction > 1.0_wp) then
      call fail_with_message('interception_fraction must lie in [0, 1]' // context)
    end if
    if (tile_parameters%canopy_air_exchange <= 0.0_wp .or. tile_parameters%canopy_ground_exchange < 0.0_wp) then
      call fail_with_message('canopy_air_exchange must be positive and canopy_ground_exchange non-negative' &
                             // context)
    end if
    if (tile_parameters%heat_capacity <= 0.0_wp .and. tile_parameters%leaf_area_index > 0.0_wp) then
      call fail_with_message('a vegetated tile needs a positive heat_capacity' // context)
    end if
  end subroutine require_valid_tile

  ! Columns: cell_id, soil_id, water_table_depth, and one tile_fraction_<tile_id> per library tile.
  function read_cell_definitions(file_path, soil_library, tile_library) result(definitions)
    character(len=*), intent(in) :: file_path
    type(SoilType), intent(in) :: soil_library(:)
    type(CanopyParameters), intent(in) :: tile_library(:)
    type(CellDefinition), allocatable :: definitions(:)

    type(CsvTable) :: table
    integer, allocatable :: cell_ids(:), soil_ids(:), tile_ids(:), library_soil_ids(:)
    real(wp), allocatable :: water_table_depth(:), fraction_column(:)
    integer :: row, tile_index, soil_index, column, prefix_length, io_status, column_tile_id

    table = read_csv_table(file_path)
    cell_ids = csv_integer_column(table, 'cell_id')
    call require_unique_ids(cell_ids, 'cell_id', file_path)
    soil_ids = csv_integer_column(table, 'soil_id')
    water_table_depth = csv_real_column(table, 'water_table_depth')
    tile_ids = tile_library%id
    library_soil_ids = soil_library%id

    ! Every tile_fraction_<id> column must refer to a known tile ID.
    prefix_length = len(tile_fraction_prefix)
    do column = 1, size(table%column_names)
      if (table%column_names(column)(1:prefix_length) == tile_fraction_prefix) then
        read(table%column_names(column)(prefix_length + 1:), *, iostat=io_status) column_tile_id
        if (io_status /= 0) call fail_with_message('bad tile fraction column name "' &
                                                   // trim(table%column_names(column)) // '" in ' // trim(file_path))
        if (find_index_by_id(tile_ids, column_tile_id) == 0) then
          call fail_with_message('column "' // trim(table%column_names(column)) // '" in ' // trim(file_path) &
                                 // ' refers to an unknown tile_id')
        end if
      end if
    end do

    allocate(definitions(number_of_rows(table)))
    do row = 1, size(definitions)
      allocate(definitions(row)%tile_fractions(size(tile_library)))
    end do
    do tile_index = 1, size(tile_library)
      fraction_column = csv_real_column(table, tile_fraction_prefix // integer_to_string(tile_ids(tile_index)))
      do row = 1, size(definitions)
        definitions(row)%tile_fractions(tile_index) = fraction_column(row)
      end do
    end do

    do row = 1, size(definitions)
      soil_index = find_index_by_id(library_soil_ids, soil_ids(row))
      if (soil_index == 0) then
        call fail_with_message('cell ' // integer_to_string(cell_ids(row)) // ' in ' // trim(file_path) &
                               // ' refers to unknown soil_id ' // integer_to_string(soil_ids(row)))
      end if
      if (any(definitions(row)%tile_fractions < 0.0_wp) &
          .or. abs(sum(definitions(row)%tile_fractions) - 1.0_wp) > fraction_sum_tolerance) then
        call fail_with_message('tile fractions of cell ' // integer_to_string(cell_ids(row)) // ' in ' &
                               // trim(file_path) // ' must be non-negative and sum to 1')
      end if
      if (water_table_depth(row) <= 0.0_wp) then
        call fail_with_message('water_table_depth of cell ' // integer_to_string(cell_ids(row)) // ' in ' &
                               // trim(file_path) // ' must be positive')
      end if
      definitions(row)%cell_id = cell_ids(row)
      definitions(row)%soil_index = soil_index
      definitions(row)%water_table_depth = water_table_depth(row)
    end do
  end function read_cell_definitions

end module parameters_mod

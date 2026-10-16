program test_parameters
  use kinds_mod, only: wp
  use types_mod, only: SoilType, CanopyParameters, CellDefinition
  use parameters_mod, only: read_soil_library, read_tile_library, read_cell_definitions, find_index_by_id
  use test_utilities_mod, only: assert_true, assert_close, write_text_file, expect_failure_of_subprocess
  implicit none

  character(len=*), parameter :: work_directory = 'test_parameters_files'
  character(len=*), parameter :: soil_file = work_directory // '/soil_types.csv'
  character(len=*), parameter :: tile_file = work_directory // '/tile_types.csv'
  character(len=*), parameter :: cells_file = work_directory // '/cells.csv'
  character(len=1), parameter :: newline = achar(10)
  character(len=*), parameter :: soil_header = 'soil_id,name,theta_sat,theta_res,van_genuchten_alpha,' &
      // 'van_genuchten_n,k_sat,specific_storage,heat_capacity_dry,conductivity_dry,conductivity_sat' // newline
  character(len=*), parameter :: tile_header = 'tile_id,name,leaf_area_index,radiation_extinction,' &
      // 'max_water_storage,interception_fraction,heat_capacity,canopy_air_exchange,canopy_ground_exchange' // newline
  ! Non-contiguous IDs on purpose: the reader must map ID to library position.
  character(len=*), parameter :: valid_soils = soil_header &
      // '7,sandy,0.43,0.045,14.5,2.68,8.25e-5,1e-4,1.3e6,0.3,2.2' // newline &
      // '3,loamy,0.43,0.078,3.6,1.56,2.89e-6,1e-4,1.2e6,0.25,1.6' // newline
  character(len=*), parameter :: valid_tiles = tile_header &
      // '5,forest,5.0,0.5,1.0,0.7,2e4,20,5' // newline &
      // '2,bare,0,0.5,0,0,0,10,10' // newline
  character(len=*), parameter :: valid_cells = 'cell_id,soil_id,water_table_depth,tile_fraction_2,tile_fraction_5' &
      // newline // '11,3,1.5,0.25,0.75' // newline // '12,7,2.0,1.0,0.0' // newline
  character(len=32) :: case_name

  if (command_argument_count() >= 1) then
    call get_command_argument(1, case_name)
    call run_failing_case(trim(case_name))
    stop
  end if

  call execute_command_line('mkdir -p ' // work_directory)
  call test_valid_files()
  call expect_failure_of_subprocess('duplicate_soil_id', 'duplicate soil_id stops the program')
  call expect_failure_of_subprocess('unphysical_soil', 'theta_res >= theta_sat stops the program')
  call expect_failure_of_subprocess('missing_soil_column', 'a missing soil parameter column stops the program')
  call expect_failure_of_subprocess('unknown_soil_id', 'a cell with an unknown soil_id stops the program')
  call expect_failure_of_subprocess('unknown_tile_column', 'a tile_fraction column with an unknown tile ID stops the program')
  call expect_failure_of_subprocess('fractions_not_one', 'tile fractions that do not sum to 1 stop the program')
  call expect_failure_of_subprocess('missing_tile_column', 'a missing tile_fraction column stops the program')

contains

  subroutine test_valid_files()
    type(SoilType), allocatable :: soil_library(:)
    type(CanopyParameters), allocatable :: tile_library(:)
    type(CellDefinition), allocatable :: definitions(:)

    call write_text_file(soil_file, '# comment' // newline // valid_soils)
    call write_text_file(tile_file, valid_tiles)
    call write_text_file(cells_file, valid_cells)

    soil_library = read_soil_library(soil_file)
    call assert_true(size(soil_library) == 2 .and. soil_library(1)%id == 7 .and. soil_library(2)%id == 3, &
                     'soil IDs are read in file order')
    call assert_close(soil_library(1)%k_sat, 8.25e-5_wp, 1.0e-12_wp, 'soil k_sat is read')
    call assert_close(soil_library(2)%van_genuchten_n, 1.56_wp, 1.0e-12_wp, 'soil van_genuchten_n is read')
    call assert_true(trim(soil_library(2)%name) == 'loamy', 'soil name is kept as a record')

    tile_library = read_tile_library(tile_file)
    call assert_true(size(tile_library) == 2 .and. tile_library(1)%id == 5 .and. tile_library(2)%id == 2, &
                     'tile IDs are read in file order')
    call assert_close(tile_library(1)%max_water_storage, 1.0_wp, 1.0e-12_wp, 'tile max_water_storage is read')

    definitions = read_cell_definitions(cells_file, soil_library, tile_library)
    call assert_true(size(definitions) == 2 .and. definitions(1)%cell_id == 11, 'cells are read')
    call assert_true(definitions(1)%soil_index == 2 .and. definitions(2)%soil_index == 1, &
                     'soil_id is mapped to the library position')
    ! tile_library(1) has ID 5 and tile_library(2) has ID 2, so fractions follow library order.
    call assert_close(definitions(1)%tile_fractions(1), 0.75_wp, 1.0e-12_wp, 'tile_fraction_5 maps to tile ID 5')
    call assert_close(definitions(1)%tile_fractions(2), 0.25_wp, 1.0e-12_wp, 'tile_fraction_2 maps to tile ID 2')
    call assert_close(definitions(2)%water_table_depth, 2.0_wp, 1.0e-12_wp, 'water table depth is read')
    call assert_true(find_index_by_id([7, 3], 3) == 2 .and. find_index_by_id([7, 3], 4) == 0, 'find_index_by_id')
  end subroutine test_valid_files

  subroutine write_valid_inputs()
    call write_text_file(soil_file, valid_soils)
    call write_text_file(tile_file, valid_tiles)
    call write_text_file(cells_file, valid_cells)
  end subroutine write_valid_inputs

  subroutine run_failing_case(name_of_case)
    character(len=*), intent(in) :: name_of_case
    type(SoilType), allocatable :: soil_library(:)
    type(CanopyParameters), allocatable :: tile_library(:)
    type(CellDefinition), allocatable :: definitions(:)
    character(len=*), parameter :: cells_header = 'cell_id,soil_id,water_table_depth,'

    call write_valid_inputs()
    select case (name_of_case)
    case ('duplicate_soil_id')
      call write_text_file(soil_file, valid_soils // '7,again,0.43,0.045,14.5,2.68,8.25e-5,1e-4,1.3e6,0.3,2.2' // newline)
    case ('unphysical_soil')
      call write_text_file(soil_file, soil_header // '1,bad,0.2,0.3,14.5,2.68,8.25e-5,1e-4,1.3e6,0.3,2.2' // newline)
    case ('missing_soil_column')
      call write_text_file(soil_file, 'soil_id,name,theta_sat' // newline // '1,x,0.4' // newline)
    case ('unknown_soil_id')
      call write_text_file(cells_file, cells_header // 'tile_fraction_2,tile_fraction_5' // newline // '1,99,1.0,0.5,0.5' // newline)
    case ('unknown_tile_column')
      call write_text_file(cells_file, cells_header // 'tile_fraction_2,tile_fraction_5,tile_fraction_9' // newline &
                           // '1,3,1.0,0.5,0.5,0.0' // newline)
    case ('fractions_not_one')
      call write_text_file(cells_file, cells_header // 'tile_fraction_2,tile_fraction_5' // newline // '1,3,1.0,0.5,0.4' // newline)
    case ('missing_tile_column')
      call write_text_file(cells_file, cells_header // 'tile_fraction_2' // newline // '1,3,1.0,1.0' // newline)
    end select

    soil_library = read_soil_library(soil_file)
    tile_library = read_tile_library(tile_file)
    definitions = read_cell_definitions(cells_file, soil_library, tile_library)
  end subroutine run_failing_case

end program test_parameters

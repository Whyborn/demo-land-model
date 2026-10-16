! Plain-text CSV output of cell diagnostics and column profiles.
module io_mod
  use kinds_mod, only: wp
  use constants_mod, only: seconds_per_day
  use types_mod, only: GridCell
  implicit none
  private

  public :: open_timeseries, write_timeseries_rows, write_profiles

contains

  subroutine open_timeseries(file_path, unit_number)
    character(len=*), intent(in) :: file_path
    integer, intent(out) :: unit_number
    open(newunit=unit_number, file=file_path, status='replace', action='write')
    write(unit_number, '(a)') 'day,cell_id,precip_kg_m2_s,throughfall_kg_m2_s,canopy_evap_kg_m2_s,' &
        // 'soil_evap_kg_m2_s,infiltration_kg_m2_s,runoff_kg_m2_s,drainage_kg_m2_s,ground_heat_W_m2,' &
        // 'canopy_store_kg_m2,soil_water_kg_m2,top_layer_theta,top_layer_temp_K,water_balance_error_kg_m2'
  end subroutine open_timeseries

  ! One row per cell.
  subroutine write_timeseries_rows(unit_number, time_seconds, cells)
    integer, intent(in) :: unit_number
    real(wp), intent(in) :: time_seconds
    type(GridCell), intent(in) :: cells(:)
    integer :: cell_index
    do cell_index = 1, size(cells)
      associate(diagnostics => cells(cell_index)%diagnostics)
        write(unit_number, '(f10.5,",",i0,13(",",es13.5))') time_seconds / seconds_per_day, cells(cell_index)%cell_id, &
            diagnostics%precipitation, diagnostics%throughfall, diagnostics%canopy_evaporation, &
            diagnostics%soil_evaporation, diagnostics%infiltration, diagnostics%surface_runoff, &
            diagnostics%drainage, diagnostics%ground_heat_flux, diagnostics%canopy_water_storage, &
            diagnostics%soil_water_storage, cells(cell_index)%column%theta(1), &
            cells(cell_index)%column%temperature(1), diagnostics%water_balance_error
      end associate
    end do
  end subroutine write_timeseries_rows

  subroutine write_profiles(file_path, cells)
    character(len=*), intent(in) :: file_path
    type(GridCell), intent(in) :: cells(:)
    integer :: unit_number, cell_index, layer
    real(wp) :: centre_depth

    open(newunit=unit_number, file=file_path, status='replace', action='write')
    write(unit_number, '(a)') 'cell_id,layer,depth_m,theta,matric_head_m,temperature_K'
    do cell_index = 1, size(cells)
      centre_depth = 0.0_wp
      do layer = 1, size(cells(cell_index)%column%theta)
        associate(column => cells(cell_index)%column)
          centre_depth = centre_depth + 0.5_wp * column%layer_thickness(layer)
          write(unit_number, '(i0,",",i0,4(",",es13.5))') cells(cell_index)%cell_id, layer, centre_depth, &
              column%theta(layer), column%matric_head(layer), column%temperature(layer)
          centre_depth = centre_depth + 0.5_wp * column%layer_thickness(layer)
        end associate
      end do
    end do
    close(unit_number)
  end subroutine write_profiles

end module io_mod

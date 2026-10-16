! Derived types describing the land surface hierarchy:
!   GridCell -> one SoilColumn (shared) + several Tiles (each with its own canopy)
! All quantities are SI. Water fluxes are kg m-2 s-1, heat fluxes W m-2.
module types_mod
  use kinds_mod, only: wp
  implicit none
  private

  public :: SoilType, SoilColumn, CanopyParameters, CanopyState, TileFluxes
  public :: CanopyStepResult, Tile, CellDiagnostics, GridCell, AtmosphericForcing, CellDefinition

  ! Static soil properties (van Genuchten hydraulics, simple thermal properties).
  type :: SoilType
    integer :: id = 0                            ! soil type ID used in input files
    character(len=16) :: name = 'unnamed'
    real(wp) :: theta_sat = 0.4_wp                ! saturated water content (m3 m-3)
    real(wp) :: theta_res = 0.05_wp               ! residual water content (m3 m-3)
    real(wp) :: van_genuchten_alpha = 2.0_wp      ! inverse air-entry head (m-1)
    real(wp) :: van_genuchten_n = 1.5_wp          ! pore-size distribution (-)
    real(wp) :: k_sat = 1.0e-6_wp                 ! saturated conductivity (m s-1)
    real(wp) :: specific_storage = 1.0e-4_wp      ! storage when saturated (m-1)
    real(wp) :: heat_capacity_dry = 1.2e6_wp      ! dry volumetric heat capacity (J m-3 K-1)
    real(wp) :: conductivity_dry = 0.25_wp        ! dry thermal conductivity (W m-1 K-1)
    real(wp) :: conductivity_sat = 1.6_wp         ! saturated thermal conductivity (W m-1 K-1)
  end type SoilType

  ! Prognostic state of one 1D soil column. Layer 1 is the top layer.
  type :: SoilColumn
    real(wp), allocatable :: layer_thickness(:)   ! m
    real(wp), allocatable :: matric_head(:)       ! m (negative when unsaturated)
    real(wp), allocatable :: theta(:)             ! volumetric water content (m3 m-3)
    real(wp), allocatable :: temperature(:)       ! K
  end type SoilColumn

  ! Prescribed canopy properties for one tile type.
  type :: CanopyParameters
    integer :: id = 0                            ! tile type ID used in input files
    character(len=16) :: name = 'unnamed'
    real(wp) :: leaf_area_index = 0.0_wp          ! m2 m-2
    real(wp) :: radiation_extinction = 0.5_wp     ! Beer's law coefficient (-)
    real(wp) :: max_water_storage = 0.0_wp        ! kg m-2
    real(wp) :: interception_fraction = 0.0_wp    ! fraction of rain over canopy that is caught (-)
    real(wp) :: heat_capacity = 0.0_wp            ! J m-2 K-1
    real(wp) :: canopy_air_exchange = 10.0_wp     ! canopy-air heat exchange (W m-2 K-1)
    real(wp) :: canopy_ground_exchange = 10.0_wp  ! canopy-ground heat exchange (W m-2 K-1)
  end type CanopyParameters

  type :: CanopyState
    real(wp) :: water_store = 0.0_wp              ! kg m-2
    real(wp) :: temperature = 283.15_wp           ! K
  end type CanopyState

  ! Fluxes leaving a tile towards the ground (or atmosphere), as rates.
  type :: TileFluxes
    real(wp) :: throughfall = 0.0_wp              ! water reaching ground (kg m-2 s-1)
    real(wp) :: canopy_evaporation = 0.0_wp       ! kg m-2 s-1
    real(wp) :: soil_evaporation_demand = 0.0_wp  ! kg m-2 s-1
    real(wp) :: ground_heat_flux = 0.0_wp         ! into the soil surface (W m-2)
  end type TileFluxes

  type :: CanopyStepResult
    type(CanopyState) :: state
    type(TileFluxes) :: fluxes
  end type CanopyStepResult

  type :: Tile
    integer :: tile_type_index = 1                ! index into the tile library
    real(wp) :: fraction = 0.0_wp                 ! area fraction of the grid cell
    type(CanopyState) :: canopy_state
    type(TileFluxes) :: fluxes
  end type Tile

  ! Cell-mean diagnostics of the last step. Rates are kg m-2 s-1 unless stated.
  type :: CellDiagnostics
    real(wp) :: precipitation = 0.0_wp
    real(wp) :: throughfall = 0.0_wp
    real(wp) :: canopy_evaporation = 0.0_wp
    real(wp) :: soil_evaporation = 0.0_wp
    real(wp) :: infiltration = 0.0_wp
    real(wp) :: surface_runoff = 0.0_wp
    real(wp) :: drainage = 0.0_wp
    real(wp) :: ground_heat_flux = 0.0_wp         ! W m-2
    real(wp) :: canopy_water_storage = 0.0_wp     ! kg m-2
    real(wp) :: soil_water_storage = 0.0_wp       ! kg m-2
    real(wp) :: water_balance_error = 0.0_wp      ! kg m-2 over the step
    logical :: solver_converged = .true.
  end type CellDiagnostics

  type :: GridCell
    integer :: cell_id = 0                        ! cell ID from the input files
    integer :: soil_type_index = 1                ! index into the soil library
    type(SoilColumn) :: column
    type(Tile), allocatable :: tiles(:)
    type(CellDiagnostics) :: diagnostics
  end type GridCell

  ! One grid cell as described in the cells input file. Indices refer to the soil and tile libraries.
  type :: CellDefinition
    integer :: cell_id = 0
    integer :: soil_index = 1
    real(wp) :: water_table_depth = 2.0_wp        ! m below the surface
    real(wp), allocatable :: tile_fractions(:)    ! one per tile library entry
  end type CellDefinition

  type :: AtmosphericForcing
    real(wp) :: precipitation = 0.0_wp            ! kg m-2 s-1
    real(wp) :: air_temperature = 288.15_wp       ! K
    real(wp) :: net_radiation = 0.0_wp            ! W m-2
    real(wp) :: potential_evaporation = 0.0_wp    ! kg m-2 s-1
  end type AtmosphericForcing
end module types_mod

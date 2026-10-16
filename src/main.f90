! Usage: land_driver [path/to/land.nml]
program land_driver
  use driver_mod, only: run_land_model
  implicit none

  character(len=256) :: namelist_path

  namelist_path = 'config/land.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, namelist_path)
  call run_land_model(trim(namelist_path))
end program land_driver

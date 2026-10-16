! Minimal CSV reading. The first data line is the header; columns are looked up by name.
! Blank lines and lines starting with '#' are ignored, so input files can carry comments
! (for example a record of what each type ID means). Fields are truncated to field_length characters.
module csv_mod
  use kinds_mod, only: wp
  use errors_mod, only: fail_with_message
  implicit none
  private

  public :: CsvTable, read_csv_table, split_fields, column_index, number_of_rows
  public :: csv_real_column, csv_integer_column, csv_string_column, integer_to_string

  integer, parameter, public :: field_length = 64
  integer, parameter :: max_line_length = 16384

  type :: CsvTable
    character(len=field_length), allocatable :: column_names(:)
    character(len=field_length), allocatable :: fields(:, :)    ! (row, column)
  end type CsvTable

contains

  pure function integer_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=:), allocatable :: text
    character(len=16) :: buffer
    write(buffer, '(i0)') value
    text = trim(buffer)
  end function integer_to_string

  pure function is_data_line(line) result(is_data)
    character(len=*), intent(in) :: line
    logical :: is_data
    character(len=len(line)) :: stripped
    stripped = adjustl(line)
    is_data = len_trim(stripped) > 0
    if (is_data) is_data = stripped(1:1) /= '#'
  end function is_data_line

  pure function split_fields(line) result(fields)
    character(len=*), intent(in) :: line
    character(len=field_length), allocatable :: fields(:)
    integer :: n_fields, field_index, start_position, comma_offset, line_length

    line_length = len_trim(line)
    n_fields = 1 + count([(line(field_index:field_index) == ',', field_index = 1, line_length)])
    allocate(fields(n_fields))
    start_position = 1
    do field_index = 1, n_fields
      comma_offset = index(line(start_position:line_length), ',')
      if (comma_offset == 0) then
        fields(field_index) = adjustl(line(start_position:line_length))
      else
        fields(field_index) = adjustl(line(start_position:start_position + comma_offset - 2))
        start_position = start_position + comma_offset
      end if
    end do
  end function split_fields

  function read_csv_table(file_path) result(table)
    character(len=*), intent(in) :: file_path
    type(CsvTable) :: table

    character(len=max_line_length) :: line
    character(len=field_length), allocatable :: line_fields(:)
    integer :: unit_number, io_status, n_data_lines, row, n_columns
    logical :: header_read

    open(newunit=unit_number, file=file_path, status='old', action='read', iostat=io_status)
    if (io_status /= 0) call fail_with_message('cannot open CSV file ' // trim(file_path))

    n_data_lines = 0
    do
      read(unit_number, '(a)', iostat=io_status) line
      if (io_status /= 0) exit
      if (is_data_line(line)) n_data_lines = n_data_lines + 1
    end do
    if (n_data_lines < 1) call fail_with_message('CSV file has no header: ' // trim(file_path))
    rewind(unit_number)

    header_read = .false.
    row = 0
    do
      read(unit_number, '(a)', iostat=io_status) line
      if (io_status /= 0) exit
      if (.not. is_data_line(line)) cycle
      line = strip_carriage_return(line)
      line_fields = split_fields(line)
      if (.not. header_read) then
        table%column_names = line_fields
        n_columns = size(line_fields)
        allocate(table%fields(n_data_lines - 1, n_columns))
        header_read = .true.
      else
        row = row + 1
        if (size(line_fields) /= n_columns) then
          call fail_with_message('wrong number of fields in data row ' // integer_to_string(row) &
                                 // ' of ' // trim(file_path))
        end if
        table%fields(row, :) = line_fields
      end if
    end do
    close(unit_number)
  end function read_csv_table

  pure function strip_carriage_return(line) result(cleaned)
    character(len=*), intent(in) :: line
    character(len=len(line)) :: cleaned
    integer :: last_position
    cleaned = line
    last_position = len_trim(cleaned)
    if (last_position > 0) then
      if (cleaned(last_position:last_position) == achar(13)) cleaned(last_position:last_position) = ' '
    end if
  end function strip_carriage_return

  pure function number_of_rows(table) result(n_rows)
    type(CsvTable), intent(in) :: table
    integer :: n_rows
    n_rows = size(table%fields, 1)
  end function number_of_rows

  ! Position of a named column, or 0 when the column is absent.
  pure function column_index(table, column_name) result(position)
    type(CsvTable), intent(in) :: table
    character(len=*), intent(in) :: column_name
    integer :: position
    integer :: candidate
    position = 0
    do candidate = 1, size(table%column_names)
      if (trim(table%column_names(candidate)) == column_name) then
        position = candidate
        return
      end if
    end do
  end function column_index

  function required_column_index(table, column_name) result(position)
    type(CsvTable), intent(in) :: table
    character(len=*), intent(in) :: column_name
    integer :: position
    position = column_index(table, column_name)
    if (position == 0) call fail_with_message('missing CSV column "' // column_name // '"')
  end function required_column_index

  function csv_string_column(table, column_name) result(values)
    type(CsvTable), intent(in) :: table
    character(len=*), intent(in) :: column_name
    character(len=field_length), allocatable :: values(:)
    values = table%fields(:, required_column_index(table, column_name))
  end function csv_string_column

  function csv_real_column(table, column_name) result(values)
    type(CsvTable), intent(in) :: table
    character(len=*), intent(in) :: column_name
    real(wp), allocatable :: values(:)
    integer :: row, position, io_status
    position = required_column_index(table, column_name)
    allocate(values(number_of_rows(table)))
    do row = 1, size(values)
      read(table%fields(row, position), *, iostat=io_status) values(row)
      if (io_status /= 0) then
        call fail_with_message('column "' // column_name // '" row ' // integer_to_string(row) &
                               // ': "' // trim(table%fields(row, position)) // '" is not a number')
      end if
    end do
  end function csv_real_column

  function csv_integer_column(table, column_name) result(values)
    type(CsvTable), intent(in) :: table
    character(len=*), intent(in) :: column_name
    integer, allocatable :: values(:)
    integer :: row, position, io_status
    position = required_column_index(table, column_name)
    allocate(values(number_of_rows(table)))
    do row = 1, size(values)
      read(table%fields(row, position), *, iostat=io_status) values(row)
      if (io_status /= 0) then
        call fail_with_message('column "' // column_name // '" row ' // integer_to_string(row) &
                               // ': "' // trim(table%fields(row, position)) // '" is not an integer')
      end if
    end do
  end function csv_integer_column

end module csv_mod

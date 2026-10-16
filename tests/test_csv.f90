program test_csv
  use kinds_mod, only: wp
  use csv_mod, only: CsvTable, read_csv_table, split_fields, column_index, number_of_rows, &
                     csv_real_column, csv_integer_column, csv_string_column
  use test_utilities_mod, only: assert_true, assert_close, write_text_file, expect_failure_of_subprocess
  implicit none

  character(len=*), parameter :: work_directory = 'test_csv_files'
  character(len=*), parameter :: valid_file = work_directory // '/valid.csv'
  character(len=1), parameter :: newline = achar(10)
  character(len=32) :: case_name

  if (command_argument_count() >= 1) then
    call get_command_argument(1, case_name)
    call run_failing_case(trim(case_name))
    stop
  end if

  call execute_command_line('mkdir -p ' // work_directory)
  call test_split_fields()
  call test_read_table()
  call expect_failure_of_subprocess('missing_column', 'asking for a missing column stops the program')
  call expect_failure_of_subprocess('not_a_number', 'a non-numeric value in a numeric column stops the program')
  call expect_failure_of_subprocess('ragged_row', 'a row with the wrong number of fields stops the program')
  call expect_failure_of_subprocess('missing_file', 'a missing file stops the program')

contains

  subroutine test_split_fields()
    character(len=64), allocatable :: fields(:)
    fields = split_fields('a, bb ,,3.5')
    call assert_true(size(fields) == 4, 'split_fields counts fields including empty ones')
    call assert_true(trim(fields(1)) == 'a' .and. trim(fields(2)) == 'bb', 'split_fields trims leading and trailing blanks')
    call assert_true(len_trim(fields(3)) == 0 .and. trim(fields(4)) == '3.5', 'split_fields keeps empty and last fields')
  end subroutine test_split_fields

  subroutine test_read_table()
    type(CsvTable) :: table
    real(wp), allocatable :: values(:)
    character(len=64), allocatable :: names(:)
    call write_text_file(valid_file, '# a comment line' // newline // newline // &
        'id,name,value' // achar(13) // newline // '10,alpha,1.5' // achar(13) // newline // &
        '# comment between rows' // newline // '20,beta,2.5e-3' // newline)
    table = read_csv_table(valid_file)
    call assert_true(number_of_rows(table) == 2, 'comments and blank lines are skipped; header is not a data row')
    call assert_true(column_index(table, 'value') == 3 .and. column_index(table, 'absent') == 0, &
                     'columns are found by header name')
    call assert_true(all(csv_integer_column(table, 'id') == [10, 20]), 'integer column is parsed (CRLF tolerated)')
    values = csv_real_column(table, 'value')
    call assert_close(values(2), 2.5e-3_wp, 1.0e-14_wp, 'real column is parsed')
    names = csv_string_column(table, 'name')
    call assert_true(trim(names(1)) == 'alpha', 'string column is returned')
  end subroutine test_read_table

  subroutine run_failing_case(name_of_case)
    character(len=*), intent(in) :: name_of_case
    type(CsvTable) :: table
    real(wp), allocatable :: values(:)

    select case (name_of_case)
    case ('missing_column')
      call write_text_file(valid_file, 'a,b' // newline // '1,2' // newline)
      table = read_csv_table(valid_file)
      values = csv_real_column(table, 'c')
      if (size(values) > 0) write(*, '(a)') 'unexpectedly succeeded'
    case ('not_a_number')
      call write_text_file(valid_file, 'a,b' // newline // '1,abc' // newline)
      table = read_csv_table(valid_file)
      values = csv_real_column(table, 'b')
      if (size(values) > 0) write(*, '(a)') 'unexpectedly succeeded'
    case ('ragged_row')
      call write_text_file(valid_file, 'a,b' // newline // '1,2,3' // newline)
      table = read_csv_table(valid_file)
    case ('missing_file')
      table = read_csv_table(work_directory // '/does_not_exist.csv')
    end select
  end subroutine run_failing_case

end program test_csv

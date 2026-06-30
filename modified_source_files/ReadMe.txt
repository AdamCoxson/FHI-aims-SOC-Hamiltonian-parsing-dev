Modified file "calculate_second_variational_soc.f90" to dump realspace spin orbit coupling Hamiltonians.

See lines 665 to 685 (pasted separately below). You can C&P this file into $your_fhiaims_dir/src/soc and overwrite the original file.
You will have to rebuild fhi-aims via cmake. Obviously, make a copy of the entire fhiaims just in case I've broken anything.



####################################################################################################
Add in block to line 665 in "calculate_second_variational_soc.f90"

  ! --- DUMP SOC MATRIX ---
  ! --- MODIFICATION BY Adam Coxson
  if (myid == 0) then
    ! Unformatted dump (smaller file, faster read in Python)
    open(unit=999, file="realspace_soc_matrix_unfmttd.out", form='unformatted', status='replace')
    ! Write the dimensions first so the reader script knows the size
    write(999) ld_soc_matrix
    write(999) 3
    ! Write the actual matrix block
    write(999) soc_matrix(1:ld_soc_matrix, 1:3)
    close(999)

    ! Alternative: Formatted dump (easier to debug, but larger file)
    open(unit=998, file="realspace_soc_matrix.out", status='replace')
    write(998, *) ld_soc_matrix, 3
    do i = 1, ld_soc_matrix
      write(998, '(3E24.15)') soc_matrix(i, 1:3)
      end do
    close(998)
  end if
  ! -----------------------

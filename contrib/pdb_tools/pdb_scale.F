c=======================================================================
c
!> \brief Scale the coordinates and box size in a PDB file
!>
!> This program scales the atomic coordinates and the box size in a 
!> PBD file. The scale factor is provided on the command line.
!> Negative scale factors are allowed but it is not clear what that
!> means if the PDB file has a CRYST1 record.
!>
!> The PDB file is read from standard input and the result written
!> to standard output.
!>
      program pdb_scale
      character(len=32) :: arg      ! command line argument
      character(len=80) :: line_in  ! line of input file
      character(len=80) :: line_out ! line for output file
      integer :: iline              ! counts lines of input file
      integer :: i                  ! generic counter
      double precision :: factor    ! scale factor
      double precision :: xyz(3)    ! coordinates or lattice parameters
      iline = 0
      if (command_argument_count().ne.1) then
        write(0,*)"ERROR: Invalid invocation of pdb_scale."
        write(0,*)"       The program expects 1 command line argument"
        write(0,*)"       that is the scale factor."
      endif
      call get_command_argument(1,arg)
      read(arg,*)factor
      do while (.true.)
        iline = iline + 1
        read(*,'(a80)',err=999,end=888)line_in
        line_out = line_in
        select case(line_in(1:6))
        case("ATOM  ","HETATM")
          read(line_in(31:54),'(3f8.3)')(xyz(i),i=1,3)
          do i = 1, 3
            xyz(i) = factor*xyz(i)
          enddo
          write(line_out(31:54),'(3f8.3)')(xyz(i),i=1,3)
        case("CRYST1")
          read(line_in(7:33),'(3f9.3)')(xyz(i),i=1,3)
          do i = 1, 3
            xyz(i) = factor*xyz(i)
          enddo
          write(line_out(7:33),'(3f9.3)')(xyz(i),i=1,3)
        case default
          ! do nothing: the input is echoed to the ouput
        end select
        write(*,'(a80)')line_out
      enddo
 999  write(0,*)"ERROR: I/O error reading line ",iline
 888  continue
      end program pdb_scale
c
c=======================================================================

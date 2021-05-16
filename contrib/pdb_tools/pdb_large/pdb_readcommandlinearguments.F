      subroutine pdb_readcommandlinearguments(pdb_input,pdb_output)
      implicit none
      character(len=*), intent(out) :: pdb_input
      character(len=*), intent(out) :: pdb_output
!
      integer :: narg ! the number of command line arguments
      integer :: iarg ! the current command line argument
!
      character(len=8) :: arg
!
      pdb_input = "-"
      pdb_output = "-"
      narg = command_argument_count()
      select case (narg)
        case (2, 4)
          iarg=1
          do while (iarg.lt.narg)
            call get_command_argument(iarg,arg)
            select case (arg)
              case ("-i")
                iarg=iarg+1
                call get_command_argument(iarg,pdb_input)
              case ("-o")
                iarg=iarg+1
                call get_command_argument(iarg,pdb_output)
              case default
                go to 10
            end select
            iarg=iarg+1
          enddo
        case default
          go to 10
      end select
      return
!
 10   write(*,*)'Usage: pdb_large -i inputfile -o outputfile'
      write(*,*)
      write(*,*)'inputfile : PDB file in either normal or large format'
      write(*,*)
      write(*,*)'outputfile: PDB file in the other format'
      stop 10
!     
      end

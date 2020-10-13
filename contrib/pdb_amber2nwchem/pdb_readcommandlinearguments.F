      subroutine pdb_readcommandlinearguments(pdb_input,pdb_output,
     +                                        pdb_dialect)
      implicit none
      character(len=*), intent(out) :: pdb_input
      character(len=*), intent(out) :: pdb_output
      character(len=*), intent(out) :: pdb_dialect
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
        case (2, 4, 6)
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
              case ("-f")
                iarg=iarg+1
                call get_command_argument(iarg,arg)
                if (arg(1:5).eq."amber") then
                  pdb_dialect = "amber"
                else if (arg(1:6).eq."charmm") then
                  pdb_dialect = "charmm"
                else
                  write(*,*)"Error: unsupported forcefield: ",arg
                  write(*,*)'Supported forcefields: "amber" or
     +"charmm"'
                  write(*,*)
                  go to 10
                endif
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
 10   write(*,*)'Usage: pdb_amber2nwchem -i inputfile -o outputfile'
      write(*,*)'                        [-f forcefield]'
      write(*,*)
      write(*,*)'inputfile : the name of the Amber PDB input file.'
      write(*,*)
      write(*,*)'outputfile: the name of the NWChem PDB output 
     +file.'
      write(*,*)
      write(*,*)'forcefield: the name of the forcefield (default
     +"amber").'
      stop 10
!     
      end

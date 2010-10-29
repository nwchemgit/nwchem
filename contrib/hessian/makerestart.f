      program hessrst
      implicit none
c
c  $Id$
c  This program helps create a fdrst file for numerical hessian restarts.  
c  There are several lines below that need to be modified so the program
c  picks up the correct hessian and fd_ddipole files and so that
c  some other info is set correctly.  Look for the comments.
c  First set the rank to be 3 times the number of atoms.
c
      integer rank
      parameter (rank = 141)
c
c     Restart file to be constructed
      CHARACTER*255 FILEATR
c     Ascii hessian (need for restart)
      CHARACTER*255 FILEHESS
c     derivative dipole file (need for restart)
      CHARACTER*255 FILEDDIPOLE

      double precision hess(rank,rank)
      double precision ddipole(3*rank)
      integer iatom_start, ixyz_start
      integer i, j, iflag_grad0, dipole_there

      logical does_it_exist
c
c
c  The following file names need to be set by the user.
c  FILEATR is the .fdrst file to be created, FILEHESS is
c  the .hessian file to be read from, and FILEDDIPOLE is
c  the .fd_ddipole file to be read from.
c
      FILEATR     = "test.fdrst"
      FILEHESS    = "test.hessian"
      FILEDDIPOLE = "test.fd_ddipole"
c
c  The following variables need to be set by the user.
c  iatom_start is the atom to start with, ixyz_start is
c  1, 2, or 3 for x, y, or z.
c
      iatom_start = 45
      ixyz_start  = 3
c
c  The user should not have to modify anything below here.
c

c
c Open the files
c
      does_it_exist = .false.
      inquire(file=FILEATR,exist=does_it_exist)
      if (does_it_exist) then
        write(6,*) 'There is already a .fdrst file.'
        stop
      else
        open(unit=69,file=FILEATR,
     &      form='unformatted',
     &      access='sequential',
     &      status='new')
      endif
c
      does_it_exist = .false.
      inquire(file=FILEHESS,exist=does_it_exist)
      if (does_it_exist) then
        open (unit=46, file=FILEHESS, form='formatted',
     &      access='sequential', status='unknown')
      else
        write(6,*) 'The hessian file does not exist.'
        stop
      endif
c
      does_it_exist = .false.
      inquire(file=FILEDDIPOLE,exist=does_it_exist)
      if (does_it_exist) then
        open (unit=67, file=FILEDDIPOLE, form='formatted',
     &      access='sequential', status='unknown')
      else
        write(6,*) 'The fd_ddipole file does not exist.'
        stop
      endif
c
c  Read hessian information
c
      rewind (unit=46)
      do i = 1, rank
        do j = 1, i
          read(46,*,end=6,err=6) hess(i,j)
          hess(j,i) = hess(i,j)
        enddo
      enddo
      close (unit=46,status='keep')
      write(6,*)
     &    '  Nuclear Hessian retrieved from ASCII file.'
c
c  Read dipole derivative information
c
      rewind (unit=67)
      read(67,*,end=6,err=6)(ddipole(i),i=1,3*rank)
      close (unit=67,status='keep')
      write(6,*)
     &   '  Dipole derivative retrieved from ASCII file.'
c
c  Write information to the fdrst file
c
      iflag_grad0 = 0
      dipole_there = 1
      write(69)iatom_start,ixyz_start,rank,iflag_grad0,dipole_there
      write(69)hess
      write(69)ddipole
      close(unit=69,status='keep')
      goto 7
c
    6 continue
      write(6,*) 'There is a problem reading a file'
    7 continue
c
      end

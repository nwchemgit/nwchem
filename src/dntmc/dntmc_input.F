      subroutine dntmc_input(rtdb)
C     $Id$
C     version 5.0  Dec. 20, 2007 full release version (lcrosby)
C     modified to include input of number of bins (nob) June 10, 2008 (lcrosby)

      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
c
      integer rtdb  ! [input]
c
      integer mprnt, procrestart
      integer node, nob
      character*255 test
      integer nspecies, mcsteps
      double precision temp, rmin, rmax
      double precision disp, alpha
      integer i, l_nmol, k_nmol, k_spec, l_spec
      double precision conv
      character*256 directory
C
C  Some defaults
C
      nspecies = 1
      temp = 273.15
      mcsteps = 100
      rmin = 1.0
      rmax = 10.0
C
C  While loop for reading
C
 10   if (.not. inp_read()) call errquit('dntmc_input: inp_read failed',
     &        0, INPUT_ERR)
      if (.not. inp_a(test))
     &     call errquit('dntmc_input: failed to read keyword', 0,
     &                   INPUT_ERR)
c
      if (inp_compare(.false.,'nspecies', test)) then
c
c  Number of species - nspecies <integer nspecies>
c
        if (.not. inp_i(nspecies)) call errquit
     &       ('dntmc_input: nspecies <nspecies>', 0, INPUT_ERR)
        if (.not. rtdb_put(rtdb, 'dntmc:nspecies', MT_INT, 1, 
     &       nspecies)) call errquit
     &     ('dntmc_input: rtdb put failed for nspecies', 0, RTDB_ERR)
c
      else if (inp_compare(.false., 'nob', test)) then
c
C Number of bins in radial distributions - nob <intger nob>
c
        if (.not. inp_i(nob)) call errquit
     &       ('dntmc_input: nob <nob>', 0, INPUT_ERR)
        if (.not. rtdb_put(rtdb, 'dntmc:nob', MT_INT, 1,
     &       nob)) call errquit
     &     ('dntmc_input: rtdb put failed for nob', 0, RTDB_ERR)
c
      else if (inp_compare(.false., 'species', test)) then
c
c Species geometry names - species <geom1, geom2,..., geomnspecies>
c
         if (.not.ma_push_get(MT_BYTE, nspecies*255, 'dntmc: species', 
     &        l_spec, k_spec)) call errquit
     &        ('dntmc_input: cannot allocate species', 0, MA_ERR)
         call dntmc_inputspecies(rtdb, nspecies, byte_mb(k_spec))
         if (.not.ma_pop_stack(l_spec)) call errquit
     &       ('dntmc_input: ma pop failed for species', 1, MA_ERR)
c
      else if (inp_compare(.false., 'nmol', test)) then
c
c  Number of molecules for each species - nmol <nmol1, nmol2, ..., nmolnspecies>
c
c        write(6,*) 'nspecies is ', nspecies
         if (.not.ma_push_get(MT_INT, nspecies, 'dntmc: nmol', l_nmol, 
     &        k_nmol)) call errquit
     &        ('dntmc_input: cannot allocate nmol', 0, MA_ERR)
         call dntmc_getnmol(rtdb, nspecies, int_mb(k_nmol))
         if (.not.ma_pop_stack(l_nmol)) call errquit
     &       ('dntmc_input: ma pop failed for nmol', 1, MA_ERR)
c
      else if (inp_compare(.false., 'temp', test)) then
c
c  Temperature for running simulation - temp <integer temp>
c
        if (.not. inp_f(temp)) call errquit
     &       ('dntmc_input: temp <temp>', 0, INPUT_ERR)
        if (.not. rtdb_put(rtdb, 'dntmc:temp', MT_DBL, 1,
     &       temp)) call errquit
     &     ('dntmc_input: rtdb put failed for temp', 0, RTDB_ERR)
c
      else if (inp_compare(.false., 'mcsteps', test)) then
c
c  Number of MC steps - mcsteps <integer mcsteps>
c
        if (.not. inp_i(mcsteps)) call errquit
     &       ('dntmc_input: mcsteps <mcsteps>', 0, INPUT_ERR)
        if (.not. rtdb_put(rtdb, 'dntmc:mcsteps', MT_INT, 1,
     &       mcsteps)) call errquit
     &     ('dntmc_input: rtdb put failed for mcsteps', 0, RTDB_ERR)
c
      else if (inp_compare(.false., 'rmin', test)) then
c
c  Minumum radius - rmin <double precision rmin>
c
        if (.not. inp_f(rmin)) call errquit
     &       ('dntmc_input: rmin <rmin>', 0, INPUT_ERR)
        if (.not. rtdb_put(rtdb, 'dntmc:rmin', MT_DBL, 1,
     &       rmin)) call errquit
     &     ('dntmc_input: rtdb put failed for rmin', 0, RTDB_ERR)
c
      else if (inp_compare(.false., 'rmax', test)) then
c
c  Maximum radius - rmax <double precision rmax>
c
        if (.not. inp_f(rmax)) call errquit
     &       ('dntmc_input: rmax <rmax>', 0, INPUT_ERR)
        if (.not. rtdb_put(rtdb, 'dntmc:rmax', MT_DBL, 1,
     &       rmax)) call errquit
     &     ('dntmc_input: rtdb put failed for rmax', 0, RTDB_ERR) 
c
      else if (inp_compare(.false.,'end', test)) then
c
c  No more keywords to check against
c
      else if (inp_compare(.false.,'norestart', test)) then
c
c  Turn off the production of restart information
c
        if (.not. rtdb_put(rtdb, 'dntmc:norestart', MT_INT, 1,
     &     0)) call errquit
     &     ('dntmc_input: rtdb put failed for norestart', 0,
     &     RTDB_ERR)
c
      else if (inp_compare(.false.,'mprnt', test)) then
c
c  MPRNT int <value> (frequency of restart write)
c
        if (.not. inp_i(mprnt)) call errquit
     &       ('dntmc_input: mprnt <mprnt>', 0, INPUT_ERR)
        if (.not. rtdb_put(rtdb, 'dntmc:mprnt', MT_INT, 1,
     &       mprnt)) call errquit
     &     ('dntmc_input: rtdb put failed for mprnt', 0, RTDB_ERR)

      else if (inp_compare(.false.,'procrestart', test)) then
c
c  procrestart int <ngroups> to build restart rtdb
c

       if (.not. inp_i(procrestart)) call errquit
     &    ('dntmc_input: procrestart <ngroups>', 0, INPUT_ERR)
       if (.not. rtdb_put(rtdb, 'dntmc:procrestart', MT_INT, 1,
     &     procrestart)) call errquit
     &     ('dntmc_input: rtdb put failed for procrestart', 0,
     &     RTDB_ERR)

      else if (inp_compare(.false.,'tdisp', test)) then
c
c  tdisp double <dist> to define translational displ.
c

       if (.not. inp_f(disp)) call errquit
     &    ('dntmc_input: tdist <displacement>', 0, INPUT_ERR)
       if (.not. rtdb_put(rtdb, 'dntmc:disp', MT_DBL, 1,
     &     disp)) call errquit
     &     ('dntmc_input: rtdb put failed for tdisp', 0,
     &     RTDB_ERR)

      else if (inp_compare(.false.,'rdisp', test)) then
c
c  rdisp double <dist> to define rotational displ.
c

       if (.not. inp_f(alpha)) call errquit
     &    ('dntmc_input: rdist <displacement>', 0, INPUT_ERR)
       if (.not. rtdb_put(rtdb, 'dntmc:alpha', MT_DBL, 1,
     &     alpha)) call errquit
     &     ('dntmc_input: rtdb put failed for rdisp', 0,
     &     RTDB_ERR)

      else if (inp_compare(.false.,'rconfig', test)) then
c
c  rweight keywork to allow rsim weighting
c
      if (.not. rtdb_put(rtdb, 'dntmc:rconf', MT_LOG, 1, .true.))
     &  call errquit('dntmc_input: rtdb put failed for rconf',0,
     &  RTDB_ERR)

      else if (inp_compare(.false.,'rsim',test)) then
c
c  rweight keywork to allow rsim weighting
c
      if (.not. rtdb_put(rtdb, 'dntmc:rconf', MT_LOG, 1, .false.))
     &  call errquit('dntmc_input: rtdb put failed for rconf',0,
     &  RTDB_ERR)
      

      else if (inp_compare(.false.,'print', test)) then
C
C  Setup print control
C
      call util_print_input(rtdb,'dntmc')

      else if (inp_compare(.false.,'noprint', test)) then
C
C  Setup print control
C
      call util_print_input(rtdb,'dntmc')


      else if (inp_compare(.false.,'convergence',test)) then
C
C  Define convergence criterion
C
       if (.not. inp_f(conv)) call errquit
     &    ('dntmc_input: convergence <limit>', 0, INPUT_ERR)
       if (.not. rtdb_put(rtdb, 'dntmc:convergence', MT_DBL, 1,
     &     conv)) call errquit
     &     ('dntmc_input: rtdb put failed for convergence', 0,
     &     RTDB_ERR)

      else if (inp_compare(.false.,'dntmc_dir',test)) then
C
C  Define directory for files
C
      if (.not. inp_a(directory)) call errquit
     &   ('dntmc_input: error reading dntmc_directory',0, INPUT_ERR)
      if (.not. rtdb_cput(rtdb,'dntmc:directory', 1, directory)) 
     &   call errquit
     &   ('dntmc_input: rtdb put failed for dntmc_dir',0,RTDB_ERR)

      else
        call errquit('dntmc_input: unrecognized directive', 0, 
     &                INPUT_ERR)
      endif
c
c     End of while loop over the input
c
      if (.not. inp_compare(.false.,'end', test)) goto 10

      return
      end
c
C Little subroutine to help with the character input
c
      subroutine dntmc_inputspecies(rtdb, nspecies, species)
c
      implicit none
c
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb   ! [input]
      integer nspecies ! [input]
      character*255 species(nspecies) ! [input/output] list of geometry names
c
      integer i
c
      do i = 1, nspecies
        if (.not. inp_a(species(i))) call errquit
     &     ('dntmc_inputspecies: species line incorrect', 0, 
     &       INPUT_ERR)
      enddo

      if (.not.rtdb_cput(rtdb, 'dntmc:species', nspecies, species))
     &   call errquit('dntmc_inputspecies: rtdb failure for species', 
     &   0, RTDB_ERR)

      return
      end
c
c Little subroutine to help with the integer input
c
      subroutine dntmc_getnmol(rtdb, nspecies, nmol)
c
      implicit none
c
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb   ! [input]
      integer nspecies ! [input]
      integer nmol(nspecies) ! [input/output]
c
      integer i
c
         do i = 1, nspecies
            if (.not. inp_i(nmol(i))) call errquit
     &         ('dntmc_getnmol: nmol <nmol1, nmol2,..,nmolnspecies>', 0,
     &         INPUT_ERR)
         enddo
         if (.not. rtdb_put(rtdb, 'dntmc:nmol', MT_INT, nspecies,
     &       nmol)) call errquit
     &       ('dntmc_getnmol: rtdb put failed for nmol', 0, RTDB_ERR)

      return
      end

      subroutine geom_input(rtdb)
      implicit none
#include "inp.fh"
#include "geom.fh"
#include "tcgmsg.fh"
      integer rtdb              ! [input]
c
      character*255 field       ! for character input
      character*255 name        ! for name of geometry
      character*12 units        ! holds units of coordinates
      integer ncenter           ! counts no. of centers as input
      integer geom              ! handle for geometry
      logical status            ! scratch for return codes
      logical oprint
      integer max_center        ! parameter for local array dimension
      parameter (max_center = 1000)
      double precision coords(3,max_center)
      double precision charge(max_center)
      character*16 tags(max_center)
c
c     read a geometry from the input deck and output it
c     to the rtdb.
c
c     current input line should begin 'geometry ...'
c
c     Cartesians only for now
c
      if (nodeid() .ne. 0) return
c
c     Check that this is indeed a geometry line
c
      call inp_set_field(0)     ! goto start of line
      if (.not. inp_a(field))
     $     call errquit('geom_input: no input present', 0)
      if (.not. inp_compare(.false., 'geom', field))
     $     call errquit('geom_input: not geometry input', 0)
c     
c     geometry [<name>] [units <units>] [print] [noprint]
c     
      units = 'atomic units'
      name  = ' '
      oprint = .true.
c
 10   if (inp_a(field)) then
         if (inp_compare(.false.,'print', field)) then
            oprint = .true.
         else if (inp_compare(.false.,'noprint',field)) then
            oprint = .false.
         else if (inp_compare(.false.,'units', field)) then
            if (.not. inp_a(units)) call errquit
     $          ('geom_input: geometry [<name>] [units <units>]', 0)
            if (units(1:2).ne.'au' .and. units(1:2).ne.'bo')
     $           call errquit('geom_input: au only!', 0)
         else
            if (name .ne. ' ') call errquit
     $          ('geom_input: geometry [<name>] [units <units>]', 0)
            name = field
         endif
         goto 10
      endif
c
      if (.not. geom_create(geom, name)) call errquit
     $     ('geom_input: geom_create failed !', 0)
c     
c     tag charge x y z
c
      ncenter = 0
 20   if (inp_read()) then
         status = inp_a(field)
         if (inp_compare(.false., 'end', field)) then
            goto 30
         else
            if ((ncenter+1) .eq. max_center) call errquit
     $           ('geom_input: too many centers?', ncenter)
            tags(ncenter+1) = field
            status = status .and. inp_f(charge(ncenter+1))
            status = status .and. inp_f(coords(1,ncenter+1))
            status = status .and. inp_f(coords(2,ncenter+1))
            status = status .and. inp_f(coords(3,ncenter+1))
            if (.not. status) call errquit
     $           ('geom_input: <tag> <charge> <x> <y> <z>', 0)
            ncenter = ncenter + 1
         endif
         goto 20
      else
         call errquit('geom_input: premature end of file', 0)
      endif
c
 30   if (.not. geom_cart_set(geom, ncenter, tags, coords, charge))
     $     call errquit('geom_input: geom_cart_set failed', 0)
c
      if (oprint) then
         if (.not. geom_print(geom))
     $        call errquit('geom_input: print failed ', 0)
      endif
c
      if (.not. geom_rtdb_store(rtdb, geom, name))
     $     call errquit('geom_input: geom_rtdb_store failed', 0)
c
      if (.not. geom_destroy(geom))
     $     call errquit('geom_input: geom_destroy failed', 0)
c
c     done
c
      end

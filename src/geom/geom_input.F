      subroutine geom_input(rtdb)
      implicit none
#include "inp.fh"
#include "geom.fh"
#include "tcgmsg.fh"
#include "rtdb.fh"
      integer rtdb              ! [input]
c
      character*255 field       ! for character input
      character*255 name        ! for name of geometry
      character*12 units        ! holds units of coordinates
      integer ncenter           ! counts no. of centers as input
      integer geom              ! handle for geometry
      integer lpr_sym           ! print sym operators or not (1 yes,2 no)
      integer nata              ! num. of unique atoms 
      logical status            ! scratch for return codes
      logical oprint
      integer max_center        ! parameter for local array dimension
      double precision scale    ! For unit conversion
      parameter (max_center = 500)
      double precision coords(3,max_center)
      double precision charge(max_center)
      character*16 tags(max_center)
      character*2 symbol
      character*16 element
      logical osymmetry
      integer atn
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
c     subdirectives
c
c     system molecule||polymer||surface||crystal
c       angles
c       lengths
c
c     symmetry group_number setting_number
c
c     
      units = 'atomic units'
      name  = ' '
      oprint = .false.
      osymmetry = .true.
c
 10   if (inp_a(field)) then
         if (inp_compare(.false.,'print', field)) then
            oprint = .true.
         else if (inp_compare(.false.,'noprint',field)) then
            oprint = .false.
         else if (inp_compare(.false.,'units', field)) then
            if (.not. inp_a(units)) call errquit
     $          ('geom_input: geometry [<name>] [units <units>]', 0)

*            if (units(1:2).ne.'au' .and. units(1:2).ne.'bo')
*     $           call errquit('geom_input: au only!', 0)

         else
            if (name .ne. ' ') call errquit
     $          ('geom_input: geometry [<name>] [units <units>]', 0)
            name = field
         endif
         goto 10
      endif
c
      if (units(1:2) .eq. 'au' .or. units(1:2) .eq. 'bo' .or.
     $     units(1:2) .eq. 'at') then
         scale = 1.0d0
      else if (units(1:2) .eq. 'an') then
         scale = 1.8897265d0
      else
         call errquit('geom_input: unknown units', 0)
      endif

      if (name .eq. ' ') name = 'geometry'
      if (.not. geom_create(geom, name)) call errquit
     $     ('geom_input: geom_create failed !', 0)

      if (scale.ne.1.0d0) write(6,1) name(1:inp_strlen(name)), scale
 1    format(/'Scaling coordinates for geometry "',a,'" by ',f10.8/)
c     
c     system ...
c     symmetry ...
c     tag charge x y z
c
      ncenter = 0
 20   if (inp_read()) then
         status = inp_a(field)
         if (inp_compare(.false., 'end', field)) then
            goto 30
         else if (inp_compare(.false.,'symmetry',field)) then
            osymmetry = .true.
            call geom_sym_input(geom)
            goto 20
         else if (inp_compare(.false.,'system', field)) then
            osymmetry = .true.
            call geom_sys_input(geom)
            goto 20
         else
            if ((ncenter+1) .eq. max_center) call errquit
     $           ('geom_input: too many centers?', ncenter)
            tags(ncenter+1) = field
            status = status .and. inp_f(charge(ncenter+1))
            status = status .and. inp_f(coords(1,ncenter+1))
            status = status .and. inp_f(coords(2,ncenter+1))
            status = status .and. inp_f(coords(3,ncenter+1))

            coords(1,ncenter+1) = coords(1,ncenter+1)*scale
            coords(2,ncenter+1) = coords(2,ncenter+1)*scale
            coords(3,ncenter+1) = coords(3,ncenter+1)*scale

            ncenter = ncenter + 1

            if (.not. status) call errquit
     $           ('geom_input: <tag> <charge> <x> <y> <z>', ncenter)
c
c     Check that charge of center is consistent with atom type
c
            if (geom_tag_to_element(tags(ncenter), symbol,
     $           element, atn)) then
               if (abs(charge(ncenter) - dble(atn)).gt.1e-8) then
                  write(6,2) tags(ncenter), symbol, element, atn,
     $                 charge(ncenter)
 2                format('geom_input: ', a, '(',a,',',a,',',i3,') ',
     $                 f6.2)
                  call util_flush(6)
                  call errquit
     $                 ('geom_input: charge-atom type mismatch',atn)
               endif
            endif
         endif
         goto 20
      else
         call errquit('geom_input: premature end of file', 0)
      endif
c
 30   if (.not. geom_cart_set(geom, ncenter, tags, coords, charge))
     $     call errquit('geom_input: geom_cart_set failed', 0)
c
***      if (osymmetry) then
c
c     Apply system and symmetry info to the list of
c     unique centers build mapping tables set up coord lists
c
        if(oprint) then
          lpr_sym=1
        else
          lpr_sym=0
        endif
        nata=ncenter
        call sym_nwc(geom,rtdb,nata,lpr_sym)
**      else 
**         if (.not. geom_charge_center(geom)) call errquit
**     $     ('geom_input: geom_charge_center failed', 0)
**      endif
c
      if (oprint) then
         if (.not. geom_print(geom))
     $        call errquit('geom_input: print failed ', 0)
      endif
c
      if (.not. geom_rtdb_store(rtdb, geom, name))
     $     call errquit('geom_input: geom_rtdb_store failed', 0)
c
      if(.not. geom_verify_coords(geom))
     &    call errquit
     &    ('geom_input: geom_verify_coords failed',911)
c
      if (.not. geom_destroy(geom))
     $     call errquit('geom_input: geom_destroy failed', 0)
c
c     done
c

c      if (.not. rtdb_print(rtdb, .true.)) call errquit('print failed',0)

      end
      subroutine geom_sys_input(geom)
      implicit none
#include "geom.fh"
#include "geomP.fh"
#include "inp.fh"
      integer geom
c
      integer nopt_system, ntokens, ind
      double precision value
      parameter (nopt_system=4, ntokens=7)
      character*10 opt(nopt_system), tokens(ntokens)
      data opt /'molecule', 'polymer', 'surface', 'crystal'/
      data tokens /'lat_a', 'lat_b', 'lat_c', 
     &     'alpha', 'beta', 'gamma', 'end'/
      character*10 test
c
      if (.not. geom_check_handle(geom,'system_input'))
     $     call errquit('system_input: geom handle invalid', 0)
c
      if (inp_a(test)) then
        if (inp_match(nopt_system, .false., test, opt, ind)) then
          isystype(geom) = ind - 1
        else
          goto 1000
        endif
      else
        goto 1000
      endif
c
c     Read new line of input
c
   20 if (.not. inp_read()) call errquit
     $     ('geom_sys_input: premature end of file', 0)
c
   10 if (.not. inp_a(test)) then
        goto 20
      else if (inp_match(ntokens,.false.,test,tokens,ind)) then
        if (ind .lt. 7) then
          if (.not. inp_f(value)) call errquit
     $         ('system:input: error reading floating number', 0)
c         
          goto (100,200,300,400,500,600) ind
c
  100     lattice_vectors(1,geom) = value
          goto 10
  200     lattice_vectors(2,geom) = value
          goto 10
  300     lattice_vectors(3,geom) = value
          goto 10
  400     lattice_angles(1,geom) = value
          goto 10
  500     lattice_angles(2,geom) = value
          goto 10
  600     lattice_angles(3,geom) = value
          goto 10
        else if(ind.eq.7) then
          goto 30
        endif
      else
        goto 1000
      endif
c
 1000 call errquit('system_input: invalid/missing system', 0)
c
c--> extra geometry stuff associated with periodic lattices
c    derived from the above information. Builds reciprocal lattice
c    vectors, a-matrix (used for fractional-> Cartesian trasforms), 
c    volume of direct space lattce and the g-matrix (metric matrix used
c    for vector algebra in oblique coordinate systems.
c
c
c--> 2d and 1d systems not active yet
c
   30 if(isystype(geom).eq.3) then
        call geom_3d(geom)
c      elseif (isystype(geom).eq.2) then
c        call geom_2d(geom)
c      elseif (isystype(geom).eq.1) then
c        call geom_1d(geom)
      endif
      return
      end
      subroutine geom_sym_input(geom)
      implicit none
#include "inp.fh"
#include "geom.fh"
#include "geomP.fh"
c
      integer geom
      character*10 test
c
      if (.not. geom_check_handle(geom,'system_input'))
     $     call errquit('system_input: geom handle invalid', 0)
c
      group_number(geom) = 1
      setting_number(geom) = 1
c
   10 if (inp_a(test)) then
        if (inp_compare(.false.,'setting',test)) then
          if (.not. inp_i(setting_number(geom))) call errquit
     $         ('symmetry_input: bad setting value', 0)
        else if (inp_compare(.false.,'group',test)) then
          if (.not. inp_i(group_number(geom))) call errquit
     $         ('symmetry_input: bad setting value', 0)
        else
          call errquit('symmetry_input: invalid option', 0)
        endif
      endif
c
c      write(6,*) ' group/setting ', group_number(geom), 
c     $     setting_number(geom)
c
      end





      subroutine geom_input(rtdb)
C     $Id: geom_input.F,v 1.26 1996-07-25 18:23:39 pg481 Exp $
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
      integer nata              ! num. of unique atoms 
      logical status            ! scratch for return codes
      logical oprint
      logical include_bqbq
      logical oprint_sym
      integer max_center        ! parameter for local array dimension
      double precision scale    ! For unit conversion
      parameter (max_center = 500)
      double precision coords(3,max_center)
      double precision charge(max_center), mass(max_center)
      character*16 tags(max_center)
      character*2 symbol
      character*16 element
      logical osymmetry, oecp(max_center)
      logical old_format_input, is_atom
      integer atn, isys
      character*16 format
      external geom_data
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
      if (.not. inp_compare(.false., 'geometry', field))
     $     call errquit('geom_input: not geometry input', 0)
c     
c     geometry [<name>] [units <units>] [print] [noprint] [bqbq]
c     
c     subdirectives
c     
c     system molecule||polymer||surface||crystal
c     angles
c     lengths
c     
c     symmetry group_number setting_number
c     
c     systype must precede symmetry
c     
c     
      old_format_input = .false.
      units = 'atomic units'
      name  = ' '
      oprint = .true.           ! Default is to print the geometry
      osymmetry = .true.
      oprint_sym = .false.
      include_bqbq = .false.
c     
 10   if (inp_a(field)) then
         if (inp_compare(.false.,'bqbq',field)) then
            include_bqbq = .true.
         else if (inp_compare(.false.,'print', field)) then
            oprint = .true.
c     
c     If the next field is a recognizable format for standard
c     geometries then process as
c     
c     print format
c     
c     currently just know about xyz
c     
            if (inp_a(format)) then
               if (inp_compare(.false.,'xyz',format)) then
                  format = 'xyz'
               else
                  call inp_prev_field()
                  format = ' '
               endif
            else
               format = ' '
            endif
         else if (inp_compare(.false.,'noprint',field)) then
            oprint = .false.
         else if (inp_compare(.false.,'units', field)) then
            if (.not. inp_a(units)) call errquit
     $           ('geom_input: geometry [<name>] [units <units>]', 0)

*     if (units(1:2).ne.'au' .and. units(1:2).ne.'bo')
*     $           call errquit('geom_input: au only!', 0)

         else
            if (name .ne. ' ') call errquit
     $           ('geom_input: geometry [<name>] [units <units>]', 0)
            name = field
         end if
         goto 10
      end if
c     
      if (name .eq. ' ') name = 'geometry'
      if (.not. geom_create(geom, name)) call errquit
     $     ('geom_input: geom_create failed !', 0)
      if (include_bqbq) status = geom_set_bqbq(geom, .true.)
c     
      if (units(1:2) .eq. 'au' .or. units(1:2) .eq. 'bo' .or.
     $     units(1:2) .eq. 'at') then
         scale = 1.0d0
         status = geom_set_user_units(geom,'a.u.')
      else if (units(1:2) .eq. 'an') then
         scale  = 1.8897265d0
         status = geom_set_user_units(geom,'angstrom')
      else
         call errquit('geom_input: unknown units', 0)
      end if
c     
c     system ...
c     symmetry ...
c     tag charge x y z  !! old format
c     tag x y z [ecp] [charge q] [mass m]
c     
      ncenter = 0
 20   if (inp_read()) then
         status = inp_a(field)
         if (inp_compare(.false., 'end', field)) then
            goto 30
         else if (inp_compare(.false.,'symmetry',field)) then
            osymmetry = .true.
            call geom_sym_input(geom, oprint_sym)
            goto 20
         else if (inp_compare(.false.,'system', field)) then
            osymmetry = .true.
            call geom_sys_input(geom)
            goto 20
         else
            if ((ncenter+1) .eq. max_center) call errquit
     $           ('geom_input: too many centers?', ncenter)
            tags(ncenter+1) = field
            oecp(ncenter+1) = .false.
            charge(ncenter+1) = 0.0d0
c     
            is_atom = geom_tag_to_element(tags(ncenter+1), symbol,
     $           element, atn)
            if ((.not. is_atom) .and. symbol.ne.'bq') 
     $        call errquit('geom_input:center is neither atom nor bq',0)
c
c.. set default mass
c
            if (.not.
     &            geom_atn_to_default_mass(atn,mass((ncenter+1))))
     &            call errquit(' geom_input: default mass failed',
     &            911)
c     
c     Old style input  tag charge x y z
c     
            status = status .and. inp_f(charge(ncenter+1))
            status = status .and. inp_f(coords(1,ncenter+1))
            status = status .and. inp_f(coords(2,ncenter+1))
            status = status .and. inp_f(coords(3,ncenter+1))
            if (status) then
               old_format_input = .true.
               if (is_atom) then
                  if (abs(charge(ncenter+1)-dble(atn)).gt.1e-8) then
                     write(6,2) tags(ncenter+1), symbol, element, atn,
     $                 charge(ncenter+1)
 2                   format('geom_input: ', a, '(',a,',',a,',',i3,') ',
     $                    f6.2)
                     call util_flush(6)
                     call errquit
     $                    ('geom_input: charge-atom type mismatch',atn)
                  end if
               end if
            else
c     
c     New style input ... <tag> <x> <y> <z> [ecp] [charge <q>] [mass <m>]
c     
               call inp_set_field(1)
               status = .true.
               status = status .and. inp_f(coords(1,ncenter+1))
               status = status .and. inp_f(coords(2,ncenter+1))
               status = status .and. inp_f(coords(3,ncenter+1))
               charge(ncenter+1) = atn
c     
 111           if (inp_a(field)) then
                  if (inp_compare(.false.,'ecp', field)) then
                     oecp(ncenter+1) = .true.
c
c     Do not use default ECP charge if user has already overridden it
c
		     if (charge(ncenter+1) .eq. 0.0d0) then
                       if (.not. geom_default_charge_with_ecp
     $                    (atn, charge(ncenter+1))) call errquit
     $                    ('geom_input: no default charge with ecp?',0)
		     endif
                     goto 111
                  else if (inp_compare(.false., 'mass', field)) then
                     status = status .and. inp_f(mass(ncenter+1))
                     goto 111
                  else if (inp_compare(.false., 'charge', field)) then
                     status = status .and. inp_f(charge(ncenter+1))
                     goto 111
                  else
                     status = .false.
                  endif
               endif
               if (.not. status) call errquit
     $              ('geom_input: <tag> <x> <y> <z> [ecp]'//
     $              ' [charge <q>] [mass <m>]', ncenter)
            endif
            ncenter = ncenter + 1
c     
         end if
         goto 20
      else
         call errquit('geom_input: premature end of file', 0)
      end if
c     
c     Apply requested conversion of units after we have
c     all of the information
c     
 30   if (.not. geom_systype_get(geom, isys)) call errquit
     $     ('geom_input: systype?', 0)
      if (isys.eq.0) then
         if (scale.ne.1.0d0) then
            write(6,1) name(1:inp_strlen(name)), scale
 1          format(/'Scaling coordinates for geometry "',a,'" by ',
     $           f10.8/)
            call dscal(3*ncenter, scale, coords, 1)
         end if
      else if (isys .gt. 0) then
         if (scale .ne. 1.0d0) call errquit
     $        ('geom_input: specified units for fractional coords', 0)
      else
         call errquit('geom_input: invalid system type', 0)
      end if
c     
      if (.not. geom_cart_set(geom, ncenter, tags, coords, charge))
     $     call errquit('geom_input: geom_cart_set failed', 0)
c
      if (.not. geom_masses_set(geom, ncenter, mass))
     &      call errquit('geom_input: geom_masses_set failed', 0)
c     
c     Apply system and symmetry info to the list of
c     unique centers build mapping tables set up coord lists
c     
      nata=ncenter
      call sym_nwc(geom,rtdb,nata,oprint_sym)
c     
c     Force exact symetry on the coordinates
c     
      call sym_geom_project(geom, 1d-3)
c     
      call geom_compute_values(geom)
c     
      if (old_format_input) write(6,77)
 77   format(' !! You have used an archaic geometry format. Use',
     $     ' instead:'/
     $     ' !!     <tag> <x> <y> <z> [ecp] [charge <q>] [mass <m>]'/)
      if (oprint) then
         if (.not. geom_print(geom))
     $        call errquit('geom_input: print failed ', 0)
         if (format .eq. 'xyz') then
            write(6,*)
            call util_print_centered(6,'XYZ format geometry',20,.true.)
            if (.not. geom_print_xyz(geom,6))
     $           call errquit('geom_input: print xyz?',0)
            write(6,*)
            call util_flush(6)
         endif
*     call sym_print_all(geom, .true., .true., .true., .true.,.true.)
      end if
c     
      if (.not. geom_rtdb_store(rtdb, geom, name))
     $     call errquit('geom_input: geom_rtdb_store failed', 0)
c     
      if(.not. geom_verify_coords(geom))
     &     call errquit
     &     ('geom_input: geom_verify_coords failed',911)
c     
      if (.not. geom_destroy(geom))
     $     call errquit('geom_input: geom_destroy failed', 0)
c     
c     done
c     

c     if (.not. rtdb_print(rtdb, .true.)) call errquit('print failed',0)

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
      logical status
      parameter (nopt_system=4, ntokens=7)
      character*10 opt(nopt_system), tokens(ntokens)
      character*10 test
      data opt /'molecule', 'polymer', 'surface', 'crystal'/
      data tokens /'lat_a', 'lat_b', 'lat_c', 
     &     'alpha', 'beta', 'gamma', 'end'/
c
      if (.not. geom_check_handle(geom,'system_input'))
     $     call errquit('system_input: geom handle invalid', 0)
c
      if (inp_a(test)) then
        if (inp_match(nopt_system, .false., test, opt, ind)) then
          isystype(geom) = ind - 1
        else
          goto 1000
        end if
      else
        goto 1000
      end if
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
        end if
      else
        goto 1000
      end if
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
        status = geom_set_user_units(geom,'fractional')
      elseif (isystype(geom).eq.2) then
        call geom_2d(geom)
      elseif (isystype(geom).eq.1) then
        call geom_1d(geom)
      end if
      return
      end
      subroutine geom_sym_input(geom, oprint_sym)
      implicit none
#include "inp.fh"
#include "geom.fh"
#include "geomP.fh"
c     
      integer geom
      character*10 test
      logical oprint_sym
c     
c     symmetry [[group] (<string group_name>|<integer group_number>)] \
c              [setting <integer setting>]
c
      logical status
c
      if (.not. geom_check_handle(geom,'geom_sym_input'))
     $     call errquit('geom_sym_input: geom handle invalid', 0)
c     
      group_number(geom) = 1
      setting_number(geom) = 1
      oprint_sym = .false.
c     
 10   if (inp_a(test)) then
         if (inp_compare(.false.,'setting',test)) then
            if (.not. inp_i(setting_number(geom))) call errquit
     $           ('geom_sym_input: bad setting value', 0)
            goto 10
         else if (inp_compare(.false.,'print',test)) then
            oprint_sym = .true.
            goto 10
         else if (inp_compare(.false.,'group',test)) then
            goto 10             ! Just ignore
         else 
            call inp_prev_field()
            if (.not. inp_i(group_number(geom))) then
               if (inp_a(test)) then
                  status = .false.
                  if (isystype(geom).eq.0) then
                     status = inp_match(46,.false.,test,sym_molgnames,
     $                    group_number(geom))
                  else
                     status = inp_match(230,.false.,test,sym_spgnames,
     $                 group_number(geom))
                  endif
                  if (.not. status)
     $                 call errquit('geom_sym_input: unknown group '//
     $                 'name/number', 0)
               endif
            endif
            goto 10
         end if
      end if
c     
      end





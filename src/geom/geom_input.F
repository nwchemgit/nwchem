      subroutine geom_input(rtdb)
C     $Id: geom_input.F,v 1.94 1999-05-25 22:06:58 d3g681 Exp $
      implicit none
#include "stdio.fh"
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
      integer ncenter_before_auto
      integer geom              ! handle for geometry
      integer nata              ! num. of unique atoms 
      logical status            ! scratch for return codes
      logical oprint
      logical ozmatrix
      logical oautoz
      logical oautosym
      logical ocenter
      logical oforce_autoz
      logical include_bqbq
      logical oprint_sym
      integer nzinput
      integer max_center        ! parameter for local array dimension
      double precision scale    ! For unit conversion
      double precision threquiv ! Threshold for autosym atom equivalence
c     .                           and for forcing symmetry
      parameter (max_center = 500)
      double precision coords(3,max_center), velocities(3,max_center)
      double precision charge(max_center), mass(max_center)
      character*16 tags(max_center)
      character*2 symbol
      character*16 element, groupname
      logical osymmetry
      logical old_format_input, is_atom
      logical found_ang2au_conv
      logical found_au2ang_conv
      double precision ang2au, au2ang
      integer atn, isys
      character*16 format
      logical dbug
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
      dbug=.false.
c     
c     Check that this is indeed a geometry line
c     
      call inp_set_field(0)     ! goto start of line
      if (.not. inp_a(field))
     $     call errquit('geom_input: no input present', 0)
      if (.not. inp_compare(.false., 'geometry', field))
     $     call errquit('geom_input: not geometry input', 0)
c     
c     geometry [<name>] [units <units>] [print] [noprint] [bqbq] \
c     [angstrom_to_au <real>] [au_to_angstrom <real>] \
c     [(no)autoz]
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
      units = 'angstroms'        ! The new default
      name  = ' '
      oprint = .true.           ! Default is to print the geometry
      ozmatrix = .false.        ! Default is cartesian input
      oautoz   = .true.         ! Default is try -autoz-
      oautosym = .false.        ! Default is don't try -autosym-
      osymmetry = .true.
      oprint_sym = .false.
      include_bqbq = .false.
      ocenter = .true.
      oforce_autoz = .false.
      nzinput=0
      threquiv = 1d-4
c     
      found_ang2au_conv = .false.
      found_au2ang_conv = .false.
c
      call dfill(3*max_center, 0.0d0, velocities, 1)
c     
 10   if (inp_a(field)) then
         if (inp_compare(.false.,'bqbq',field)) then
            include_bqbq = .true.
         else if ((inp_compare(.false.,'angstrom_to_au',field)).or.
     &           (inp_compare(.false.,'ang2au',field))) then
            found_ang2au_conv = .true.
            if (.not.inp_f(ang2au)) call errquit
     &           ('error parsing real value of angstrom_to_au',911)
         else if ((inp_compare(.false.,'au_to_angstrom',field)).or.
     &           (inp_compare(.false.,'au2ang',field))) then
            found_au2ang_conv = .true.
            if (.not.inp_f(au2ang)) call errquit
     &           ('error parsing real value of angstrom_to_au',911)
         else if (inp_compare(.false.,'autoz',field)) then
            oautoz       = .true.  
            oforce_autoz = .true.
         else if (inp_compare(.false.,'noautoz',field)) then
            oautoz       = .false.
         else if (inp_compare(.false.,'autosym',field)) then
            oautosym     = .true.  
            status = inp_f(threquiv)
         else if (inp_compare(.false.,'noautosym',field)) then
            oautosym     = .false.
         else if (inp_compare(.false.,'center', field)) then
            ocenter = .true.
         else if (inp_compare(.false.,'nocenter', field)) then
            ocenter = .false.
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
            endif
         else if (inp_compare(.false.,'noprint',field)) then
            oprint = .false.
         else if (inp_compare(.false.,'units', field)) then
            if (.not. inp_a(units)) call errquit
     $           ('geom_input: geometry [<name>] [units <units>]', 0)
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
      if (found_au2ang_conv.and.found_ang2au_conv) call errquit
     &     ('geom_input: user specified both au->ang and ang->au '//
     &     'conversion factors.  This is not allowed', 911)
      if (found_au2ang_conv) then
         if (.not.geom_set_au2ang(geom,au2ang)) call errquit
     &        ('geom_input: error setting au2ang conversion unit',911)
      elseif (found_ang2au_conv) then
         if (.not.geom_set_ang2au(geom,ang2au)) call errquit
     &        ('geom_input: error setting ang2au conversion unit',911)
      endif
c     
c     ----- check units -----
c     
      call inp_lcase(units)
c     
      if (units(1:2) .eq. 'au' .or. units(1:2) .eq. 'bo' .or.
     &     units(1:2) .eq. 'at') then
         units = 'a.u.'
      else if (units(1:2) .eq. 'an') then
         units = 'angstroms'
      else if (units(1:2).eq.'nm' .or. units(1:2).eq.'na') then
         units = 'nanometer'
      else  if (units(1:2).eq.'pm' .or. units(1:2).eq.'pi') then
         units = 'picometer'
      else
         call errquit('geom_input: unknown units', 0)
      end if
c     
      if (.not. geom_set_user_units(geom,units))
     $     call errquit('geom_input: failed setting user units',0)
      if (.not. geom_get_user_scale(geom,scale))
     $     call errquit('geom_input: failed getting user scale',0)
c     
c     ----- coordinate input -----
c     
c     system ...
c     symmetry ...
c     tag charge x y z  !! old format rak 9/96 removed
c     tag x y z [charge q] [mass m]
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
            call geom_sys_input(geom,scale)
            goto 20
         else if (inp_compare(.false.,'zcoord', field)) then
            call geom_autoz_input(geom)
            goto 20
         else if (inp_compare(.false.,'zmt', field).or.
     $           inp_compare(.false.,'zmat', field).or.
     $           inp_compare(.false.,'zmatr', field).or.
     $           inp_compare(.false.,'zmatri', field).or.
     $           inp_compare(.false.,'zmatrix', field)) then
            ozmatrix=.true.   
            oautoz = .false.
            call geom_zmt_input(geom,coords,tags,charge,mass,ncenter,
     $           nzinput,scale,oprint)
            oautoz=.true.
            goto 20
         else
            if ((ncenter+1) .eq. max_center) call errquit
     $           ('geom_input: too many centers?', ncenter)
            tags(ncenter+1) = field
            charge(ncenter+1) = 0.0d0
c     
            is_atom = geom_tag_to_element(tags(ncenter+1), symbol,
     $           element, atn)
            if ((.not. is_atom) .and. symbol.ne.'bq') 
     $           call errquit
     &           ('geom_input:center is neither atom nor bq',0)
c     
c..   set default mass
c     
            if (.not.
     &           geom_atn_to_default_mass(atn,mass((ncenter+1))))
     &           call errquit(' geom_input: default mass failed',
     &           911)
c     
c     New style input ... <tag> <x> <y> <z> [vx vvy vz] \
c                         [charge <q>] [mass <m>]
c     
            call inp_set_field(1)
            status = .true.
            status = status .and. inp_f(coords(1,ncenter+1))
            status = status .and. inp_f(coords(2,ncenter+1))
            status = status .and. inp_f(coords(3,ncenter+1))
            if (inp_f(velocities(1,ncenter+1))) then
               status = status .and. inp_f(velocities(2,ncenter+1))
               status = status .and. inp_f(velocities(3,ncenter+1))
            endif
            charge(ncenter+1) = atn
c     
 111        if (inp_a(field)) then
               if (inp_compare(.false., 'mass', field)) then
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
     $           ('geom_input: <tag> <x> <y> <z> '//
     $           ' [charge <q>] [mass <m>]', ncenter)
            ncenter = ncenter + 1
c     
         end if
         goto 20
      else
         call errquit('geom_input: premature end of file', 0)
      end if
c     
c     For molecules, apply requested conversion of units after we have
c     all of the information
c     
 30   if (.not. geom_systype_get(geom, isys)) call errquit
     $     ('geom_input: systype?', 0)
      if (isys.eq.0) then
         if (scale.ne.1.0d0) then
            if (oprint)
     $           write(luout,1) name(1:inp_strlen(name)), scale, 
     $           1.0d0/scale
 1          format(/' Scaling coordinates for geometry "',a,'" by ',
     $           f12.9/' (inverse scale = ',f12.9,')',/)
            call dscal(3*ncenter, scale, coords, 1)
         end if
      else if (isys .lt. 0) then
         call errquit('geom_input: invalid system type', 0)
      end if
c     
      call geom_check_input_quants
     &    (tags,charge,ncenter,1.0d-12,'charge')
      call geom_check_input_quants
     &    (tags,mass,ncenter,1.0d-8,'mass')
c
      if (.not. geom_cart_set(geom, ncenter, tags, coords, charge))
     $     call errquit('geom_input: geom_cart_set failed', 0)
c
      if (.not. geom_vel_set(geom, velocities))
     $     call errquit('geom_input: geom_vel_set failed', 0)
c     
      if (.not. geom_masses_set(geom, ncenter, mass))
     &     call errquit('geom_input: geom_masses_set failed', 0)
c
      if(oautosym) then
c
c     ----- detect symmetry -----
c
         if (.not. geom_cart_get(geom, ncenter, tags, coords, charge))
     $        call errquit('geom_input: geom_cart_get failed', 0)
         ncenter_before_auto = ncenter
         call geom_auto_sym(rtdb,geom,coords,charge,tags,ncenter,
     $        threquiv,groupname)
         if (geom_group_set(geom,groupname)) then
            if (.not.geom_cart_set(geom, ncenter, tags, coords, charge))
     $           call errquit('geom_input: geom_cart_set failed', 0)
         else
            write(6,*) ' autosym detected uknown group ', groupname
            call errquit('geom_input: autosym: invalid group',0)
         endif
      endif
c
c     Apply system and symmetry info to the list of
c     unique centers build mapping tables set up coord lists
c     
      nata=ncenter
      call sym_nwc(geom,rtdb,nata,oprint_sym,scale)
c
c     Check that if we used autosym that we ended up with the 
c     same no. of atoms ... if we don't then autosym and nwchemsym
c     don't agree on the orientation of point group elements
c
      if (oautosym) then
         if (.not. geom_ncent(geom,ncenter)) call errquit
     $        ('geom_input: geom_cent?',0)
         if (ncenter_before_auto .ne. ncenter) call errquit
     $        ('geom_input: autosym bug : too many atoms',ncenter)
      endif
c
      if (oprint_sym.and.(isys.eq.0)) call sym_print_all(geom,
     $     .true.,.true.,.true.,.true.,.true.)
c     
c     Center origin at the center of charge if requested
c     
      if (ocenter .and. isys.eq.0) then
         if (.not. geom_charge_center(geom)) call errquit
     $        ('geom_input: failed centering geometry',0)
      endif
c     
c     Force exact symetry on the coordinates
c     
      call sym_geom_project(geom, threquiv)
c     
      if(oautoz) then
         if (.not. geom_cart_get(geom, ncenter, tags, coords, charge))
     $        call errquit('geom_input: geom_cart_get failed', 0)
c     
c     ----- call auto-z-matrix generation -----
c     
         call geom_zmt_autoz(rtdb,geom,coords,charge,tags,ncenter,
     $        nzinput,oforce_autoz,oprint)
      else
c     
c     ----- no auto-z -----
c     
         if (.not. geom_cart_get(geom, ncenter, tags, coords, charge))
     $        call errquit('geom_input: geom_cart_get failed', 0)
         call geom_no_autoz(geom,coords,charge,tags,ncenter,nzinput)
      endif
c
      if (oprint) then
         if (.not. geom_print(geom))
     $        call errquit('geom_input: print failed ', 0)
         if (format .eq. 'xyz') then
            write(luout,*)
            call util_print_centered
     &           (luout,'XYZ format geometry',20,.true.)
            if (.not. geom_print_xyz(geom,6))
     $           call errquit('geom_input: print xyz?',0)
            write(luout,*)
            call util_flush(luout)
         endif
         if(dbug) then
            call sym_print_all(geom,.true.,.true.,.true.,.true.,.true.)
         endif
         if (.not.geom_print_distances(geom))
     &         call errquit('geom_input: print_distances failed ',911)
         if (.not.geom_print_angles(geom))
     &         call errquit('geom_input: print_angles failed ',911)
      end if
c     
      if (.not. geom_rtdb_store(rtdb, geom, name))
     $     call errquit('geom_input: geom_rtdb_store failed', 0)
c     
      if(.not. geom_verify_coords(geom))
     &     call errquit
     &     ('geom_input: geom_verify_coords failed',911)
c
c
c
c     
      if (.not. geom_destroy(geom))
     $     call errquit('geom_input: geom_destroy failed', 0)
c     
c     done

c     if (.not. rtdb_print(rtdb, .true.)) call errquit('print failed',0)

      end
      subroutine geom_sys_input(geom,scale)
      implicit none
#include "geom.fh"
#include "geomP.fh"
#include "inp.fh"
      integer geom
c
      integer nopt_system, ntokens, ind
c
      double precision scale
c
      double precision value
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
        call geom_3d(geom,scale)
c        status = geom_set_user_units(geom,'fractional')
      elseif (isystype(geom).eq.2) then
        call geom_2d(geom,scale)
      elseif (isystype(geom).eq.1) then
        call geom_1d(geom,scale)
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
c     symemtry [[group] (<string group_name>|<integer group_number>)] \
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
      logical function geom_autoz_info_get(geom,
     $     xcvr_scaling, xijbond, xijkang, xijklto, xijklop, 
     $     xijklnb)
      implicit none
#include "geomP.fh"
      integer geom
      integer xijbond(2,max_zcoord) 
      integer xijkang(3,max_zcoord)
      integer xijklto(4,max_zcoord)
      integer xijklop(4,max_zcoord)
      integer xijklnb(4,max_zcoord)
      double precision xcvr_scaling
c
      logical geom_check_handle
      external geom_check_handle
      integer i
c
      geom_autoz_info_get = geom_check_handle(geom,
     $     'geom_autoz_info_get')
      if (.not. geom_autoz_info_get) return
c
***      write(6,*) '  cvr_scaling(geom) ', cvr_scaling(geom)
      if (cvr_scaling(geom) .eq. 0d0) then
         geom_autoz_info_get = .false.
         return
      endif
c
      xcvr_scaling = cvr_scaling(geom)
      do i = 1, max_zcoord
         xijbond(1,i) = ijbond(1,i,geom)
         xijbond(2,i) = ijbond(2,i,geom)
         xijkang(1,i) = ijkang(1,i,geom)
         xijkang(2,i) = ijkang(2,i,geom)
         xijkang(3,i) = ijkang(3,i,geom)
         xijklto(1,i) = ijklto(1,i,geom)
         xijklto(2,i) = ijklto(2,i,geom)
         xijklto(3,i) = ijklto(3,i,geom)
         xijklto(4,i) = ijklto(4,i,geom)
         xijklop(1,i) = ijklop(1,i,geom)
         xijklop(2,i) = ijklop(2,i,geom)
         xijklop(3,i) = ijklop(3,i,geom)
         xijklop(4,i) = ijklop(4,i,geom)
         xijklnb(1,i) = ijklnb(1,i,geom)
         xijklnb(2,i) = ijklnb(2,i,geom)
         xijklnb(3,i) = ijklnb(3,i,geom)
         xijklnb(4,i) = ijklnb(4,i,geom)
      enddo
c     
      end
      subroutine geom_autoz_input(geom)
      implicit double precision (a-h,o-z)
#include "geomP.fh"
      integer geom
c
      call hnd_autoz_input(
     $     ijbond(1,1,geom),
     $     ijkang(1,1,geom),
     $     ijklto(1,1,geom),
     1     ijklop(1,1,geom),
     $     ijklnb(1,1,geom),
     $     max_zcoord,
     2     cvr_scaling(geom))
c     
      return
      end
      subroutine hnd_autoz_input(ijbond,ijkang,ijklto,
     1                           ijklop,ijklnb,max_zcoord,
     2                           cvr_scaling)
      implicit none
#include "inp.fh"
c
      double precision cvr_scaling
      double precision cvr_factor
      parameter (cvr_factor=1.20d+00)
      character*30 field
      integer ijbond(2,*), max_zcoord
      integer ijkang(3,*)
      integer ijklto(4,*)
      integer ijklop(4,*)
      integer ijklnb(4,*)
c
      integer nbond, nangl, ntors, nopla, nlinb, i
c
      cvr_scaling=cvr_factor
      nbond=0
      nangl=0
      ntors=0
      nopla=0
      nlinb=0
      do i=1,max_zcoord
         ijbond(1,i)=0
         ijbond(2,i)=0
         ijkang(1,i)=0
         ijkang(2,i)=0
         ijkang(3,i)=0
         ijklto(1,i)=0
         ijklto(2,i)=0
         ijklto(3,i)=0
         ijklto(4,i)=0
         ijklop(1,i)=0
         ijklop(2,i)=0
         ijklop(3,i)=0
         ijklop(4,i)=0
         ijklnb(1,i)=0
         ijklnb(2,i)=0
         ijklnb(3,i)=0
         ijklnb(4,i)=0
      enddo
c
c     ----- read values from input -----
c
 10   if(.not.inp_read()) 
     $     call errquit('hnd_autoz_input: unexpected eof',911)

      if (.not.inp_a(field)) 
     $     call errquit('hnd_autoz_input: read of field failed',911)
c
      if (inp_compare(.false.,'zend',field).or.   
     1     inp_compare(.false.,'end',field)) then
         return
      else if(inp_compare(.false.,'cvr_scaling',field)) then
         if(.not. inp_f(cvr_scaling)) call errquit
     1        ('hnd_autoz_input:  reading cvr_scaling',911)
      else if (inp_compare(.false.,'ijbond',field)) then
         nbond = nbond + 1
         if(.not. inp_i(ijbond(1,nbond))) call errquit
     $        ('hnd_autoz_input: reading ijbond',1)
         if(.not. inp_i(ijbond(2,nbond))) call errquit
     $        ('hnd_autoz_input: reading ijbond',2)
      else if(inp_compare(.false.,'ijkang',field)) then
         nangl = nangl + 1
         if(.not. inp_i(ijkang(1,nangl))) call errquit
     $        ('hnd_autoz_input: reading ijkang',1)
         if(.not. inp_i(ijkang(2,nangl))) call errquit
     $        ('hnd_autoz_input: reading ijkang',2)
         if(.not. inp_i(ijkang(3,nangl))) call errquit
     $        ('hnd_autoz_input: reading ijkang',3)
      else if(inp_compare(.false.,'ijklto',field)) then
         ntors=ntors+1
         if(.not. inp_i(ijklto(1,ntors))) call errquit
     $        ('hnd_autoz_input: reading ijklto',1)
         if(.not. inp_i(ijklto(2,ntors))) call errquit
     $        ('hnd_autoz_input: reading ijklto',2)
         if(.not. inp_i(ijklto(3,ntors))) call errquit
     $        ('hnd_autoz_input: reading ijklto',3)
         if(.not. inp_i(ijklto(4,ntors))) call errquit
     $        ('hnd_autoz_input: reading ijklto',4)
      else if(inp_compare(.false.,'ijklop',field)) then
         nopla=nopla+1
         if(.not. inp_i(ijklop(1,nopla))) call errquit
     $        ('hnd_autoz_input: reading ijklop',1)
         if(.not. inp_i(ijklop(2,nopla))) call errquit
     $        ('hnd_autoz_input: reading ijklop',2)
         if(.not. inp_i(ijklop(3,nopla))) call errquit
     $        ('hnd_autoz_input: reading ijklop',3)
         if(.not. inp_i(ijklop(4,nopla))) call errquit
     $        ('hnd_autoz_input: reading ijklop',4)
      else if(inp_compare(.false.,'ijklnb',field)) then
         nlinb=nlinb+1
         if(.not. inp_i(ijklnb(1,nlinb))) call errquit
     $        ('hnd_autoz_input: reading ijklnb',1)
         if(.not. inp_i(ijklnb(2,nlinb))) call errquit
     $        ('hnd_autoz_input: reading ijklnb',2)
         if(.not. inp_i(ijklnb(3,nlinb))) call errquit
     $        ('hnd_autoz_input: reading ijklnb',3)
         if(.not. inp_i(ijklnb(4,nlinb))) call errquit
     $        ('hnd_autoz_input: reading ijklnb',4)
      else
         call errquit('hnd_autoz_input: unkown directive',0)
      endif
      goto 10                   ! Read next line
c
      end
      subroutine  geom_zmt_input(geom,
     1                           coords,tags,charge,mass,ncenter,
     2                           nzinput,units,oprint)
      implicit none
#include "stdio.fh"
#include "geom.fh"
      integer      mxizmt
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxizmt=6000)
      parameter    (mxzmat=1500)
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxcoor=1500)

      character*8 zvarname(mxcoor)
      double precision zvarsign(mxcoor)

      integer      i
      integer      nzinput
      integer      nizmat
      integer      izmat
      integer      nzmat
      integer      nzvar
      integer      nvar
      integer      izfrz
      integer      nzfrz
      integer      nat
      integer      nuc
      integer      ir
      integer      iw
      integer      ncenter
      integer      icenter
      integer      atn
      character*16 tags
      character*16 element
      character*16 atmlab 
      character*2  symbol
      logical      is_atom
      integer      geom
      logical      dbug 
      logical      some 
      logical      oprint
      double precision zmat
      double precision c
      double precision zan
      double precision units
      common/hnd_iofile/ir,iw
      common/hnd_molnuc/nuc(mxatom)
      common/hnd_mollab/atmlab(mxatom)
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_zmtdat/zmat(mxcoor),izmat(mxizmt),nizmat
      common/hnd_zmtfrz/izfrz(mxcoor),nzfrz
      double precision coords(3,*)
      double precision charge(*)
      double precision   mass(*)
      dimension tags(*)
c
      dbug=.false.
      some=.true.
      some=some.or.dbug
c
      ir=LuIn
      iw=LuOUT 
      if(dbug) then
         write(iw,9999)
      endif
c
c     ----- call routine to read z-matrix data -----
c
      call geom_zmt_geo(coords,charge,tags,ncenter,
     1                  nizmat,izmat,nzvar,izfrz,nzfrz,
     2                  units,zvarname,zvarsign)
c
c     -----    set parameters for -nwchem-     -----
c
      do icenter=1,ncenter
              is_atom = geom_tag_to_element(tags(icenter),
     1                                      symbol,element,atn)
         atmlab(icenter)=tags(icenter)
              is_atom = geom_atn_to_default_mass(atn,     
     1                                      mass(icenter))     
         if(dbug) then
            write(iw,9998) icenter,tags(icenter),charge(icenter),
     1      coords(1,icenter),coords(2,icenter),coords(3,icenter)
         endif
      enddo
c
c     ----- store all zmatrix data for nwchem -----
c
      if (.not. geom_zmt_put_data(geom, nizmat, izmat,
     $     nzfrz, izfrz, nzvar, zvarname, zvarsign))
     $     call errquit('geom_zmt_input: bad geom handle',geom)

      if(dbug) then
         write(iw,9994)
         write(iw,9997)
         write(iw,9996) (izmat(i),i=1,nizmat)
         write(iw,9995)
         write(iw,9990) nzfrz
         if(nzfrz.gt.0) then
            write(iw,9992)
            write(iw,9991)
            write(iw,9996) (izfrz(i),i=1,nzfrz)
            write(iw,9995)
         endif
      endif
c
      nat=ncenter
      do i=1,ncenter
         nuc(i)=nint(charge(i))
         c(1,i)=coords(1,i)*units
         c(2,i)=coords(2,i)*units
         c(3,i)=coords(3,i)*units
      enddo
      if(dbug) then
         do i=1,nat
            write(iw,9993) i,nuc(i),c(1,i),c(2,i),c(3,i)
         enddo
      endif
c
c     ----- print -zmat- data -----
c
      if (oprint)
     $     call geom_prt_autoz(zmat,nzvar,izmat,nizmat,c)
c
c     ----- set number of inputed -z- coordinates -----
c
      nzinput=nzvar
      nzmat=nizmat
      nvar =nzvar
c
      return
 9999 format(' in -geom_zmt_input- ')
 9998 format(' --- ',i5,1x,a16,1x,f4.0,1x,3f12.6)
 9997 format(' $zmat izmat = ')
 9996 format(12(i5,','))
 9995 format(' $end ')
 9994 format(' storing -izmat- via - geom_zmt_put_izmat - ')
 9993 format(1x,i5,i5,3f10.6)
 9992 format(' storing -izfrz- via - geom-zmt-put-izfrz - ')
 9991 format(' $opt  izfrz = ')
 9990 format(' nzfrz = ',i5)
      end
      logical function geom_zmt_opn_fil(
     1                 opn_bmat,opn_binv,
     2                 opn_zsym,opn_zind)
      implicit none
c
      character*255 zmt_bmat_fil
      character*255 zmt_binv_fil
      character*255 zmt_zsym_fil
      character*255 zmt_zind_fil
      common/hnd_zmt_fil/zmt_bmat_fil,zmt_binv_fil,
     1                   zmt_zsym_fil,zmt_zind_fil
c
      integer     ir
      integer     iw
      logical     dbug
      common/hnd_iofile/ir,iw
c
      logical opn_bmat
      logical opn_binv
      logical opn_zsym
      logical opn_zind
      integer ift2
      integer ift3
      integer ift4
      integer ift5
      character*8 errmsg(3)
      data errmsg /'program ','stop in ','zmt_opn '/
      data    ift2 /32/
      data    ift3 /33/
      data    ift4 /34/
      data    ift5 /35/
c
      dbug=.false.
c
      if(opn_bmat) then
         call util_file_name('zmt.bmat',.true.,.true.,zmt_bmat_fil)
         open(unit=ift2,file=zmt_bmat_fil,form='unformatted',
     1        access='sequential',status='unknown',err=911)
         if(dbug) then
            write(iw,*) 'zmt_bmat_fil = ',zmt_bmat_fil
         endif
      endif
      if(opn_binv) then
         call util_file_name('zmt.binv',.true.,.true.,zmt_binv_fil)
         open(unit=ift3,file=zmt_binv_fil,form='unformatted',
     1        access='sequential',status='unknown',err=911)
         if(dbug) then
            write(iw,*) 'zmt_binv_fil = ',zmt_binv_fil
         endif
      endif
      if(opn_zind) then
         call util_file_name('zmt.zind',.true.,.true.,zmt_zind_fil)
         open(unit=ift4,file=zmt_zind_fil,form='unformatted',
     1        access='sequential',status='unknown',err=911)
         if(dbug) then
            write(iw,*) 'zmt_zind_fil = ',zmt_zind_fil
         endif
      endif
      if(opn_zsym) then
         call util_file_name('zmt.zsym',.true.,.true.,zmt_zsym_fil)
         open(unit=ift5,file=zmt_zsym_fil,form='unformatted',
     1        access='sequential',status='unknown',err=911)
         if(dbug) then
            write(iw,*) 'zmt_zsym_fil = ',zmt_zsym_fil
         endif
      endif
         geom_zmt_opn_fil=.true.
      return
c
  911 continue
      if(dbug) then
         write(iw,*) 'zmt_bmat_fil = ',zmt_bmat_fil
         write(iw,*) 'zmt_binv_fil = ',zmt_binv_fil
         write(iw,*) 'zmt_zsym_fil = ',zmt_zsym_fil
         write(iw,*) 'zmt_zind_fil = ',zmt_zind_fil
      endif
      call hnd_hnderr(3,errmsg)
      geom_zmt_opn_fil=.false.
      return
      end
      logical function geom_zmt_cls_fil(
     1                 opn_bmat,opn_binv,
     2                 opn_zsym,opn_zind)
      implicit none
c
      character*255 zmt_bmat_fil
      character*255 zmt_binv_fil
      character*255 zmt_zsym_fil
      character*255 zmt_zind_fil
      common/hnd_zmt_fil/zmt_bmat_fil,zmt_binv_fil,
     1                   zmt_zsym_fil,zmt_zind_fil
c
      integer     ir
      integer     iw
      logical     dbug
      common/hnd_iofile/ir,iw
c
      logical opn_bmat
      logical opn_binv
      logical opn_zind
      logical opn_zsym
      integer ift2
      integer ift3
      integer ift4
      integer ift5
      character*8 errmsg(3)
      data errmsg /'program ','stop in ','zmt_cls '/
      data    ift2 /32/
      data    ift3 /33/
      data    ift4 /34/
      data    ift5 /35/
c
      dbug=.false.
c
      close(unit=ift2,status='keep',err=911)
      close(unit=ift3,status='keep',err=911)
      close(unit=ift4,status='keep',err=911)
      close(unit=ift5,status='keep',err=911)
      geom_zmt_cls_fil=.true.
      return
c
  911 continue
      if(dbug) then
         write(iw,*) 'zmt_bmat_fil = ',zmt_bmat_fil
         write(iw,*) 'zmt_binv_fil = ',zmt_binv_fil
         write(iw,*) 'zmt_zsym_fil = ',zmt_zsym_fil
         write(iw,*) 'zmt_zind_fil = ',zmt_zind_fil
      endif
      call hnd_hnderr(3,errmsg)
      geom_zmt_cls_fil=.false.
      return
      end
      logical function geom_zmt_put_bmat(bmat,ncart,nzvar,mcart)
      implicit none
c
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
c                   mxcoor=max(mxcart,mxzmat)
      integer nzvar
      integer ncart
      integer mcart
      integer ir
      integer iw
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision bmat(mcart,nzvar)
c
      integer ift2
      data    ift2 /32/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      if(out) then
         write(iw,*) ' _put_bmat ,ncart,nzvar = ',ncart,nzvar
         if(dbug) then
            write(iw,*) ' _put_bmat '
            call hnd_prsq(bmat,nzvar,ncart,mcart)
         endif
      endif
      if(mcart.gt.mxcart) then
         write(iw,*) ' in geom_zmt_put_bmat, mcart.gt.mxcart. stop '
         call errquit('geom_input: mcart too big',0)
      endif
      if(nzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_put_bmat, nzvar.gt.mxzmat. stop '
         call errquit('geom_input: nzvar too big',0)
      endif
c
      rewind ift2
      write(ift2) bmat
c
      geom_zmt_put_bmat=.true.
      return
      end
      logical function geom_zmt_get_bmat(bmat,ncart,nzvar,mcart)
      implicit none
c     
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
c                   mxcoor=max(mxcart,mxzmat)
      integer nzvar
      integer ncart
      integer mcart
      integer ir
      integer iw
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision bmat(mcart,nzvar)
c
      integer ift2
      data    ift2 /32/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      if(out) then
         write(iw,*) ' _get_bmat ,ncart,nzvar = ',ncart,nzvar
      endif
c
      if(mcart.gt.mxcart) then
         write(iw,*) ' in geom_zmt_get_bmat, mcart.gt.mxcart. stop '
         call errquit('geom_input: mcart too big',0)
      endif
      if(nzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_get_bmat, nzvar.gt.mxzmat. stop '
         call errquit('geom_input: nzvar too big',0)
      endif
c
      rewind ift2
      read(ift2) bmat
c
      if(dbug) then
         write(iw,*) ' _get_bmat '
         call hnd_prsq(bmat,nzvar,ncart,mcart)
      endif
c
      geom_zmt_get_bmat=.true.
      return
      end
      logical function geom_zmt_put_zsym(zsym,nzvar,mzvar)
      implicit none
c
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
c                   mxcoor=max(mxcart,mxzmat)
      integer nzvar
      integer mzvar
      integer ir
      integer iw
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision zsym(mzvar,mzvar)
c
      integer ift5
      data    ift5 /35/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      if(out) then
         write(iw,*) ' _put_zsym ,mzvar,nzvar = ',mzvar,nzvar
         if(dbug) then
            write(iw,*) ' _put_zsym '
            call hnd_prsq(zsym,nzvar,nzvar,mzvar)
         endif
      endif
c
      if(mzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_put_zsym, mzvar.gt.mxzmat. stop '
         call errquit('geom_input: mzvar too big',0)
      endif
      if(nzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_put_zsym, nzvar.gt.mxzmat. stop '
         call errquit('geom_input: nzvar too big',0)
      endif
c
      rewind ift5
      write(ift5) zsym
c
      geom_zmt_put_zsym=.true.
      return
      end
      logical function geom_zmt_get_zsym(zsym,nzvar,mzvar)
      implicit none
c     
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
c                   mxcoor=max(mxcart,mxzmat)
      integer nzvar
      integer mzvar
      integer ir
      integer iw
      logical dbug
      logical out 
      common/hnd_iofile/ir,iw
      double precision zsym(mzvar,mzvar)
c
      integer ift5
      data    ift5 /35/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      if(out) then
         write(iw,*) ' _get_zsym ,mzvar,nzvar = ',mzvar,nzvar
      endif
c
      if(mzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_get_zsym, mzvar.gt.mxzmat. stop '
         call errquit('geom_input: mzvar too big',0)
      endif
      if(nzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_get_zsym, nzvar.gt.mxzmat. stop '
         call errquit('geom_input: nzvar too big',0)
      endif
c
      rewind ift5
      read(ift5) zsym
c
      if(dbug) then
         write(iw,*) ' _get_zsym'
         call hnd_prsq(zsym,nzvar,nzvar,mzvar)
      endif
c
      geom_zmt_get_zsym=.true.
      return
      end
      logical function geom_zmt_put_zind(zind,nzvar,mzvar)
      implicit none
c
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
c                   mxcoor=max(mxcart,mxzmat)
      integer nzvar
      integer mzvar
      integer ir
      integer iw
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision zind(mzvar,mzvar)
c
      integer ift4
      data    ift4 /34/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      if(out) then
         write(iw,*) ' _put_zind ,mzvar,nzvar = ',mzvar,nzvar
         if(dbug) then
            write(iw,*) ' _put_zind '
            call hnd_prsq(zind,nzvar,nzvar,mzvar)
         endif
      endif
c
      if(mzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_put_zind, mzvar.gt.mxzmat. stop '
         call errquit('geom_input: mzvar too big',0)
      endif
      if(nzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_put_zind, nzvar.gt.mxzmat. stop '
         call errquit('geom_input: nzvar too big',0)
      endif
c
      rewind ift4
      write(ift4) zind
c
      geom_zmt_put_zind=.true.
      return
      end
      logical function geom_zmt_get_zind(zind,nzvar,mzvar)
      implicit none
c     
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
c                   mxcoor=max(mxcart,mxzmat)
      integer nzvar
      integer mzvar
      integer ir
      integer iw
      logical dbug
      logical out 
      common/hnd_iofile/ir,iw
      double precision zind(mzvar,mzvar)
c
      integer ift4
      data    ift4 /34/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      if(out) then
         write(iw,*) ' _get_zind ,mzvar,nzvar = ',mzvar,nzvar
      endif
c
      if(mzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_put_zind, mzvar.gt.mxzmat. stop '
         call errquit('geom_input: mzvar too big',0)
      endif
      if(nzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_put_zind, nzvar.gt.mxzmat. stop '
         call errquit('geom_input: nzvar too big',0)
      endif
c
      rewind ift4
      read(ift4) zind
c
      if(dbug) then
         write(iw,*) ' _get_zind'
         call hnd_prsq(zind,nzvar,nzvar,mzvar)
      endif
c
      geom_zmt_get_zind=.true.
      return
      end
      logical function geom_zmt_put_binv(binv,ncart,nzvar,mcart)
      implicit none
c     
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
c                   mxcoor=max(mxcart,mxzmat)
      integer nzvar
      integer ncart
      integer mcart
      integer ir
      integer iw
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision binv(mcart,nzvar)
c
      integer ift3
      data    ift3 /33/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      if(out) then
         write(iw,*) ' _get_zind ,ncart,nzvar = ',ncart,nzvar
         if(dbug) then
            write(iw,*) ' _put_binv '
            call hnd_prsq(binv,nzvar,ncart,mcart)
         endif
      endif
c
      if(mcart.gt.mxcart) then
         write(iw,*) ' in geom_zmt_put_binv, mcart.gt.mxcart. stop '
         call errquit('geom_input: mcart too big',0)
      endif
      if(nzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_put_binv, nzvar.gt.mxzmat. stop '
         call errquit('geom_input: nzvar too big',0)
      endif
c
      rewind ift3
      write(ift3) binv
c
      geom_zmt_put_binv=.true.
      return
      end
      logical function geom_zmt_get_binv(binv,ncart,nzvar,mcart)
      implicit none
c     
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
c                   mxcoor=max(mxcart,mxzmat)
      integer nzvar
      integer ncart
      integer mcart
      integer ir
      integer iw
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision binv(mcart,nzvar)
c
      integer ift3
      data    ift3 /33/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      if(out) then
         write(iw,*) ' _get_binv ,ncart,nzvar = ',ncart,nzvar
      endif
c
      if(mcart.gt.mxcart) then
         write(iw,*) ' in geom_zmt_put_binv, mcart.gt.mxcart. stop '
         call errquit('geom_input: mcart too big',0)
      endif
      if(nzvar.gt.mxzmat) then
         write(iw,*) ' in geom_zmt_put_binv, nzvar.gt.mxzmat. stop '
         call errquit('geom_input: nzvar too big',0)
      endif
c
      rewind ift3
      read(ift3) binv
c
      if(dbug) then
         write(iw,*) ' _get_binv '
         call hnd_prsq(binv,nzvar,ncart,mcart)
      endif
c
      geom_zmt_get_binv=.true.
      return
      end
      subroutine  geom_no_autoz(geom,coords,charge,tags,ncenter,
     1                          nzinput)
      implicit none
#include "stdio.fh"
      integer      mxizmt
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxizmt=6000)
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
c                   mxcoor=max(mxcart,mxzmat)
      integer      ir
      integer      iw
      integer      i
      integer      nzinput
      integer      nizmat
      integer      izmat
      integer      nzmat
      integer      nzvar
      integer      nvar 
      integer      ncenter
      integer      nuc
      integer      nat
      integer      geom
      logical      dbug 
      logical      some 
      double precision coords
      double precision charge
      double precision zmat  
      double precision c
      double precision zan
      character*16 tags
      character*16 atmlab
      character*8  errmsg
      common/hnd_iofile/ir,iw
      common/hnd_molnuc/nuc(mxatom)
      common/hnd_mollab/atmlab(mxatom)
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_zmtdat/zmat(mxcoor),izmat(mxizmt),nizmat
      dimension coords(3,*)
      dimension charge(  *)
      dimension   tags(  *)
      dimension errmsg(3)
      data errmsg /'program ','stop in ','-auto_z-'/
c
      some=.false.
      dbug=.false.
c
      ir=LuIN 
      iw=LuOut 
      if(some) then
         write(iw,9999)
      endif
c
      nat=ncenter
      do i=1,ncenter
         atmlab(i)=tags(i)
         nuc(i)=nint(charge(i))
         c(1,i)=coords(1,i)
         c(2,i)=coords(2,i)
         c(3,i)=coords(3,i)
      enddo
      if(dbug) then
         do i=1,nat    
            write(iw,9998) i,nuc(i),c(1,i),c(2,i),c(3,i)
         enddo
      endif
c
c     ----- no z-matrix data -----
c
      nzmat=0
      nzvar=0
      nvar =3*nat
c
      return
 9999 format(' in -geom_no_autoz- ')
 9998 format(1x,i5,i5,3f10.6)
      end
      subroutine  geom_zmt_autoz(rtdb,geom,
     1                           coords,charge,tags,ncenter,
     2                           nzinput,oforce_autoz, oprint)
      implicit none
#include "stdio.fh"
      integer      rtdb
      integer      mxizmt
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxizmt=6000)
      parameter    (mxzmat=1500)
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxcoor=1500)
      integer      ir
      integer      iw
      integer      i
      integer      nzinput
      integer      nizmat
      integer      izmat
      integer      nzmat
      integer      nzvar
      integer      nvar 
      integer      ncenter
      integer      nuc
      integer      nat
      integer      geom
      logical      oforce_autoz
      logical      zdone
      logical      dbug 
      logical      some 
      logical     oprint
      double precision coords(3,*)
      double precision charge(*)
      double precision zmat  
      double precision c
      double precision zan
      character*16 tags(*)
      character*16 atmlab
      character*8  errmsg(3)
      common/hnd_iofile/ir,iw
      common/hnd_molnuc/nuc(mxatom)
      common/hnd_mollab/atmlab(mxatom)
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_zmtdat/zmat(mxcoor),izmat(mxizmt),nizmat
      character*1 zvarname(mxcoor)
      double precision zvarsign(mxcoor)
      logical geom_zmt_put_data
      external geom_zmt_put_data
      data errmsg /'program ','stop in ','-auto_z-'/
c
      if(.not.oforce_autoz) then
         if(nzinput.eq.(3*ncenter-6)) then
            return
         endif
      endif
c
      dbug=.false.
      some=.false.
      some=some.or.dbug
c
      ir=LuIN 
      iw=LuOut 
      if(some) then
         write(iw,9999)
      endif
c
c     ----- call auto-zmat -----
c
      nat=ncenter
      do i=1,ncenter
         atmlab(i)=tags(i)
         nuc(i)=nint(charge(i))
         c(1,i)=coords(1,i)
         c(2,i)=coords(2,i)
         c(3,i)=coords(3,i)
      enddo
      if(dbug) then
         do i=1,nat    
            write(iw,9998) i,nuc(i),c(1,i),c(2,i),c(3,i)
         enddo
      endif
c
      call hnd_zmtmak(zdone,geom,oprint)
c
      if(.not.zdone) then
         write(iw,9997)
         nzmat=0
      endif
c
c     ----- save -izmat- and -zmat- -----
c
      do i = 1, nzvar
         zvarname(i) = ' '
         zvarsign(i) = 1d0
      enddo
      if (.not. geom_zmt_put_data(geom, nizmat, izmat,
     $     0, 0, nzvar, zvarname, zvarsign))
     $     call errquit('geom_zmt_input: bad geom handle',geom)
c
c     ----- print -zmat- data -----
c
      if (oprint) 
     $     call geom_prt_autoz(zmat,nzvar,izmat,nizmat,c)
c
      if (.not. zdone) call errquit
     $     ('geom_input: autoz failed ... adjust AUTOZ options'//
     $     ' or use NOAUTOZ',0)
c
      return
 9999 format(' in -geom_zmt_autoz- ')
 9998 format(1x,i5,i5,3f10.6)
 9997 format(' -auto_z- has failed to create a good z-matrix.',/,
     1       ' program will print the z-matrix, but will not use it.')
 9996 format(' $zmat izmat = ')
 9995 format(12(i5,','))
 9994 format(' $end ')
 9993 format(' storing -izmat- via - geom_zmt_put_izmat - ')
      end
      subroutine geom_prt_autoz(zmat,nzvar,izmat,nzmat,c)
      implicit double precision (a-h,o-z)
#include "stdio.fh"
c
c     ----- get the -zmat- values ----- 
C
      parameter   (mxatom=500)
      parameter   (mxcart=3*mxatom)
      parameter   (na=10)
      LOGICAL     DBUG
      LOGICAL     OUT
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      DIMENSION B(MXCART,1)
      DIMENSION ZMAT(*),IZMAT(*)
      DIMENSION A(3,NA)
      DIMENSION C(3,*)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-prt:z -'/
      DATA PI2    /6.28318530717958D+00/
      DATA DEGREE /360.0D+00/
      DATA BOHR   /0.52917715D+00/
      DATA ZERO   /0.0D+00/
      DATA FOUR   /4.0D+00/
C
      DBUG  =.FALSE. 
      OUT   =.FALSE. 
      OUT   =OUT.OR.DBUG
C
      PIO2  =PI2/FOUR
      RTOD  =PI2/DEGREE
C
      IF(OUT) THEN
         WRITE(IW,9997) 
         WRITE(IW,9996) (IZMAT(I),I=1,NZMAT)
      ENDIF
C
C     ----- ZERO OUT THE B MATRIX -----
C
      NCART=MXCART
      DO J = 1,1     
         DO I = 1,NCART 
            B(I,J) = ZERO
         ENDDO
      ENDDO
C
      IADD=1
      I   =0
      II  =0
      DO IZVAR=1,NZVAR
      I    =I +1            
      II   =II+1
      I    =   1            
      ITYPE=IZMAT(IADD)
      IF(DBUG) THEN
         WRITE(IW,9999) IZVAR,IADD,ITYPE,NZMAT,NZVAR
      ENDIF
      IF(ITYPE.EQ.1) THEN
C
C     ----- BOND STRETCH -----
C
         CALL HND_BSTR(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             C,B,NCART,BOHR)
         ZMAT(II) = EQVAL
         IADD = IADD + 3
      ELSEIF(ITYPE.EQ.2 .or. itype.eq.7) THEN
C
C     ----- ANGLE BEND -----
C
         CALL HND_BEND(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),
     2             C,B,NCART,RTOD,BOHR)
         ZMAT(II) = EQVAL
         IADD = IADD + 4
         if (itype .eq. 7) iadd = iadd + 1 ! Additional orientation
      ELSEIF(ITYPE.EQ.3) THEN
C
C     ----- TORSION -----
C
         CALL HND_TORS(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,BOHR)
         ZMAT(II) = EQVAL
         IADD = IADD + 5
      ELSEIF(ITYPE.EQ.4) THEN
C
C     ----- OUT OF PLANE ANGLE BEND -----
C
         CALL HND_OPLA(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,PIO2,BOHR)
         ZMAT(II) = EQVAL
         IADD = IADD + 5
      ELSEIF(ITYPE.EQ.5) THEN
C
C     ----- LINEAR ANGLE BEND -----
C
         CALL HND_LIBE(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,A,NA,BOHR)
         ZMAT(II  ) = EQVAL
         ZMAT(II+1) = EQVAL
         I=I+1
         IADD=IADD+5
      ELSEIF(ITYPE.EQ.6) THEN
C
C     ----- DIHEDRAL ANGLE BETWEEN TWO PLANES SHARING ONE ATOM -----
C
         CALL HND_DIHPLA(EQVAL,I,IZMAT(IADD+1),C,B,NCART,RTOD)
         ZMAT(II)=EQVAL
         IADD=IADD+6
C
      ELSE
         WRITE(IW,9995) IZVAR,IADD,ITYPE
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF    
C
      ENDDO
C
C     ----- PRINT -----
C
      CALL HND_PRTZMT(NZVAR,IZMAT,ZMAT,DUM,0)
C
      RETURN
 9999 FORMAT(' IN -prt-autoz- IZVAR,IADD,ITYPE,NZMAT,NZVAR = ',5I5)
 9997 FORMAT(' IN -prt-autoz- IZMAT = ')
 9996 FORMAT(12I5)
 9995 FORMAT(' ERROR IN -prt-autoz- , IZVAR,IADD,ITYPE = ',3I5)
      END
      subroutine geom_zmt_icoord(rtdb,geom)
      implicit none
#include "rtdb.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer      rtdb
      integer      geom
      integer      mxizmt
      integer      mxzmat
      integer      mxatom
      integer      mxcart
      integer      mxcoor
      parameter    (mxizmt=6000)
      parameter    (mxatom=500)
      parameter    (mxcart=3*mxatom)
      parameter    (mxzmat=1500)
      parameter    (mxcoor=1500)
      integer      ir
      integer      iw
      integer      i
      integer      nprint
      integer      list    
      integer      ifzmat
      integer      natzmt
      integer      nizmat
      integer      izmat
      integer      nzmat
      integer      nzvar
      integer      nvar
      integer      nuc
      integer      nat
      logical      geom_zmt_get_izmat
      logical      geom_zmt_get_nizmat
      logical      dbug
      logical      some
      double precision zmat
      double precision c
      double precision zan
      character*8  errmsg
      common/hnd_output/nprint
      common/hnd_listng/list
      common/hnd_iofile/ir,iw
      common/hnd_molnuc/nuc(mxatom)
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      common/hnd_zmtdat/zmat(mxcoor),izmat(mxizmt),nizmat
      common/hnd_zmtpar/nzmat,nzvar,nvar
      common/hnd_zmtopt/ifzmat,natzmt
      dimension errmsg(3)
      data errmsg /'program ','stop in ','-icoord-'/
c
      dbug=.false.
      some=.false.
      some=some.or.dbug
c
      nprint=0
      list  =0
      ir    =LuIn 
      iw    =LuOut 
      if(some) then
         write(iw,9999)
      endif
c
c     ----- read in -izmat- from -nwchem- -----
c
      do i=1,mxizmt
         izmat(i)=0
      enddo
c
*      status = rtdb_parallel(.false.)
*      if(ga_nodeid().eq.0) then
         if(.not.geom_zmt_get_nizmat(geom,nizmat))
     1   call errquit('geom_input: geom_zmt_get_nizmat failed',0)
         if(.not.geom_zmt_get_izmat(geom,izmat,nizmat))
     1   call errquit('geom_input: geom_zmt_get_izmat failed',0)
*      endif
*      call ga_brdcst(18,nizmat,mitob(1     ),0)
*      call ga_brdcst(19,izmat ,mitob(mxizmt),0)

      if(dbug) then
         write(iw,9998)
         write(iw,9997)
         write(iw,9996) (izmat(i),i=1,nizmat)
         write(iw,9995)
      endif
*      call ga_sync()
*      status = rtdb_parallel(.true.)
c
c     ----- call -icoord- -----
c
      ifzmat=1
      call hnd_icoord(geom)
c
c     RJH.  Since AUTOZ does not presently handle linear
c     bends it might not have generated a complete set
c     of coordinates.  Check that the no. of linearly independent
c     variables (NVAR) is at least NCART-6.  
c
c     When this fails and what else this breaks I don't yet know
c
**      write(6,*) ' NVAR NAT ', nvar, nat
      if (nvar .lt. 3*nat-6) call errquit('geom_zmt_icoord:'//
     $     ' insufficient linearly independent internals ...'//
     $     ' adjust AUTOZ options or use NOAUTOZ', nvar)
      return
 9999 format(' in -geom_zmt_icoord ')
 9998 format(' reading -izmat- via -geom_zmt_get_izmat- ')
 9997 format(' $zmat izmat = ')
 9996 format(12(i5,','))
 9995 format(' $end ')
      end
      subroutine geom_zmat_to_cart(s,nzvar,c)
      implicit none
#include "global.fh"
      integer nzvar
      double precision s(nzvar), c(3,*)
c
c     Return the cartesians corresponding to the given internals.
c
c     The internals are in units of angstrom and degrees.
c     The cartesians will be returned in units of bohr.
c     
      double precision zmat
      integer izmat, nat
      integer mxcoor, mxizmt, nizmat
      integer mxatom
      double precision cc, zan
      parameter (mxatom=500)
      parameter (mxcoor = 1500)
      parameter (mxizmt=6000)
      common/hnd_zmtdat/zmat(mxcoor),izmat(mxizmt),nizmat
      common/hnd_molxyz/cc(3,mxatom),zan(mxatom),nat
c
      integer iadd, itype, izvar, i, j, k, l, ind, iat, orient
      double precision rij, aijk, aijkl, angs_to_au, pi, deg_to_rad
      double precision ekj(3), ekl(3), ejkl(3), alpha, beta, gamma
      double precision dkj, dkl, djk, djl, dklkj, aijl, djkl
      double precision ejk(3), ejl(3), d
      double precision ddot
      external ddot
c
      angs_to_au = 1d0/0.52917715D+00
      pi  = 4d0*atan(1d0)
      deg_to_rad = pi/180d0
c     
      iadd=1
      izvar = 1
c
      c(1,1) = 0d0
      c(2,1) = 0d0
      c(3,1) = 0d0
      if (nat .ge. 2) then
         itype = izmat(iadd)
         i = izmat(iadd+1)
         j = izmat(iadd+2)
**         write(6,*) izvar, itype, i, j, s(izvar)
         c(1,2) = 0d0
         c(2,2) = 0d0
         c(3,2) = s(izvar)*angs_to_au
         izvar = izvar + 1
         iadd  = iadd + 3
      endif
      if (nat .ge. 3) then
         itype = izmat(iadd)
         i = izmat(iadd+1)
         j = izmat(iadd+2)
**         write(6,*) izvar, itype, i, j, s(izvar)
         rij = s(izvar)*angs_to_au
         iadd = iadd + 3
         izvar = izvar + 1
c
         itype = izmat(iadd)
         k = izmat(iadd+3)
**         write(6,*) izvar, itype, i, j, k, s(izvar)
         if (s(izvar).lt.0d0 .or. s(izvar).gt.180d0) goto 1000
         aijk = s(izvar)*deg_to_rad
         iadd = iadd + 4
         izvar = izvar + 1
c         
         c(1,3) = rij*sin(pi-aijk)
         c(2,3) = 0d0 
         c(3,3) = c(3,j) + sign(1d0,c(3,j)-c(3,k))*rij*cos(pi-aijk)
      endif
c
      do iat = 4, nat
         itype = izmat(iadd)
         i = izmat(iadd+1)
         j = izmat(iadd+2)
**         write(6,*) izvar, itype, i, j, s(izvar)
         rij = s(izvar)*angs_to_au
         iadd = iadd + 3
         izvar = izvar + 1
c
         itype = izmat(iadd)
         k = izmat(iadd+3)
**         write(6,*) izvar, itype, i, j, k, s(izvar)
         if (s(izvar).lt.0d0 .or. s(izvar).gt.180d0) goto 1000
         aijk = s(izvar)*deg_to_rad
         iadd = iadd + 4
         izvar = izvar + 1
c
         itype = izmat(iadd)
         if (itype.eq.3) then
            l = izmat(iadd+4)
**         write(6,*) izvar, itype, i, j, k, l, s(izvar)
            aijkl = s(izvar)*deg_to_rad
         else if (itype.eq.7) then
            l = izmat(iadd+3)
            orient= izmat(iadd+4)
            aijl  = s(izvar)*deg_to_rad
**            write(6,*) izvar, itype, i, j, k, l, orient, s(izvar)
         else
            call errquit('geom_zmat_to_cart: unknown 3rd type',itype)
         endif
         iadd = iadd + 5
         izvar = izvar + 1
c
         if (itype .eq. 3) then
c     
c     Form orthogonal left-handed coord system
c     ekj = (k->j)
c     ekl = (k->l) orthoged to ekj
c     ejkl= (ekj ^ ekl)
c     
            do ind = 1, 3
               ekj(ind) = c(ind,j) - c(ind,k)
               ekl(ind) = c(ind,l) - c(ind,k)
            enddo
            dkj = 1d0/sqrt(ddot(3, ekj, 1, ekj, 1))
            call dscal(3, dkj, ekj, 1)
            dklkj = ddot(3, ekl, 1, ekj, 1)
            do ind = 1, 3
               ekl(ind) = ekl(ind) - dklkj*ekj(ind)
            enddo
            dkl = 1d0/sqrt(ddot(3, ekl, 1, ekl, 1))
            call dscal(3, dkl, ekl, 1)
            call cross_product(ekj, ekl, ejkl)
c     
c     Unit vector from j to i is then written
c     eji = alpha*ekj + beta*ekl + gamma*ejkl
c     alpha = cos(180-aijk) = -cos(aijk)
c     beta  = cos(aijkl)*|sin(180-aijk)| = cos(aijkl)*|sin(aijk)|
c     |gamma| = from condition alpha**2 + beta**2 + gamma**2 = 1
c     with sign(gamma) = sign(aijkl) if -180<=aijkl<=180.
c     
            alpha = -cos(aijk)
            beta = cos(aijkl)*abs(sin(aijk))
            gamma = sqrt(max(0d0,1d0 - alpha**2 - beta**2))
            gamma = sign(gamma,aijkl)
c     
            do ind = 1, 3
               c(ind,i) = c(ind,j) + rij*(
     $              alpha*ekj(ind) + beta*ekl(ind) + gamma*ejkl(ind))
            enddo
         else if (itype .eq. 7) then
c     
c     Form unit vectors in the non-orthogonal directions
c     ejl = (j->l)
c     ejk = (j->k)
c     ejkl= (ejk^ejl)
c            
            do ind = 1, 3
               ejk(ind) = c(ind,k) - c(ind,j)
               ejl(ind) = c(ind,l) - c(ind,j)
            enddo
            call cross_product(ejk, ejl, ejkl)
            djk  = 1d0/sqrt(ddot(3, ejk, 1, ejk, 1))
            djl  = 1d0/sqrt(ddot(3, ejl, 1, ejl, 1))
            djkl = 1d0/sqrt(ddot(3, ejkl,1, ejkl,1))
            call dscal(3, djk, ejk, 1)
            call dscal(3, djl, ejl, 1)
            call dscal(3, djkl,ejkl,1)
c
c     unit vector from j to i is then
c     eji = alpha*ejl + beta*ejk + gamma*ejkl
c     cos(aijl) = alpha + beta*(ejk.ejl) = alpha   + beta*d
c     cos(aijk) = alpha*(ejk.ekl) + beta = alpha*d + beta
c
c     alpha = (cos(aijl) - d*cos(aijk)) / (1-d*d)
c     beta  = (cos(aijk) - d*cos(aijl)) / (1-d*d)
c     |gamma| = from condition |eji|=1d0
c     with sign of gamma from orientation
c
            d = ddot(3, ejk, 1, ejl, 1)
            alpha = (cos(aijl) - d*cos(aijk)) / (1d0-d*d)
            beta  = (cos(aijk) - d*cos(aijl)) / (1d0-d*d)
            gamma = sqrt(max(0d0,1d0 - alpha**2 - beta**2
     $           - 2d0*d*alpha*beta))
            if (orient .lt. 0) gamma = -gamma
c
            do ind = 1, 3
               c(ind,i) = c(ind,j) + rij*(
     $              alpha*ejl(ind) + beta*ejk(ind) + gamma*ejkl(ind))
            enddo
         endif
      enddo
         
c$$$
c$$$
c$$$         itype=izmat(iadd)
c$$$         if(itype.eq.1) then
c$$$c     
c$$$c     ----- bond stretch -----
c$$$c     
c$$$            i = izmat(iadd+1)
c$$$            j = izmat(iadd+2)
c$$$            write(6,*) izvar, ' bond ', i, j
c$$$            rij = 
c$$$            iadd = iadd + 3
c$$$         elseif(itype.eq.2) then
c$$$c     
c$$$c     ----- angle bend -----
c$$$c     
c$$$            i = izmat(iadd+1)
c$$$            j = izmat(iadd+2)   ! Central atom
c$$$            k = izmat(iadd+3)
c$$$            write(6,*) izvar
c$$$            iadd = iadd + 4
c$$$         elseif(itype.eq.3) then
c$$$c     
c$$$c     ----- torsion -----
c$$$c     
c$$$            i = izmat(iadd+1)
c$$$            j = izmat(iadd+2)
c$$$            k = izmat(iadd+3)
c$$$            l = izmat(iadd+4)
c$$$            iadd = iadd + 5
c$$$         elseif(itype.eq.4) then
c$$$c     
c$$$c     ----- out of plane angle bend -----
c$$$c     
c$$$            call errquit('ztoc: opla',0)
c$$$            iadd = iadd + 5
c$$$         elseif(itype.eq.5) then
c$$$c     
c$$$c     ----- linear angle bend -----
c$$$c     
c$$$            call errquit('ztoc: libe',0)
c$$$            iadd=iadd+5
c$$$         elseif(itype.eq.6) then
c$$$c     
c$$$c     ----- dihedral angle between two planes sharing one atom -----
c$$$c     
c$$$            call errquit('ztoc:twop',0)
c$$$            iadd=iadd+6
c$$$c     
c$$$         else
c$$$            call errquit('geom_zmat_to cart: ?', itype)
c$$$         endif    
c$$$c
c$$$      enddo
c
      return
c
c     Jump here if an angle is out of range ... try to help the user
c
 1000 if (ga_nodeid() .eq. 0) then
         write(6,3) izvar, s(izvar)
 3       format(//' Error ... internal coordinate number',i5,' is out',
     $        ' of range.'/
     $        ' Bond-angle = ',f12.6,' but we require 0<angle<180.'//
     $        ' You have several options: ',/
     $        ' 1. Disable use of internals for taking the step '/
     $        '    (they will still be used to generate the step)'//
     $        '    driver; noqstep; end'//
     $        '    and restart at the current geometry.  However, '/
     $        '    constants will not be exactly enforced, and',/
     $        '    constraints will not work at all.',
     $        ' 2. Manually adjust the current internal coordinates',
     $        ' to give '/
     $        '    a better guess that won''t step thru a linear bond.'/
     $        ' 3. Respecify the internal coordinates so that linear'/
     $        '    bond angles will be avoided.  Note that the range',/
     $        '    of dihedral angles is not restricted.'/)
         call util_flush(6)
      endif
      call ga_sync()
      call errquit('geom_zmat_to_cart: invalid bond angle',0)
c     
      end
      subroutine geom_zmat_ico_scale(s,nzvar,bscale,ascale)
      implicit none
c
      integer nzvar
      double precision s(nzvar)
      double precision bscale, ascale
c     
c     Scale the internal coordinates in s 
c     
      double precision zmat
      integer izmat
      integer mxcoor, mxizmt, nizmat
      parameter (mxcoor = 1500)
      parameter (mxizmt=6000)
      common/hnd_zmtdat/zmat(mxcoor),izmat(mxizmt),nizmat
c
      integer iadd, itype, izvar
c     
      iadd=1
c
*      write(6,*) ' ICO SCALE ', nzvar, bscale, ascale
c
      do izvar=1,nzvar
         itype=izmat(iadd)
*         write(6,*) izvar, ' is ', itype
         if(itype.eq.1) then
c     
c     ----- bond stretch -----
c     
            iadd = iadd + 3
         elseif(itype.eq.2 .or. itype.eq.7) then
c     
c     ----- angle bend -----
c     
            iadd = iadd + 4
            if (itype .eq. 7) iadd = iadd + 1
         elseif(itype.eq.3) then
c     
c     ----- torsion -----
c     
            iadd = iadd + 5
         elseif(itype.eq.4) then
c     
c     ----- out of plane angle bend -----
c     
            iadd = iadd + 5
         elseif(itype.eq.5) then
c     
c     ----- linear angle bend -----
c     
            iadd=iadd+5
         elseif(itype.eq.6) then
c     
c     ----- dihedral angle between two planes sharing one atom -----
c     
            iadd=iadd+6
c     
         else
            call errquit('geom_zmat_ico_scale: ?', itype)
         endif    
c
         if (itype .eq. 1) then
            s(izvar) = s(izvar) * bscale
         else
            s(izvar) = s(izvar) * ascale
         endif
c     
      enddo
c     
      end
      SUBROUTINE GEOM_ZMT_GEO(COORDS,CHARGE,TAGS,NCENTER,
     1                        IZ,IZMAT,NZMOD,ICFRZ,NCFRZ,
     2                        UNITS,zvarname,zvarsign)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -----  PARAMETERS DEFINING MAXIMUM VALUES -----
C
C       MAXGEO   = MAXIMUM NUMBER OF ATOMS
C       MAXWRD   = MAXIMUM NUMBER OF WORDS ON A LINE
C       MAXVAR   = MAXIMUM NUMBER OF VARIABLES IN -ZMT-
C
C     ----- VARIABLES FOR INPUTTING OF -Z- MATRIX -----
C
C         FLGZMT  =  TYPE OF VARIABLE  0=BLANK
C                                      1=FLOATING
C                                      2=INTEGER
C                                      3=ALPHANUMERIC
C         NUMZMT  =  NUMBER OF WORDS ON A LINE
C         ZMTCHR  =  NUMBER OF CHARACTERS IN EACH WORD
C         PRSZMT  =  ARRAY OF EXTERNAL WORDS FROM EACH LINE
C
C     ----- -Z- MATRIX DATA -----
C
C      ZVAL        =  DISTANCE, ANGLE AND TORSION ANGLE VALUE
C       ZMT        =  I, J, K, L INDICES
C
      PARAMETER     (MXATOM=500)
      PARAMETER     (MXCOOR=1500)
      PARAMETER     (MAXGEO=MXATOM+1,MAXWRD=40,MAXVAR=256)
      PARAMETER     (MAXPRM=100)
      PARAMETER     (MXIZMT=6000)
      PARAMETER     (MAXLST=10+1)
      CHARACTER*16  TAGS
      LOGICAL       GEOM_LST_PUT_COORD
      EXTERNAL      GEOM_LST_PUT_COORD
      LOGICAL       GEOM_LST_GET_COORD
      EXTERNAL      GEOM_LST_GET_COORD
      LOGICAL       STATUS
      LOGICAL       DBUG
      LOGICAL       OUT 
      LOGICAL       LST
      LOGICAL       NEWFIL
      CHARACTER*255 STRING
      INTEGER       FLGWRD
      INTEGER       FLGZMT
      INTEGER       FLGVAR
      INTEGER       ZMTCHR
      INTEGER       VARCHR
      CHARACTER*80  PRSWRD
      CHARACTER*80  PRSZMT
      CHARACTER*80  PRSVAR
      CHARACTER*80  WORD
      CHARACTER*80  TITLE
      INTEGER       ZMT
      CHARACTER*8   CHREND
      LOGICAL       IZFRZ
      LOGICAL       CART0
      LOGICAL       CART
      LOGICAL       FRZVAL
      LOGICAL       FRZVAR
      LOGICAL       GHOST
      CHARACTER*8   ATNAME
      CHARACTER*2   ATLBL1,ATLBL2,ATLBL3,ATLBL4
      CHARACTER*8   DUMMY
      CHARACTER*2   GH
      CHARACTER*1   BLK
      CHARACTER*80  BLNK80
      CHARACTER*80  GEOFIL
      CHARACTER*17  BIOSYM
      CHARACTER*5   BIOEND
      CHARACTER*5   CHAR5
      CHARACTER*12  CHAR12
      CHARACTER*2   CHAR2
      CHARACTER*2   SYMBOL
      CHARACTER*8   ATMNAM
      CHARACTER*8   BLNK8
      CHARACTER*5   BLNK5
      CHARACTER*2   BLNK2
      CHARACTER*1   UNDERS
      CHARACTER*1   DASH
      CHARACTER*1   DIGIT(10)
      CHARACTER*8   ERRMSG
      CHARACTER*8   WRDBAS
      CHARACTER*8   WRDEND
      CHARACTER*8   WRDXYZ
      CHARACTER*8   WRDGEO
      CHARACTER*8   WRDZMT
      CHARACTER*8   WRDIZM
      CHARACTER*8   WRDOPT
      CHARACTER*8   WRDFRZ
      CHARACTER*10  WR1VAR
      CHARACTER*10  WR2VAR
      CHARACTER*10  WR1CON
      CHARACTER*10  WR2CON
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_FREERD/PRSWRD(40),NUMCHR(40),FLGWRD(40),NUMWRD
      COMMON/HND_XYZGEO/XX(MAXGEO),YY(MAXGEO),ZZ(MAXGEO),
     1                  ATNAME(MAXGEO),ATNUM(MAXGEO),CART(MAXGEO)
      DIMENSION COORDS(3,*)
      DIMENSION CHARGE(  *)
      DIMENSION   TAGS(  *)
      DIMENSION  IZMAT(  *)   
      DIMENSION  ICFRZ(  *)    
      DIMENSION  IZFRZ(MXCOOR)
      DIMENSION  CART0(MAXGEO)
      DIMENSION  XXLST(MAXGEO,     2)
      DIMENSION  YYLST(MAXGEO,     2)
      DIMENSION  ZZLST(MAXGEO,     2)
      DIMENSION  GHOST(       MAXGEO)
      DIMENSION NUMZMT(       MAXGEO)
      DIMENSION PRSZMT(MAXWRD,MAXGEO)
      DIMENSION FLGZMT(MAXWRD,MAXGEO)
      DIMENSION ZMTCHR(MAXWRD,MAXGEO)
      DIMENSION FRZVAL(     3,MAXGEO)
      DIMENSION   ZVAL(     3,MAXGEO)
      DIMENSION    ZMT(     5,MAXGEO)
      DIMENSION   ZLST(     3,MAXGEO)
      DIMENSION   ZSTP(     3,MAXGEO)
      DIMENSION FRZVAR(       MAXVAR)
      DIMENSION NUMVAR(       MAXVAR)
      DIMENSION PRSVAR(MAXWRD,MAXVAR)
      DIMENSION FLGVAR(MAXWRD,MAXVAR)
      DIMENSION VARCHR(MAXWRD,MAXVAR)
      character*8 zvarname(*) ! Name of Z matrix variables
      double precision zvarsign(*) ! Attached signs
      DIMENSION ATLBL1(105),ATLBL2(105),ATLBL3(105),ATLBL4(105)
      DIMENSION BLK(80)
      DIMENSION GH(4)
      DIMENSION ERRMSG(3)
      EQUIVALENCE (BLNK80,BLK(1))
      EQUIVALENCE (CHREND,WRDEND)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- ZGEO -'/
      DATA DIGIT  /'1' ,'2' ,'3' ,'4' ,'5' ,'6' ,'7' ,'8' ,'9' ,'0' /
      DATA UNDERS /'_'/
      DATA DASH   /'-'/
      DATA WRDBAS /' $BASIS '/
      DATA WRDEND /'zend    '/
      DATA WRDXYZ /' $XYZ   '/
      DATA WRDGEO /' $GEO   '/
      DATA WRDZMT /' $ZMT   '/
      DATA WRDIZM /' IZMAT ='/
      DATA WRDOPT /' $OPTZ  '/
      DATA WRDFRZ /' ICFRZ ='/
      DATA WR1VAR /'VARIABLES:'/
      DATA WR2VAR /'variables:'/
      DATA WR1CON /'CONSTANTS:'/
      DATA WR2CON /'constants:'/
      DATA DUMMY  /'   DUMMY'/
      DATA BLK    /80*' '/
      DATA BLNK8  /'        '/
      DATA BLNK5  /'     '/
      DATA BLNK2  /'  '/
      DATA GH     /'GH','Gh','gH','gh'/
      DATA ZERO   /0.0D+00/
      DATA ATLBL1 /'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',
     1             'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA',
     2             'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     3             'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR',
     4             'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN',
     5             'SB','TE','I ','XE','CS','BA','LA','CE','PR','ND',
     6             'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB',
     7             'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG',
     8             'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH',
     9             'PA','U ','NP','PU','AM','CM','BK','CF','ES','FM',
     1             'MD','NO','LR','X ','BQ'/
      DATA ATLBL2 /'h ','he','li','be','b ','c ','n ','o ','f ','ne',
     1             'na','mg','al','si','p ','s ','cl','ar','k ','ca',
     2             'sc','ti','v ','cr','mn','fe','co','ni','cu','zn',
     3             'ga','ge','as','se','br','kr','rb','sr','y ','zr',
     4             'nb','mo','tc','ru','rh','pd','ag','cd','in','sn',
     5             'sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     6             'pm','sm','eu','gd','tb','dy','ho','er','tm','yb',
     7             'lu','hf','ta','w ','re','os','ir','pt','au','hg',
     8             'tl','pb','bi','po','at','rn','fr','ra','ac','th',
     9             'pa','u ','np','pu','am','cm','bk','cf','es','fm',
     1             'md','no','lr','x ','bq'/
      DATA ATLBL3 /'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     1             'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     2             'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     3             'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     4             'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     5             'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     6             'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     7             'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     8             'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     9             'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     1             'Md','No','Lr','X ','Bq'/
      DATA ATLBL4 /'h ','hE','lI','bE','b ','c ','n ','o ','f ','nE',
     1             'nA','mG','aL','sI','p ','s ','cL','aR','k ','cA',
     2             'sC','tI','v ','cR','mN','fE','cO','nI','cU','zN',
     3             'gA','gE','aS','sE','bR','kR','rB','sR','y ','zR',
     4             'nB','mO','tC','rU','rH','pD','aG','cD','iN','sN',
     5             'sB','tE','i ','xE','cS','bA','lA','cE','pR','nD',
     6             'pM','sM','eU','gD','tB','dY','hO','eR','tM','yB',
     7             'lU','hF','tA','w ','rE','oS','iR','pT','aU','hG',
     8             'tL','pB','bI','pO','aT','rN','fR','rA','aC','tH',
     9             'pA','u ','nP','pU','aM','cM','bK','cF','eS','fM',
     1             'mD','nO','lR','x ','bQ'/
C
      DATA BIOSYM /'!BIOSYM archive 3'/
      DATA BIOEND /'end  '/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
      LST =.FALSE.
c
      title = ' '
C
C     ----- PROGRAM TO READ THE -Z- MATRIX IN FREE-FORMAT INPUT -----
C
      DO I=1,MAXGEO
         CART(I)  = .FALSE.
         CART0(I) = .FALSE.
         DO J=1,3
            FRZVAL(J,I) =.FALSE.
              ZVAL(J,I) = ZERO
              ZLST(J,I) = ZERO
         ENDDO
         DO J=1,5
            ZMT(J,I) = 0
         ENDDO
      ENDDO
      DO I=1,MAXVAR
         FRZVAR(I)=.FALSE.
      ENDDO
      DO I=1,MXCOOR
         IZFRZ(I)=.FALSE.
         ICFRZ(I)=0
      ENDDO
C
C     ----- CHECK FOR ALTERNATE INPUT FILE FOR COORDINATES -----
C
      NEWFIL=.FALSE.
      DO I=1,80
                          GEOFIL(I:I) =  ' '
         NEWFIL=NEWFIL.OR.GEOFIL(I:I).NE.' '
      ENDDO
      IF(NEWFIL) THEN
         IRSAV=IR
         IR=5
         CALL HND_GEOCLS(IR)
         CALL HND_GEOOPN(IR,GEOFIL)
         REWIND IR
         READ(IR,9999) WORD
C
C     ----- CHECK FOR BIOSYM INPUT -----
C
         IF(WORD(1:17).EQ.BIOSYM(1:17)) THEN
            READ(IR,9999)
            READ(IR,9999) TITLE
            READ(IR,9999)
C
            IAT = 0
    5       IAT = IAT + 1
            CALL HND_RDFREE(IR,STRING,IERR)
            IF(IERR.NE.0) THEN
               WRITE(IW,7775)
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            IF(NUMWRD.GT.1) THEN
               DO I=1,NUMWRD
                  PRSZMT(I,IAT)=PRSWRD(I)
                  FLGZMT(I,IAT)=FLGWRD(I)
                  ZMTCHR(I,IAT)=NUMCHR(I)-2
                  IF(DBUG) THEN
                    WRITE(IW,9992) I,NUMWRD,FLGWRD(I),NUMCHR(I),
     1                                                PRSWRD(I)
                  ENDIF
               ENDDO
               IF(FLGWRD(1).NE.3.OR.FLGWRD(2).NE.1.OR.
     1            FLGWRD(3).NE.1.OR.FLGWRD(4).NE.1    ) THEN
                  WRITE(IW,7774)
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
               CHAR5 =BLNK5
               SYMBOL=BLNK2
               ATMNAM=BLNK8
c              READ(PRSWRD(1),*) CHAR5
c              READ(PRSWRD(2),*) XXIAT
c              READ(PRSWRD(3),*) YYIAT
c              READ(PRSWRD(4),*) ZZIAT
c              READ(PRSWRD(5),*) CHAR12
c              READ(PRSWRD(6),*) CHAR2
c              READ(PRSWRD(7),*) CHAR2
c              READ(PRSWRD(8),*) SYMBOL
               call hnd_dparsc(prswrd(1),numchr(1),char5 , 5)
               call hnd_dparsr(prswrd(2),numchr(2),xxiat    )
               call hnd_dparsr(prswrd(3),numchr(3),yyiat    )
               call hnd_dparsr(prswrd(4),numchr(4),zziat    )
               call hnd_dparsc(prswrd(5),numchr(5),char12,12)
               call hnd_dparsc(prswrd(6),numchr(6),char2 , 2)
               call hnd_dparsc(prswrd(7),numchr(7),char2 , 2)
               call hnd_dparsc(prswrd(8),numchr(8),symbol, 2)
               ATMNAM(1:ZMTCHR(8,IAT))=SYMBOL(1:ZMTCHR(8,IAT))
               ATMNAM(ZMTCHR(8,IAT)+1:ZMTCHR(8,IAT)              +1)='_'
               ATMNAM(ZMTCHR(8,IAT)+2:ZMTCHR(8,IAT)+ZMTCHR(1,IAT)+1)=
     1                                  CHAR5(1:ZMTCHR(1,IAT))
               ZMTCHR(1,IAT)=ZMTCHR(1,IAT)+ZMTCHR(8,IAT)+1
               ATNAME(IAT)=ATMNAM
                   XX(IAT)=XXIAT
                   YY(IAT)=YYIAT
                   ZZ(IAT)=ZZIAT
               NUMZMT(IAT)=4
               NUMWRD     =4
               GO TO 5
            ELSEIF(NUMWRD.EQ.1) THEN
               IF(FLGWRD(1).NE.3.OR.NUMCHR(1).NE.5.OR.
     1            PRSWRD(1)(2:4).NE.BIOEND(1:3)) THEN
                  IF(DBUG) THEN
                    WRITE(IW,9992) I,NUMWRD,FLGWRD(I),NUMCHR(I),
     1                                                PRSWRD(I)
                  ENDIF
                  WRITE(IW,7774)
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
            ELSE
               WRITE(IW,7774)
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            NAT = IAT - 1
C
         ELSE
            WRITE(IW,7779) WORD
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(DBUG) THEN
            WRITE(IW,7777)
            WRITE(IW,9999) TITLE
            DO IAT=1,NAT
               WRITE(IW,7776) ATNAME(IAT)(1:8),XX(IAT),YY(IAT),ZZ(IAT)
            ENDDO
         ENDIF
         DBUG=.FALSE.
         CALL HND_GEOCLS(IR)
         IR=IRSAV
         GO TO 200
      ENDIF
C
C     ----- THIS IS THE -HONDO- INPUT -$GEO- -----
C
      IAT=0
      NAT=0
      IVAR=0
      NVAR=0
C
C     ----- IN THE -NWCHEM- CODE , WE ARE POSITIONED CORRECTLY -----
C
C     ----- READ THE LINES FOR ALL THE ATOMS   -----
C           THE END OF THE DEFINITIONS OF THE
C           ATOMS IS DETECTED VIA A BLANK LINE
C
      IAT = 0
  110 IAT = IAT + 1
      CALL HND_RDFREE(IR,STRING,IERR)
      IF(IERR.NE.0) THEN
         WRITE(IW,8886)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF((NUMWRD.GT.1).OR.
     1   (NUMWRD.EQ.1) .AND. .NOT. 
     2     ( (PRSWRD(1)(2:4).EQ.WR1VAR(1:3)) .OR.
     3       (PRSWRD(1)(2:4).EQ.WR2VAR(1:3)) .OR.
     4       (PRSWRD(1)(2:4).EQ.WR1CON(1:3)) .OR.
     5       (PRSWRD(1)(2:4).EQ.WR2CON(1:3)) )    ) THEN
         DO I=1,NUMWRD
            PRSZMT(I,IAT)=PRSWRD(I)
            FLGZMT(I,IAT)=FLGWRD(I)
            ZMTCHR(I,IAT)=NUMCHR(I)
            IF(DBUG) THEN
               WRITE(IW,9992) I,NUMWRD,FLGWRD(I),NUMCHR(I),PRSWRD(I)
            ENDIF
         ENDDO
         NUMZMT(  IAT)=NUMWRD
         IF(NUMWRD.GT.1) THEN
            GO TO 110
         ELSE
            IF((NUMCHR(1).EQ.6.AND.PRSWRD(1)(2:5).EQ.CHREND(1:4)).OR.
     1         (NUMCHR(1).EQ.5.AND.PRSWRD(1)(2:4).EQ.CHREND(2:4))) THEN
               NAT = IAT -1
               GO TO 140
            ELSE
               GO TO 110
            ENDIF
         ENDIF
      ENDIF
      NAT = IAT - 1
C
C     ----- READ THE LINES FOR ALL THE VARIABLES   -----
C           THE END OF THE DEFINITIONS OF THE
C           VARIABLES IS DETECTED VIA A BLANK LINE
C
      IVAR = 0
  120 IVAR = IVAR + 1
      CALL HND_RDFREE(IR,STRING,IERR)
      IF(IERR.NE.0) THEN
         WRITE(IW,8886)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF((NUMWRD.GT.1).OR.
     1   (NUMWRD.EQ.1) .AND. .NOT. (
     2       (PRSWRD(1)(2:4).EQ.WR1CON(1:3)) .OR.
     3       (PRSWRD(1)(2:4).EQ.WR2CON(1:3)) )    ) THEN
         DO I=1,NUMWRD
            PRSVAR(I,IVAR)=PRSWRD(I)
            FLGVAR(I,IVAR)=FLGWRD(I)
            VARCHR(I,IVAR)=NUMCHR(I)
            IF(DBUG) THEN
               WRITE(IW,9992) I,NUMWRD,FLGWRD(I),NUMCHR(I),PRSWRD(I)
            ENDIF
         ENDDO
         NUMVAR(  IVAR)=NUMWRD
         FRZVAR(  IVAR)=.FALSE.
         IF(NUMWRD.GT.1) THEN
            IF(NUMWRD.GT.2) THEN
               LST=.TRUE.
               IF(NUMWRD.GT.3) THEN
                  NUMWRD=3
               ENDIF
            ENDIF
            GO TO 120
         ELSE
            IF((NUMCHR(1).EQ.6.AND.PRSWRD(1)(2:5).EQ.CHREND(1:4)).OR.
     1         (NUMCHR(1).EQ.5.AND.PRSWRD(1)(2:4).EQ.CHREND(2:4))) THEN
               NVAR = IVAR -1
               GO TO 140
            ELSE
               GO TO 120
            ENDIF
         ENDIF
      ENDIF
      NVAR = IVAR - 1
C
C     ----- NOW READ VARIABLES THAT GET MARKED 'FROZEN' -----
C              ONLY AFTER THE 'VARIABLE' VARIABLES.
C
      IF(NVAR.GT.0) THEN
         IVAR = NVAR
  130    IVAR = IVAR + 1
         CALL HND_RDFREE(IR,STRING,IERR)
         IF(IERR.NE.0) THEN
            WRITE(IW,8886)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(NUMWRD.GT.1) THEN
            DO I=1,NUMWRD
               PRSVAR(I,IVAR)=PRSWRD(I)
               FLGVAR(I,IVAR)=FLGWRD(I)
               VARCHR(I,IVAR)=NUMCHR(I)
               IF(DBUG) THEN
                  WRITE(IW,9992) I,NUMWRD,FLGWRD(I),NUMCHR(I),PRSWRD(I)
               ENDIF
            ENDDO
               NUMVAR(  IVAR)=NUMWRD
               FRZVAR(  IVAR)=.TRUE.
            GO TO 130
         ELSEIF(NUMWRD.EQ.1) THEN
            IF((NUMCHR(1).EQ.6.AND.PRSWRD(1)(2:5).EQ.CHREND(1:4)).OR.
     1         (NUMCHR(1).EQ.5.AND.PRSWRD(1)(2:4).EQ.CHREND(2:4))) THEN
            ELSE
               WRITE(IW,8880)
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
         ENDIF
         NVAR = IVAR - 1
      ENDIF
C
  140 CONTINUE
C
C     ----- TAKE CARE OF -GHOST- ATOMS BY CHECKING THE LAST -----
C           PIECE OF DATA FOR EACH ATOM.
C
      DO IAT=1,NAT
         GHOST(IAT)=.FALSE.
         IWRD=NUMZMT(IAT)
         IFLG=FLGZMT(IWRD,IAT)
         ILEN=ZMTCHR(IWRD,IAT)
         WORD=PRSZMT(IWRD,IAT)
         IF(IFLG.EQ.3) THEN
            IF(ILEN.GE.4) THEN
               IF((WORD(2:3).EQ.GH(1)(1:2)).OR.
     1            (WORD(2:3).EQ.GH(2)(1:2)).OR.
     2            (WORD(2:3).EQ.GH(3)(1:2)).OR.
     3            (WORD(2:3).EQ.GH(4)(1:2))    ) THEN
                  GHOST(IAT)=.TRUE.
                  NUMZMT(IAT)=NUMZMT(IAT)-1
               ENDIF
            ENDIF
         ENDIF
      ENDDO
C
C     -----  NOW WE HAVE ALL THE VARIABLES AND THE -Z- MATRIX  -----
C            WITH SOME OF THE ELEMENTS EXPRESSED AS VARIABLES.
C           WE NOW NEED TO SUBSTITUTE VALUES FOR THE VARIABLES.
C
      IZ   =0
      NZMOD=0
      DO IAT=1,NAT
         IF(DBUG) THEN
            WRITE(IW,9982) IAT,NUMZMT(IAT)
            WRITE(IW,9981) (FLGZMT(I,IAT),I=1,NUMZMT(IAT))
            WRITE(IW,9980) (ZMTCHR(I,IAT),I=1,NUMZMT(IAT))
            WRITE(IW,9979) (PRSZMT(I,IAT),I=1,NUMZMT(IAT))
            WRITE(IW,9978) GHOST(IAT)
         ENDIF
         CALL HND_ZDAT(IAT,NAT,NVAR,ZVAL,ZLST,
     1             ZMT,NUMZMT,PRSZMT,FLGZMT,ZMTCHR,
     2             NUMVAR,PRSVAR,VARCHR,FRZVAR,FRZVAL,LST,
     3             IZMAT,IZ,IZFRZ,MXIZMT,NZMOD,DBUG,
     $        zvarname, zvarsign)
      ENDDO
c
      NCFRZ=0
      IF(NZMOD.GT.0) THEN
         DO IZMOD=1,NZMOD
            IF(IZFRZ(IZMOD)) THEN
               NCFRZ=NCFRZ+1
               ICFRZ(NCFRZ)=IZMOD
            ENDIF
         ENDDO
      ENDIF
C
C     ----- IF SOME ATOMS ARE GIVEN IN CARTESIAN COORDINATES -----
C           THEN THE -$ZMAT- DATA, AUTOMATICALLY CREATED
C           WILL BE IGNORED .
C
      DO IAT=1,NAT
         IF(CART(IAT)) THEN
            IZ   =0
            NZMOD=0
            NCFRZ=0
         ENDIF
      ENDDO
C
C     ----- WE HAVE THE -Z- MATRIX WITH THE INDICES IN -ZMT- -----
C           AND THE VALUES IN -ZVAL-
C
      IF(DBUG) THEN
         WRITE(IW,9998) TITLE
         WRITE(IW,9997)
         DO IAT=1,NAT
            IF(CART(IAT)) THEN
               WRITE(IW,9986) IAT,ATNAME(IAT),XX(IAT),YY(IAT),ZZ(IAT)
            ELSE
               WRITE(IW,9974) IAT,ATNAME(IAT)
               WRITE(IW,9996) IAT,(   ZMT(J,IAT),J=1,4),
     1                            (  ZVAL(J,IAT),J=1,3),
     2                            (FRZVAL(J,IAT),J=1,3)
            ENDIF
         ENDDO
         IF(LST) THEN
            WRITE(IW,9998) TITLE
            WRITE(IW,9997)
            DO IAT=1,NAT
               IF(CART(IAT)) THEN
                  WRITE(IW,9986) IAT,ATNAME(IAT),XX(IAT),YY(IAT),ZZ(IAT)
               ELSE
                  WRITE(IW,9996) IAT,(   ZMT(J,IAT),J=1,4),
     1                               (  ZLST(J,IAT),J=1,3),
     2                               (FRZVAL(J,IAT),J=1,3)
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      IF(OUT) THEN
         IF(IZ.GT.0) THEN
            WRITE(IW,9988)
            WRITE(IW,9991) WRDZMT,NZMOD,WRDIZM
            WRITE(IW,9990) (IZMAT(I),I=1,IZ)
            WRITE(IW,9989) WRDEND
         ENDIF
         IF(NCFRZ.GT.0) THEN
            WRITE(IW,9988)
            WRITE(IW,9987) WRDOPT,WRDFRZ
            WRITE(IW,9990) (ICFRZ(I),I=1,NCFRZ)
            WRITE(IW,9989) WRDEND
         ENDIF
      ENDIF
C
C     ----- CALCULATE CARTESIAN COORDINATES -----
C
      IF(LST) THEN
         DO IAT=1,NAT
            CART0(IAT)=CART(IAT)
         ENDDO
      ENDIF
      CALL HND_ZXYZ(NAT,ZMT,ZVAL)
      IF(OUT.OR.DBUG) THEN
         WRITE(IW,9995) TITLE
         WRITE(IW,9994)
         DO IAT=1,NAT
            WRITE(IW,9993) ATNAME(IAT),XX(IAT),YY(IAT),ZZ(IAT)
         ENDDO
      ENDIF
C
      IF(LST) THEN
         DO IAT=1,NAT
            DO J=1,3
               ZSTP(J,IAT)=(ZLST(J,IAT)-ZVAL(J,IAT))/DBLE(MAXLST-1)
               ZVAL(J,IAT)= ZVAL(J,IAT)-ZSTP(J,IAT)
            ENDDO
         ENDDO
         DO ILST=1,MAXLST
            IF(DBUG) THEN
               WRITE(IW,9976) ILST
            ENDIF
            DO IAT=1,NAT
               DO J=1,3
                  ZVAL(J,IAT)=ZVAL(J,IAT)+ZSTP(J,IAT)
               ENDDO
               CART(IAT)=CART0(IAT)
            ENDDO
            CALL HND_ZXYZ(NAT,ZMT,ZVAL)
            IF(DBUG) THEN
               WRITE(IW,9995) TITLE
               WRITE(IW,9994)
               DO IAT=1,NAT
                  WRITE(IW,9993) ATNAME(IAT),XX(IAT),YY(IAT),ZZ(IAT)
               ENDDO
            ENDIF
            DO IAT=1,NAT
               XXLST(IAT,   2)=XX(IAT)*UNITS
               YYLST(IAT,   2)=YY(IAT)*UNITS
               ZZLST(IAT,   2)=ZZ(IAT)*UNITS
            ENDDO
            STATUS=GEOM_LST_PUT_COORD(XXLST(1,2),
     1                                YYLST(1,2),
     2                                ZZLST(1,2),NAT)
            IF(.NOT.STATUS) THEN
               CALL ERRQUIT('GEOM_LST_PUT_COORD : PUT FAILED',911)
            ENDIF
            IF(ILST.EQ.1) THEN
               DO IAT=1,NAT
                  XXLST(IAT,   1)=XX(IAT)
                  YYLST(IAT,   1)=YY(IAT)
                  ZZLST(IAT,   1)=ZZ(IAT)
               ENDDO
            ENDIF
         ENDDO
         DO IAT=1,NAT
            XX(IAT)=XXLST(IAT,1)
            YY(IAT)=YYLST(IAT,1)
            ZZ(IAT)=ZZLST(IAT,1)
         ENDDO
      ENDIF
C
  200 CONTINUE
C
C     ----- EXTRACT ATOMIC NUMBER FROM CHEMICAL SYMBOL -----
C
      DO IAT=1,NAT
         ATNUM(IAT)=ZERO
         NCHR=0
  210    NCHR=NCHR+1
         IF(NCHR.GT.ZMTCHR(1,IAT)) GO TO 220
C
         IF(ATNAME(IAT)(NCHR:NCHR).EQ.UNDERS.OR.
     1      ATNAME(IAT)(NCHR:NCHR).EQ.DASH      ) GO TO 220
         DO IDIGIT=1,10
            IF(ICHAR(ATNAME(IAT)(NCHR:NCHR)).EQ.ICHAR(DIGIT(IDIGIT)))
     1                                            GO TO 220
         ENDDO
         GO TO 210
  220    CONTINUE
         NCHR=NCHR-1
         IF(NCHR.EQ.1) THEN
            DO ISYMBL=1,105
               IF((ATLBL1(ISYMBL)(1:1).EQ.ATNAME(IAT)(1:1).AND.
     1             ATLBL1(ISYMBL)(2:2).EQ.' '                  ).OR.
     2            (ATLBL2(ISYMBL)(1:1).EQ.ATNAME(IAT)(1:1).AND.
     3             ATLBL2(ISYMBL)(2:2).EQ.' '                  ).OR.
     4            (ATLBL3(ISYMBL)(1:1).EQ.ATNAME(IAT)(1:1).AND.
     5             ATLBL3(ISYMBL)(2:2).EQ.' '                  ).OR.
     6            (ATLBL4(ISYMBL)(1:1).EQ.ATNAME(IAT)(1:1).AND.
     7             ATLBL4(ISYMBL)(2:2).EQ.' '                  )) THEN
                  ATNUM(IAT)=DBLE(ISYMBL)
                  IF(GHOST(IAT)) THEN
                     ATNUM(IAT)=-ATNUM(IAT)
                  ENDIF
               ENDIF
            ENDDO
         ELSEIF(NCHR.EQ.2) THEN
            DO ISYMBL=1,105
               IF((ATLBL1(ISYMBL)(1:1).EQ.ATNAME(IAT)(1:1).AND.
     1             ATLBL1(ISYMBL)(2:2).EQ.ATNAME(IAT)(2:2)     ).OR.
     2            (ATLBL2(ISYMBL)(1:1).EQ.ATNAME(IAT)(1:1).AND.
     3             ATLBL2(ISYMBL)(2:2).EQ.ATNAME(IAT)(2:2)     ).OR.
     4            (ATLBL3(ISYMBL)(1:1).EQ.ATNAME(IAT)(1:1).AND.
     5             ATLBL3(ISYMBL)(2:2).EQ.ATNAME(IAT)(2:2)     ).OR.
     6            (ATLBL4(ISYMBL)(1:1).EQ.ATNAME(IAT)(1:1).AND.
     7             ATLBL4(ISYMBL)(2:2).EQ.ATNAME(IAT)(2:2)     )) THEN
                  ATNUM(IAT)=DBLE(ISYMBL)
                  IF(GHOST(IAT)) THEN
                     ATNUM(IAT)=-ATNUM(IAT)
                  ENDIF
               ENDIF
            ENDDO
         ELSE
            WRITE(IW,8883) IAT,NCHR
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(ATNUM(IAT).EQ.ZERO) THEN
            WRITE(IW,8879) IAT
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(ATNUM(IAT).GT.DBLE(103)) THEN
            ATNUM(IAT)=ZERO
         ENDIF
      ENDDO
C
C     ----- -LST- GEOMETRIES IF MADE UP -----
C
      IF(LST.AND.OUT) THEN
         WRITE(IW,9975) MAXLST
         DO ILST=1,MAXLST
            STATUS=GEOM_LST_GET_COORD(XXLST(1,2),
     1                                YYLST(1,2),ZZLST(1,2),NAT,ILST)
            IF(.NOT.STATUS) THEN
               CALL ERRQUIT('GEOM_LST_GET_COORD : GET FAILED',911)
            ENDIF
            WRITE(IW,8888) WRDXYZ
            DO IAT=1,NAT
               WRITE(IW,8884) ATNAME(IAT),ATNUM(IAT),
     1                         XXLST(IAT,   2),YYLST(IAT,   2),
     2                         ZZLST(IAT,   2)
            ENDDO
            WRITE(IW,8888) WRDEND
         ENDDO
      ENDIF
C
C     ----- CREATE DATA FOR -NWCHEM- -----
C
      NCENTER=NAT
      DO ICENTER=1,NCENTER
         COORDS(1,ICENTER)=   XX(ICENTER)
         COORDS(2,ICENTER)=   YY(ICENTER)
         COORDS(3,ICENTER)=   ZZ(ICENTER)
         CHARGE(  ICENTER)=ATNUM(ICENTER)
           TAGS(  ICENTER)(1:16)=' '
           TAGS(  ICENTER)(1: 2)=ATNAME(ICENTER)(1:2)
        IF(TAGS(  ICENTER)(2: 2).EQ.'_') THEN         
           TAGS(  ICENTER)(2: 2)=' '                  
        ENDIF
      ENDDO
C
      RETURN
 9999 FORMAT(A80)
 9998 FORMAT(' -Z- MATRIX DEFINITION FOR MOLECULE:',/,1X,A80,/)
 9997 FORMAT(' ATOM # ',' I ','   J   ',' K ','   L   ','   DIST    ','
     1       ANGLE   ','      TORSION   ',/)
 9996 FORMAT(1X,I4,3X,I3,2X,I3,2X,I3,2X,I3,4X,F8.4,7X,F8.3,7X,F8.3,
     1                                           4X,L1,3X,L1,3X,L1)
 9995 FORMAT(/,' CARTESIAN COORDINATES FOR MOLECULE:',/,1X,A80,/)
 9994 FORMAT(' NAME  ',13X,'X',18X,'Y',18X,'Z',//)
 9993 FORMAT(1X,A5,3X,3(F15.5,4X))
 9992 FORMAT(4I4,4X,A40)
 9991 FORMAT(A8,' NZMOD =',I4,1H,,A8)
 9990 FORMAT(12(I4,1H,))
 9989 FORMAT(A8)
 9988 FORMAT(/)
 9987 FORMAT(2A8)
 9986 FORMAT(1X,I5,1X,A5,3X,3(F15.5,4X))
 9985 FORMAT(' INCORRECT -$BAS- INPUT DATA. -NBASIS- AND -NAT-',
     1       ' MUST BE EQUAL. STOP')
 9984 FORMAT(' INCORRECT -$BAS- INPUT DATA. -NBASIS- = ',I4)
 9983 FORMAT(' THIS BASIS SET SPECIFICATION OF -GLOBAL- AND -INTERNAL-',
     1       ' MAY NOT BE USED WHEN -DUMMY- ATOMS ARE PRESENT.',/,
     2       ' USE AN ATOM-BY-ATOM SPECIFICATION IN -$BAS-',
     3       ' INCLUDING THE -DUMMY- SPECIFICATION. STOP')
 9982 FORMAT(' CALLING -HND_ZDAT- WITH -IAT, NUMZMT(IAT)- = ',2I5)
 9981 FORMAT(' -FLGZMT( ,IAT) = ',10I4)
 9980 FORMAT(' -ZMTCHR( ,IAT) = ',10I4)
 9979 FORMAT(' -PRSZMT( ,IAT) = --- ',A80)
 9978 FORMAT(' GHOST ATOM ? = ',L4)
 9977 FORMAT(' SOMETHING IS NOT RIGHT WITH -$BAS- :',/,
     1       ' AN ATOM WITH NON-DUMMY ATOMIC NUMBER IS GIVEN AN',
     2       ' EMPTY -DUMMY- BASIS SET. STOP. -IAT- = ',I5)
 9976 FORMAT(' ----- -LST- POINT NO. = ',I3,' -----')
 9975 FORMAT(' $PES     NPES =',I4,', IUNIT =1, $END')
 9974 FORMAT(' IAT=',I5,' ATNAME=',2X,A8)
 8889 FORMAT(/,10X,12(1H-),/,10X,'-$GEO- INPUT',/,10X,12(1H-))
 8888 FORMAT(A8)
 8887 FORMAT(' NO DATA GROUP -$GEO- FOUND AS ALTERNATE INPUT.     ')
 8886 FORMAT(' END-OF-FILE ENCOUNTERED WHILE READING -$GEO- . STOP')
 8885 FORMAT(A80)
 8884 FORMAT(A8,F4.0,3F15.7)
 8883 FORMAT(' ATOMIC SYMBOL FOR IAT = ',I4,' HAS MORE THAN 2',
     1       ' CHARACTERS. NCHR = ',I3,' STOP.',/,
     2       ' THE FIRST ONE OR TWO CHARACTERS OF THE ATOM NAME',
     3       ' MUST FORM THE ATOMIC SYMBOL.',/,' BEYOND THAT THE',
     4       ' OTHER CHARACTERS MUST BE DIGITS.')
 8882 FORMAT(' ATOM -IAT- = ',I4,' -',A8,'- HAS BASIS SET # ',I4,
     1       ' -',A8,'- .')
 8881 FORMAT(' -EXTNAL- = ',L4,' -INTNAL- = ',L4,' -GLOBAL- = ',L4,
     1     /,I5,' ATOMS AND',I5,' INDIVIDUAL BASIS SET(S) LISTED.',/,
     2          ' THERE OUGHT TO BE AT LEAST AS MANY BASIS SETS',
     3          ' AS ATOMS. IF NOT, STOP .')
 8880 FORMAT(' THE END OF THE -$GEO- DATA GROUP IS INCORRECTLY',
     1       ' INDICATED. STOP',/,' IT SHOULD BE A BLANK CARD OR',
     2       ' A - ZEND - CARD.')
 8879 FORMAT(' THE ATOMIC SYMBOL FOR IAT = ',I3,' IS NOT LEGAL. STOP.',
     1     /,' THE ATOMIC SYMBOL IS EXTRACTED FROM THE FIRST TWO',
     2       ' CHARACTERS OF THE ATOM NAME.')
 7779 FORMAT(' THE ALTERNATE FILE FOR COORDINATES IS NOT A',
     1       ' -BIOSYM.CAR- FILE. THE FIRST LINE IS = ',/,
     2       1X,3H---,A80,3H---,/,1X,'STOP.')
 7778 FORMAT(A5,3F15.9)
 7777 FORMAT(1X,'ATOMIC COORDINATES READ FROM -BIOSYM.CAR- FILE =',/,
     1       1X,'------------------------------------------------')
 7776 FORMAT(1X,A8,3F15.9)
 7775 FORMAT(' END-OF-FILE ENCOUNTERED WHILE READING -BIOSYM.CAR- .',
     1       ' STOP .')
 7774 FORMAT(' ERROR WHILE READING -BIOSYM.CAR- . STOP .')
      END
      logical function geom_lst_put_coord(x,y,z,n)        
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
c
      integer n
      integer i
      integer ir
      integer iw
      integer ift
      integer irc
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision x(*),y(*),z(*) 
      character*255 lst_coord_fil
      data ift /32/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      call util_file_name('lst.coord',
     1     .false.,.false.,lst_coord_fil)
      open(unit=ift, file=lst_coord_fil, form='unformatted',
     1     access='sequential', status='unknown', err=911)
c
      if(out) then
         write(iw,*) 'lst_put_coord = ',lst_coord_fil
         write(iw,*) 'n = ',n
      endif
      if(dbug) then
         call hnd_prsq(x,1,n,n)            
         call hnd_prsq(y,1,n,n)            
         call hnd_prsq(z,1,n,n)            
      endif
c
      rewind ift
      irc=0
   10 read(ift,end=20,err=20)
         irc=irc+1
         go to 10
   20 continue
      if(out) then
         write(iw,*) irc,' records found on -lst_coord_fil- '
      endif
      rewind ift
      if(irc.gt.0) then
         do i=1,irc
            read(ift)
         enddo
      endif
c
      write(ift) (x(i),i=1,n),
     &           (y(i),i=1,n),                  
     &           (z(i),i=1,n)                   
c
      close(ift,status='keep')
c
      geom_lst_put_coord=.true.
      return
c
  911 call errquit('geom_lst_put_coord : open failed',0)
c
      end
      logical function geom_lst_get_coord(x,y,z,n,irc)      
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
c
      integer n
      integer i
      integer ir
      integer iw
      integer ift
      integer irc
      logical dbug
      logical out
      common/hnd_iofile/ir,iw
      double precision x(*),y(*),z(*) 
      character*255 lst_coord_fil
      data ift /32/
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
c
      call util_file_name('lst.coord',
     1     .false.,.false.,lst_coord_fil)
      open(unit=ift, file=lst_coord_fil, form='unformatted',
     1     access='sequential', status='unknown', err=911)
      if(out) then
         write(iw,*) 'lst_coord_fil = ', lst_coord_fil
         write(iw,*) 'n = ',n
         write(iw,*) 'irc = ',irc
      endif
c
      rewind ift
      if(irc.gt.1) then
         do i=1,irc-1
            read(ift)
         enddo
      endif
      read(ift,end=910,err=910) (x(i),i=1,n),
     &                          (y(i),i=1,n),
     &                          (z(i),i=1,n)
c
      close(ift,status='keep')
c
      if(dbug) then
         call hnd_prsq(x,1,n,n)         
         call hnd_prsq(y,1,n,n)         
         call hnd_prsq(z,1,n,n)         
      endif
c
      geom_lst_get_coord=.true.
      return
c
c     ----- unable to complete the read ---
c
 910  close(ift,status='keep')
      if(out) then
         write(iw,*) 'unable to complete the read in geom_lst_get_coord'
      endif
      geom_lst_get_coord=.false.
      return
c
c     ----- unable to open ; stop ... -----
c
 911  call errquit('geom_lst_get_coord : open failed',0)
      return
c
      end
      subroutine geom_check_input_quants(
     &    tags,vector,ncenter,diff_thresh,identity)
      implicit none
#include "stdio.fh"
#include "inp.fh"
      integer ncenter
      character*(*) identity
      character*16 tags(ncenter)
      double precision vector(ncenter)
      double precision diff_thresh
c
      integer i,j,k, icent, jcent
      logical status
      double precision diff
c
      status = .true.
      do i = 2,ncenter
        do j = 1,(i-1)
          if (tags(i).eq.tags(j))
     &        status = status .and.
     &        (abs(vector(i)-vector(j)).lt.diff_thresh)
          if (.not.status) then
            icent = i
            jcent = j
            diff = abs(vector(i)-vector(j))
            goto 00010
          endif
        enddo
      enddo
      return
00010 continue
c
      k = inp_strlen(identity)
      write(luout,*)
     &    ' geom_input: error mis-match on quantities of ',
     &    identity(1:k)
      write(luout,*)
     &    '           : tags must be different for different values',
     &    ' of "',identity(1:k),'"'
      write(luout,00001)
c
      do i = 1,ncenter
        j = inp_strlen(tags(i))
        write(luout,00002)i,tags(i)(1:j),identity(1:k),vector(i)
      enddo
      write(luout,*)
      write(luout,00003)
     &    icent,tags(icent),identity(1:k),vector(icent),
     &    jcent,tags(jcent),identity(1:k),vector(jcent),
     &    diff,diff_thresh
      write(luout,00001)
      call errquit('geom_check_input_quants: fatal error',911)
c
00001 format(1x,80('-'),/,/)
00002 format(1x,'center',i4,1x,a16,' has a ',a,' of ',f10.4)
00003 format(1x,'detected error occured with the input for:',/,
     &    1x,'center',i4,1x,a16,' has a ',a,' of ',f20.13,/,
     &    1x,'center',i4,1x,a16,' has a ',a,' of ',f20.13,/,
     &    1x,' differ by ',31x,f20.13,/,
     &    1x,' which is larger than the threshold:',6x,f20.13,/,/)
      end
      subroutine  geom_auto_sym(rtdb,geom,
     1                          coords,charge,tags,ncenter,threquiv,
     $     group)
      implicit none
#include "stdio.fh"
      integer      rtdb
      integer      geom
      integer      mxatom
      parameter    (mxatom=500)
      integer      ir
      integer      iw
      integer      i
      integer      ncenter
      integer      nuc
      integer      nat
      logical      dbug 
      logical      some 
      double precision threquiv
      double precision coords
      double precision charge
      double precision c
      double precision zan
      character*16 tags
      character*16 atmlab
      character*8  errmsg,groupname
      character*(*) group
      common/hnd_iofile/ir,iw
      common/hnd_molnuc/nuc(mxatom)
      common/hnd_mollab/atmlab(mxatom)
      common/hnd_molxyz/c(3,mxatom),zan(mxatom),nat
      dimension coords(3,*)
      dimension charge(  *)
      dimension   tags(  *)
      dimension errmsg(3)
      logical odone
      data errmsg /'program ','stop in ','-autsym-'/
c
      dbug=.false.
      some=.false.
      some=some.or.dbug
c
      ir=LuIN 
      iw=LuOut 
      if(some) then
         write(iw,9999)
         if(dbug) then
            do i=1,nat    
               write(iw,9998) i,charge(i),
     1                        coords(1,i),coords(2,i),coords(3,i),
     2                        tags(i)
            enddo
         endif
      endif
c
      nat=ncenter
      do i=1,ncenter
         atmlab(i)=tags(i)
         zan(i)=     charge(i) 
         nuc(i)=nint(charge(i))
         c(1,i)=coords(1,i)
         c(2,i)=coords(2,i)
         c(3,i)=coords(3,i)
      enddo
      if(dbug) then
         do i=1,nat    
            write(iw,9997) i,nuc(i),c(1,i),c(2,i),c(3,i),zan(i)
         enddo
      endif
c
      odone=.false.
c
c     ----- call auto_sym -----
c
      call hnd_autsym(odone,rtdb,threquiv,groupname)
c
      group = groupname
c
      do i=1,ncenter
         coords(1,i)=c(1,i)
         coords(2,i)=c(2,i)
         coords(3,i)=c(3,i)
      enddo
c
      return
 9999 format(' in -geom_auto_sym- ')
 9998 format(1x,i5,f5.1,3f10.6,2x,a16)
 9997 format(1x,i5,i5,3f10.6,f10.3)
c
      end
      SUBROUTINE HND_AUTSYM(ODONE,RTDB,THREQUIV,groupname)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER     RTDB
      LOGICAL     ODONE
      PARAMETER   (MXATOM=500)
      PARAMETER   (MXSYM =120)
      CHARACTER*8 ERRMSG,groupname
      CHARACTER*1 XYZNAM
      CHARACTER*1 ABCNAM
      LOGICAL     SOME
      LOGICAL     OUT
      LOGICAL     DBUG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_SYMTRY/INVT(MXSYM),NT,NTMAX,NTWD,NOSYM
      DIMENSION   CI(3),AI(3,3)
      DIMENSION   C1(3,MXATOM)
      DIMENSION   C2(3,MXATOM)
      DIMENSION   C3(3,MXATOM)
      DIMENSION   AXS(3,3),EIG(3)
      DIMENSION   RT(3,3)
      DIMENSION   TR(3)
      DIMENSION   XYZNAM(3),ABCNAM(3)
      DIMENSION   ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-AUTSYM-'/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA XYZNAM /'X','Y','Z'/
      DATA ABCNAM /'A','B','C'/
C
      DBUG=.FALSE.
      OUT =.FALSE. 
      OUT =OUT.OR.DBUG
      SOME=.FALSE.  
      SOME=SOME.OR.OUT
C
      NT=1
      NTMAX = 1                 ! Never seems to be set
      IF(NT.GT.1.OR.NTMAX.GT.1) THEN
         ODONE=.TRUE.
         RETURN
      ELSE
         ODONE=.FALSE.
      ENDIF
C
C     ----- AUTOMATIC DETECTION OF SYMMETRY -----
C
      IF(SOME) THEN
         WRITE(IW,9997)
         IF(OUT) THEN
            WRITE(IW,*) 'NT,NTMAX,NTWD,NOSYM = ',
     1                   NT,NTMAX,NTWD,NOSYM
         ENDIF
      ENDIF
C
C     ----- GET MOMENTS OF INERTIA -----                
C
      CALL HND_MOLINR(C,C1,NAT,CI,AI,TR)
C
C     ----- GET PRINCIPAL AXES -----
C
      CALL HND_MOLAXS(AI,AXS,EIG,3,3,3)
      IF(OUT) THEN
         WRITE(IW,9999)
         CALL HND_PRSQ(AXS,3,3,3)
      ENDIF
C
C     ----- IDENTIFY SYMMETRY OPERATIONS -----
C
      CALL HND_MOLOPS(C,C1,C2,C3,NAT,AXS,EIG,RT,TR,ODONE,THREQUIV,
     $     groupname)
C
C     ----- CREATE ATOM AND SHELL MAPPINGS -----
C
      CALL HND_MOLMAP(C,C3,NAT)
C
      IF(SOME) THEN
         IF(OUT) THEN
            WRITE(IW,*) 'NT,NTMAX,NTWD,NOSYM = ',
     1                   NT,NTMAX,NTWD,NOSYM
         ENDIF
         WRITE(IW,9996)
      ENDIF
C
C     ----- SYMMETRIZE CARTESIAN COORDINATES -----
C
      CALL HND_SYMATM
c
      do i = 1, nat
         c(1,i) = c3(1,i)
         c(2,i) = c3(2,i)
         c(3,i) = c3(3,i)
      enddo
C
      RETURN
 9999 FORMAT(/,1X,'PRINCIPAL AXES OF INERTIA',/,1X,25(1H-))
 9998 FORMAT(/,1X,'ANGLES ASSOCIATED WITH AXES',/,1X,27(1H-))
 9997 FORMAT(/,10X,7(1H-),                             
     1       /,10X,'autosym', 
     2       /,10X,7(1H-))
 9996 FORMAT(/)
      END
      SUBROUTINE HND_MOLOPS(C,C1,C2,C3,NAT,AXS,EIG,RT,TR,ODONE,
     $     THREQUIV,groupname)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- ROUTINE DETECTS MOLECULAR SYMMETRY OPERATIONS -----
C           CODE ADAPTED FROM A.M.CHAKA's ORIGINAL CODE.
c
c     RJH ... added parameter THREQUIV to provide user control
c     .       for noisy geometries, and groupname to return value
C
      LOGICAL     ODONE
      PARAMETER   (MXORDR=24)
      PARAMETER   (MXATOM=500)
      CHARACTER*8 ERRMSG,groupname
      LOGICAL     MUCH
      LOGICAL     DBUG
      LOGICAL     OUT
      LOGICAL     SOME
      LOGICAL     ONLY
      LOGICAL     GOTX,GOTY,GOTZ
      LOGICAL     NOSYM
      LOGICAL     PROPER
      LOGICAL     PRPAXS
      LOGICAL     IMPROP
      LOGICAL     IMPAXS
      LOGICAL     SYMINV
      LOGICAL     INVERS
      LOGICAL     SYMC2X
      LOGICAL     SYMC2Y
      LOGICAL     MIRRYZ
      LOGICAL     MIRRZX
      LOGICAL     MIRRXY
      LOGICAL     MIRROR
      LOGICAL     ATOMIC
      LOGICAL     LINEAR
      LOGICAL     DEGNR2
      LOGICAL     DEGNR3
      LOGICAL     CUBIC
      LOGICAL     C2ROT
      LOGICAL     C2AXS
      LOGICAL     C4ROT
      LOGICAL     C4AXS
      LOGICAL     S4ROT
      LOGICAL     S4AXS
      LOGICAL     GRPOH
      LOGICAL     GRPTH
      LOGICAL     GRPTD
      LOGICAL     GRPT
      LOGICAL     GRPO
      INTEGER     AXORDR 
      COMPLEX*16  QX,QY,QZ
      COMPLEX*16  QZERO
      COMPLEX*16  QDUMX,QDUMY,QDUMZ
      COMPLEX*16  QDUMI,QDUMJ
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLNUC/NUC(MXATOM)
      DIMENSION IEQU(MXATOM)
      DIMENSION   QX(MXORDR),  QY(MXORDR),  QZ(MXORDR)
      DIMENSION XMAG(MXORDR),YMAG(MXORDR),ZMAG(MXORDR)
      DIMENSION PROPER(3)
      DIMENSION IMPROP(3)
      DIMENSION AXORDR(3)
      DIMENSION NEWAXS(3)
      DIMENSION NUORDR(3)
      DIMENSION AXS(3,3),EIG(3)
      DIMENSION  RT(3,3), TR(3)
      DIMENSION PRM(3,3)
      DIMENSION  C(3,*)
      DIMENSION C1(3,*)
      DIMENSION C2(3,*)
      DIMENSION C3(3,*)
      DIMENSION CM(3),AM(3,3),TM(3)
      DIMENSION AXM(3,3)
      DIMENSION C2AXS(3)
      DIMENSION C4AXS(3)
      DIMENSION S4AXS(3)
      DIMENSION ERRMSG(3)
      EQUIVALENCE (DEGNR3,CUBIC)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-MOLOPS-'/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
      DATA FOUR   /4.0D+00/
      DATA TENM02 /1.0D-02/
      DATA TENM04 /1.0D-04/
      DATA TENM05 /1.0D-05/
      DATA TEN06  /1.0D+06/
      DATA EPS    /1.01D+00/
C
      QZERO=cmplx(ZERO,ZERO)
         PI=    FOUR*ATAN(ONE)
      TWOPI=TWO*FOUR*ATAN(ONE)
      DEGREE=360.00D+00/TWOPI
c
****      THREQUIV = TENM04
C
      MUCH=.FALSE.
      DBUG=.FALSE. 
      DBUG=DBUG.OR.MUCH
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      SOME=.FALSE.
      SOME=SOME.OR.OUT
      ONLY=.FALSE.
      ONLY=ONLY.OR.SOME
C
C     ----- INERTIA AXES AND DATA -----
C
      IF(OUT) THEN
         WRITE(IW,9986)
         WRITE(IW,9978) (I,EIG(I),I=1,3)
         CALL HND_PRSQ(AXS,3,3,3)
      ENDIF
C
C     ----- CHECK FOR ATOM OR LINEAR MOLECULE -----
C
      NZER=0
      DO I=1,3
         IF(ABS(EIG(I)).LT.TENM04) NZER=NZER+1
      ENDDO
      IF(NZER.EQ.0) THEN
         ATOMIC=.FALSE.
         LINEAR=.FALSE.
      ELSEIF(NZER.EQ.1) THEN
         ATOMIC=.FALSE.
         LINEAR=.TRUE.
      ELSEIF(NZER.GT.1) THEN
         ATOMIC=.TRUE.
         LINEAR=.FALSE.
      ENDIF
      IF(SOME) THEN
         IF(ATOMIC) WRITE(IW,*) 'THIS IS AN ATOM'
         IF(LINEAR) WRITE(IW,*) 'THIS IS A LINEAR MOLECULE'
      ENDIF
      IF(ATOMIC) THEN
         RETURN
      ENDIF
C 
C     ----- CHECK FOR DEGENERACY -----
C
      DEGNR2=(ABS(EIG(2)-EIG(1)).LT.TENM04).OR.
     1       (ABS(EIG(2)-EIG(3)).LT.TENM04)
      DEGNR3=(ABS(EIG(2)-EIG(1)).LT.TENM04).AND.
     1       (ABS(EIG(2)-EIG(3)).LT.TENM04)
C
      IF(DEGNR3) THEN
C
C     ----- CUBIC GROUPS -----
C
         DEGNR2=.FALSE.
         IF(SOME) THEN
            WRITE(IW,*) '3-FOLD DEGENERATE MOMENTS ... CUBIC GROUP'
         ENDIF
         GRPTH=.FALSE.
         GRPOH=.FALSE.
         GRPTD=.FALSE.
         GRPT =.FALSE.
         GRPO =.FALSE.
C
C     ----- DETERMINE PRESENCE OF CENTER OF INVERSION -----
C                   Th AND Oh GROUPS
C
         SYMINV=.TRUE.
         DO IAT=1,NAT
            INVERS=.FALSE.
            DO JAT=1,NAT
               IF(NUC(JAT).EQ.NUC(IAT).AND..NOT.INVERS) THEN
                  DX=C1(1,IAT)+C1(1,JAT)
                  DY=C1(2,IAT)+C1(2,JAT)
                  DZ=C1(3,IAT)+C1(3,JAT)
                  INVERS=(ABS(DX).LT.THREQUIV).AND.
     1                   (ABS(DY).LT.THREQUIV).AND.
     2                   (ABS(DZ).LT.THREQUIV)
                  IF(DBUG) THEN
                     WRITE(IW,*) ' IAT,JAT,INVERS = ',
     1                             IAT,JAT,INVERS,DX,DY,DZ
                  ENDIF
               ENDIF
            ENDDO
            SYMINV=SYMINV.AND.INVERS
            IF(DBUG) THEN
               WRITE(IW,*) ' IAT,INVERS,SYMINV = ',IAT,INVERS,SYMINV
            ENDIF
         ENDDO
         IF(SOME) THEN
            WRITE(IW,9993) SYMINV
            IF(SYMINV) THEN
               WRITE(IW,*) ' Th OR Oh GROUP '
            ELSE
               WRITE(IW,*) ' T  OR Td OR O GROUP '
            ENDIF
         ENDIF
C
C     ----- DETERMINE ATOM EQUIVALENCIES -----
C
         DO IAT=1,NAT
            IEQU(IAT)=IAT
         ENDDO
         DO IAT=1,NAT
            IF(IEQU(IAT).EQ.IAT) THEN
               DUM=C1(1,IAT)**2+C1(2,IAT)**2+C1(3,IAT)**2
               DISTI=SQRT(DUM)
               IF(DISTI.GT.THREQUIV) THEN
                  DO JAT=1,NAT
                     IF(NUC(JAT).EQ.NUC(IAT).AND.JAT.NE.IAT) THEN
                        DUM=C1(1,JAT)**2+C1(2,JAT)**2+C1(3,JAT)**2
                        DISTJ=SQRT(DUM)
                        IF(ABS(DISTJ-DISTI).LT.THREQUIV) THEN
                           IEQU(JAT)=IAT
                        ENDIF
                     ENDIF
                  ENDDO
               ELSE
                  IEQU(IAT)=0
               ENDIF
            ENDIF
         ENDDO
         IF(OUT) THEN
           WRITE(IW,*) 'IEQU = ',(IEQU(JAT),JAT=1,NAT)
         ENDIF
C
C     ----- LOOK FOR 'GOOD' AXES BY DISTORTING THE MOLECULE -----
C           C2 AXES FOR -T - GROUP 
C           C2 AXES FOR -Th- GROUP 
C           S4 AXES FOR -Td- GROUP
C           C4 AXES FOR -O - GROUP
C           C4 AXES FOR -Oh- GROUP
C
         NAXM=0
         DO IAT=1,NAT
            IF(IEQU(IAT).EQ.IAT) THEN
               DO JAT=1,NAT
                  IF(JAT.NE.IAT.AND.IEQU(JAT).EQ.IAT) THEN
                     DO KAT=1,NAT
                        DO I=1,3
                           C2(I,KAT)=C1(I,KAT)
                        ENDDO
                     ENDDO
                     DO I=1,3
                        C2(I,IAT)=C1(I,IAT)*EPS
                        C2(I,JAT)=C1(I,JAT)*EPS
                     ENDDO
                     CALL HND_MOLINR(C2,C3,NAT,CM,AM,TM)
                     CALL HND_MOLAXS(AM,AXS,EIG,3,3,3)
                     IF(OUT) THEN
                        WRITE(IW,*) 'PRINCIPAL AXES OF ',
     1                              'DISTORTED MOLECULE ',
     2                              'IAT, JAT = ',IAT,JAT
                        CALL HND_PREV(AXS,EIG,3,3,3)
                     ENDIF
                     DO IIAT=1,NAT
                        DO J=1,3
                           DUM=ZERO
                           DO I=1,3
                              DUM=DUM+C1(I,IIAT)*AXS(I,J)
                           ENDDO
                           C2(J,IIAT)=DUM
                        ENDDO
                     ENDDO
                     IF(OUT) THEN
                        WRITE(IW,*) 'OLD COORDINATES'
                        DO IIAT=1,NAT
                           WRITE(IW,9989) IIAT,(C1(I,IIAT),I=1,3)
                        ENDDO
                        WRITE(IW,*) 'NEW COORDINATES'
                        DO IIAT=1,NAT
                           WRITE(IW,9989) IIAT,(C2(I,IIAT),I=1,3)
                        ENDDO
                     ENDIF
C
C     ----- CHECK THE AXES -----
C
                     DO IAXIS=1,3
                        IF(IAXIS.EQ.1) THEN
                           JAXIS=2
                           KAXIS=3
                        ELSEIF(IAXIS.EQ.2) THEN
                           JAXIS=3
                           KAXIS=1
                        ELSEIF(IAXIS.EQ.3) THEN
                           JAXIS=1
                           KAXIS=2
                        ENDIF
                        IF(DBUG) THEN
                           WRITE(IW,*) 'AXIS CHECK FOR -IAXIS- = ',IAXIS
                        ENDIF
                        C2AXS(IAXIS)=.TRUE.
                        C4AXS(IAXIS)=.TRUE.
                        S4AXS(IAXIS)=.TRUE.
                        DO IIAT=1,NAT
                           C2ROT=.FALSE.
                           JJC2 =0
                           DO JJAT=1,NAT
                              IF(.NOT.C2ROT) THEN
                                 D1=C2(IAXIS,IIAT)-C2(IAXIS,JJAT)                  
                                 D2=C2(JAXIS,IIAT)+C2(JAXIS,JJAT)                  
                                 D3=C2(KAXIS,IIAT)+C2(KAXIS,JJAT)                  
                                 C2ROT=(ABS(D1).LT.THREQUIV).AND.
     1                                 (ABS(D2).LT.THREQUIV).AND.
     2                                 (ABS(D3).LT.THREQUIV)
                                 IF(C2ROT) JJC2=JJAT
                              ENDIF
                           ENDDO
                           C4ROT=.FALSE.
                           JJC4 =0
                           DO JJAT=1,NAT
                              IF(.NOT.C4ROT) THEN
                                 D1=C2(IAXIS,IIAT)-C2(IAXIS,JJAT)
                                 D2=C2(JAXIS,IIAT)-C2(KAXIS,JJAT)
                                 D3=C2(KAXIS,IIAT)+C2(JAXIS,JJAT)
                                 C4ROT=(ABS(D1).LT.THREQUIV).AND.
     1                                 (ABS(D2).LT.THREQUIV).AND.
     2                                 (ABS(D3).LT.THREQUIV)
                                 IF(C4ROT) JJC4=JJAT
                              ENDIF
                           ENDDO
                           S4ROT=.FALSE.
                           JJS4 =0
                           DO JJAT=1,NAT
                              IF(.NOT.S4ROT) THEN
                                 D1=C2(IAXIS,IIAT)+C2(IAXIS,JJAT)
                                 D2=C2(JAXIS,IIAT)-C2(KAXIS,JJAT)
                                 D3=C2(KAXIS,IIAT)+C2(JAXIS,JJAT)
                                 S4ROT=(ABS(D1).LT.THREQUIV).AND.
     1                                 (ABS(D2).LT.THREQUIV).AND.
     2                                 (ABS(D3).LT.THREQUIV)
                                 IF(S4ROT) JJS4=JJAT
                              ENDIF
                           ENDDO
                           IF(DBUG) THEN
                              WRITE(IW,*) 'IIAT,JJC2,JJC4,JJS4 = ',
     1                                     IIAT,JJC2,JJC4,JJS4
                           ENDIF
                           C2AXS(IAXIS)=C2AXS(IAXIS).AND.C2ROT
                           C4AXS(IAXIS)=C4AXS(IAXIS).AND.C4ROT
                           S4AXS(IAXIS)=S4AXS(IAXIS).AND.S4ROT
                        ENDDO
                     ENDDO
                     IF(OUT) THEN
                        WRITE(IW,*) 'C2 AXES CHECK',(C2AXS(I),I=1,3)
                        WRITE(IW,*) 'C4 AXES CHECK',(C4AXS(I),I=1,3)
                        WRITE(IW,*) 'S4 AXES CHECK',(S4AXS(I),I=1,3)
                     ENDIF
C
C     ----- WE HAVE AXES THAT MAY BE C2 AND/OR C4 AND/OR S4 AXES -----
C           C2 AXES FOR -T - GROUP 
C           S4 AXES FOR -Td- GROUP
C           C4 AXES FOR -O - GROUP
C           C2 AXES FOR -Th- GROUP 
C           C4 AXES FOR -Oh- GROUP
C
                     IF(SYMINV) THEN
                        GRPOH=.FALSE.
                        DO IAXIS=1,3
                           GRPOH=GRPOH.OR.C4AXS(IAXIS)
                        ENDDO
                        GRPTH=.NOT.GRPOH
                        IF(OUT) THEN
                           WRITE(IW,*) 'TH , OH = ',GRPTH,GRPOH
                        ENDIF
                        DO IAXIS=1,3
                           IF((GRPTH.AND.
     1                         C2AXS(IAXIS)                 ).OR. 
     2                        (GRPOH.AND.
     3                         C2AXS(IAXIS).AND.C4AXS(IAXIS))) THEN
                              NAXM=NAXM+1
                              IF(NAXM.LE.3) THEN
                                 DO I=1,3
                                    AXM(I,NAXM)=AXS(I,IAXIS)
                                 ENDDO
                              ENDIF
                           ENDIF
                        ENDDO
                     ELSE
                        GRPTD=.FALSE.
                        GRPO =.FALSE.
                        DO IAXIS=1,3
                           GRPTD=GRPTD.OR.S4AXS(IAXIS)
                           GRPO =GRPO .OR.C4AXS(IAXIS)
                        ENDDO
                        GRPT =.NOT.GRPTD.AND..NOT.GRPO
                        IF(OUT) THEN
                           WRITE(IW,*) 'T , TD , O = ',
     1                                  GRPT,GRPTD,GRPO
                        ENDIF
                        DO IAXIS=1,3
                           IF((GRPO .AND.
     1                         (C2AXS(IAXIS).AND.C4AXS(IAXIS))).OR.
     2                        (GRPTD.AND.
     3                         (C2AXS(IAXIS).AND.S4AXS(IAXIS))).OR.
     4                        (GRPT .AND.
     5                         (C2AXS(IAXIS)                 ))) THEN
                              NAXM=NAXM+1
                              IF(NAXM.LE.3) THEN
                                 DO I=1,3
                                    AXM(I,NAXM)=AXS(I,IAXIS)
                                 ENDDO
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDIF
                     IF(OUT.AND.NAXM.GT.0) THEN
                        CALL HND_PRSQ(AXM,NAXM,3,3)
                     ENDIF
                     IF(NAXM.GE.2) GO TO 10
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
   10    CONTINUE
         DO IAXM=1,2
            DO I=1,3
               AXS(I,IAXM)=AXM(I,IAXM)
            ENDDO
         ENDDO
         AXS(1,3)=AXS(2,1)*AXS(3,2)-AXS(2,2)*AXS(3,1)
         AXS(2,3)=AXS(3,1)*AXS(1,2)-AXS(3,2)*AXS(1,1)
         AXS(3,3)=AXS(1,1)*AXS(2,2)-AXS(1,2)*AXS(2,1)
         IF(SOME) THEN
            WRITE(IW,*) ' NEW FRAME OF CUBIC GROUP '
            CALL HND_PRSQ(AXS,3,3,3)
         ENDIF
      ENDIF
C
      IF(DEGNR2) THEN
C
C     ----- IF DEGENERACY, THEN MAKE ONE CHOICE AND -----
C           EXPRESS COORDINATES IN THAT NEW FRAME.
C
         IF(SOME) THEN
            WRITE(IW,*) '2-FOLD DEGENERATE MOMENTS OF INERTIA'
         ENDIF
         IF(ABS(EIG(2)-EIG(1)).LT.TENM05) KAXIS=3
         IF(ABS(EIG(2)-EIG(3)).LT.TENM05) KAXIS=1
         IF(OUT) THEN
            WRITE(IW,*) 'KAXIS = ',KAXIS
         ENDIF
         DO IAT=1,NAT
            DUM=C1(1,IAT)*C1(1,IAT)+C1(2,IAT)*C1(2,IAT)+    
     1                              C1(3,IAT)*C1(3,IAT)      
            DUM=SQRT(DUM)
            IF(DUM.GT.TENM05) THEN
               IF(OUT) THEN
                  WRITE(IW,*) 
     1         '-DUM- NEW FRAME DEFINED FROM ATOM -IAT- ',IAT,DUM
               ENDIF
               DIRX=C1(1,IAT)/DUM
               DIRY=C1(2,IAT)/DUM
               DIRZ=C1(3,IAT)/DUM
               TMP=DIRX*AXS(1,KAXIS)+DIRY*AXS(2,KAXIS)+
     1                               DIRZ*AXS(3,KAXIS)
               DIRX=DIRX-TMP*AXS(1,KAXIS)
               DIRY=DIRY-TMP*AXS(2,KAXIS)
               DIRZ=DIRZ-TMP*AXS(3,KAXIS)
               TMP=DIRX*DIRX+DIRY*DIRY+DIRZ*DIRZ
               TMP=SQRT(TMP)
               IF(TMP.GT.TENM05) THEN
                  IF(OUT) THEN
                     WRITE(IW,*) 
     1         '-TMP- NEW FRAME DEFINED FROM ATOM -IAT- ',IAT,TMP
                  ENDIF
                  JAXIS=2
                  AXS(1,JAXIS)=DIRX/TMP
                  AXS(2,JAXIS)=DIRY/TMP
                  AXS(3,JAXIS)=DIRZ/TMP
                  DIRX=AXS(2,JAXIS)*AXS(3,KAXIS)-
     1                 AXS(2,KAXIS)*AXS(3,JAXIS)
                  DIRY=AXS(3,JAXIS)*AXS(1,KAXIS)-
     1                 AXS(3,KAXIS)*AXS(1,JAXIS)
                  DIRZ=AXS(1,JAXIS)*AXS(2,KAXIS)-
     1                 AXS(1,KAXIS)*AXS(2,JAXIS)
                  IF(KAXIS.EQ.3) THEN
                     IAXIS=1
                     AXS(1,IAXIS)=DIRX
                     AXS(2,IAXIS)=DIRY
                     AXS(3,IAXIS)=DIRZ
                  ELSEIF(KAXIS.EQ.1) THEN
                     IAXIS=3
                     AXS(1,IAXIS)=-DIRX
                     AXS(2,IAXIS)=-DIRY
                     AXS(3,IAXIS)=-DIRZ
                  ENDIF
                  GO TO 20
               ENDIF
            ENDIF
         ENDDO
   20    CONTINUE
         IF(OUT) THEN
            WRITE(IW,*) 'NEW FRAME = '
            CALL HND_PRSQ(AXS,3,3,3)
         ENDIF
      ENDIF
C
C     ----- NO DEGENERACY -----
C
      IF(.NOT.DEGNR2.AND..NOT.DEGNR3) THEN
         IF(SOME) THEN
            WRITE(IW,*) 'MOMENTS OF INERTIA ARE NOT DEGENERATE'
         ENDIF
      ENDIF 
C
   30 CONTINUE
      IF(OUT) THEN
         WRITE(IW,*) '    FRAME = '
         CALL HND_PRSQ(AXS,3,3,3)
      ENDIF
C
C     ----- TRANSFORM COORDINATES TO PRINCIPAL AXES -----
C
      DO IAT=1,NAT
         DO J=1,3
            DUM=ZERO
            DO I=1,3
               DUM=DUM+C1(I,IAT)*AXS(I,J)
            ENDDO
            C2(J,IAT)=DUM
         ENDDO
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9985)
         DO IAT=1,NAT
            WRITE(IW,9984) IAT,(C(I,IAT),I=1,3),(C1(I,IAT),I=1,3),
     1                                          (C2(I,IAT),I=1,3)
         ENDDO
      ENDIF
C
C     ----- DETERMINE ORDER OF PROPER/IMPROPER PRINCIPAL AXES -----  
C
      DO IORDR=1,MXORDR
         QX(IORDR)=QZERO 
         QY(IORDR)=QZERO 
         QZ(IORDR)=QZERO 
         DO IAT=1,NAT
            QDUMX=(CMPLX(C2(2,IAT),C2(3,IAT)))**IORDR
            QDUMY=(CMPLX(C2(1,IAT),C2(3,IAT)))**IORDR
            QDUMZ=(CMPLX(C2(1,IAT),C2(2,IAT)))**IORDR
            QX(IORDR)=QX(IORDR)+QDUMX
            QY(IORDR)=QY(IORDR)+QDUMY
            QZ(IORDR)=QZ(IORDR)+QDUMZ
            IF(MUCH) THEN
               WRITE(IW,9996) IORDR,IAT
               WRITE(IW,9995) QDUMX,QDUMY,QDUMZ
            ENDIF
         ENDDO
         XMAG(IORDR)=ABS(QX(IORDR))     
         YMAG(IORDR)=ABS(QY(IORDR))
         ZMAG(IORDR)=ABS(QZ(IORDR))
         IF(MUCH) THEN
            WRITE(IW,9998) QX(IORDR),QY(IORDR),QZ(IORDR)
            WRITE(IW,9999) IORDR,XMAG(IORDR),YMAG(IORDR),ZMAG(IORDR)
         ENDIF
      ENDDO
C
      AXORDR(1)=1
      AXORDR(2)=1
      AXORDR(3)=1
      GOTX=.FALSE.
      GOTY=.FALSE.
      GOTZ=.FALSE.
      DO IORDR=1,MXORDR
         IF((XMAG(IORDR).GT.TENM02).AND.(.NOT.GOTX)) THEN
            AXORDR(1)=IORDR
            GOTX=.TRUE.
         ENDIF
         IF((YMAG(IORDR).GT.TENM02).AND.(.NOT.GOTY)) THEN
            AXORDR(2)=IORDR
            GOTY=.TRUE.
         ENDIF
         IF((ZMAG(IORDR).GT.TENM02).AND.(.NOT.GOTZ)) THEN
            AXORDR(3)=IORDR
            GOTZ=.TRUE.
         ENDIF
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9997) (AXORDR(I),I=1,3)
      ENDIF
C
C     ----- PROPER OR IMPROPER ROTATION ?  -----
C
      DO IAXIS=1,3
         IF(IAXIS.EQ.1) THEN
            I=2
            J=3
            K=1
         ELSEIF(IAXIS.EQ.2) THEN
            I=3
            J=1
            K=2
         ELSEIF(IAXIS.EQ.3) THEN
            I=1
            J=2
            K=3
         ENDIF
         THETA=TWOPI/DBLE(AXORDR(IAXIS))
         IF(OUT) THEN
            WRITE(IW,9979) IAXIS,AXORDR(IAXIS)
         ENDIF
         PROPER(IAXIS)=.TRUE.
         DO IAT=1,NAT
            PRPAXS=.FALSE.
            DO JAT=1,NAT
               IF(NUC(JAT).EQ.NUC(IAT).AND..NOT.PRPAXS) THEN
                  QDUMI=CMPLX(C2(I,IAT),C2(J,IAT))
                  RR=ABS(QDUMI)
                  QDUMI=CMPLX(COS(THETA),SIN(THETA))*QDUMI
                  QDUMJ=CMPLX(C2(I,JAT),C2(J,JAT))
                  DD=ABS(QDUMI-QDUMJ)
                  DZ=C2(K,IAT)-C2(K,JAT)
                  PRPAXS=(ABS(DD).LT.THREQUIV).AND.
     1                   (ABS(DZ).LT.THREQUIV)
                  IF(DBUG) THEN
                     WRITE(IW,*) ' IAT,JAT,PRPAXS = ',
     1                             IAT,JAT,PRPAXS,RR,DD,DZ
                  ENDIF
               ENDIF
            ENDDO
            PROPER(IAXIS)=PROPER(IAXIS).AND.PRPAXS
            IF(DBUG) THEN
               WRITE(IW,*) ' IAT,PRPAXS,PROPER(IAXIS) = ',
     1                       IAT,PRPAXS,PROPER(IAXIS)
            ENDIF
         ENDDO
         IF(SOME) THEN
            WRITE(IW,9982) IAXIS,PROPER(IAXIS)
         ENDIF
         IMPROP(IAXIS)=.TRUE.
         DO IAT=1,NAT
            IMPAXS=.FALSE.
            DO JAT=1,NAT
               IF(NUC(JAT).EQ.NUC(IAT).AND..NOT.IMPAXS) THEN
                  QDUMI=CMPLX(C2(I,IAT),C2(J,IAT))
                  RR=ABS(QDUMI)
                  QDUMI=CMPLX(COS(THETA),SIN(THETA))*QDUMI
                  QDUMJ=CMPLX(C2(I,JAT),C2(J,JAT))
                  DD=ABS(QDUMI-QDUMJ)
                  DZ=C2(K,IAT)+C2(K,JAT)
                  IMPAXS=(ABS(DD).LT.THREQUIV).AND.
     1                   (ABS(DZ).LT.THREQUIV)
                  IF(DBUG) THEN
                     WRITE(IW,*) ' IAT,JAT,IMPAXS = ',
     1                             IAT,JAT,IMPAXS,RR,DD,DZ
                  ENDIF
               ENDIF
            ENDDO
            IMPROP(IAXIS)=IMPROP(IAXIS).AND.IMPAXS
            IF(DBUG) THEN
               WRITE(IW,*) ' IAT,IMPAXS,IMPROP(IAXIS) = ',
     1                       IAT,IMPAXS,IMPROP(IAXIS)
            ENDIF
         ENDDO
         IF(SOME) THEN
            WRITE(IW,9981) IAXIS,IMPROP(IAXIS)
         ENDIF
      ENDDO
C
C     ----- ALIGN AXIS OF HIGHEST ORDER ALONG -Z- -----
C
      IF(AXORDR(2).GT.AXORDR(1)) THEN
         KAXIS=2
         JAXIS=1
         IAXIS=3
      ELSE
         KAXIS=1
         JAXIS=3
         IAXIS=2
      ENDIF
      IF(AXORDR(3).GT.AXORDR(KAXIS)) THEN
         KAXIS=3
         JAXIS=2
         IAXIS=1
      ENDIF
      PRPAXS=PROPER(1).OR.PROPER(2).OR.PROPER(3)
      IF(PRPAXS) THEN
         IF(.NOT.PROPER(KAXIS)) THEN
            IF(AXORDR(JAXIS).GE.AXORDR(IAXIS)) THEN
               KAXIS=JAXIS
               IF(.NOT.PROPER(KAXIS)) THEN
                  KAXIS=IAXIS
               ENDIF
            ELSE
               KAXIS=IAXIS
               IF(.NOT.PROPER(KAXIS)) THEN
                  KAXIS=JAXIS
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF(SOME) THEN
         WRITE(IW,9994) AXORDR(KAXIS)
      ENDIF
C
      IF(KAXIS.EQ.3) THEN
         NEWAXS(1)=1
         NEWAXS(2)=2
         NEWAXS(3)=3
      ELSEIF(KAXIS.EQ.2) THEN
         NEWAXS(1)=3
         NEWAXS(2)=1
         NEWAXS(3)=2
      ELSEIF(KAXIS.EQ.1) THEN
         NEWAXS(1)=2
         NEWAXS(2)=3
         NEWAXS(3)=1
      ENDIF
      DO K=1,3
         NUORDR(K)=AXORDR(NEWAXS(K))
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9987) KAXIS
         WRITE(IW,9988) (NEWAXS(K),K=1,3)
         WRITE(IW,9980) (NUORDR(K),K=1,3)
      ENDIF
C
C     ----- C'2 AXES PERPENDICULAR TO PRINCIPAL AXIS -----
C
      SYMC2X=MOD(AXORDR(NEWAXS(1)),2).EQ.0.AND.
     1           PROPER(NEWAXS(1))
      SYMC2Y=MOD(AXORDR(NEWAXS(2)),2).EQ.0.AND.
     1           PROPER(NEWAXS(1))
      IF(SOME) THEN
         WRITE(IW,9992) SYMC2X,SYMC2Y
      ENDIF
C
C     ----- ASSEMBLE PERMUTATION MATRIX -----
C
      DO J=1,3
         DO I=1,3
            PRM(I,J)=ZERO
         ENDDO
      ENDDO
      PRM(NEWAXS(1),1)=ONE
      PRM(NEWAXS(2),2)=ONE
      PRM(NEWAXS(3),3)=ONE
      IF(DBUG) THEN
         WRITE(IW,*) 'PERMUTATION MATRIX'
         CALL HND_PRSQ(PRM,3,3,3)
      ENDIF
C
C     ----- FIRST WAY OF GETTING -C3- FROM -C2- -----
C
      DO IAT=1,NAT
         DO J=1,3
            DUM=ZERO
            DO I=1,3   
               DUM=DUM+C2(I,IAT)*PRM(I,J) 
            ENDDO
            C3(J,IAT)=DUM
         ENDDO
      ENDDO
      IF(DBUG) THEN
         WRITE(IW,9990)
         DO IAT=1,NAT
            WRITE(IW,9989) IAT,(C3(I,IAT),I=1,3)
         ENDDO
      ENDIF
C
C     ----- SECOND WAY OF GETTING -C3- FROM -C2- -----
C
      DO IAT=1,NAT
         XX=C2(NEWAXS(1),IAT)
         YY=C2(NEWAXS(2),IAT)
         ZZ=C2(NEWAXS(3),IAT)
         C3(1,IAT)=XX
         C3(2,IAT)=YY
         C3(3,IAT)=ZZ
      ENDDO
      IF(DBUG) THEN
         WRITE(IW,9990)
         DO IAT=1,NAT
            WRITE(IW,9989) IAT,(C3(I,IAT),I=1,3)
         ENDDO
      ENDIF
C
C     ----- ASSEMBLE FINAL ROTATION MATRIX -----
C
      DO J=1,3
         DO I=1,3
            DUM=ZERO
            DO K=1,3
               DUM=DUM+AXS(I,K)*PRM(K,J)
            ENDDO
            RT(I,J)=DUM
         ENDDO
      ENDDO
C
C     ----- NOW GET -C3- FROM -C1- -----
C
      DO IAT=1,NAT
         DO J=1,3
            DUM=ZERO
            DO I=1,3
               DUM=DUM+C1(I,IAT)*RT(I,J)
            ENDDO
            C3(J,IAT)=DUM
         ENDDO
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9990)
         DO IAT=1,NAT
            WRITE(IW,9989) IAT,(C3(I,IAT),I=1,3)
         ENDDO
      ENDIF
C
C     ----- CHECK : AXES OF ORDER HIGHER THAN 1 MUST  -----
C           BE EITHER PROPER OR IMPROPER. IF ALL AXES
C           ARE NEITHER PROPER NOR IMPROPER, THEN THE
C           ALGORITHM FAILED, AND IT SHOULD BE A 
C           - NO SYMMETRY - CASE .
C
      NOSYM=.FALSE.
      DO IAXIS=1,3
         NOSYM=NOSYM.AND.
     1        (.NOT.PROPER(IAXIS).AND..NOT.IMPROP(IAXIS))
      ENDDO
      IF(NOSYM) THEN
         AXORDR(KAXIS)=1
         DO J=1,3
            DO I=1,3
               RT(I,J)=ZERO
            ENDDO
            RT(J,J)=ONE
            TR(J  )=ZERO
         ENDDO
         DO IAT=1,NAT
            DO I=1,3
               C3(I,IAT)=C(I,IAT)
               C2(I,IAT)=C(I,IAT)
               C1(I,IAT)=C(I,IAT)
            ENDDO
         ENDDO
      ENDIF
C
C     ----- DETERMINE PRESENCE OF CENTER OF INVERSION -----
C
      SYMINV=.TRUE. 
      DO IAT=1,NAT
         INVERS=.FALSE.
         DO JAT=1,NAT
            IF(NUC(JAT).EQ.NUC(IAT).AND..NOT.INVERS) THEN
               DX=C3(1,IAT)+C3(1,JAT)
               DY=C3(2,IAT)+C3(2,JAT)
               DZ=C3(3,IAT)+C3(3,JAT)
               INVERS=(ABS(DX).LT.THREQUIV).AND.
     1                (ABS(DY).LT.THREQUIV).AND.
     2                (ABS(DZ).LT.THREQUIV)
               IF(DBUG) THEN
                  WRITE(IW,*) ' IAT,JAT,INVERS = ',
     1                          IAT,JAT,INVERS,DX,DY,DZ
               ENDIF
            ENDIF
         ENDDO
         SYMINV=SYMINV.AND.INVERS
         IF(DBUG) THEN
            WRITE(IW,*) ' IAT,INVERS,SYMINV = ',IAT,INVERS,SYMINV
         ENDIF
      ENDDO
      IF(SOME) THEN
         WRITE(IW,9993) SYMINV
      ENDIF
C
C     ----- MIRROR PLANES PERPENDICULAR TO PRINCIPAL AXES -----
C
      MIRRYZ=.TRUE. 
      DO IAT=1,NAT
         MIRROR=.FALSE.
         DO JAT=1,NAT
            IF(NUC(JAT).EQ.NUC(IAT).AND..NOT.MIRROR) THEN
               DX=C3(1,IAT)+C3(1,JAT)
               DY=C3(2,IAT)-C3(2,JAT)
               DZ=C3(3,IAT)-C3(3,JAT)
               MIRROR=(ABS(DX).LT.THREQUIV).AND.
     1                (ABS(DY).LT.THREQUIV).AND.
     2                (ABS(DZ).LT.THREQUIV)
               IF(DBUG) THEN
                  WRITE(IW,*) ' IAT,JAT,MIRRYZ = ',
     1                          IAT,JAT,MIRROR,DX,DY,DZ
               ENDIF
            ENDIF
         ENDDO
         MIRRYZ=MIRRYZ.AND.MIRROR
         IF(OUT) THEN
            WRITE(IW,*) ' IAT,MIRROR,MIRRYZ = ',IAT,MIRROR,MIRRYZ
         ENDIF
      ENDDO
      MIRRZX=.TRUE.
      DO IAT=1,NAT
         MIRROR=.FALSE.
         DO JAT=1,NAT
            IF(NUC(JAT).EQ.NUC(IAT).AND..NOT.MIRROR) THEN
               DX=C3(1,IAT)-C3(1,JAT)
               DY=C3(2,IAT)+C3(2,JAT)
               DZ=C3(3,IAT)-C3(3,JAT)
               MIRROR=(ABS(DX).LT.THREQUIV).AND.
     1                (ABS(DY).LT.THREQUIV).AND.
     2                (ABS(DZ).LT.THREQUIV)
               IF(DBUG) THEN
                  WRITE(IW,*) ' IAT,JAT,MIRRZX = ',
     1                          IAT,JAT,MIRROR,DX,DY,DZ
               ENDIF
            ENDIF
         ENDDO
         MIRRZX=MIRRZX.AND.MIRROR
         IF(OUT) THEN
            WRITE(IW,*) ' IAT,MIRROR,MIRRZX = ',IAT,MIRROR,MIRRZX
         ENDIF
      ENDDO
      MIRRXY=.TRUE.
      DO IAT=1,NAT
         MIRROR=.FALSE.
         DO JAT=1,NAT
            IF(NUC(JAT).EQ.NUC(IAT).AND..NOT.MIRROR) THEN
               DX=C3(1,IAT)-C3(1,JAT)
               DY=C3(2,IAT)-C3(2,JAT)
               DZ=C3(3,IAT)+C3(3,JAT)
               MIRROR=(ABS(DX).LT.THREQUIV).AND.
     1                (ABS(DY).LT.THREQUIV).AND.
     2                (ABS(DZ).LT.THREQUIV)
               IF(DBUG) THEN
                  WRITE(IW,*) ' IAT,JAT,MIRRXY = ',
     1                          IAT,JAT,MIRROR,DX,DY,DZ
               ENDIF
            ENDIF
         ENDDO
         MIRRXY=MIRRXY.AND.MIRROR
         IF(OUT) THEN
            WRITE(IW,*) ' IAT,MIRROR,MIRRXY = ',IAT,MIRROR,MIRRXY
         ENDIF
      ENDDO
      IF(SOME) THEN
         WRITE(IW,9991) MIRRYZ,MIRRZX,MIRRXY
      ENDIF         
C
C     ----- SET MOLECULAR POINT GROUP SYMBOL -----
C
      CALL HND_MOLSMB(MIRRYZ,MIRRZX,MIRRXY,
     1            SYMC2X,SYMC2Y,SYMINV,
     2            PROPER(NEWAXS(1)),PROPER(NEWAXS(2)),
     3            PROPER(NEWAXS(3)),
     4            CUBIC,GRPOH,GRPTH,GRPTD,GRPT,GRPO,
     5            AXORDR(KAXIS),TR,RT,ODONE,groupname)
C
      IF(OUT) THEN
         WRITE(IW,9983)
         DO IAT=1,NAT
            WRITE(IW,9984) IAT,(C(I,IAT),I=1,3),(C1(I,IAT),I=1,3),
     1                                          (C3(I,IAT),I=1,3)
         ENDDO
      ENDIF
      IF(ONLY) THEN
         WRITE(IW,9977)
         DO IAT=1,NAT
            WRITE(IW,9989) IAT,(C3(I,IAT),I=1,3)
         ENDDO
      ENDIF
C
      RETURN
 9999 FORMAT(' IORDR= ',I3,' XMAG= ',G12.4,' YMAG= ',G12.4,
     1                                     ' ZMAG= ',G12.4)
 9998 FORMAT(' ABOUT X-AXIS = ',G12.4,' + i ',G12.4,/,
     1       ' ABOUT Y-AXIS = ',G12.4,' + i ',G12.4,/,
     2       ' ABOUT Z-AXIS = ',G12.4,' + i ',G12.4)
 9997 FORMAT(' X-AXIS ORDER = ',I3,' Y-AXIS ORDER = ',I3,
     1                             ' Z-AXIS ORDER = ',I3)
 9996 FORMAT(' IORDR,IAT = ',2I4)
 9995 FORMAT(2G12.4)
 9994 FORMAT(' HIGHEST ORDER OF PROPER/IMPROPER ROTATION = ',I2)
 9993 FORMAT(' INVERSION CENTER = ',L4)
 9992 FORMAT('           C2 -X- = ',L4,/,
     1       '           C2 -Y- = ',L4) 
 9991 FORMAT(' MIRROR PLANE -YZ-= ',L4,/,
     1       ' MIRROR PLANE -ZX-= ',L4,/,
     2       ' MIRROR PLANE -XY-= ',L4)
 9990 FORMAT(' COORDINATES AFTER ALIGNEMENT OF AXES ')
 9989 FORMAT(1X,I5,3F12.8)
 9988 FORMAT('           PERMUTATION IS = ',I1,1X,I1,1X,I1)
 9987 FORMAT(' AXIS OF HIGHEST ORDER IS = ',I1)
 9986 FORMAT(' INERTIA AXES AND DATA ')
 9985 FORMAT(/,17X,'ORIGINAL -XYZ-',22X,'MASS-CENTERED -XYZ-',
     1         17X,'ROTATED -XYZ-')
 9984 FORMAT(1X,I5,9F12.8)
 9983 FORMAT(/,17X,'ORIGINAL -XYZ-',22X,'MASS-CENTERED -XYZ-',
     1         17X,' FINAL  -XYZ-')
 9982 FORMAT('   PROPER PRINCIPAL AXIS ',I1,' = ',L4)
 9981 FORMAT(' IMPROPER PRINCIPAL AXIS ',I1,' = ',L4)
 9980 FORMAT(' ORDER OF PERMUTED AXES = ',3I3)
 9979 FORMAT(' PROPER/IMPROPER CHECK FOR -IAXIS- AND -ORDER- = ',2I3)
 9978 FORMAT(' I EIG = ',I5,F20.8)
 9977 FORMAT(' coordinates in local symmetry frame ',/)
      END
      SUBROUTINE HND_MOLSMB(MIRRYZ,MIRRZX,MIRRXY,
     1                     SYMC2X,SYMC2Y,SYMINV,
     2                     PROPRX,PROPRY,PROPRZ,
     3                     CUBIC,GRPOH,GRPTH,GRPTD,GRPT,GRPO,
     4                     AXORDR,TR,RT,ODONE,groupname)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL  ODONE
      LOGICAL  DBUG
      LOGICAL  SOME
      LOGICAL  PUNCH
      CHARACTER*8 WRDSYM
      CHARACTER*8 WRDEND
      CHARACTER*8 GRPSAV
      CHARACTER*8 GROUP,groupname
      CHARACTER*8 BLANK
      CHARACTER*1 B
      CHARACTER*1 C
      CHARACTER*1 S
      CHARACTER*1 D
      CHARACTER*1 N
      CHARACTER*1 V
      CHARACTER*1 H
      CHARACTER*2 C1
      CHARACTER*2 CI
      CHARACTER*2 CS
      CHARACTER*2 T 
      CHARACTER*2 TH
      CHARACTER*2 TD
      CHARACTER*2 O 
      CHARACTER*2 OH
      INTEGER  AXORDR
      LOGICAL  CUBIC
      LOGICAL  GRPOH
      LOGICAL  GRPTH
      LOGICAL  GRPTD
      LOGICAL  GRPT
      LOGICAL  GRPO
      LOGICAL  PROPRX
      LOGICAL  PROPRY
      LOGICAL  PROPRZ
      LOGICAL  SYMINV
      LOGICAL  SYMC2X
      LOGICAL  SYMC2Y
      LOGICAL  MIRRYZ
      LOGICAL  MIRRZX
      LOGICAL  MIRRXY
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_SYMMOL/COMPLX,IGROUP,NAXIS,LINEAR
      COMMON/HND_SYMNAM/GROUP
      DIMENSION TR(3),RT(3,3)
      DATA BLANK  /'        '/
      DATA B      /' '/
      DATA S      /'S'/
      DATA C      /'C'/
      DATA D      /'D'/
      DATA N      /'N'/
      DATA V      /'V'/
      DATA H      /'H'/
      DATA C1     /'C1'/
      DATA CI     /'CI'/
      DATA CS     /'CS'/
      DATA T      /'T '/
      DATA TH     /'TH'/
      DATA TD     /'TD'/
      DATA O      /'O '/
      DATA OH     /'OH'/
      DATA WRDSYM /' $SYM   '/
      DATA WRDEND /' $END   '/
C
      DBUG=.FALSE.
      SOME=.TRUE. 
      SOME=SOME.OR.DBUG
C
      PUNCH=.FALSE.
C
      IF(DBUG) THEN
         WRITE(IW,*) 'IN -MOLSMB- '
         WRITE(IW,*) 'MIRRYZ = ',MIRRYZ
         WRITE(IW,*) 'MIRRZX = ',MIRRZX
         WRITE(IW,*) 'MIRRXY = ',MIRRXY
         WRITE(IW,*) 'SYMC2X = ',SYMC2X
         WRITE(IW,*) 'SYMC2Y = ',SYMC2Y
         WRITE(IW,*) 'SYMINV = ',SYMINV
         WRITE(IW,*) 'PROPRX = ',PROPRX
         WRITE(IW,*) 'PROPRY = ',PROPRY
         WRITE(IW,*) 'PROPRZ = ',PROPRZ
         WRITE(IW,*) 'CUBIC  = ',CUBIC 
         WRITE(IW,*) 'GRPOH  = ',GRPOH 
         WRITE(IW,*) 'GRPTH  = ',GRPTH 
         WRITE(IW,*) 'GRPTD  = ',GRPTD 
         WRITE(IW,*) 'GRPT   = ',GRPO  
         WRITE(IW,*) 'GRPO   = ',GRPO  
         WRITE(IW,*) 'AXORDR = ',AXORDR
      ENDIF
C
      NAXSAV=NAXIS
      GRPSAV=GROUP
      IGRSAV=IGROUP
      LINSAV=LINEAR
      CMXSAV=COMPLX
C
C     ----- SET GROUP ORDER -----
C
      NAXIS=AXORDR
C
C     ----- SET GROUP SYMBOL -----
C
      GROUP=BLANK
      IF(NAXIS.EQ.1) THEN
         GROUP(1:2)=C1(1:2)
         IF(SYMINV) THEN
            GROUP(1:2)=CI
         ENDIF
         IF(MIRRYZ.OR.MIRRZX.OR.MIRRXY) THEN
            GROUP(1:2)=CS
         ENDIF
      ELSE
         IF(.NOT.CUBIC) THEN
            IF(SYMC2X.OR.SYMC2Y) THEN
               GROUP(1:1)=D
               GROUP(2:2)=N
               IF(MIRRXY) THEN
                  GROUP(3:3)=H
               ELSE
                  IF(MIRRYZ.OR.MIRRZX) THEN
                     GROUP(3:3)=B
                  ELSE
                     GROUP(3:3)=D
                     NAXIS=NAXIS/2
                  ENDIF
               ENDIF
            ELSE
               IF(PROPRZ) THEN
                  GROUP(1:1)=C
                  GROUP(2:2)=N
                  IF(MIRRXY) THEN
                     GROUP(3:3)=H
                  ELSE
                     IF(MIRRYZ.OR.MIRRZX) THEN
                        GROUP(3:3)=V
                     ELSE
                        GROUP(3:3)=B
                     ENDIF
                  ENDIF
               ELSE
                  GROUP(1:1)=S
                  GROUP(2:2)=N
               ENDIF
            ENDIF
         ELSE
C
C     ----- CUBIC GROUPS -----
C
            IF(GRPOH) THEN
               GROUP(1:2)=OH
            ELSEIF(GRPTH) THEN
               GROUP(1:2)=TH
            ELSEIF(GRPTD) THEN
               GROUP(1:2)=TD
            ELSEIF(GRPT ) THEN
               GROUP(1:2)=T
            ELSEIF(GRPO ) THEN
               GROUP(1:2)=O
            ENDIF
         ENDIF
      ENDIF
C
      if (group(2:2).eq.'N') then
         groupname = ' '
         write(groupname,'(a1,i1,a6)') group(1:1),naxis,group(3:8)
      else
         groupname = group
      endif
c
      IF(SOME) THEN
         WRITE(IW,9999) groupname,
     1                   (TR(I  ),I=1,3),
     2                  ((RT(I,J),J=1,3),I=1,3)
      ENDIF
c
      IF(PUNCH) THEN
         WRITE(IW,9998) WRDSYM
         WRITE(IW,9997) GROUP(1:3),NAXIS
         WRITE(IW,9998) WRDEND
         WRITE(IW,*) 'TRANSLATION'
         WRITE(IW,9996) (TR(I),I=1,3)
         WRITE(IW,*) 'ROTATION   '
         WRITE(IW,9996) ((RT(I,J),J=1,3),I=1,3)
      ENDIF
C
C     ----- GET OPERATOR MATRICES -----
C
      IF(GROUP(1:2).EQ.CS) THEN
         XPT1=TR(1)
         YPT1=TR(2)
         ZPT1=TR(3)
         IF(MIRRYZ) THEN
            XPT2=XPT1+RT(1,2)
            YPT2=YPT1+RT(2,2)
            ZPT2=ZPT1+RT(3,2)
            XPT3=XPT1+RT(1,3)
            YPT3=YPT1+RT(2,3)
            ZPT3=ZPT1+RT(3,3)
         ELSEIF(MIRRZX) THEN
            XPT2=XPT1+RT(1,3)
            YPT2=YPT1+RT(2,3)
            ZPT2=ZPT1+RT(3,3)
            XPT3=XPT1+RT(1,1)
            YPT3=YPT1+RT(2,1)
            ZPT3=ZPT1+RT(3,1)
         ELSEIF(MIRRXY) THEN
            XPT2=XPT1+RT(1,1)
            YPT2=YPT1+RT(2,1)
            ZPT2=ZPT1+RT(3,1)
            XPT3=XPT1+RT(1,2)
            YPT3=YPT1+RT(2,2)
            ZPT3=ZPT1+RT(3,2)
         ENDIF
      ELSE
         IF(MIRRZX) THEN
            XPT0=TR(1)
            YPT0=TR(2)
            ZPT0=TR(3)
            XPT1=XPT0+RT(1,3)
            YPT1=YPT0+RT(2,3)
            ZPT1=ZPT0+RT(3,3)
            XPT2=XPT0+RT(1,1)
            YPT2=YPT0+RT(2,1)
            ZPT2=ZPT0+RT(3,1)
         ELSEIF(MIRRYZ) THEN
            XPT0=TR(1)
            YPT0=TR(2)
            ZPT0=TR(3)
            XPT1=XPT0+RT(1,3)
            YPT1=YPT0+RT(2,3)
            ZPT1=ZPT0+RT(3,3)
            XPT2=XPT0+RT(1,2)
            YPT2=YPT0+RT(2,2)
            ZPT2=ZPT0+RT(3,2)
         ELSE
            XPT0=TR(1)
            YPT0=TR(2)
            ZPT0=TR(3)
            XPT1=XPT0+RT(1,3)
            YPT1=YPT0+RT(2,3)
            ZPT1=ZPT0+RT(3,3)
            XPT2=XPT0+RT(1,1)
            YPT2=YPT0+RT(2,1)
            ZPT2=ZPT0+RT(3,1)
         ENDIF
      ENDIF
      CALL HND_MOLGRP(XPT0,YPT0,ZPT0,XPT1,YPT1,ZPT1,
     1                XPT2,YPT2,ZPT2,XPT3,YPT3,ZPT3)
C
      ODONE=.TRUE.
C
      RETURN
 9999 FORMAT(
     $       ' ',a3,' symmetry detected'/
     1       ' Molecular frame translation',3F10.5,/,
     2       ' Molecular frame rotation   ',3F10.5,/,
     3       '                            ',3F10.5,/,
     4       '                            ',3F10.5)
 9998 FORMAT(A8)
 9997 FORMAT(1X,A3,2X,I5)
 9996 FORMAT(3F12.8)
      END
      SUBROUTINE HND_MOLAXS(A,VEC,EIG,NVEC,N,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- ROUTINE TO SUBSTITUTE DIAGIV FOR DIAGONALIZATION -----
C           OF SYMMETRIC 3X3 MATRIX A  IN TRIANGULAR FORM
C
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION A(3,3),VEC(3,3),EIG(3)
      DIMENSION AA(6),IA(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-MOLAXS-'/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA CONV   /1.0D-10/
      DATA MAXIT  /50/
C
      IF(N.EQ.3.AND.NDIM.EQ.3) GO TO 10
         WRITE(IW,9999)
         CALL HND_HNDERR(3,ERRMSG)
   10 CONTINUE
C
      DO 30 I=1,3
         DO 20 J=1,3
            VEC(J,I)=ZERO
   20    CONTINUE
         VEC(I,I)=ONE
   30 CONTINUE
      AA(1)=A(1,1)
      AA(2)=A(2,1)
      AA(3)=A(2,2)
      AA(4)=A(3,1)
      AA(5)=A(3,2)
      AA(6)=A(3,3)
      IA(1)=0
      IA(2)=1
      IA(3)=3

      CALL HND_DIAGIV(AA,VEC,EIG,IA,N,N,N)
C
C     ----- CHECK FOR RIGHT HANDEDNESS, CORRECT IF NOT -----
C
      TEST =   VEC(1,3)*( VEC(2,1)*VEC(3,2) - VEC(3,1)*VEC(2,2) )
     1       + VEC(2,3)*( VEC(3,1)*VEC(1,2) - VEC(1,1)*VEC(3,2) )
     2       + VEC(3,3)*( VEC(1,1)*VEC(2,2) - VEC(2,1)*VEC(1,2) )
      IF(TEST.GT.ZERO) RETURN
      IF( ABS(EIG(1)-EIG(2)).GT.CONV) GO TO 60
         T = EIG(1)
         EIG(1) = EIG(2)
         EIG(2) = T
         DO 50 I=1,3
            T = VEC(I,1)
            VEC(I,1) = VEC(I,2)
            VEC(I,2) = T
   50    CONTINUE
         RETURN
   60 IF( ABS(EIG(2)-EIG(3)).GT.CONV) GO TO 80
         T = EIG(2)
         EIG(2) = EIG(3)
         EIG(3) = T
         DO 70 I=1,3
            T = VEC(I,2)
            VEC(I,2) = VEC(I,3)
            VEC(I,3) = T
   70    CONTINUE
         RETURN
   80 DO 90 I=1,3
         VEC(I,3) = - VEC(I,3)
   90 CONTINUE
      RETURN
 9999 FORMAT(/,' -DIAAXS- DIAGONALIZATION ONLY SET UP FOR 3X3 MATRIX')
      END
      SUBROUTINE HND_MOLINR(C,C1,NAT,CI,AI,TR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- CENTER AND MOMENTS OF INERTIA ----- 
C
      PARAMETER (MXATOM=500)
      PARAMETER (ZERO=0.0D+00)
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      DIMENSION AM(MXATOM)
      DIMENSION C(3,*),C1(3,*)
      DIMENSION CI(3),AI(3,3)
      DIMENSION TR(3)
C
      OUT=.FALSE.
C
C     -----                GET ATOMIC MASSES                -----
C     ----- CALCULATE CENTRE OF MASS AND MOMENTS OF INERTIA -----
C
      CALL HND_MOLAMS(AM,OUT)
C
      DO J=1,3
         CI(J)=ZERO
         DO I=1,3
            AI(J,I)=ZERO
         ENDDO
      ENDDO
      AMASS=ZERO
      DO IAT=1,NAT
         AMASS=AMASS+AM(IAT)
         DO I=1,3
            CI(I)=CI(I)+AM(IAT)*C(I,IAT)
         ENDDO
      ENDDO
      DO I=1,3
         CI(I)=CI(I)/AMASS
         TR(I)=CI(I)
      ENDDO    
      DO IAT=1,NAT
         DO I=1,3
            C1(I,IAT)=C(I,IAT)-TR(I) 
         ENDDO    
      ENDDO
C
      DO IAT=1,NAT
         WT=AM(IAT)
         X =C1(1,IAT)
         Y =C1(2,IAT)
         Z =C1(3,IAT)
         AI(1,1)=AI(1,1)+WT*(Y*Y+Z*Z)
         AI(2,1)=AI(2,1)-WT* X*Y
         AI(1,2)=AI(2,1)
         AI(3,1)=AI(3,1)-WT* X*Z
         AI(1,3)=AI(3,1)
         AI(2,2)=AI(2,2)+WT*(X*X+Z*Z)
         AI(3,2)=AI(3,2)-WT* Y*Z
         AI(2,3)=AI(3,2)
         AI(3,3)=AI(3,3)+WT*(X*X+Y*Y)
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9999)
         WRITE(IW,9998) (CI(I),I=1,3)
         WRITE(IW,9997)
         CALL HND_PRSQ(AI,3,3,3)
      ENDIF
C
      RETURN
 9999 FORMAT(/,1X,'CENTER OF MASS',/,1X,14(1H-))
 9998 FORMAT(' X = ',F12.8,' Y = ',F12.8,' Z = ',F12.8)
 9997 FORMAT(/,1X,'MOMENTS OF INERTIA',/,1X,18(1H-))
      END
      SUBROUTINE HND_MOLAMS(BMASS,SOME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER   (MXATOM=500)
      LOGICAL     SOME
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      DIMENSION   AMASS(MXATOM)
      DIMENSION   BMASS(*)
      DIMENSION   AMS(105)
      DIMENSION   ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-MOLAMS-'/
      DATA (AMS(I),I=1,54)  /
     1   1.007825D+00,  4.002600D+00,  7.016000D+00,  9.012180D+00,
     2  11.009310D+00, 12.000000D+00, 14.003070D+00, 15.994910D+00,
     3  18.998400D+00, 19.992440D+00, 22.989800D+00, 23.985040D+00,
     4  26.981530D+00, 27.976930D+00, 30.973760D+00, 31.972070D+00,
     5  34.968850D+00, 39.948000D+00, 38.963710D+00, 39.962590D+00,
     6  44.955920D+00, 47.900000D+00, 50.944000D+00, 51.940500D+00,
     7  54.938100D+00, 55.934900D+00, 58.933200D+00, 57.935300D+00,
     8  62.929800D+00, 63.929100D+00, 68.925700D+00, 73.921900D+00,
     9  74.921600D+00, 79.916500D+00, 78.918300D+00, 83.911500D+00,
     1  84.911700D+00, 87.905600D+00, 89.905400D+00, 89.904300D+00,
     2  92.906000D+00, 97.905500D+00, 97.000000D+00,101.903700D+00,
     3 102.904800D+00,105.903200D+00,106.904100D+00,113.903600D+00,
     4 114.904100D+00,119.902200D+00,120.903800D+00,129.906700D+00,
     5 126.904400D+00,131.904200D+00/
      DATA (AMS(I),I=55,105)  /
     1   132.9054D+00,137.9052D+00,138.9063D+00,139.9054D+00,
     2   140.9076D+00,141.9077D+00,144.9127D+00,151.9197D+00,
     3   152.9212D+00,157.9241D+00,158.9253D+00,163.9292D+00,
     4   164.9303D+00,165.9303D+00,168.9342D+00,173.9389D+00,
     5   174.9408D+00,179.9465D+00,180.9480D+00,183.9509D+00,
     6   186.9557D+00,191.9615D+00,192.9629D+00,194.9648D+00,
     7   196.9665D+00,201.9706D+00,204.9744D+00,207.9766D+00,
     8   208.9804D+00,208.9824D+00,209.9871D+00,222.0176D+00,
     9   223.0197D+00,226.0254D+00,227.0278D+00,232.0381D+00,
     1   231.0359D+00,238.0508D+00,237.0482D+00,244.0642D+00,
     2   243.0614D+00,247.0703D+00,247.0703D+00,251.0796D+00,
     3   252.0829D+00,257.0751D+00,258.0986D+00,259.1009D+00,
     4   260.1053D+00,  0.0000D+00,  0.0000D+00/
      DATA ZERO /0.0D+00/
C
      DO IAT=1,NAT
         NUCZ=NUC(IAT)
         IF(NUCZ.GT.103) THEN
            WRITE(IW,9996)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(NUCZ.EQ.0) ZNUC=ZERO
         IF(NUCZ.GT.0) ZNUC=AMS(NUCZ)
         AMASS(IAT)=ZNUC
      ENDDO
      DO IAT=1,NAT
         BMASS(IAT)= ABS(AMASS(IAT))
      ENDDO
      IF(SOME) THEN
         WRITE(IW,9998)
         DO IAT=1,NAT
            WRITE(IW,9997) IAT,BMASS(IAT)
         ENDDO
      ENDIF
      RETURN
 9998 FORMAT(/,10X,21(1H-),/,10X,'ATOMIC WEIGHTS (A.U.)',
     1       /,10X,21(1H-),/)
 9997 FORMAT(I5,5X,F15.5)
 9996 FORMAT(' NO ATOMIC MASS TABULATED FOR  -NUCZ.GT.103-  . STOP')
      END
      SUBROUTINE HND_MOLGRP(XPT0,YPT0,ZPT0,XPT1,YPT1,ZPT1,         
     1                      XPT2,YPT2,ZPT2,XPT3,YPT3,ZPT3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER     (MXSYM =120)
      PARAMETER     (MXSYMT=120*9)
      PARAMETER     (MXIODA=255)
      CHARACTER*8   ERRMSG
      LOGICAL       OUT
      LOGICAL       ALIGN
      CHARACTER*8   GRPCHR
      CHARACTER*8   DRCCHR
      CHARACTER*8   CHRGRP
      CHARACTER*8   CHRDIR
      CHARACTER*8   CHRBLK
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_SYMTRY/INVT(MXSYM),NT,NTMAX,NTWD,NOSYM
      COMMON/HND_SYMTRF/XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
      COMMON/HND_FRAME/U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
      COMMON/HND_SYMMOL/COMPLX,IGROUP,NAXIS,LINEAR
      COMMON/HND_SYMNAM/GROUP
      COMMON/HND_SYMMAT/T(MXSYMT)
      COMMON/HND_ORIENT/ALIGN
      DIMENSION     GRPCHR(19)
      DIMENSION     DRCCHR(2)
      DIMENSION     GRP(19)
      DIMENSION     DRC(2)
      DIMENSION     ERRMSG(3)
      EQUIVALENCE   (CHRGRP,GROUP )
      EQUIVALENCE   (CHRDIR,DIRECT)
      EQUIVALENCE   (CHRBLK,BLANK)
      EQUIVALENCE   (GRP(1),GRPCHR(1))
      EQUIVALENCE   (DRC(1),DRCCHR(1)) 
      DATA ERRMSG   /'PROGRAM ','STOP IN ','- PTGRP-'/
      DATA GRPCHR   /'C1      ','CS      ','CI      ',
     1               'CN      ','S2N     ','CNH     ',
     2               'CNV     ','DN      ','DNH     ',
     3               'DND     ','CINFV   ','DINFH   ',
     4               'T       ','TH      ','TD      ',
     5               'O       ','OH      ',
     6               'I       ','IH      '/
      DATA CHRBLK   /'        '/
      DATA DRCCHR   /'NORMAL  ','PARALLEL'/
      DATA TOL      /1.0D-10/
      DATA PI2      /6.28318530717958D+00/
      DATA ZERO     /0.0D+00/
      DATA PT5      /0.5D+00/
      DATA ONE      /1.0D+00/
      DATA THREE    /3.0D+00/
C
      OUT=.FALSE.
C
C     ----- GROUP INFO -----
C
      LINEAR=0
C
      IGROUP=20
      DO I=1,19
         IF(GROUP.EQ.GRP(I)) IGROUP=I
      ENDDO
      IF(NAXIS.LE.0) NAXIS=1
      IF(OUT) THEN
         WRITE(IW,9991) GROUP,NAXIS
         IF(LINEAR.NE.0) THEN
            WRITE(IW,9978)
         ENDIF
      ENDIF
      IF(IGROUP.GT.19) THEN
         WRITE(IW,9994)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(IGROUP.EQ.18.OR.IGROUP.EQ.19) THEN
         WRITE(IW,9998)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(IGROUP.LE.3                 ) GO TO 200
      IF(IGROUP.EQ.11.OR.IGROUP.EQ.12) GO TO 200
C
C     ----- DEFINE LOCAL FRAME
C     ----- READ IN PRINCIPAL AXIS   ( 1 CARD )
C     ----- READ IN X-LOCAL AXIS     ( 1 CARD )
C     ----- DEFAULT OPTION _ LOCAL FRAME IDENTICAL TO MASTER FRAME
C
      X0=XPT0
      Y0=YPT0
      Z0=ZPT0
      X1=XPT1
      Y1=YPT1
      Z1=ZPT1
C
      RHO= SQRT((X1-X0)**2+(Y1-Y0)**2+(Z1-Z0)**2)
      IF(RHO.LE.TOL) THEN
         X0=ZERO
         Y0=ZERO
         Z0=ZERO
         X1=ZERO
         Y1=ZERO
         Y2=ZERO
         Z2=ZERO
         Z1=ONE
         X2=ONE
         DIRECT=DRC(2)
         RHO=ONE
      ELSE
         X2=XPT2
         Y2=YPT2
         Z2=ZPT2
         DIRECT=BLANK
      ENDIF
      IF(DIRECT.NE.DRC(1)) DIRECT=DRC(2)
      ALIGN=X0.EQ.ZERO.AND.Y0.EQ.ZERO.AND.Z0.EQ.ZERO.AND.
     1      X1.EQ.ZERO.AND.Y1.EQ.ZERO               .AND.
     2                     Y2.EQ.ZERO.AND.Z2.EQ.ZERO.AND.
     3      DIRECT.EQ.DRC(2)
      W1=(X1-X0)/RHO
      W2=(Y1-Y0)/RHO
      W3=(Z1-Z0)/RHO
      WW=W1*W1+W2*W2+W3*W3
      X02=X2-X0
      Y02=Y2-Y0
      Z02=Z2-Z0
      RHO=(W1*X02+W2*Y02+W3*Z02)/WW
      DUM=RHO*W1
      X0=X0+DUM
      X02=X02-DUM
      DUM=RHO*W2
      Y0=Y0+DUM
      Y02=Y02-DUM
      DUM=RHO*W3
      Z0=Z0+DUM
      Z02=Z02-DUM
      UU=(X02*X02+Y02*Y02+Z02*Z02)
      U= SQRT(UU)
      U1=X02/U
      U2=Y02/U
      U3=Z02/U
      V3=W1*U2-W2*U1
      V2=W3*U1-W1*U3
      V1=W2*U3-W3*U2
      IF(DIRECT.NE.DRC(2)) THEN
         DUM=U1
         U1=V1
         V1=-DUM
         DUM=U2
         U2=V2
         V2=-DUM
         DUM=U3
         U3=V3
         V3=-DUM
      ENDIF
      IF(OUT) THEN
         WRITE(IW,9990) X0,Y0,Z0,U1,V1,W1,U2,V2,W2,U3,V3,W3
      ENDIF
      IF(IGROUP.GE.13) GO TO 200
C
C     ----- ROTATION ABOUT PRINCIPAL AXIS
C
      NN=0
      N=NAXIS
      ALPHA=ZERO
      ALPH=PI2/ DBLE(N)
   10 NN=NN+1
      IF(NN.GT.N) GO TO 20
      COSA=COS(ALPHA)
      SINA=SIN(ALPHA)
      I=9*(NN-1)
      T(I+1)=COSA
      T(I+5)=COSA
      T(I+2)=-SINA
      T(I+4)=SINA
      T(I+3)=ZERO
      T(I+6)=ZERO
      T(I+7)=ZERO
      T(I+8)=ZERO
      T(I+9)=ONE
      ALPHA=ALPHA+ALPH
      GO TO 10
C
C     ----- END OF GROUP 4
C
   20 NT=N
      II=9*NT
      IF(OUT) THEN
         WRITE(IW,9989)
         N1=1
         N2=NAXIS
         CALL HND_PRTSYM(N1,N2)
      ENDIF
      IF(IGROUP.EQ.4                ) GO TO 1000
      IF(IGROUP.EQ.5                ) GO TO 500
      IF(IGROUP.EQ.7                ) GO TO 115
      IF(IGROUP.NE.6.AND.IGROUP.NE.9) GO TO 55
C
C     ----- SIGMA-H PLANE  EQUATION (Z=0) IN LOCAL FRAME
C
      NN=0
   30 NN=NN+1
      IF(NN.GT.NT) GO TO 50
C
C     ----- GROUP 6 0R 9
C
      I=II+9*(NN-1)
      DO 40 J=1,8
   40 T(I+J)=T(I+J-II)
      T(I+9)=-T(I+9-II)
      GO TO 30
   50 NT=NT+NT
      II=9*NT
      IF(OUT) THEN
         WRITE(IW,9988)
         N1=N2+1
         N2=N2+NAXIS
         CALL HND_PRTSYM(N1,N2)
      ENDIF
C
C     ----- END OF GROUP 6
C
      IF(IGROUP.EQ.6) GO TO 1000
C
C     ----- ONE CP2 AXIS IS THE X-AXIS OF THE LOCAL FRAME
C     ----- GROUP 8 , 9 ,10
C
   55 CONTINUE
      NN=0
   60 NN=NN+1
      IF(NN.GT.NT) GO TO 70
      I=II+9*(NN-1)
      T(I+1)=T(I+1-II)
      T(I+2)=-T(I+2-II)
      T(I+3)=-T(I+3-II)
      T(I+4)=T(I+4-II)
      T(I+5)=-T(I+5-II)
      T(I+6)=-T(I+6-II)
      T(I+7)=T(I+7-II)
      T(I+8)=-T(I+8-II)
      T(I+9)=-T(I+9-II)
      GO TO 60
   70 NT=NT+NT
      II=9*NT
      IF(OUT) THEN
         WRITE(IW,9987)
         N1=N2+1
         N2=N2+NAXIS
         CALL HND_PRTSYM(N1,N2)
         IF(IGROUP.EQ.9) THEN
            WRITE(IW,9981)
            N1=N2+1
            N2=N2+NAXIS
            CALL HND_PRTSYM(N1,N2)
         ENDIF
      ENDIF
C
C     ----- END OF GROUP 8 AND 9
C
      IF(IGROUP.EQ.8.OR.IGROUP.EQ.9) GO TO 1000
C
C     ----- DND GROUP . EQUATION OF PLANE SIGMA-D IS _
C     ----- SIN(ALPH/4)*X-COS(ALPH/4)*Y=0
C     ----- THE X-AXIS IS THE CP2 AXIS.
C
C     ----- GROUP 10
C
      BETA=PT5*ALPH
      COSA=COS(BETA)
      SINA=SIN(BETA)
      NN=0
  100 NN=NN+1
      IF(NN.GT.NT) GO TO 110
      I=II+9*(NN-1)
      T(I+1)=COSA*T(I+1-II) + SINA*T(I+2-II)
      T(I+2)=SINA*T(I+1-II) - COSA*T(I+2-II)
      T(I+3)=     T(I+3-II)
      T(I+4)=COSA*T(I+4-II) + SINA*T(I+5-II)
      T(I+5)=SINA*T(I+4-II) - COSA*T(I+5-II)
      T(I+6)=     T(I+6-II)
      T(I+7)=COSA*T(I+7-II) + SINA*T(I+8-II)
      T(I+8)=SINA*T(I+7-II) - COSA*T(I+8-II)
      T(I+9)=     T(I+9-II)
      GO TO 100
  110 NT=NT+NT
      II=9*NT
      IF(OUT) THEN
         WRITE(IW,9986)
         N1=N2+1
         N2=N2+NAXIS
         CALL HND_PRTSYM(N1,N2)
         WRITE(IW,9980)
         N1=N2+1
         N2=N2+NAXIS
         CALL HND_PRTSYM(N1,N2)
      ENDIF
C
C     ----- END OF GROUP 10
C
      GO TO 1000
C
C     ----- GROUP 7
C     ----- SIGMA-V IS THE (X-Z) PLANE OF LOCAL FRAME
C
  115 CONTINUE
      NN=0
  120 NN=NN+1
      IF(NN.GT.NT) GO TO 130
      I=II+9*(NN-1)
      T(I+1)=T(I+1-II)
      T(I+2)=-T(I+2-II)
      T(I+3)=T(I+3-II)
      T(I+4)=T(I+4-II)
      T(I+5)=-T(I+5-II)
      T(I+6)=T(I+6-II)
      T(I+7)=T(I+7-II)
      T(I+8)=-T(I+8-II)
      T(I+9)=T(I+9-II)
      GO TO 120
  130 NT=NT+NT
      II=9*NT
C
C     ----- END OF GROUP 7
C
      IF(OUT) THEN
         WRITE(IW,9985)
         N1=N2+1
         N2=N2+NAXIS
         CALL HND_PRTSYM(N1,N2)
      ENDIF
      GO TO 1000
C
  200 CONTINUE
      T(1)=ONE
      T(5)=ONE
      T(9)=ONE
      T(2)=ZERO
      T(3)=ZERO
      T(4)=ZERO
      T(6)=ZERO
      T(7)=ZERO
      T(8)=ZERO
      IF(IGROUP.EQ.1                 ) GO TO 210
      IF(IGROUP.EQ.2                 ) GO TO 250
      IF(IGROUP.EQ.3                 ) GO TO 300
      IF(IGROUP.EQ.11.OR.IGROUP.EQ.12) GO TO 400
      GO TO 600
  210 NT=1
      ALIGN=.TRUE.
      X0=ZERO
      Y0=ZERO
      Z0=ZERO
      U1=ONE
      V2=ONE
      W3=ONE
      U2=ZERO
      U3=ZERO
      V1=ZERO
      V3=ZERO
      W1=ZERO
      W2=ZERO
      GO TO 1000
C
C     ----- CS SYMMETRY GROUP
C     ----- THE 3 POINTS 1,2,3 DEFINE SIGMA-H PLANE
C
  250 CONTINUE
C
      X1=XPT1
      Y1=YPT1
      Z1=ZPT1
      X2=XPT2
      Y2=YPT2
      Z2=ZPT2
      RHO=(X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2
      IF(RHO.LE.TOL) THEN
C
C     ----- DEFAULT OPTION _ PLANE IS THE (X,Y) PLANE
C
         X1=ZERO
         Y1=ZERO
         Z1=ZERO
         Y2=ZERO
         Z2=ZERO
         X3=ZERO
         Z3=ZERO
         X2=ONE
         Y3=ONE
      ELSE
         X3=XPT3
         Y3=YPT3
         Z3=ZPT3
      ENDIF
      ALIGN=X1.EQ.ZERO.AND.Y1.EQ.ZERO.AND.Z1.EQ.ZERO.AND.
     1                     Y2.EQ.ZERO.AND.Z2.EQ.ZERO.AND.
     2      X3.EQ.ZERO.AND.               Z3.EQ.ZERO
      NT=2
      W1=(Y2-Y1)*(Z3-Z1)-(Y3-Y1)*(Z2-Z1)
      W2=(Z2-Z1)*(X3-X1)-(Z3-Z1)*(X2-X1)
      W3=(X2-X1)*(Y3-Y1)-(X3-X1)*(Y2-Y1)
      RHO= SQRT(W1*W1+W2*W2+W3*W3)
      W1=W1/RHO
      W2=W2/RHO
      W3=W3/RHO
      U1=X2-X1
      U2=Y2-Y1
      U3=Z2-Z1
      RHO= SQRT(U1*U1+U2*U2+U3*U3)
      U1=U1/RHO
      U2=U2/RHO
      U3=U3/RHO
      V1=W2*U3-W3*U2
      V2=W3*U1-W1*U3
      V3=W1*U2-W2*U1
      X0=X1
      Y0=Y1
      Z0=Z1
      T(10)=ONE
      T(14)=ONE
      T(18)=-ONE
      T(11)=ZERO
      T(12)=ZERO
      T(13)=ZERO
      T(15)=ZERO
      T(16)=ZERO
      T(17)=ZERO
      IF(OUT) THEN
         WRITE(IW,9983) W1,W2,W3
         WRITE(IW,9982) U1,V1,W1,U2,V2,W2,U3,V3,W3
      ENDIF
      GO TO 1000
C
C     ----- CI SYMMETRY GROUP
C     ----- CENTER OF INVERSION IS (X0,Y0,Z0)
C
  300 CONTINUE
      X0=XPT0
      Y0=YPT0
      Z0=ZPT0
C
      ALIGN=X0.EQ.ZERO.AND.Y0.EQ.ZERO.AND.Z0.EQ.ZERO
      IF(OUT) THEN
         WRITE(IW,9984) X0,Y0,Z0
      ENDIF
      T(10)=-ONE
      T(14)=-ONE
      T(18)=-ONE
      T(11)=ZERO
      T(12)=ZERO
      T(13)=ZERO
      T(15)=ZERO
      T(16)=ZERO
      T(17)=ZERO
      NT=2
      U1=ONE
      V2=ONE
      W3=ONE
      U2=ZERO
      U3=ZERO
      V1=ZERO
      V3=ZERO
      W1=ZERO
      W2=ZERO
      GO TO 1000
  400 WRITE(IW,9995)
      CALL HND_HNDERR(3,ERRMSG)
  500 NN=0
      BETA=PT5*ALPH
      COSB=COS(BETA)
      SINB=SIN(BETA)
  510 NN=NN+1
      IF(NN.GT.NT) GO TO 520
C
C     ----- S2N GROUP
C     ----- THE PLANE OF SYMMETRY FOR THE IMPROPER ROTATION
C     ----- IS THE (X,Y) PLANE OF THE LOCAL FRAME
C
      I=II+9*(NN-1)
      T(I+1)= COSB*T(I+1-II)+SINB*T(I+2-II)
      T(I+2)=-SINB*T(I+1-II)+COSB*T(I+2-II)
      T(I+3)=     -T(I+3-II)
      T(I+4)= COSB*T(I+4-II)+SINB*T(I+5-II)
      T(I+5)=-SINB*T(I+4-II)+COSB*T(I+5-II)
      T(I+6)=     -T(I+6-II)
      T(I+7)= COSB*T(I+7-II)+SINB*T(I+8-II)
      T(I+8)=-SINB*T(I+7-II)+COSB*T(I+8-II)
      T(I+9)=     -T(I+9-II)
      GO TO 510
  520 NT=NT+NT
      II=9*NT
      IF(OUT) THEN
         WRITE(IW,9979)
         N1=N2+1
         N2=N2+NAXIS
         CALL HND_PRTSYM(N1,N2)
      ENDIF
      GO TO 1000
C
C     ----- T GROUP AND OTHERS CONTAINING A SUBGROUP T _
C     ----- LOCAL X,Y,Z ARE THE C2 AXES
C
  600 DO 610 I=10,36
  610    T(I)=ZERO
      T(10)=ONE
      T(23)=ONE
      T(36)=ONE
      T(14)=-ONE
      T(18)=-ONE
      T(19)=-ONE
      T(27)=-ONE
      T(28)=-ONE
      T(32)=-ONE
      DO 620 II=5,12
         I=9*(II-1)
         J=9*(II-5)
         T(I+1)=T(J+7)
         T(I+2)=T(J+8)
         T(I+3)=T(J+9)
         T(I+4)=T(J+1)
         T(I+5)=T(J+2)
         T(I+6)=T(J+3)
         T(I+7)=T(J+4)
         T(I+8)=T(J+5)
  620    T(I+9)=T(J+6)
      NT=12
      IF(IGROUP.EQ.13) GO TO 1000
      IF(IGROUP.EQ.14) GO TO 650
      IF(IGROUP.EQ.15) GO TO 680
      GO TO 670
C
C     ----- TH GROUP
C     ----- EXPAND GROUP BY TAKING DIRECT PRODUCT WITH CI
C
  650 I=9*NT
      DO 660 J=1,I
  660    T(J+I)=-T(J)
      NT=NT+NT
      GO TO 1000
C
C     ----- OCTAHEDRAL GROUP IS DIRECT PRODUCT OF T AND A C4 ROTATION
C     ----- ABOUT Z AXIS
C
  670 SIGN=-ONE
      GO TO 685
C
C     ----- TD GROUP IS DIRECT PRODUCT OF T AND A REFLECTION IN A
C     ----- PLANE ( EQUATION OF THE PLANE   X=Y  )
C
  680 SIGN=ONE
  685 DO 690 II=13,24
         I=9*(II-1)
         J=9*(II-13)
         T(I+1)=T(J+4)*SIGN
         T(I+2)=T(J+5)*SIGN
         T(I+3)=T(J+6)*SIGN
         T(I+4)=T(J+1)
         T(I+5)=T(J+2)
         T(I+6)=T(J+3)
         T(I+7)=T(J+7)
         T(I+8)=T(J+8)
  690    T(I+9)=T(J+9)
      NT=24
      IF(IGROUP.NE.17) GO TO 1000
C
C     ----- OH GROUP IS DIRECT PRODUCT OF O AND CI
C
      I=9*NT
      DO 700 J=1,I
  700    T(J+I)=-T(J)
      NT=48
 1000 CONTINUE
C
C     ----- FIND THE INVERSE TRANSFORMATIONS
C
      DO 1002 ITR=1,NT
         NN=9*(ITR-1)
         DO 1001 IT=1,NT
            II=9*(IT-1)
            TEST= T(NN+1)*T(II+1)+T(NN+2)*T(II+4)+T(NN+3)*T(II+7)
     1           +T(NN+4)*T(II+2)+T(NN+5)*T(II+5)+T(NN+6)*T(II+8)
     1           +T(NN+7)*T(II+3)+T(NN+8)*T(II+6)+T(NN+9)*T(II+9)
     1           -THREE
            IF( ABS(TEST).GT.TOL) GO TO 1001
            INVT(ITR)=IT
            GO TO 1002
 1001       CONTINUE
 1002    CONTINUE
C
C     ----- GENERATE TRANSFORMATION MATRICES -----
C           FOR  P, D, F, G  BASIS FUNCTIONS
C
      CALL HND_MATSYM
C
      NTMAX=NT
      IF(NTMAX.EQ.1) THEN
         NOSYM=1
      ELSE
         NOSYM=0
      ENDIF
      IF(OUT) THEN
         WRITE(IW,*) 'NT,NTMAX,NOSYM = ',NT,NTMAX,NOSYM
      ENDIF
C
      RETURN
 9999 FORMAT(A5,2I5)
 9998 FORMAT(' ARE YOU KIDDING... YOU SHOULD GIVE UP...')
 9997 FORMAT(9F10.5)
 9996 FORMAT(3F10.5,A8)
 9995 FORMAT(' LINEAR MOLECULE , POINT GROUP IS CINFV OR DINFH ',/,
     1       ' PLEASE USE GROUP C4V OR D4H...')
 9994 FORMAT(' ILLEGAL POINT GROUP INPUT. STOP..')
 9991 FORMAT(/,' THE POINT GROUP OF THE MOLECULE IS ...',A8,/,
     1         ' THE ORDER OF THE PRINCIPAL AXIS IS ...',I5)
 9990 FORMAT(/,' THE ORIGIN OF THE LOCAL FRAME IS AT X =  ',F10.5,
     1         ' Y = ',F10.5,' Z = ',F10.5,/,
     2         ' DIRECTOR COSINES OF THE NEW AXES ',/,
     3         34X,3(5X,F10.5),/,34X,
     4         3(5X,F10.5),/,34X,3(5X,F10.5))
 9989 FORMAT(' ROTATIONS ABOUT PRINCIPAL AXIS')
 9988 FORMAT(' SIGMA-H FOLLOWED BY ROTATIONS')
 9987 FORMAT(' C2 FOLLOWED BY ROTATIONS ')
 9986 FORMAT(' SIGMA-D FOLLOWED BY ROTATIONS')
 9985 FORMAT(' SIGMA-V FOLLOWED BY ROTATIONS')
 9984 FORMAT(/,10X,' CENTER OF SYMMETRY AT X = ',F10.5,' Y = ',F10.5,
     1             ' Z = ',F10.5)
 9983 FORMAT(/,' PLANE OF SYMMETRY DEFINED BY ITS NORMAL U = ',F10.5,
     1         ' V = ',F10.5,' W = ',F10.5)
 9982 FORMAT(/,10X,3F15.9,/,10X,3F15.9,/,10X,3F15.9)
 9981 FORMAT(' C2 FOLLOWED BY SIGMA-H FOLLOWED BY ROTATIONS')
 9980 FORMAT(' SIGMA-D FOLLOWED BY C2 FOLLOWED BY ROTATIONS')
 9979 FORMAT(' S2N ROTATION FOLLOWED BY ROTATIONS')
 9978 FORMAT(' THE MOLECULE IS LINEAR ')
      END
      SUBROUTINE HND_MOLMAP(C,C3,NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
      SUBROUTINE HND_SYMATM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
      SUBROUTINE HND_MATSYM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
      SUBROUTINE HND_PRTSYM(N1,N2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXSYM =120)
      PARAMETER (MXSYMT=120*9)
      PARAMETER (MXIODA=255)
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_DAFILE/IDAF,NAV,IODA(MXIODA)
      COMMON/HND_SYMTRY/INVT(MXSYM),NT,NTMAX,NTWD,NOSYM
      COMMON/HND_SYMMAT/T(MXSYMT)
      DIMENSION NN(MXSYM)
      IMAX=N1-1
  100 IMIN=IMAX+1
      IMAX=IMAX+4
      IF(IMAX.GT.N2) IMAX=N2
      NJ=9*N1-8
      DO 200 J=1,3
      NI=0
      DO 150 I=IMIN,IMAX
      NN(I)=NJ+NI
  150 NI=NI+9
      WRITE(IW,1000) (T(NN(I)),T(NN(I)+1),T(NN(I)+2),I=IMIN,IMAX)
  200 NJ=NJ+3
      WRITE(IW,1001)
      IF(IMAX.LT.N2) GO TO 100
 1000 FORMAT(4X,4(3F10.5,2H '))
 1001 FORMAT(/)
      RETURN
      END

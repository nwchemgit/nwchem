      subroutine geom_input(rtdb)
C     $Id: geom_input.F,v 1.87 1999-04-08 19:21:35 d3g681 Exp $
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
      oautoz   = .false.        ! Default is don't try -autoz-
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
            call geom_autoz_input(rtdb)
            goto 20
         else if (inp_compare(.false.,'zmt', field).or.
     $           inp_compare(.false.,'zmat', field).or.
     $           inp_compare(.false.,'zmatr', field).or.
     $           inp_compare(.false.,'zmatri', field).or.
     $           inp_compare(.false.,'zmatrix', field)) then
            ozmatrix=.true.   
            call geom_zmt_input(geom,coords,tags,charge,mass,ncenter,
     $           nzinput,scale)
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
c     geom_compute_values is called in the very end of 
c     sym_geom_project. Current call of geom_compute_values 
c     is commented out. 
c
c     call geom_compute_values(geom)
c     
      if(oautoz) then
         if (.not. geom_cart_get(geom, ncenter, tags, coords, charge))
     $        call errquit('geom_input: geom_cart_get failed', 0)
c     
c     ----- call auto-z-matrix generation -----
c     
         call geom_zmt_autoz(rtdb,geom,coords,charge,tags,ncenter,
     $        nzinput,oforce_autoz)
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
      if (.not. geom_destroy(geom))
     $     call errquit('geom_input: geom_destroy failed', 0)
c     
c     done
c     

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
      subroutine geom_autoz_input(rtdb)
      implicit double precision (a-h,o-z)
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer  rtdb
      logical  status
      logical  zcoord_input
c
      parameter (max_zcoord=16)
      integer ir
      integer iw
      common/hnd_iofile/ir,iw
      dimension ijbond(2,max_zcoord)
      dimension ijkang(3,max_zcoord)
      dimension ijklto(4,max_zcoord)
      dimension ijklop(4,max_zcoord)
      dimension ijklnb(4,max_zcoord)
c
      ir=LuIn
      iw=LuOUT
c
      if(ga_nodeid().eq.0) then
c
         call hnd_autoz_input(ijbond,ijkang,ijklto,
     1                        ijklop,ijklnb,max_zcoord,
     2                        cvr_scaling)
c
         zcoord_input=.true.
         status=rtdb_put(rtdb,'autoz:zcoord',mt_log,1,zcoord_input)
         status=status.or.
     1          rtdb_put(rtdb,'autoz:cvr_scaling',mt_dbl,1,  
     2                                                 cvr_scaling)
         status=status.or.
     1          rtdb_put(rtdb,'autoz:ijbond',mt_int,2*max_zcoord,  
     2                                                    ijbond)
         status=status.or.
     1          rtdb_put(rtdb,'autoz:ijkang',mt_int,3*max_zcoord,  
     2                                                    ijkang)
         status=status.or.
     1          rtdb_put(rtdb,'autoz:ijklto',mt_int,4*max_zcoord,  
     2                                                    ijklto)
         status=status.or.
     1          rtdb_put(rtdb,'autoz:ijklop',mt_int,4*max_zcoord,  
     2                                                    ijklop)
         status=status.or.
     1          rtdb_put(rtdb,'autoz:ijklnb',mt_int,4*max_zcoord,  
     2                                                    ijklnb)
         if(.not.status) then
            call errquit('zcoord_input : rtdb_put failed',911)
         endif
      endif
c
      return
      end
      subroutine hnd_autoz_input(ijbond,ijkang,ijklto,
     1                           ijklop,ijklnb,max_zcoord,
     2                           cvr_scaling)
      implicit double precision (a-h,o-z)
#include "inp.fh"
c
      parameter (mxatom=500)
      parameter (mxcart=3*mxatom)
      parameter (mxzmat=1500)
      parameter (mxcoor=1500)
c                mxcoor=max(mxcart,mxzmat)
      parameter (cvr_factor=1.20d+00)
      character*255 field
      logical       ready
      logical       dbug
      logical       out 
      common/hnd_iofile/ir,iw
      dimension ijbond(2,*)
      dimension ijkang(3,*)
      dimension ijklto(4,*)
      dimension ijklop(4,*)
      dimension ijklnb(4,*)
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
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
      call inp_set_field(0)
      ifield=inp_cur_field()
      if(dbug) then
         write(iw,*) 'in hnd_autoz_input, ifield = ',ifield
      endif
      ready=.false.
c
   10 continue
      if(.not.inp_a(field)) then
c
         if(.not.inp_read()) then
            write(iw,*)  'unexpected eof while reading ',
     1                   'hnd_autoz_input data'
            call errquit('hnd_autoz_input: unexpected eof',911)
         else
            if(dbug) then
               write(iw,*) 'new line read in hnd_autoz_input'
            endif
            go to 10
         endif
c
      else
         if(dbug) then
            write(iw,*) 'in hnd_autoz_input, field = ',field
         endif
         if(.not.ready) then
            if(inp_compare(.false.,'zcoord',field)) then
               ready=.true.
               go to 10
            else
               write(iw,*) 'un-expected value of field'
               call errquit('hnd_autoz_input: misposition',911)
            endif
         endif
c
         if(inp_compare(.false.,'zend',field).or.   
     1      inp_compare(.false.,'end',field)     ) then
            go to 100
         else
            if(inp_compare(.false.,'cvr_scaling',field)) then
   20          continue
               if(inp_f(dum)) then
                  cvr_scaling=dum
                  go to 10
               else
                  if(.not.inp_read()) then
                     write(iw,*)  'unexpected eof while reading ',
     1                            'cvr_scaling input'
                     call errquit
     1                    ('hnd_autoz_input: unexpected eof',911)
                  else
                     go to 20
                  endif
               endif
            elseif(inp_compare(.false.,'ijbond',field)) then
   30          continue
               if(inp_i(idum)) then
   40             continue
                  if(inp_i(jdum)) then
                     nbond=nbond+1
                     ijbond(1,nbond)=idum
                     ijbond(2,nbond)=jdum
                     go to 30
                  else
                     if(.not.inp_read()) then
                        write(iw,*)  'unexpected eof while reading ',
     1                               'ijbond data'
                        call errquit
     1                       ('hnd_autoz_input: unexpected eof',911)
                     else
                        go to 40
                     endif
                  endif
               else
                  go to 10
               endif
            elseif(inp_compare(.false.,'ijkang',field)) then
  130          continue
               if(inp_i(idum)) then
  140             continue
                  if(inp_i(jdum)) then
  150                continue
                     if(inp_i(kdum)) then
                        nangl=nangl+1
                        ijkang(1,nangl)=idum
                        ijkang(2,nangl)=jdum
                        ijkang(3,nangl)=kdum
                        go to 130
                     else
                        if(.not.inp_read()) then
                           write(iw,*)  'unexpected eof while reading ',
     1                                  'ijkang data'
                           call errquit
     1                          ('hnd_autoz_input: unexpected eof',911)
                        else
                           go to 150
                        endif
                     endif
                  else
                     if(.not.inp_read()) then
                        write(iw,*)  'unexpected eof while reading ',
     1                               'ijkang data'
                        call errquit
     1                       ('hnd_autoz_input: unexpected eof',911)
                     else
                        go to 140
                     endif
                  endif
               else 
                  go to 10
               endif
            elseif(inp_compare(.false.,'ijklto',field)) then
  230          continue
               if(inp_i(idum)) then
  240             continue
                  if(inp_i(jdum)) then
  250                continue
                     if(inp_i(kdum)) then
  260                   continue
                        if(inp_i(ldum)) then
                           ntors=ntors+1
                           ijklto(1,ntors)=idum
                           ijklto(2,ntors)=jdum
                           ijklto(3,ntors)=kdum
                           ijklto(4,ntors)=ldum
                           go to 230
                        else
                           if(.not.inp_read()) then
                           write(iw,*)  'unexpected eof while reading ',
     1                                  'ijklto data'
                           call errquit
     1                          ('hnd_autoz_input: unexpected eof',911)
                           else
                              go to 260
                           endif
                        endif
                     else
                        if(.not.inp_read()) then
                           write(iw,*)  'unexpected eof while reading ',
     1                                  'ijklto data'
                           call errquit
     1                          ('hnd_autoz_input: unexpected eof',911)
                        else
                           go to 250
                        endif
                     endif
                  else
                     if(.not.inp_read()) then
                        write(iw,*)  'unexpected eof while reading ',
     1                               'ijklto data'
                        call errquit
     1                       ('hnd_autoz_input: unexpected eof',911)
                     else
                        go to 240
                     endif
                  endif
               else
                  go to 10
               endif
            elseif(inp_compare(.false.,'ijklop',field)) then
  330          continue
               if(inp_i(idum)) then
  340             continue
                  if(inp_i(jdum)) then
  350                continue
                     if(inp_i(kdum)) then
  360                   continue
                        if(inp_i(ldum)) then
                           nopla=nopla+1
                           ijklop(1,nopla)=idum
                           ijklop(2,nopla)=jdum
                           ijklop(3,nopla)=kdum
                           ijklop(4,nopla)=ldum
                           go to 330
                        else
                           if(.not.inp_read()) then
                           write(iw,*)  'unexpected eof while reading ',
     1                                  'ijklop data'
                           call errquit
     1                          ('hnd_autoz_input: unexpected eof',911)
                           else
                              go to 360
                           endif
                        endif
                     else
                        if(.not.inp_read()) then
                           write(iw,*)  'unexpected eof while reading ',
     1                                  'ijklop data'
                           call errquit
     1                          ('hnd_autoz_input: unexpected eof',911)
                        else
                           go to 350
                        endif
                     endif
                  else
                     if(.not.inp_read()) then
                        write(iw,*)  'unexpected eof while reading ',
     1                               'ijklop data'
                        call errquit
     1                       ('hnd_autoz_input: unexpected eof',911)
                     else
                        go to 340
                     endif
                  endif
               else
                  go to 10
               endif
            elseif(inp_compare(.false.,'ijklnb',field)) then
  430          continue
               if(inp_i(idum)) then
  440             continue
                  if(inp_i(jdum)) then
  450                continue
                     if(inp_i(kdum)) then
  460                   continue
                        if(inp_i(ldum)) then
                           nlinb=nlinb+1
                           ijklnb(1,nlinb)=idum
                           ijklnb(2,nlinb)=jdum
                           ijklnb(3,nlinb)=kdum
                           ijklnb(4,nlinb)=ldum
                           go to 430
                        else
                           if(.not.inp_read()) then
                           write(iw,*)  'unexpected eof while reading ',
     1                                  'ijklnb data'
                           call errquit
     1                          ('hnd_autoz_input: unexpected eof',911)
                           else
                              go to 460
                           endif
                        endif
                     else
                        if(.not.inp_read()) then
                           write(iw,*)  'unexpected eof while reading ',
     1                                  'ijklnb data'
                           call errquit
     1                          ('hnd_autoz_input: unexpected eof',911)
                        else
                           go to 450
                        endif
                     endif
                  else
                     if(.not.inp_read()) then
                        write(iw,*)  'unexpected eof while reading ',
     1                               'ijklnb data'
                        call errquit
     1                       ('hnd_autoz_input: unexpected eof',911)
                     else
                        go to 440
                     endif
                  endif
               else
                  go to 10
               endif
            endif
         endif
c
      endif
c
  100 continue
      if(out) then
         write(iw,*) 'zcoord input :'
         write(iw,*) 'cvr_scaling', cvr_scaling
         if(nbond.gt.0) then
            do mbond=1,nbond
               write(iw,*) 'ijbond',
     1                      ijbond(1,mbond),ijbond(2,mbond) 
            enddo
         endif
         if(nangl.gt.0) then
            do mangl=1,nangl
               write(iw,*) 'ijkang',
     1                      ijkang(1,mangl),ijkang(2,mangl), 
     2                      ijkang(3,mangl)
            enddo
         endif
         if(ntors.gt.0) then
            do mtors=1,ntors
               write(iw,*) 'ijklto',
     1                      ijklto(1,mtors),ijklto(2,mtors), 
     2                      ijklto(3,mtors),ijklto(4,mtors)
            enddo
         endif
         if(nopla.gt.0) then
            do mopla=1,nopla
               write(iw,*) 'ijklop',
     1                      ijklop(1,mopla),ijklop(2,mopla), 
     2                      ijklop(3,mopla),ijklop(4,mopla)
            enddo
         endif
         if(nlinb.gt.0) then
            do mlinb=1,nlinb
               write(iw,*) 'ijklnb',
     1                      ijklnb(1,mlinb),ijklnb(2,mlinb), 
     2                      ijklnb(3,mlinb),ijklnb(4,mlinb)
            enddo
         endif
      endif
c
      return
      end
      subroutine  geom_zmt_input(geom,
     1                           coords,tags,charge,mass,ncenter,
     2                           nzinput,units)
      implicit none
#include "stdio.fh"
#include "rtdb.fh"
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
      logical      geom_tag_to_element
      logical      geom_atn_to_default_mass
      logical      geom_zmt_put_izmat          
      logical      geom_zmt_put_nizmat          
      logical      geom_usr_put_izmat          
      logical      geom_usr_put_nizmat          
      logical      geom_usr_put_izfrz          
      logical      geom_usr_put_nzfrz           
      logical      dbug 
      logical      some 
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
      external     geom_tag_to_element
      external     geom_atn_to_default_mass
      external     geom_zmt_put_izmat          
      external     geom_zmt_put_nizmat          
      external     geom_usr_put_izmat          
      external     geom_usr_put_nizmat          
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
     2                  units)
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
c     ----- store -izmat- for -nwchem- -----
c
      if(.not.geom_zmt_put_nizmat(geom,nizmat))
     1   call errquit('geom_input: geom_zmt_put_nizmat failed', 0)
      if(.not.geom_zmt_put_izmat(geom,izmat,nizmat))
     1   call errquit('geom_input: geom_zmt_put_izmat failed', 0)
      if(.not.geom_usr_put_nizmat(geom,nizmat))
     1   call errquit('geom_input: geom_usr_put_nizmat failed', 0)
      if(.not.geom_usr_put_izmat(geom,izmat,nizmat))
     1   call errquit('geom_input: geom_usr_put_izmat failed', 0)
      if(.not.geom_usr_put_nzfrz(geom,nzfrz))
     1   call errquit('geom_input: geom_usr_put_nzfrz failed', 0)
      if(.not.geom_usr_put_izfrz(geom,izfrz,nzfrz))
     1   call errquit('geom_input: geom_usr_put_izfrz failed', 0)
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
      call geom_prt_autoz(zmat,nzvar,izmat,nizmat,c)
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
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geomP.fh"
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
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geomP.fh"
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
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geomP.fh"
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
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geomP.fh"
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
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geomP.fh"
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
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geomP.fh"
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
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geomP.fh"
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
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geomP.fh"
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
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geomP.fh"
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
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geomP.fh"
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
     2                           nzinput,oforce_autoz)
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
      logical      geom_zmt_put_izmat          
      logical      geom_zmt_put_nizmat          
      logical      oforce_autoz
      logical      zdone
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
      call hnd_zmtmak(zdone,rtdb)
c
      if(.not.zdone) then
         write(iw,9997)
         nzmat=0
      endif
c
c     ----- save -izmat- and -zmat- -----
c
      if(.not.geom_zmt_put_nizmat(geom,nizmat))
     1   call errquit('geom_input: geom_zmt_put_nizmat failed', 0)
      if(.not.geom_zmt_put_izmat(geom,izmat,nizmat))
     1   call errquit('geom_input: geom_zmt_put_izmat failed', 0)
      if(dbug) then
         write(iw,9993)
         write(iw,9996)
         write(iw,9995) (izmat(i),i=1,nizmat)
         write(iw,9994)
      endif
c
c     ----- print -zmat- data -----
c
      call geom_prt_autoz(zmat,nzvar,izmat,nizmat,c)
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
#include "rtdb.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "mafdecls.fh"
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
      ELSEIF(ITYPE.EQ.2) THEN
C
C     ----- ANGLE BEND -----
C
         CALL HND_BEND(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),
     2             C,B,NCART,RTOD,BOHR)
         ZMAT(II) = EQVAL
         IADD = IADD + 4
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
#include "msgtypesf.h"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      logical      status
      integer      rtdb
      integer      geom
      integer      mitob
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
      status = rtdb_parallel(.false.)
      if(ga_nodeid().eq.0) then
         if(.not.geom_zmt_get_nizmat(geom,nizmat))
     1   call errquit('geom_input: geom_zmt_get_nizmat failed',0)
         if(.not.geom_zmt_get_izmat(geom,izmat,nizmat))
     1   call errquit('geom_input: geom_zmt_get_izmat failed',0)
         call ga_brdcst(18+msgint,nizmat,mitob(1     ),0)
         call ga_brdcst(19+msgint,izmat ,mitob(mxizmt),0)
      else
         call ga_brdcst(18+msgint,nizmat,mitob(1     ),0)
         call ga_brdcst(19+msgint,izmat ,mitob(mxizmt),0)
      endif
      if(dbug) then
         write(iw,9998)
         write(iw,9997)
         write(iw,9996) (izmat(i),i=1,nizmat)
         write(iw,9995)
      endif
      call ga_sync()
      status = rtdb_parallel(.true.)
c
c     ----- call -icoord- -----
c
      ifzmat=1
      call hnd_icoord(rtdb,geom)
      return
 9999 format(' in -geom_zmt_icoord ')
 9998 format(' reading -izmat- via -geom_zmt_get_izmat- ')
 9997 format(' $zmat izmat = ')
 9996 format(12(i5,','))
 9995 format(' $end ')
      end
      SUBROUTINE HND_ICOORD(RTDB,GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
c
      logical status
      integer rtdb
      integer geom
      logical  geom_zmt_opn_fil
C
C     ---                 CONSTRUCT -Z- MATRIX                     ---
C     ---             CALCULATE ATOM-ATOM DISTANCES                ---
C     --- -NZMAT-  =  NUMBER OF ELEMENTS OF -IZMAT- GIVEN IN $ZMAT ---
C     --- -NZVAR-  =  NUMBER OF INTERNAL COORDINATES FROM -IZMAT-  ---
C     --- -NVAR-   =  (3N-6) OR (3N-5) OR -NZMOD-  .LE. -NZVAR-    ---
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      LOGICAL ZDONE
      LOGICAL DBUG
      LOGICAL SOME
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_MOLLAB/ANAM(MXATOM),BNAM(MXATOM)
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOPT/IFZMAT,NATZMT          
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION R(MXATOM)
      DATA ONE   /1.0D+00/
      DATA UNIT  /0.52917715D+00/
C
      DBUG=.FALSE.
      SOME=.FALSE.
      SOME=SOME.OR.DBUG
      IF(DBUG) THEN
         WRITE(IW,*) 'IN HND_ICOORD'
         WRITE(IW,*) 'IFZMAT',IFZMAT
      ENDIF
C
C     ----- in -nwchem- , some MUST BE true -----
C
      SOME=.TRUE.
C
      IF(IFZMAT.NE.0) THEN
C
C     ----- LOOK FOR -Z- MATRIX DATA AND READ IT IN -----
C
         IF(SOME) THEN
            CALL HND_RDZMAT(NAT)
         ENDIF
      ENDIF
C
      IF(DBUG) THEN
         WRITE(IW,*) 'IFZMAT',IFZMAT
         WRITE(IW,*) 'NZMAT ',NZMAT
      ENDIF
      IF(IFZMAT.NE.0.AND.NZMAT.GT.0) THEN
C
C     -----    THERE IS A -Z- MATRIX    -----
C     ----- CALCULATE -B- AND -B**(-1)- -----
C
         STATUS=GEOM_ZMT_OPN_FIL(.TRUE.,.TRUE.,.FALSE.,.FALSE.)
         IF(.NOT.STATUS) THEN
            CALL ERRQUIT('GEOM_ZMT_OPN_FIL?',911)
         ENDIF
         CALL HND_BANDBI(GEOM,NAT,C)
         CALL HND_ZMTPRT(GEOM)
      ELSEIF(IFZMAT.NE.0.AND.NZMAT.EQ.0) THEN
C
C     -----       THERE IS NO -Z- MATRIX     -----
C     ----- TRY TO CREAT ONE ; IF SUCCESSFUL -----
C     -----    CALCULATE -B- AND -B**(-1)-   -----
C
         CALL HND_ZMTMAK(ZDONE,RTDB)
         IF(ZDONE) THEN
            STATUS=GEOM_ZMT_OPN_FIL(.TRUE.,.TRUE.,.FALSE.,.FALSE.)
            IF(.NOT.STATUS) THEN
               CALL ERRQUIT('GEOM_ZMT_OPN_FIL?',911)
            ENDIF
            CALL HND_BANDBI(GEOM,NAT,C)
            CALL HND_ZMTPRT(GEOM)
         ELSE
C
C     -----  DEFINE ROTATIONS+TRANSLATIONS  -----
C           BY MEANS OF HND_ECKARTs CONDITIONS
C
            NCART=3*NAT
            NZVAR=3*NAT
            NVAR =3*NAT
            NZMAT=0
            NCOOR=MAX(NCART,NZVAR)
            CALL HND_TRNROT(NAT,C,NCART,NCOOR)
C
C     ----- PRINT INTERNUCLEAR DISTANCES -----
C
            IF(SOME) THEN
               IPASS=1
               FAC=UNIT
               WRITE(IW,9995)
   10          CONTINUE
               JMAX=0
   20          JMIN=JMAX+1
               JMAX=JMAX+10
               IF(JMAX.GT.NAT) JMAX=NAT
               WRITE(IW,9999)
               WRITE(IW,9998) (ANAM(J),BNAM(J),J=JMIN,JMAX)
               WRITE(IW,9999)
               DO I=1,NAT
                  DO J=JMIN,JMAX
                     RR=(C(1,I)-C(1,J))**2+(C(2,I)-C(2,J))**2
     1                                    +(C(3,I)-C(3,J))**2
                     R(J)= SQRT(RR)*FAC
                  ENDDO
                  WRITE(IW,9997) I,ANAM(I),BNAM(I),(R(J),J=JMIN,JMAX)
               ENDDO
               IF(JMAX.LT.NAT) GO TO 20
               IPASS=2
               IF(IPASS.EQ.2) RETURN
               FAC=ONE
               WRITE(IW,9996)
               GO TO 10
            ENDIF
C
         ENDIF
      ELSEIF(IFZMAT.EQ.0) THEN
C
C     -----  DEFINE ROTATIONS+TRANSLATIONS  -----
C           BY MEANS OF HND_ECKARTs CONDITIONS
C
         NCART=3*NAT
         NZVAR=3*NAT
         NVAR =3*NAT
         NCOOR=MAX(NCART,NZVAR)
         CALL HND_TRNROT(NAT,C,NCART,NCOOR)
C
C     ----- PRINT INTERNUCLEAR DISTANCES -----
C
         IF(SOME) THEN
            IPASS=1
            FAC=UNIT
            WRITE(IW,9995)
  110       CONTINUE
            JMAX=0
  120       JMIN=JMAX+1
            JMAX=JMAX+10
            IF(JMAX.GT.NAT) JMAX=NAT
            WRITE(IW,9999)
            WRITE(IW,9998) (ANAM(J),BNAM(J),J=JMIN,JMAX)
            WRITE(IW,9999)
            DO I=1,NAT
               DO J=JMIN,JMAX
                  RR=(C(1,I)-C(1,J))**2+(C(2,I)-C(2,J))**2
     1                                 +(C(3,I)-C(3,J))**2
                  R(J)= SQRT(RR)*FAC
               ENDDO
               WRITE(IW,9997) I,ANAM(I),BNAM(I),(R(J),J=JMIN,JMAX)
            ENDDO
            IF(JMAX.LT.NAT) GO TO 120
            IPASS=2
            IF(IPASS.EQ.2) RETURN
            FAC=ONE
            WRITE(IW,9996)
            GO TO 110
         ENDIF
C
      ENDIF
      RETURN
 9999 FORMAT(/)
 9998 FORMAT(15X,10(1X,A8,A2))
 9997 FORMAT(1X,I3,1X,A8,A2,10(1X,F8.4,2X))
 9995 FORMAT(/,10X,'INTERNUCLEAR DISTANCES (ANGS.)',/,10X,30(1H-))
 9996 FORMAT(/,10X,'INTERNUCLEAR DISTANCES ( A.U.)',/,10X,30(1H-))
      END
      SUBROUTINE HND_TRNROT(NAT,C,NCART,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "geomP.fh"
C
      LOGICAL OUT 
      LOGICAL DBUG
      COMMON/HND_IOFILE/IR,IW
C
      DBUG=.TRUE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
C
C     ----- MEMORY POINTERS -----
C
C     I10=1
C     I20=I10+NCART*NCART
C     I30=I20+NCART*3     
C     I40=I30+NCART*3     
C     I50=I40+NDIM*NDIM  
C
      IF(.NOT.MA_PUSH_GET(MT_DBL,NCART*NCART,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_TRNROT, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NCART*3    ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_TRNROT, MALLOC OF I20  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NCART*3    ,'MEM I30 ',I_I30,I30))
     & CALL ERRQUIT('HND_TRNROT, MALLOC OF I30  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NDIM*NDIM  ,'MEM I40 ',I_I40,I40))
     & CALL ERRQUIT('HND_TRNROT, MALLOC OF I40  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I50 ',I_I50,I50))
     & CALL ERRQUIT('HND_TRNROT, MALLOC OF I50  FAILED',911)
C
C     ----- CREATE HND_ECKARTs PROJECTOR AND WRITE ON -94- -----
C
      CALL HND_ECKART(NAT,C,DBL_MB(I10),DBL_MB(I20),
     &                      DBL_MB(I30),DBL_MB(I20),NCART,DBUG)
C
      DO J=1,NCART
         DO I=1,NCART
            DBL_MB(I+(J-1)*NDIM+I40-1)=DBL_MB(I+(J-1)*NCART+I10-1)
         ENDDO
      ENDDO
C
      NDUM=NDIM*NDIM    
      IF(OUT) THEN
         WRITE(IW,9999)
         CALL HND_PRSQ(DBL_MB(I40),NCART,NCART,NDIM)
      ENDIF
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I50))
     & CALL ERRQUIT('HND_TRNROT, MA_POP_STACK OF I50  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I40))
     & CALL ERRQUIT('HND_TRNROT, MA_POP_STACK OF I40  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I30))
     & CALL ERRQUIT('HND_TRNROT, MA_POP_STACK OF I30  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_TRNROT, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_TRNROT, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
 9999 FORMAT(' IN -HND_TRNROT- PROJECTOR = ')
      END
      SUBROUTINE HND_ECKART(NAT,C,P,T,R,TR,NCART,OUT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,1)
      DIMENSION G(3)
      DIMENSION P(NCART,*),T(3,NAT,*),R(3,NAT,*),TR(NCART,*)
      DATA ZERO /0.0D+00/
      DATA ONE  /1.0D+00/
      DATA TOL  /1.0D-09/ 
C
      DO JTR=1,6
         DO ICART=1,NCART
            TR(ICART,JTR)=ZERO
         ENDDO
      ENDDO
C
C     ----- DEFINE CENTER OF MASS ( UNIT MASS ) -----
C
      DO IXYZ=1,3
         DUM=ZERO
         DO IAT=1,NAT
            DUM=DUM+C(IXYZ,IAT)
         ENDDO
         DUM=DUM/DBLE(NAT)
         G(IXYZ)=DUM
      ENDDO
C
C     ----- DEFINE TRANSLATIONS -----
C
      DO IXYZ=1,3
         DO IAT=1,NAT
            T(IXYZ,IAT,IXYZ)=ONE
         ENDDO
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9997)
         CALL HND_PRSQ(TR,3,NCART,NCART)
      ENDIF
C
C     ----- DEFINE ROTATIONS -----
C
      DO IXYZ=1,3
         IF(IXYZ.EQ.1) THEN
            DO IAT=1,NAT
               R(2,IAT,IXYZ)=-C(3,IAT)+G(3)
               R(3,IAT,IXYZ)= C(2,IAT)-G(2)
            ENDDO
         ELSEIF(IXYZ.EQ.2) THEN
            DO IAT=1,NAT
               R(1,IAT,IXYZ)= C(3,IAT)-G(3)
               R(3,IAT,IXYZ)=-C(1,IAT)+G(1)
            ENDDO
         ELSEIF(IXYZ.EQ.3) THEN
            DO IAT=1,NAT
               R(1,IAT,IXYZ)=-C(2,IAT)+G(2)
               R(2,IAT,IXYZ)= C(1,IAT)-G(1)
            ENDDO
         ENDIF
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9996)
         CALL HND_PRSQ(TR,6,NCART,NCART)
      ENDIF
C
C     ----- ORTHONORMALIZE -----
C
      NTR=0
      DO ITR=1,6
         DUM=dnrm2(NCART,TR(1,ITR),1)
         IF(ABS(DUM).GT.TOL) THEN
            NTR=NTR+1
            DUM=ONE/DUM
            CALL dscal(NCART,DUM,TR(1,ITR),1)
            IF(ITR.LT.6) THEN
               DO JTR=ITR+1,6
                  DUM= -ddot(NCART,    TR(1,ITR),1,TR(1,JTR),1)
                  CALL daxpy(NCART,DUM,TR(1,ITR),1,TR(1,JTR),1)
               ENDDO
            ENDIF
         ELSE
            DO ICART=1,NCART
               TR(ICART,ITR)=ZERO
            ENDDO
         ENDIF
      ENDDO
      IF(OUT.AND.NTR.NE.6) THEN
         WRITE(IW,9999) NTR
      ENDIF
      IF(OUT) THEN
         WRITE(IW,9995)
         CALL HND_PRSQ(TR,6,NCART,NCART)
      ENDIF
C
C     ----- CREATE PROJECTOR ONTO INTERNAL SPACE -----
C
      DO JCART=1,NCART
         DO ICART=1,NCART
            P(ICART,JCART)=ZERO
         ENDDO
         P(JCART,JCART)=ONE
         DO ICART=1,NCART
            DUM=ZERO
            DO ITR=1,6
               DUM=DUM+TR(ICART,ITR)*TR(JCART,ITR)
            ENDDO
            P(ICART,JCART)=P(ICART,JCART)-DUM
         ENDDO
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9998)
         CALL HND_PRSQ(P,NCART,NCART,NCART)
      ENDIF
C
      RETURN
 9999 FORMAT(/,' IN -HND_ECKART- , ONLY ',I1,' TRAN+ROT WERE FOUND .')
 9998 FORMAT(/,' HND_ECKART PROJECTOR',/,' ----------------')
 9997 FORMAT(/,' TRANSLATIONS',/,' ------------')
 9996 FORMAT(/,' TRANS. + ROT.',/,' -------------')
 9995 FORMAT(/,' TRANS. + ROT. (ORTHONORM)',/,' ',25(1h-))
      END
      SUBROUTINE HND_RDZMAT(NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- READ IN Z-MATRIX INPUT -----
C
      PARAMETER (    NA=10 )
      PARAMETER (MXNIJS=500)
      PARAMETER (MXIZMT=6000)
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)
      PARAMETER (MXCOOR=1500)
      LOGICAL OUT
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOPT/IFZMAT,NATZMT          
      COMMON/HND_ZMTDAT/ZMAT(MXCOOR),IZMAT(MXIZMT),NIZMAT
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION SIJ(MXNIJS),IJS(2,MXNIJS)
      DIMENSION A(3,NA)
      DIMENSION ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-RDZMAT-'/
      DATA LINEAR  /0/
      DATA NZMOD   /0/
      DATA NATURAL /2/
      DATA NPFLG   /0/
      DATA ZERO    /0.0D+00/
C
C     NAMELIST /ZMAT/ IZMAT,IJS,SIJ,A,LINEAR,NZMOD,NATURAL,NPFLG
C     NAMELIST /ZMT/  IZMAT,IJS,SIJ,A,LINEAR,NZMOD,NATURAL,NPFLG
C
      NZMAT=0
      NZVAR=3*NAT
      NVAR =3*NAT
      IF(NIZMAT.LE.0) THEN
         DO I=1,MXIZMT
            IZMAT(I)=0
         ENDDO
      ENDIF
      DO I=1,MXNIJS
         IJS(1,I)=0
         IJS(2,I)=0
         SIJ(I)=ZERO
      ENDDO
      DO J=1,NA
         DO I=1,3
            A(I,J)=ZERO
         ENDDO
      ENDDO
      NATZMT=-1
C
C     ----- READ NAMELIST -$ZMAT- OR NAMELIST -$ZMT- -----
C           BOTH PROVIDE THE SAME INPUT DATA. ONLY
C           $ZMAT HAS PRIORITY OVER $ZMT IF BOTH ARE
C           GIVEN IN INPUT DATA.
C
C     REWIND IR
C     READ(IR,ZMAT,END=100,ERR=100)
C     GO TO 110
C 100 CONTINUE
C     REWIND IR
C     READ(IR, ZMT,END=140,ERR=140)
C 110 CONTINUE
C
      IF(NPFLG.NE.1.AND.NPFLG.NE.2) THEN
         NPFLG=0
      ENDIF
      IPRZMT=NPFLG
      IF(NATURAL.NE.1.AND.NATURAL.NE.2) THEN
         NATURAL=2
      ENDIF
      IF(NATURAL.GT.0) THEN
         NATZMT=NATURAL
      ENDIF
C
      OUT=.FALSE.
      OUT=OUT.OR.NPRINT.EQ.1
      OUT=OUT.OR.IPRZMT.GE.1
      IF(OUT) THEN
         WRITE(IW,*) 'IN HND_RDZMAT'
      ENDIF
C
C     ----- INCLUDE SYMMETRICALLY EQUIVALENT -ZMAT- IF NEEDED -----
C
      CALL HND_ZMTSYM
C
C     ----- SCAN AND ANALYZE IZMAT -----
C
      NZMAT = 0
      DO I = 1,MXIZMT
         IF (IZMAT(I) .NE. 0) THEN
            NZMAT = NZMAT + 1
         ENDIF
      ENDDO
      IF(NZMAT.EQ.0) RETURN
C
      IF (OUT) WRITE (IW,9048)
      ICOUNT = 0
      NCOUNT = 0
  120 CONTINUE
      ICOUNT = ICOUNT + 1
      IF(ICOUNT.GT.NZMAT) GO TO 130
C
C     ----- CHECK FOR A VALID TYPE -----
C
C      1 ..... STRETCH
C      2 ..... BEND
C      3 ..... TORSION (DIHEDRAL ANGLE BETWEEN PLANES HAVING
C                       TWO ATOMS IN COMMON )
C      4 ..... OUT OF PLANE BEND
C      5 ..... LINEAR BEND
C      6 ..... DIHEDRAL ANGLE BETWEEN PLANES HAVING ONE ATOM
C              IN COMMON
C
      IF (IZMAT(ICOUNT) .LE. 0 .OR. IZMAT(ICOUNT) .GE. 7) THEN
         WRITE (IW,9008)
         WRITE (IW,9068) IZMAT(ICOUNT)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF (IZMAT(ICOUNT) .EQ. 1) INC = 2
      IF (IZMAT(ICOUNT) .EQ. 2) INC = 3
      IF (IZMAT(ICOUNT) .EQ. 3) INC = 4
      IF (IZMAT(ICOUNT) .EQ. 4) INC = 4
      IF (IZMAT(ICOUNT) .EQ. 5) INC = 4
      IF (IZMAT(ICOUNT) .EQ. 6) INC = 5
C
C     ----- MAKE SURE THE ATOM NUMBERS ARE VALID -----
C
      JDUM = INC
      IF (IZMAT(ICOUNT) .EQ. 5) JDUM = INC - 1
      DO J = 1,JDUM
         IF (IZMAT(ICOUNT+J) .GT. NAT) THEN
            WRITE (IW,9088) IZMAT(ICOUNT+J)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      IF (IZMAT(ICOUNT).EQ.5.AND.IZMAT(ICOUNT+INC).GT.NA) THEN
         WRITE (IW,9088) IZMAT(ICOUNT+INC)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF (IZMAT(ICOUNT).NE.5) NCOUNT=NCOUNT+1
      IF (IZMAT(ICOUNT).EQ.5) NCOUNT=NCOUNT+2
C
      IF (OUT) THEN
         WRITE (IW,9108) IZMAT(ICOUNT),(IZMAT(ICOUNT+J),J = 1,INC)
      ENDIF
      ICOUNT = ICOUNT + INC
      GO TO 120
  130 CONTINUE
C
C     ----- SAVE IZMAT ON DIRECT ACCESS -----
C     -----   -IZMAT-  ON -96-          -----
C     -----   -A-      ON -90-          -----
C
C
C     ----- SET -NVAR- AND -NZVAR- -----
C
      IF(LINEAR.EQ.0) NVAR=3*NAT-6
      IF(LINEAR.NE.0) NVAR=3*NAT-5
      IF( NZMOD.NE.0) NVAR=NZMOD
      NZVAR=NCOUNT
      IF(NZVAR.LT.NVAR) THEN
         WRITE(IW,9028)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(NZVAR.GT.MXZMAT) THEN
         WRITE(IW,9038)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- READ MATRIX FOR SYMMETRIC INTERNAL COORDINATES -----
C                 ( ALSO FOR LINEAR DEPENDENCIES )
C
      CALL HND_ZMTLIN(NZVAR,MXNIJS,IJS,SIJ)
C
C     ----- FINAL ADJUSTMENT OF -NATZMT- IF NEEDED -----
C
      IF(NATZMT.LT.0) THEN
         NATZMT=2
      ENDIF
      RETURN
C
  140 CONTINUE
C
C     ----- FINAL ADJUSTMENT OF -NATZMT- IF NEEDED -----
C
      IF(OUT) THEN
         WRITE(IW,9018)
      ENDIF
      IF(NATZMT.LT.0) THEN
         NATZMT=2
      ENDIF
      RETURN
 9108 FORMAT(5X,6I10)
 9088 FORMAT(' INVALID ATOM NUMBER ',I5,' . STOP')
 9068 FORMAT(' INVALID ZMAT TYPE ',I5,' . STOP')
 9048 FORMAT(/,10X,'INPUT Z-MATRIX ELEMENTS',
     1       /,5X,6X,'TYPE',9X,'I',9X,'J',9X,'K',9X,'L',9X,'M')
 9038 FORMAT(' NUMBER OF INTERNAL COORDINATES IS GREATER THAN ',
     1       ' -MXZMAT- . STOP ' )
 9028 FORMAT(' NUMBER OF INTERNAL COORDINATES IS LESS THAN (3N-6)',
     1       ' OR -NZMOD- .STOP')
 9018 FORMAT(' NO -$ZMAT- OR -$ZMT- NAMELISTS FOUND; PROCEED ...')
 9008 FORMAT(' ERROR IN $ZMAT INPUT . STOP')
      END
      SUBROUTINE HND_ZMTSYM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      RETURN
      END
      SUBROUTINE HND_ZMTLIN(NZVAR,MXNIJS,IJS,SIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "geomP.fh"
c
      logical status
      logical geom_zmt_put_zsym
      logical geom_zmt_opn_fil
C
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (  LBUF=9  )
      LOGICAL SOME
      LOGICAL OUT
      LOGICAL DBUG
      LOGICAL SYM
      CHARACTER*8 ERRMSG
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTOUT/IPRZMT
      COMMON/HND_ZMTOPT/IFZMAT,NATZMT          
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLNUC/NUC(MXATOM)
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      DIMENSION IJS(2,1),SIJ(1)
      DIMENSION IBUF(LBUF),BUF(LBUF)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTLIN-'/
      DATA TENM1  /1.0D-01/
      DATA TOL    /1.0D-08/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
C
      DBUG=.FALSE.
      OUT =.FALSE. 
      OUT =OUT.OR.DBUG
      OUT =OUT.OR.NPRINT.EQ.1 
      OUT =OUT.OR.IPRZMT.GE.1 
      SOME=       NPRINT.NE.-5
      SOME=SOME.OR.OUT
      IF(OUT) THEN
         WRITE(IW,*) 'IN HND_ZMTLIN'
      ENDIF
C
      NDIM=MAX(3*NAT,NZVAR)
C
C     ----- MEMORY POINTERS -----
C
      IF(.NOT.MA_PUSH_GET(MT_DBL,NDIM *NDIM ,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_ZMTLIN, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_ZMTLIN, MALLOC OF I20  FAILED',911)
C
      DO IJ=1,NDIM *NDIM 
         DBL_MB(IJ+I10-1)=ZERO
      ENDDO
      DO I=1,NZVAR        
         II=I+NDIM*(I-1)
         DBL_MB(II+I10-1)=ONE
      ENDDO
      IF(DBUG) THEN
         CALL HND_PRSQ(DBL_MB(I10),NZVAR,NZVAR,NDIM)
      ENDIF
C
C     ----- CHECK TO SEE IF WE HAVE A -LIN- MATRIX -----
C
      SYM=.FALSE.
      NIJS=0
      DO N=1,MXNIJS
         I=IJS(1,N)
         J=IJS(2,N)
         S=SIJ(  N)
         IF(      I) 100,50,10
   10    IF(NZVAR-I) 100,20,20
   20    IF(      J) 100,50,30
   30    IF(NZVAR-J) 100,40,40
   40    CONTINUE
         NIJS=NIJS+1
         DO K=1,NZVAR
            KJ=K+NDIM*(J-1)
            DBL_MB(KJ+I10-1)=ZERO
         ENDDO
         SYM=.TRUE.
   50    CONTINUE
      ENDDO
C
C     ----- SET AND PRINT -LIN- MATRIX -----
C
      IF(SYM) THEN
         MAXJ=0
         DO N=1,NIJS
            I=IJS(1,N)
            J=IJS(2,N)
            S=SIJ(  N)
            IJ=I+NDIM*(J-1)
            DBL_MB(IJ+I10-1)=S
            IF(J.GT.MAXJ) THEN
               MAXJ=J
            ENDIF
         ENDDO   
      ENDIF
      IF(DBUG) THEN
         WRITE(IW,9993)
         CALL HND_PRSQ(DBL_MB(I10),NZVAR,NZVAR,NDIM)
      ENDIF
C
      IF(.NOT.SYM) THEN
C
C     ----- NO SYMMETRY COORDINATES GIVEN -----
C
         DO I=1,NZVAR
            II=I+NDIM*(I-1)
            DBL_MB(II+I10-1)=ONE
         ENDDO
      ELSE
C
C     ----- SYMMETRY COORDINATES ARE GIVEN -----
C
         DO J=1,MAXJ 
C
C     -----           NORMALIZE            -----
C
            DUM=ZERO
            DO K=1,NZVAR
               KJ=K+NDIM*(J-1)
               DUM=DUM+DBL_MB(KJ+I10-1)**2
            ENDDO
            IF(DUM.LT.TOL) THEN
               WRITE(IW,9998) J
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            DUM=ONE/ SQRT(DUM)
            DO K=1,NZVAR
               KJ=K+NDIM*(J-1)
               DBL_MB(KJ+I10-1)=DBL_MB(KJ+I10-1)*DUM
            ENDDO
C
C     -----         ORTHONORMALIZE         -----
C
            IF(J.LT.MAXJ) THEN
               DO I=J+1,NZVAR
                  DUM=ZERO
                  DO K=1,NZVAR
                     KI=K+NDIM*(I-1)
                     KJ=K+NDIM*(J-1)
                     DUM=DUM+DBL_MB(KJ+I10-1)*DBL_MB(KI+I10-1)
                  ENDDO
                  DO K=1,NZVAR
                     KI=K+NDIM*(I-1)
                     KJ=K+NDIM*(J-1)
                     DBL_MB(KI+I10-1)=DBL_MB(KI+I10-1)
     &                               -DBL_MB(KJ+I10-1)*DUM
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         IF(DBUG) THEN
            WRITE(IW,9993)
            CALL HND_PRSQ(DBL_MB(I10),NZVAR,NZVAR,NDIM)
         ENDIF
C
C     -----    COMPLETE THE -LIN- MATRIX   -----
C      
         IF(MAXJ.LT.NZVAR) THEN
            I=0 
            DO J=MAXJ+1,NZVAR
   60          I=I+1
               IF(I.GT.NZVAR) THEN
                  WRITE(IW,9992) J
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
               DO K=1,NZVAR
                  KJ=K+(J-1)*NDIM 
                  DBL_MB(KJ+I10-1)=ZERO
               ENDDO
               IJ=I+(J-1)*NDIM
               DBL_MB(IJ+I10-1)=ONE 
C
C     -----         ORTHOGONALIZE        -----
C
               DO L=1,J-1
                  DUM=ZERO
                  DO K=1,NZVAR
                     KJ=K+(J-1)*NDIM 
                     KL=K+(L-1)*NDIM 
                     DUM=DUM+DBL_MB(KL+I10-1)*DBL_MB(KJ+I10-1)
                  ENDDO
                  DO K=1,NZVAR
                     KJ=K+(J-1)*NDIM 
                     KL=K+(L-1)*NDIM 
                     DBL_MB(KJ+I10-1)=DBL_MB(KJ+I10-1)
     &                               -DBL_MB(KL+I10-1)*DUM
                  ENDDO
               ENDDO
C
C     -----            NORMALIZE         -----
C
               DUM=ZERO
               DO K=1,NZVAR
                  KJ=K+(J-1)*NDIM 
                  DUM=DUM+DBL_MB(KJ+I10-1)*DBL_MB(KJ+I10-1)
               ENDDO
               IF(DUM.GT.TOL) THEN
                  DUM=ONE/SQRT(DUM)
                  DO K=1,NZVAR
                     KJ=K+(J-1)*NDIM 
                     DBL_MB(KJ+I10-1)=DBL_MB(KJ+I10-1)*DUM            
                  ENDDO
               ELSE
                  GO TO 60
               ENDIF
            ENDDO
         ENDIF
         IF(OUT) THEN
            WRITE(IW,9993)
            IF(DBUG) THEN
               CALL HND_PRSQ(DBL_MB(I10),NZVAR,NZVAR,NDIM)
            ENDIF
         ENDIF
C
C     ----- PRINT -LIN- MATRIX -----
C
         IF(SOME) THEN
            WRITE(IW,9997)
            DO J=1,NZVAR
               NBUF=1
               N=0
               DO I=1,NZVAR
                  IJ=I+NDIM*(J-1)
                  DUM=DBL_MB(IJ+I10-1)
                  IF( ABS(DUM).GE.TENM1) THEN
                     N=N+1
                     IBUF(N)=I
                      BUF(N)=DUM
                      IF(N.GE.LBUF) THEN
                   IF(NBUF.EQ.1) WRITE(IW,9995) J,(BUF(M),IBUF(M),M=1,N)
                   IF(NBUF.GT.1) WRITE(IW,9994)   (BUF(M),IBUF(M),M=1,N)
                         NBUF=NBUF+1
                         N=0
                      ENDIF
                  ENDIF
               ENDDO
               IF(N.GT.0) THEN
                  IF(NBUF.EQ.1) WRITE(IW,9995) J,(BUF(M),IBUF(M),M=1,N)
                  IF(NBUF.GT.1) WRITE(IW,9994)   (BUF(M),IBUF(M),M=1,N)
               ENDIF
            ENDDO
C
         ENDIF
C
      ENDIF
      GO TO 1000
C
C     ----- ERRORS -----
C
  100 CONTINUE
      WRITE(IW,9999) I,J,S
      CALL HND_HNDERR(3,ERRMSG)
C
 1000 CONTINUE
C
C     ----- SAVE SYMMETRY MATRIX -----
C           -S-       AT  -92-
C
      STATUS=GEOM_ZMT_OPN_FIL(.FALSE.,.FALSE.,.TRUE.,.TRUE.)
      IF(.NOT.STATUS) THEN
         CALL ERRQUIT('GEOM_ZMT_OPN_FIL?',911)
      ENDIF
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF(.NOT.GEOM_ZMT_PUT_ZSYM(DBL_MB(I10),NZVAR,NDIM))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_ZIND FAILED',0)
      CALL GA_SYNC
      STATUS=RTDB_PARALLEL(.TRUE.)
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_ZMTLIN, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_ZMTLIN, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
 9999 FORMAT(' INPUT ERROR IN -ZMTLIN- . I,J,S = ',2I5,F10.5)
 9998 FORMAT(' NORMALIZATION ERROR IN -ZMTLIN- FOR VECTOR J = ',I5)
 9997 FORMAT(/,10X,'SYMMETRIC INTERNAL COORDINATES',/,10X,30(1H-))
 9995 FORMAT(' COORD. = ',I3,9(1X,F5.1,2H (,I3,1H)))
 9994 FORMAT(13X,9(1X,F5.1,2H (,I3,1H)))
 9993 FORMAT(/,10X,'INT. COORD. SYM. MATRIX',/,10X,23(1H-))
 9992 FORMAT(' IN -ZMTLIN- THE ',I4,'-TH INDEPENDENT VECTOR',
     1       ' COULD NOT BE FOUND. STOP')
      END
      SUBROUTINE HND_BANDBI(GEOM,NAT,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "geomP.fh"
      INTEGER GEOM
C
C     ----- CALCULATE -B- AND -BI- MATRICES -----
C
      LOGICAL DBUG,OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION C(3,1)
C
C     ----- RETURN IF NO INTERNAL COORDINATES -----
C
      IF(NZMAT.EQ.0) RETURN
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      OUT =OUT.OR.NPRINT.EQ.1
      OUT =OUT.OR.IPRZMT.GE.1
      IF(OUT) THEN
         WRITE(IW,9999) NZMAT,NZVAR,NVAR,NAT
      ENDIF
C
      NCART=3*NAT
      NX1  =NCART
      NQ1  =NZVAR
      ND1  =MAX(NX1,NQ1)
      NZ1  =NZMAT
C
C     ----- MEMORY POINTERS -----
C
C     I10=1
C     I20=I10+NX1*NQ1
C     I30=I20+    NQ1
C     I40=I30+ABS(NZ1)  
C     LAST1=I40-1
      IF(.NOT.MA_PUSH_GET(MT_DBL,NX1*NQ1    ,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NQ1        ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I20  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ABS(NZ1)   ,'MEM I30 ',I_I30,I30))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I30  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I40 ',I_I40,I40))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I40  FAILED',911)
C
C     I11=I10
C     I21=I11+NX1*NQ1
C     I31=I21+NX1*NQ1
C     I41=I31+ND1*ND1
C     I51=I41+ND1*ND1
C     I61=I51+ND1
C     I71=I61+ND1*ND1
C     I81=I71+ND1*ND1
C     I91=I81+ND1
C     I01=I91+ND1
C     LAST2=I01-1
C     LAST=MAX0(LAST1,LAST2)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NX1*NQ1    ,'MEM I11 ',I_I11,I11))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I11  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NX1*NQ1    ,'MEM I21 ',I_I21,I21))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I21  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1*ND1    ,'MEM I31 ',I_I31,I31))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I31  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1*ND1    ,'MEM I41 ',I_I41,I41))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I41  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1        ,'MEM I51 ',I_I51,I51))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I51  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1*ND1    ,'MEM I61 ',I_I61,I61))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I61  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1*ND1    ,'MEM I71 ',I_I71,I71))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I71  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1        ,'MEM I81 ',I_I81,I81))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I81  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ND1        ,'MEM I91 ',I_I91,I91))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I91  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I01 ',I_I01,I01))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I01  FAILED',911)
C
C     ----- CALCULATE -B- MATRIX -----
C
      CALL HND_BMAT(GEOM,NZVAR,NCART,
     1              DBL_MB(I10),DBL_MB(I20),DBL_MB(I30),
     2              NZMAT,C)
      IF(OUT) THEN
         WRITE(IW,9998)
      ENDIF
C
C     ----- CALCULATE THE B INVERSE MATRIX -----
C
      CALL HND_BINVR(NVAR,NZVAR,NCART,DBL_MB(I11),DBL_MB(I21),
     1        DBL_MB(I31),DBL_MB(I41),DBL_MB(I51),DBL_MB(I61),
     2        DBL_MB(I71),DBL_MB(I81),DBL_MB(I91),ND1)
      IF(OUT) THEN
         WRITE(IW,9997)
      ENDIF
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I01))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I01  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I91))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I91  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I81))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I81  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I71))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I71  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I61))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I61  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I51))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I51  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I41))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I41  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I31))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I31  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I21))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I21  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I11))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I11  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I40))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I40  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I30))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I30  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
 9999 FORMAT(' IN -BANDBI- NZMAT,NZVAR,NVAR,NAT = ',4I5)
 9998 FORMAT(' IN -BANDBI- AFTER -BMAT- , BEFORE -BINVR- ')
 9997 FORMAT(' IN -BANDBI- AFTER -BINVR- ')
      END
      SUBROUTINE HND_BMAT(GEOM,NZVAR,NCART,B,ZMAT,IZMAT,NZMAT,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "rtdb.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      logical status
      integer geom
C
C     ----- CONSTRUCT THE B MATRIX -----
C
      PARAMETER (    NA=10 )
      LOGICAL     DBUG
      LOGICAL     OUT
      LOGICAL     GEOM_ZMT_GET_IZMAT
      LOGICAL     GEOM_ZMT_GET_NIZMAT
      LOGICAL     GEOM_ZMT_PUT_ZMAT 
      LOGICAL     GEOM_ZMT_PUT_BMAT 
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION B(NCART,*)
      DIMENSION ZMAT(*),IZMAT(*)
      DIMENSION A(3,NA)
      DIMENSION C(3,*)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- BMAT -'/
      DATA PI2    /6.28318530717958D+00/
      DATA DEGREE /360.0D+00/
      DATA BOHR   /0.52917715D+00/
      DATA ZERO   /0.0D+00/
      DATA FOUR   /4.0D+00/
C
      DBUG  =.FALSE.
      OUT   =.FALSE. 
      OUT   =OUT.OR.DBUG
      OUT   =OUT.OR.NPRINT.EQ.1
      OUT   =OUT.OR.IPRZMT.GE.1
C
      PIO2  =PI2/FOUR
      RTOD  =PI2/DEGREE
C
C     ----- READ -IZMAT- -----
C          -A-     AT -90-
C          -IZMAT- AT -96-
C
      NDUM=3*NA
C
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF( GA_NODEID().EQ.0) THEN
         IF(.NOT.GEOM_ZMT_GET_NIZMAT(GEOM,NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_NIZMAT FAILED',0)
         IF(.NOT.GEOM_ZMT_GET_IZMAT(GEOM,IZMAT,NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_IZMAT FAILED',0)
         CALL GA_BRDCST(20+MSGINT,NZMAT,MITOB(1     ),0)
         CALL GA_BRDCST(21+MSGINT,IZMAT,MITOB(NZMAT ),0)
      ELSE 
         CALL GA_BRDCST(20+MSGINT,NZMAT,MITOB(1     ),0)
         CALL GA_BRDCST(21+MSGINT,IZMAT,MITOB(NZMAT ),0)
      ENDIF
      CALL GA_SYNC()
      STATUS=RTDB_PARALLEL(.TRUE.)
C
      IF(OUT) THEN
         WRITE(IW,9997) 
         WRITE(IW,9996) (IZMAT(I),I=1,NZMAT)
      ENDIF
C
C     ----- ZERO OUT THE B MATRIX -----
C
      DO J = 1,NZVAR
         DO I = 1,NCART
            B(I,J) = ZERO
         ENDDO
      ENDDO
C
      IADD=1
      I   =0
      DO IZVAR=1,NZVAR
      I    =I+1            
      ITYPE=IZMAT(IADD)
      IF(OUT) THEN
         WRITE(IW,9999) IZVAR,IADD,ITYPE,NZMAT,NZVAR,NCART
      ENDIF
      IF(ITYPE.EQ.1) THEN
C
C     ----- BOND STRETCH -----
C
         CALL HND_BSTR(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             C,B,NCART,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 3
      ELSEIF(ITYPE.EQ.2) THEN
C
C     ----- ANGLE BEND -----
C
         CALL HND_BEND(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),
     2             C,B,NCART,RTOD,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 4
      ELSEIF(ITYPE.EQ.3) THEN
C
C     ----- TORSION -----
C
         CALL HND_TORS(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 5
      ELSEIF(ITYPE.EQ.4) THEN
C
C     ----- OUT OF PLANE ANGLE BEND -----
C
         CALL HND_OPLA(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,PIO2,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 5
      ELSEIF(ITYPE.EQ.5) THEN
C
C     ----- LINEAR ANGLE BEND -----
C
         CALL HND_LIBE(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,A,NA,BOHR)
         ZMAT(I  ) = EQVAL
         ZMAT(I+1) = EQVAL
         I=I+1
         IADD=IADD+5
      ELSEIF(ITYPE.EQ.6) THEN
C
C     ----- DIHEDRAL ANGLE BETWEEN TWO PLANES SHARING ONE ATOM -----
C
         CALL HND_DIHPLA(EQVAL,I,IZMAT(IADD+1),C,B,NCART,RTOD)
         ZMAT(I)=EQVAL
         IADD=IADD+6
C
      ELSE
         WRITE(IW,9995) IZVAR,IADD,ITYPE
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF    
C
      ENDDO
C
      IF(OUT) THEN
         CALL HND_PRTZMT(NZVAR,IZMAT,ZMAT,DUM,0)
         WRITE(IW,9998)
         IF(DBUG) THEN
            CALL HND_PRSQ(B,NZVAR,NCART,NCART)
         ENDIF
      ENDIF
C
C     ----- SAVE -B- AND -Z- ON -IDAF- -----
C          -A-                                 AT -90-
C          -B-  (pristine)                     AT -91-
C          -S-  (pristine)                     AT -92-
C          -BINV(T)= B * ( B(T)*B )**(-1) )(T) AT -93-
C                                              AT -94-
C          -ZMAT-                              AT -95-
C          -IZMAT-                             AT -96-
C
      NDUM=NCART*NZVAR
C
C     ----- STORE FOR -NWCHEM- -----
C
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF(.NOT.GEOM_ZMT_PUT_ZMAT(GEOM,ZMAT,NZVAR))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_ZMAT FAILED',0)
      IF(.NOT.GEOM_ZMT_PUT_BMAT(B   ,NCART,NZVAR,NCART))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_BMAT FAILED',0)
      CALL GA_SYNC()
      STATUS=RTDB_PARALLEL(.TRUE.)
C
      RETURN
 9999 FORMAT(' IN -BMAT- IZVAR,IADD,ITYPE,NZMAT,NZVAR,NCART = ',
     1       6I5)
 9998 FORMAT(/,' -B- MATRIX (ROWS = CARTESIAN COORDINATE.',
     1         ' COLUMNS = INTERNAL COORDINATES.)')
 9997 FORMAT(' IN -BMAT- IZMAT = ')
 9996 FORMAT(12I5)
 9995 FORMAT(' ERROR IN -BMAT- , IZVAR,IADD,ITYPE = ',3I5)
      END
      SUBROUTINE HND_BINVR(NVAR,NZVAR,NCART,B,BINV,
     1                 G,V,E,H,S,T,IA,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- CONSTRUCT THE INVERSE OF THE -B- MATRIX -----
C
C     B            NCART *NZVAR
C     BINV         NCART *NZVAR
C     G            NZVAR *NZVAR
C     V            NZVAR *NZVAR
C     E            NZVAR
C     H       MAX( NZVAR *NZVAR , NCART*NCART )
C     S            NZVAR *NZVAR
C     T            NZVAR
C     IA      MAX( NZVAR        , NCART       )
C
      LOGICAL     OUT
      LOGICAL     DBUG
      LOGICAL     NATURL
      LOGICAL     SYMTRC
      LOGICAL     INDEP
      LOGICAL     GEOM_ZMT_PUT_BINV 
      LOGICAL     GEOM_ZMT_GET_BMAT
      LOGICAL     GEOM_ZMT_PUT_ZIND
      LOGICAL     GEOM_ZMT_PUT_ZSYM
      LOGICAL     GEOM_ZMT_GET_ZSYM
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTOUT/IPRZMT
      COMMON/HND_ZMTOPT/IFZMAT,NATZMT          
      DIMENSION B(NCART,*),BINV(NCART,*)
      DIMENSION G(*),V(NDIM,*),E(*),IA(*)
      DIMENSION      S(NDIM,*),T(*)
      DIMENSION      H(NDIM,*)     
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- BINVR-'/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TM08   /1.0D-08/
C
      DBUG=.FALSE.
      DBUG=DBUG.OR.IPRZMT.GE.2
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      OUT =OUT.OR.NPRINT.EQ.1
      OUT =OUT.OR.IPRZMT.EQ.1
C
      DO I=1,NDIM           
         IA(I)=(I*(I-1))/2
      ENDDO
C
C     --- NATURL .. FINAL BASIS = NATURAL   INTERNAL        COORD. ---
C     --- SYMTRC .. FINAL BASIS = SYMMETRIC INTERNAL        COORD. ---
C     --- INDEP  .. FINAL BASIS = (NON-REDUNDANT+REDUNDANT) COORD. ---
C
      NATURL=NATZMT.EQ.0
      SYMTRC=NATZMT.EQ.1
      INDEP =NATZMT.EQ.2              
C
      IF(OUT) THEN
         WRITE(IW,9984) NZVAR,NVAR,NCART,NDIM
         WRITE(IW,9999) NATURL,SYMTRC,INDEP,NATZMT
         IF(NATURL) WRITE(IW,9998)
         IF(SYMTRC) WRITE(IW,9997)
         IF(INDEP ) WRITE(IW,9996)
      ENDIF
C
C     ----- READ -B- MATRIX FOR NATURAL INTERNAL   COORDINATES -----
C     -----  AND -S- MATRIX FOR                                -----
C     -----  TRANSFORMATION  TO SYMMETRIC INTERNAL COORDINATES -----
C
      IF(.NOT.GEOM_ZMT_GET_BMAT(B   ,NCART,NZVAR,NCART))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_BMAT FAILED',0)
      IF(OUT) THEN
         WRITE(IW,9995)
         IF(DBUG) THEN
            CALL HND_PRSQ(B,NZVAR,NCART,NCART)
         ENDIF
      ENDIF
C
C     -----   TRANSFORM -B- TO BASIS OF    -----
C     ----- SYMMETRIC INTERNAL COORDINATES -----
C
      IF(SYMTRC.OR.INDEP) THEN
         IF(.NOT.GEOM_ZMT_GET_ZSYM(S,NZVAR,NDIM))            
     1      CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_ZSYM FAILED',0)
         DO I=1,NCART
            DO J=1,NZVAR
               DUM=ZERO
               DO K=1,NZVAR
                  DUM=DUM+B(I,K)*S(K,J)
               ENDDO
               T(J)=DUM
            ENDDO
            DO J=1,NZVAR
               B(I,J)=T(J)
            ENDDO
         ENDDO
         IF(OUT) THEN
            WRITE(IW,9994)
            WRITE(IW,9993)
            IF(DBUG) THEN
               CALL HND_PRSQ(B,NZVAR,NCART,NCART)
            ENDIF
         ENDIF
      ENDIF
C
C          -A-                                 AT -90-
C          -B-  (maybe ... B tilde = B * S )   AT -91-
C          -S-  (pristine)                     AT -92-
C          -BINV(T)= B * ( B(T)*B )**(-1) )(T) AT -93-
C                                              AT -94-
C          -ZMAT-                              AT -95-
C          -IZMAT-                             AT -96-
C         NOTE THAT AT THIS POINT WE MAY BE WORKING WITH 
C           B(tilde) WHICH IS THE -S- TRANSFORMED OF B
C                      B(tilde) = B * S
C
C     -----      CONSTRUCT THE -G- MATRIX       -----
C                      G = B(T) * B
C             ACTUALLY WORK WITH ( - B(T) * B )
C           ( EASY WAY TO RE-ORDER THE VECTORS )
C
      DO I=1,NZVAR
         DO J=1,I
            IJ=IA(I)+J
            DUM=ZERO
            DO K=1,NCART
               DUM=DUM+B(K,I)*B(K,J)
            ENDDO
            G(IJ)=-DUM
         ENDDO
      ENDDO
      IF(DBUG) THEN
         WRITE(IW,9992)
         CALL HND_PRTR(G,NZVAR)
      ENDIF
C
C     ----- GET THE ( NON-REDUNDANT + REDUNDANT ) COORDINATES -----
C
      CALL HND_DIAGIV(G,V,E,IA,NZVAR,NZVAR,NDIM) 
C
C     -----  WE WORKED WITH ( - B(T) * B ) TO GET THE CORRECT ORDER  -----
C           FOR THE VECTORS. WE NEED TO CHANGE THE SIGN OF THE E's
C
      DO I=1,NZVAR
         E(I)=-E(I)
      ENDDO
C
C     -----    IF WE ARE TO WORK IN THE SPACE OF     -----
C     -----  (NON-REDUNDANT+REDUNDANT) COORDINATES   -----
C     -----   WE MUST MAINTAIN A PHASE CONSISTENCY   -----
C     -----   ... IF NOT ... WE ARE IN TROUBLES ...  -----
C     -----  SOME (PATHOLOGICAL) CASES WILL NOT WORK -----
C
      DO J=1,NZVAR 
         IMAX=0
         VMAX=ZERO
         DO I=1,NZVAR
            IF(ABS(V(I,J)).GT.VMAX) THEN
               VMAX=ABS(V(I,J))
               IMAX=I
            ENDIF
         ENDDO
         IF(V(IMAX,J).LT.ZERO) THEN
            DO I=1,NZVAR
               V(I,J)=-V(I,J)
            ENDDO
         ENDIF
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,9991)
         WRITE(IW,9990)
         CALL HND_PREV(V,E,NZVAR,NZVAR,NDIM) 
      ENDIF
C
C     ----- CHECK FOR SINGULARITY -----
C
      NINDEP=0
      NREDUN=0
      DO I=1,NZVAR
         IF( ABS(E(I)).LE.TM08) THEN
            E(I)  =ZERO
            NREDUN=NREDUN+1
         ELSE
            E(I)  =ONE/E(I)
            NINDEP=NINDEP+1
         ENDIF
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9989) NINDEP,NREDUN,NZVAR,NVAR
      ENDIF
C
C     IF NREDUN.EQ.0 , WE CAN WORK WITH NATURAL COORDINATES -----
C
      IF(NREDUN.EQ.0) THEN
         NATZMT=0
         NATURL=NATZMT.EQ.0
         SYMTRC=NATZMT.EQ.1
         INDEP =NATZMT.EQ.2              
         IF(OUT) THEN
            WRITE(IW,9979)
         ENDIF
      ENDIF
C
      IF(NATURL.OR.SYMTRC) THEN
         DO J=1,NZVAR
            DO I=1,NZVAR
               S(I,J)=ZERO
            ENDDO
            S(J,J)=ONE
         ENDDO
         IF(.NOT.GEOM_ZMT_PUT_ZIND(S,NZVAR,NDIM))            
     1      CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_ZIND FAILED',0)
      ENDIF
C
C     ----- IF WORKING BASIS IS (NON-REDUNDANT+REDUNDANT) -----
C     -----         TRANSFORM -B- TO THAT BASIS           -----
C
      IF(INDEP) THEN
         IF(.NOT.GEOM_ZMT_GET_ZSYM(S,NZVAR,NDIM))            
     1      CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_ZSYM FAILED',0)
         DO I=1,NZVAR
            DO J=1,NZVAR
               DUM=ZERO
               DO K=1,NZVAR
                  DUM=DUM+S(I,K)*V(K,J)
               ENDDO
               T(J)=DUM
            ENDDO
            DO J=1,NZVAR
               S(I,J)=T(J)
            ENDDO
         ENDDO
         IF(.NOT.GEOM_ZMT_PUT_ZSYM(S,NZVAR,NDIM))            
     1      CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_ZSYM FAILED',0)
         IF(OUT) THEN
            WRITE(IW,9980)
            IF(DBUG) THEN
               CALL HND_PRSQ(S,NZVAR,NZVAR,NDIM)
            ENDIF
         ENDIF
C
         DO J=NINDEP+1,NZVAR
            DO I=1,NZVAR
               V(I,J)=ZERO
            ENDDO
         ENDDO
C
         DO I=1,NCART
            DO J=1,NZVAR
               DUM=ZERO
               DO K=1,NZVAR
                  DUM=DUM+B(I,K)*V(K,J)
               ENDDO
               T(J)=DUM
            ENDDO
            DO J=1,NZVAR
               B(I,J)=T(J)
            ENDDO
         ENDDO
         IF(OUT) THEN
            WRITE(IW,9985)
            IF(DBUG) THEN
               CALL HND_PRSQ(B,NZVAR,NCART,NCART)
            ENDIF
         ENDIF
C
         DO J=1,NINDEP
            DO I=1,NZVAR
               V(I,J)=ZERO
            ENDDO
            V(J,J)=ONE 
         ENDDO
      ENDIF
C
C     ----- CALCULATE  K * K(T) PROJECTOR -----
C
      DO I=1,NZVAR
         DO J=1,NZVAR
            DUM=ZERO
            DO K=1,NINDEP
               DUM=DUM+V(I,K)*V(J,K)
            ENDDO
            H(I,J)=DUM
         ENDDO
      ENDDO
      NDUM=NDIM*NDIM  
      IF(OUT) THEN
         IF(DBUG) THEN
            WRITE(IW,9981)
            CALL HND_PRSQ(H,NZVAR,NZVAR,NDIM) 
         ENDIF
      ENDIF
C
C     ----- SET  ... NVAR = NZVAR ... IF WE WORK -----
C           WITH THE SYMMETRIC INTERNAL SPACE OR
C           WITH THE NATURAL   INTERNAL SPACE
C
      IF(NATURL.OR.SYMTRC) THEN
         NVAR=NZVAR
      ELSE
         IF(OUT) THEN
            IF(NINDEP.NE.NVAR) THEN
               WRITE(IW,9988) NZVAR,NINDEP,NREDUN,NVAR
            ENDIF
            IF(NREDUN.NE.NZVAR-NVAR) THEN
               WRITE(IW,9988) NZVAR,NINDEP,NREDUN,NVAR
            ENDIF
         ENDIF
      ENDIF
C
C     ----- GENERALIZED INVERSE OF -G- IN THE BASIS OF -----
C     ----- NATURAL/SYMMETRIC INTERNAL COORDINATES  OR -----
C     ----- (NON-REDUNDANT+REDUNDANT)  COORDINATES     -----
C
      DO I=1,NZVAR
         DO J=1,NZVAR
            DUM=ZERO
            DO K=1,NZVAR
               DUM=DUM+V(I,K)*E(K)*V(J,K)
            ENDDO
            H(I,J)=DUM
         ENDDO
      ENDDO
      IF(DBUG) THEN
         WRITE(IW,9987)
         CALL HND_PRSQ(H,NZVAR,NZVAR,NDIM) 
      ENDIF
C
C     -----   BINV(T) = B * ( ( B(T)*B ) **(-1) )(T)  -----
C
      DO I=1,NCART
         DO J=1,NZVAR
            DUM=ZERO
            DO K=1,NZVAR
               DUM=DUM+B(I,K)*H(K,J)
            ENDDO
            BINV(I,J)=DUM
         ENDDO
      ENDDO
      IF(OUT) THEN
         WRITE(IW,9982)
         WRITE(IW,9986)
         IF(DBUG) THEN
            CALL HND_PRSQ(BINV,NZVAR,NCART,NCART)
         ENDIF
      ENDIF
      NDUM=NCART*NVAR
C
C     ----- STORE FOR -NWCHEM- -----
C
      IF(.NOT.GEOM_ZMT_PUT_BINV(BINV,NCART,NZVAR,NCART))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_PUT_BINV FAILED',0)
C
C     -----        DEBUG CHECK OF INVERSE             -----
C     -----         H = Bi * B = K * K(T)             -----
C     -----         WHERE K IS THE MATRIX             -----
C     -----    OF (NON-REDUNDANT) COORDINATES         -----
C     -----          H IS NOT EQUAL TO   1            -----
C     ----- UNLESS THERE ARE NO REDUNDANT COORDINATES -----
C
      IF(DBUG) THEN
         NDUM= NCART*NZVAR
         DO I=1,NZVAR
            DO J=1,NZVAR
               DUM=ZERO
               DO K=1,NCART
                  DUM=DUM+BINV(K,I)*B(K,J)
               ENDDO
               H(I,J)=DUM
            ENDDO
         ENDDO
         WRITE(IW,9983)
         CALL HND_PRSQ(H,NZVAR,NZVAR,NDIM) 
      ENDIF
C
C     ----- NEXT TIME AROUND WE WILL WORK WITH -----
C           WHAT ARE NOW SYM.INTERN.COORD.
C
      NATZMT=1
C
C          -A-                                 AT -90-
C          -B-  (B tilde)                      AT -91-
C          -S-  (pristine)                     AT -92-
C          -BINV(T)= B * ( B(T)*B )**(-1) )(T) AT -93-
C          -P = K * K(T) PROJECTOR             AT -94-
C          -ZMAT-                              AT -95-
C          -IZMAT-                             AT -96-
C     -NZMAT-   = N.OF. -IZMAT- ELEMENTS GIVEN IN -$ZMAT-
C     -NZVAR-   = N.OF. INTERNAL COORDINATES DEFINED IN -IZMAT-
C     -NVAR-    = (3N-6), OR (3N-5), OR -NZMOD- , .LE. NZVAR
C
      RETURN
 9999 FORMAT(/,' IN -BINVR- -NATURL-, -SYMTRC- , -INDEP- , -NATZMT- = ',
     1       3L5,I4)
 9998 FORMAT(/,' FINAL BASIS WILL BE =    NATURAL   INTERNAL    ')
 9997 FORMAT(/,' FINAL BASIS WILL BE =    SYMMETRIC INTERNAL    ')
 9996 FORMAT(/,' FINAL BASIS WILL BE = (NON-REDUNDANT+REDUNDANT)')
 9995 FORMAT(/,' -B - IN NATURAL   INTERNAL        COORDINATES ')
 9994 FORMAT(/,' -B - IN NATURAL   INTERNAL        COORDINATES ',
     1                                             'TRANSFORMED',
     2          '     TO SYMMETRIC INTERNAL        COORDINATES ')
 9993 FORMAT(/,' -B - IN SYMMETRIC INTERNAL        COORDINATES ')
 9992 FORMAT(/,' -G - IN SYMMETRIC INTERNAL        COORDINATES ')
 9991 FORMAT(/,'         (NON-REDUNDANT+REDUNDANT) COORDINATES = ',
     1            'EIGENVECTORS OF -G- ')
 9990 FORMAT(/,' -V - OF (NON-REDUNDANT+REDUNDANT) COORDINATES ')
 9989 FORMAT(/,' NUM. OF LIN.INDEP. COORD. (NINDEP) = ',I5,
     1       /,' NUM. OF REDUNDANT  COORD. (NREDUN) = ',I5,
     2       /,'                            -NZVAR- = ',i5,
     3       /,'                            -NVAR-  = ',I5)
 9988 FORMAT(/,' -NINDEP- SHOULD BE EQUAL TO -NVAR-, BUT IT IS NOT.',
     1         ' STOP. NZVAR,NINDEP,NREDUN,NVAR = ',4I5)
 9987 FORMAT(/,' - G**(-1) - ') 
 9986 FORMAT(/,' -Bi- IN ..... CHOSEN .....        COORDINATES ')
 9985 FORMAT(/,' -B - IN (NON-REDUNDANT+REDUNDANT) COORDINATES ')
 9984 FORMAT(/,' IN -BINVR- NZVAR,NVAR,NCART,NDIM = ',4I5)
 9983 FORMAT(/,' CHECK OF ( Bi * B = K * K(T) )')
 9982 FORMAT(/,' -Bi(T) = B * ( G**(-1) )(T)')
 9981 FORMAT(/,' K * K(T) PROJECTOR ')         
 9980 FORMAT(/,' -ZIND- FROM NATURAL TO (INDEP+DEP) COORDINATES')
 9979 FORMAT(/,' COORDINATE BASIS CHANGED BACK TO -NATURAL-')
      END
      SUBROUTINE HND_BANDZ(GEOM,NAT,C,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "geomP.fh"
C
C     ----- CALCULATE -B- AND -BI- MATRICES -----
C
      LOGICAL OUT
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION C(3,*),Z(*)
C
C     ----- RETURN IF NO INTERNAL COORDINATES -----
C
      IF(NZMAT.EQ.0) RETURN
C
      OUT=.FALSE.
      OUT=OUT.OR.NPRINT.EQ.1
      OUT=OUT.OR.IPRZMT.GE.1
      IF(OUT) THEN
         WRITE(IW,9999) NZMAT,NZVAR,NVAR,NAT
      ENDIF
C
      NCART=3*NAT
      NX1  =NCART
      NQ1  =NZVAR
      ND1  =MAX(NX1,NQ1)
      NZ1  =NZMAT
C
C     ----- MEMORY POINTERS -----
C
C     I10=1
C     I20=I10+NX1*NQ1
C     I30=I20+    NQ1
C     I40=I30+ABS(NZ1)  
C     LAST1=I40-1
C
      IF(.NOT.MA_PUSH_GET(MT_DBL,NX1*NQ1    ,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,NQ1        ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I20  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ABS(NZ1)   ,'MEM I30 ',I_I30,I30))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I30  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I40 ',I_I40,I40))
     & CALL ERRQUIT('HND_BANDBI, MALLOC OF I40  FAILED',911)
C
C     ----- CALCULATE -B- MATRIX -----
C
      CALL HND_BZ(GEOM,NZVAR,NCART,
     1            DBL_MB(I10),DBL_MB(I20),DBL_MB(I30),
     2            NZMAT,C)
      DO I=1,NZVAR
         Z(I)=DBL_MB(I+I20-1)
      ENDDO
C
      IF(OUT) THEN
         WRITE(IW,9998)
      ENDIF
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I40))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I40  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I30))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I30  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_BANDBI, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
 9999 FORMAT(' IN -BANDZ- NZMAT,NZVAR,NVAR,NAT = ',4I5)
 9998 FORMAT(' IN -BANDZ- AFTER -BZ- ')
      END
      SUBROUTINE HND_BZ(GEOM,NZVAR,NCART,B,ZMAT,IZMAT,NZMAT,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer geom
C
C     ----- CONSTRUCT THE B MATRIX -----
C
      PARAMETER (    NA=10 )
      LOGICAL     OUT
      LOGICAL     GEOM_ZMT_GET_IZMAT
      LOGICAL     GEOM_ZMT_GET_NIZMAT
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_OUTPUT/NPRINT
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION B(NCART,*)
      DIMENSION ZMAT(*),IZMAT(*)
      DIMENSION A(3,NA)
      DIMENSION C(3,*)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- BMAT -'/
      DATA PI2    /6.28318530717958D+00/
      DATA DEGREE /360.0D+00/
      DATA BOHR   /0.52917715D+00/
      DATA ZERO   /0.0D+00/
      DATA FOUR   /4.0D+00/
C
      OUT   =.FALSE. 
      OUT   =OUT.OR.NPRINT.EQ.1
      OUT   =OUT.OR.IPRZMT.GE.1
C
      PIO2  =PI2/FOUR
      RTOD  =PI2/DEGREE
C
C     ----- READ -IZMAT- -----
C          -A-     AT -90-
C          -IZMAT- AT -96-
C
      NDUM=3*NA
      IF(.NOT.GEOM_ZMT_GET_NIZMAT(GEOM,NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_NIZMAT FAILED',0)
      IF(.NOT.GEOM_ZMT_GET_IZMAT(GEOM,IZMAT,NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_IZMAT FAILED',0)
      IF(OUT) THEN
         WRITE(IW,9997) 
         WRITE(IW,9996) (IZMAT(I),I=1,NZMAT)
      ENDIF
C
C     ----- ZERO OUT THE B MATRIX -----
C
      DO J = 1,NZVAR
         DO I = 1,NCART
            B(I,J) = ZERO
         ENDDO
      ENDDO
C
      IADD=1
      I   =0
      DO IZVAR=1,NZVAR
      I    =I+1            
      ITYPE=IZMAT(IADD)
      IF(OUT) THEN
         WRITE(IW,9999) IZVAR,IADD,ITYPE,NZMAT,NZVAR,NCART
      ENDIF
      IF(ITYPE.EQ.1) THEN
C
C     ----- BOND STRETCH -----
C
         CALL HND_BSTR(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             C,B,NCART,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 3
      ELSEIF(ITYPE.EQ.2) THEN
C
C     ----- ANGLE BEND -----
C
         CALL HND_BEND(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),
     2             C,B,NCART,RTOD,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 4
      ELSEIF(ITYPE.EQ.3) THEN
C
C     ----- TORSION -----
C
         CALL HND_TORS(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 5
      ELSEIF(ITYPE.EQ.4) THEN
C
C     ----- OUT OF PLANE ANGLE BEND -----
C
         CALL HND_OPLA(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,PIO2,BOHR)
         ZMAT(I) = EQVAL
         IADD = IADD + 5
      ELSEIF(ITYPE.EQ.5) THEN
C
C     ----- LINEAR ANGLE BEND -----
C
         CALL HND_LIBE(EQVAL,I,IZMAT(IADD+1),IZMAT(IADD+2),
     1             IZMAT(IADD+3),IZMAT(IADD+4),
     2             C,B,NCART,RTOD,A,NA,BOHR)
         ZMAT(I  ) = EQVAL
         ZMAT(I+1) = EQVAL
         I=I+1
         IADD=IADD+5
      ELSEIF(ITYPE.EQ.6) THEN
C
C     ----- DIHEDRAL ANGLE BETWEEN TWO PLANES SHARING ONE ATOM -----
C
         CALL HND_DIHPLA(EQVAL,I,IZMAT(IADD+1),C,B,NCART,RTOD)
         ZMAT(I)=EQVAL
         IADD=IADD+6
C
      ELSE
         WRITE(IW,9995) IZVAR,IADD,ITYPE
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF    
C
      ENDDO
C
      IF(OUT) THEN
         CALL HND_PRTZMT(NZVAR,IZMAT,ZMAT,DUM,0)
         WRITE(IW,9998)
         CALL HND_PRSQ(B,NZVAR,NCART,NCART)
      ENDIF
C
C     ----- SAVE -B- AND -Z- ON -IDAF- -----
C          -A-                                 AT -90-
C          -B-  (pristine)                     AT -91-
C          -S-  (pristine)                     AT -92-
C          -BINV(T)= B * ( B(T)*B )**(-1) )(T) AT -93-
C                                              AT -94-
C          -ZMAT-                              AT -95-
C          -IZMAT-                             AT -96-
C
      NDUM=NCART*NZVAR
C
      RETURN
 9999 FORMAT(' IN -BZ- IZVAR,IADD,ITYPE,NZMAT,NZVAR,NCART = ',
     1       6I5)
 9998 FORMAT(/,' -B- MATRIX (ROWS = CARTESIAN COORDINATE.',
     1         ' COLUMNS = INTERNAL COORDINATES.)')
 9997 FORMAT(' IN -BZ- IZMAT = ')
 9996 FORMAT(12I5)
 9995 FORMAT(' ERROR IN -BZ- , IZVAR,IADD,ITYPE = ',3I5)
      END
      SUBROUTINE HND_ZMTPRT(GEOM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "geomP.fh"
C
C     ----- PRINT THE CURRENT Z MATRIX -----
C
      LOGICAL STATUS
      INTEGER GEOM
      LOGICAL GEOM_ZMT_GET_IZMAT
      LOGICAL GEOM_ZMT_GET_NIZMAT
      LOGICAL GEOM_ZMT_GET_ZMAT
      LOGICAL DBUG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOUT/IPRZMT
c
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999) NZMAT,NZVAR,NVAR
      ENDIF
C
      IF(NZMAT.EQ.0) RETURN
C
C     ----- MEMORY POINTERS -----
C
C     I10 = 1
C     I20 = I10 + NZVAR
C     I30 = I20 + ABS(NZMAT)
C     LAST = I30 - 1
      IF(.NOT.MA_PUSH_GET(MT_DBL,NZVAR      ,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_ZMTPRT, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ABS(NZMAT) ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_ZMTPRT, MALLOC OF I20  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I30 ',I_I30,I30))
     & CALL ERRQUIT('HND_ZMTPRT, MALLOC OF I30  FAILED',911)
C
C     ----- GET IZMAT AND ZMAT FROM DIRECT ACCESS -----
C
      STATUS=RTDB_PARALLEL(.FALSE.)
      IF( GA_NODEID().EQ.0) THEN 
         IF(.NOT.GEOM_ZMT_GET_NIZMAT(GEOM,NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_NIZMAT FAILED',0)
         IF(.NOT.GEOM_ZMT_GET_IZMAT(GEOM,DBL_MB(I20),NZMAT))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_IZMAT FAILED',0)
         IF(.NOT.GEOM_ZMT_GET_ZMAT(GEOM,DBL_MB(I10),NZVAR))
     1   CALL ERRQUIT('GEOM_INPUT: GEOM_ZMT_GET_ZMAT FAILED',0)
         CALL GA_BRDCST(22+MSGDBL,DBL_MB(I20),MDTOB(NZMAT),0)
         CALL GA_BRDCST(23+MSGDBL,DBL_MB(I10),MDTOB(NZVAR),0)
      ELSE
         CALL GA_BRDCST(22+MSGDBL,DBL_MB(I20),MDTOB(NZMAT),0)
         CALL GA_BRDCST(23+MSGDBL,DBL_MB(I10),MDTOB(NZVAR),0)
      ENDIF
      CALL GA_SYNC()
      STATUS=RTDB_PARALLEL(.TRUE.)
C
C     ----- PRINT THE Z MATRIX -----
C
      IF( GA_NODEID().EQ.0.OR.DBUG) THEN 
         CALL HND_PRTZMT(NZVAR,DBL_MB(I20),DBL_MB(I10),DBL_MB(I30),0)
      ENDIF
      CALL GA_SYNC()
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I30))
     & CALL ERRQUIT('HND_ZMTPRT, MA_POP_STACK OF I30  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_ZMTPRT, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_ZMTPRT, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
 9999 FORMAT(' IN -ZMTPRT- NZMAT, NZVAR, NVAR = ',3i5)
      END
      SUBROUTINE HND_ZMTPRG(GRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "geomP.fh"
C
C     ----- PRINT THE CURRENT Z MATRIX AND ITs GRADIENT -----
C
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTOUT/IPRZMT
      DIMENSION GRAD(1)
      DATA ZERO /0.0D+00/
C
      IF(NZMAT.EQ.0) RETURN
C
C     ----- MEMORY POINTERS -----
C
C     I10 = 1
C     I20 = I10 + NZVAR
C     I30 = I20 + ABS(NZMAT)
C     LAST = I30 - 1
      IF(.NOT.MA_PUSH_GET(MT_DBL,NZVAR      ,'MEM I10 ',I_I10,I10))
     & CALL ERRQUIT('HND_ZMTPRG, MALLOC OF I10  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,ABS(NZMAT) ,'MEM I20 ',I_I20,I20))
     & CALL ERRQUIT('HND_ZMTPRG, MALLOC OF I20  FAILED',911)
      IF(.NOT.MA_PUSH_GET(MT_DBL,1          ,'MEM I30 ',I_I30,I30))
     & CALL ERRQUIT('HND_ZMTPRG, MALLOC OF I30  FAILED',911)
C
C     ----- GET IZMAT AND ZMAT FROM DIRECT ACCESS -----
C
C
C     ----- PRINT THE Z MATRIX -----
C
      IF(NZVAR.GT.NVAR) THEN
         DO I=NVAR+1,NZVAR
            GRAD(I)=ZERO
         ENDDO
      ENDIF
      CALL HND_PRTZMT(NZVAR,DBL_MB(I20),DBL_MB(I10),GRAD,1)
C
C     ----- RELEASE MEMORY POINTERS -----
C
      IF(.NOT.MA_POP_STACK(I_I30))
     & CALL ERRQUIT('HND_ZMTPRG, MA_POP_STACK OF I30  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I20))
     & CALL ERRQUIT('HND_ZMTPRG, MA_POP_STACK OF I20  FAILED',911)
      IF(.NOT.MA_POP_STACK(I_I10))
     & CALL ERRQUIT('HND_ZMTPRG, MA_POP_STACK OF I10  FAILED',911)
C
      RETURN
      END
      SUBROUTINE HND_PRTZMT(NZVAR,IZMAT,ZMAT,ZGRAD,IPRGRD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- PRINT THE Z-MATRIX AND ITs GRADIENT -----
C
C     ----- IPRGRD = 0  PRINT JUST THE Z MATRIX
C                  = 1  PRINT THE Z MATRIX AND ITs GRADIENT
C
      CHARACTER*8 RTYPE
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_ZMTOPT/IFZMAT,NATZMT          
      DIMENSION IZMAT(1),ZMAT(1),ZGRAD(1)
      DIMENSION RTYPE(6)
      DATA RTYPE/'STRETCH ','BEND    ','TORSION ','PLA.BEND',
     1           'LIN.BEND','DIH.PLA '/
      IF (IPRGRD.EQ.0) WRITE (IW,9008)
      IF (IPRGRD.EQ.1) WRITE (IW,9028)
      IADD = 1
      I = 0
   80 I = I + 1
      IF(I.GT.NZVAR) GO TO 220
      ITYPE = IZMAT(IADD)
      GO TO (100,120,140,160,180,190),ITYPE
C
C     ----- BOND STRETCH -----
C
  100 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9048) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9048) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),ZMAT(I),ZGRAD(I)
      IADD = IADD + 3
      GO TO 200
C
C     ----- ANGLE BEND -----
C
  120 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9068) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9068) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),ZMAT(I),ZGRAD(I)
      IADD = IADD + 4
      GO TO 200
C
C     ----- TORSION -----
C
  140 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I),ZGRAD(I)
      IADD = IADD + 5
      GO TO 200
C
C     ----- OUT OF PLANE ANGLE BEND -----
C
  160 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I),ZGRAD(I)
      IADD = IADD + 5
      GO TO 200
C
C     ----- LINEAR ANGLE BEND -----
C
  180 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9088) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),ZMAT(I),ZGRAD(I)
      I = I + 1
      IADD = IADD + 5
      GO TO 200
C
C     ----- DIHEDRAL ANGLE BETWEEN PLANES WITH ONE COMMON ATOM -----
C
  190 CONTINUE
      IF (IPRGRD.EQ.0) WRITE (IW,9108) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),IZMAT(IADD+5),
     2    ZMAT(I)
      IF (IPRGRD.EQ.1) WRITE (IW,9108) I,RTYPE(ITYPE),IZMAT(IADD+1),
     1    IZMAT(IADD+2),IZMAT(IADD+3),IZMAT(IADD+4),IZMAT(IADD+5),
     2    ZMAT(I),ZGRAD(I)
      IADD = IADD + 6
C
  200 CONTINUE
      GO TO 80
  220 CONTINUE
C
      IF(IPRGRD.EQ.1) THEN
         IF(NATZMT.EQ.0) WRITE(IW,9128)
         IF(NATZMT.EQ.1) WRITE(IW,9148)
         IF(NATZMT.EQ.2) WRITE(IW,9168)
      ENDIF
      RETURN
 9008 FORMAT(/,22X,20(1H-),/,22X,20Hinternal coordinates,/,22X,20(1H-),
     1       /,2X,3Hno.,2X,8H  type  ,
     2         2X,6H -i-  ,6H -j-  ,6H -k-  ,6H -l-  ,6H -m-  ,
     3         5X,7Hvalue  ,
     3       /,2X,57(1H-))
 9028 FORMAT(/,22X,20(1H-),/,22X,20Hinternal coordinates,/,22X,20(1H-),
     1       /,2X,3Hno.,2X,8H  type  ,
     2         2X,6H -i-  ,6H -j-  ,6H -k-  ,6H -l-  ,6H -m-  ,
     3         5X,7Hvalue  ,4X,8Hgradient,
     4       /,2X,69(1H-))
 9048 FORMAT(1X,I4,2X,A8,2X,2(I3,3X),18X,F12.5,F12.7)
 9068 FORMAT(1X,I4,2X,A8,2X,3(I3,3X),12X,F12.5,F12.7)
 9088 FORMAT(1X,I4,2X,A8,2X,4(I3,3X), 6X,F12.5,F12.7)
 9108 FORMAT(1X,I4,2X,A8,2X,5(I3,3X),    F12.5,F12.7)
 9128 FORMAT(7X,'                  gradients are for natural ',
     1          'internal coordinates')
 9148 FORMAT(7X,'                gradients are for symmetric ',                 
     1          'internal coordinates')
 9168 FORMAT(7X,'gradients are for (non-redundant+redundant) ',
     1          'internal coordinates')
      END
      SUBROUTINE HND_BSTR(EQVAL,NOINT,I,J,C,B,NCART,BOHR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR A
C          BOND STRETCH AS DEFINED BY WILSON (SEE WDC P.55) -----
C
      DIMENSION C(3,1),B(NCART,1)
      DIMENSION RIJ(3)
      DATA ZERO /0.0D+00/
      DIJSQ = ZERO
      DO 100 M = 1,3
      RIJ(M) = C(M,J)-C(M,I)
  100 DIJSQ = DIJSQ+RIJ(M)*RIJ(M)
      DO 120 M = 1,3
      NOCOL1 = 3*(I-1)+M
      NOCOL2 = 3*(J-1)+M
      B(NOCOL1,NOINT) = -RIJ(M)/ SQRT(DIJSQ)
  120 B(NOCOL2,NOINT) = -B(NOCOL1,NOINT)
      EQVAL = SQRT(DIJSQ)*BOHR
      RETURN
      END
      SUBROUTINE HND_BEND(EQVAL,NOINT,I,J,K,C,B,NCART,RTOD,BOHR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS OF A
C          VALENCE ANGLE BENDING COORDINATE AS DEFINED BY WILSON.
C          SEE WDC P. 56 -----
C
C     -----I AND K ARE THE NUMBERS OF THE END ATOMS.  J IS THE
C          NUMBER OF THE CENTRAL ATOM -----
C
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,1),B(NCART,1)
      DIMENSION RJI(3),RJK(3),EJI(3),EJK(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','-BEND  -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
C
      DJISQ = ZERO
      DJKSQ = ZERO
      DO 100 M = 1,3
      RJI(M) = C(M,I)-C(M,J)
      RJK(M) = C(M,K)-C(M,J)
      DJISQ = DJISQ+RJI(M)*RJI(M)
  100 DJKSQ = DJKSQ+RJK(M)*RJK(M)
      DJI =  SQRT(DJISQ)
      DJK =  SQRT(DJKSQ)
      DOTJ = ZERO
      DO 120 M = 1,3
      EJI(M) = RJI(M)/DJI
      EJK(M) = RJK(M)/DJK
  120 DOTJ = DOTJ+EJI(M)*EJK(M)
      IF (ONE- ABS(DOTJ)) 200,200,140
  140 SINJ = SQRT(ONE-DOTJ*DOTJ)
      DJI=DJI*BOHR
      DJK=DJK*BOHR
      DO 160 M = 1,3
      NOCOL1 = 3*(I-1)+M
      B(NOCOL1,NOINT) = ((DOTJ*EJI(M)-EJK(M)))/(DJI*SINJ)
      NOCOL2 = 3*(K-1)+M
      B(NOCOL2,NOINT) = ((DOTJ*EJK(M)-EJI(M)))/(DJK*SINJ)
      NOCOL3 = 3*(J-1)+M
  160 B(NOCOL3,NOINT) = -B(NOCOL1,NOINT)-B(NOCOL2,NOINT)
      DOT = ZERO
      DO 180 M = 1,3
  180 DOT = DOT+EJI(M)*EJK(M)
      EQVAL = ACOS(DOT)/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
      GO TO 220
  200 CONTINUE
      WRITE(IW,9008)
      CALL HND_HNDERR(3,ERRMSG)
  220 CONTINUE
      RETURN
 9008 FORMAT(' ROUNDOFF ERROR IN -BEND- STOP ')
      END
      SUBROUTINE HND_TORS(EQVAL,NOINT,I,J,K,L,C,B,NCART,RTOD,BOHR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR THE
C           TORSION AS DEFINED BY WILSON.  SEE WDC P60.
C           NOTE THAT THE SIGN OF THE SECOND TERM ON THE RIGHT HAND
C           SIDE IN EQUATION (22) IS INCORRECT AS PRINTED ON P60 OF
C           WDC.
C
C     ----- THE ORDER OF ATOMS IS I--J--K--L. IF THE OBSERVER LOOKS
C           ALONG THE VECTOR J-->K WITH J NEARER THE OBSERVER, THEN
C           THE CLOCKWISE ROTATION OF J-->I WHICH SUPERPOSES J--> I
C           WITH K-->L IS A NEGATIVE TORSION ANGLE.
C
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,1),B(NCART,1)
      DIMENSION RIJ(3),RJK(3),RKL(3)
      DIMENSION EIJ(3),EJK(3),EKL(3)
      DIMENSION CR1(3),CR2(3),D(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','-TORS  -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOLRD     /1.0001D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
C
C     ----- CALCULATING UNIT VECTORS ALONG I-->J, J-->K AND K-->L.
C
      DIJSQ = ZERO
      DJKSQ = ZERO
      DKLSQ = ZERO
      DO 120 M = 1,3
         RIJ(M) = C(M,J)-C(M,I)
         DIJSQ  = DIJSQ+RIJ(M)*RIJ(M)
         RJK(M) = C(M,K)-C(M,J)
         DJKSQ  = DJKSQ+RJK(M)*RJK(M)
         RKL(M) = C(M,L)-C(M,K)
         DKLSQ  = DKLSQ+RKL(M)*RKL(M)
  120 CONTINUE
      DIJ =  SQRT(DIJSQ)
      DJK =  SQRT(DJKSQ)
      DKL =  SQRT(DKLSQ)
      DO 180 M = 1,3
         EIJ(M) = RIJ(M)/DIJ
         EJK(M) = RJK(M)/DJK
         EKL(M) = RKL(M)/DKL
  180 CONTINUE
C
C                   -->   --->   --->     -->   --->   --->
C     ----- FORMING CR1 = I--J X J--K AND CR2 = J--K X K--L
C
      CR1(1) = EIJ(2)*EJK(3)-EIJ(3)*EJK(2)
      CR1(2) = EIJ(3)*EJK(1)-EIJ(1)*EJK(3)
      CR1(3) = EIJ(1)*EJK(2)-EIJ(2)*EJK(1)
      CR2(1) = EJK(2)*EKL(3)-EJK(3)*EKL(2)
      CR2(2) = EJK(3)*EKL(1)-EJK(1)*EKL(3)
      CR2(3) = EJK(1)*EKL(2)-EJK(2)*EKL(1)
C
C     ----- CALCULATING COS AND SIN OF ANGLES I-J-K AND J-K-L.
C
      DOTPJ = ZERO
      DOTPK = ZERO
      DO 220 M = 1,3
         DOTPJ = DOTPJ-EIJ(M)*EJK(M)
         DOTPK = DOTPK-EJK(M)*EKL(M)
  220 CONTINUE
      IF (ONE- ABS(DOTPJ)) 320,340,240
  240 IF (ONE- ABS(DOTPK)) 360,380,260
  260 SINPJ =  SQRT(ONE-DOTPJ*DOTPJ)
      SINPK =  SQRT(ONE-DOTPK*DOTPK)
C
C     ----- CALCULATING THE B MATRIX ELEMENTS.
C
      DIJ=DIJ*BOHR
      DJK=DJK*BOHR
      DKL=DKL*BOHR
      DO 280 M = 1,3
         SMI = -CR1(M)/(DIJ*SINPJ*SINPJ)
         NOCOL1 = 3*(I-1)+M
         B(NOCOL1,NOINT) = SMI
         F1 = (CR1(M)*(DJK-DIJ*DOTPJ))/(DJK*DIJ*SINPJ*SINPJ)
         F2 = (DOTPK*CR2(M))/(DJK*SINPK*SINPK)
         SMJ = F1-F2
         NOCOL2 = 3*(J-1)+M
         B(NOCOL2,NOINT) = SMJ
         SML = CR2(M)/(DKL*SINPK*SINPK)
         NOCOL3 = 3*(L-1)+M
         B(NOCOL3,NOINT) = SML
         NOCOL4 = 3*(K-1)+M
         B(NOCOL4,NOINT) = -(SMI+SMJ+SML)
  280 CONTINUE
C
C     ----- CALCULATING THE MAGNITUDE OF THE TORSION ANGLE, WHICH IS
C           BETWEEN 0.0 AND 180.0 DEGREES.
C
      DOT = ZERO
      DO 300 M = 1,3
         DOT = DOT+(CR1(M)*CR2(M))/(SINPJ*SINPK)
  300 CONTINUE
      IF( ABS(DOT) .GT. TOLRD) GO TO 400
      IF( ABS(DOT) .GT. ONE  ) DOT = SIGN(ONE,DOT)
      EQVAL = ACOS(DOT)/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
C
C     ----- CALCULATING THE TORSION ANGLE WITH THE CONVENTIONAL SIGN.
C
C                       -PI .LE. TORSION ANGLE .LE. +PI
C
C               -->  -->   -->
C           HERE D = CR1 X CR2
C
      D(1) = CR1(2)*CR2(3)-CR1(3)*CR2(2)
      D(2) = CR1(3)*CR2(1)-CR1(1)*CR2(3)
      D(3) = CR1(1)*CR2(2)-CR1(2)*CR2(1)
      EQVAL =   SIGN(EQVAL,(D(1)*EJK(1)+D(2)*EJK(2)+D(3)*EJK(3)))
      RETURN
C
C     ----- ERROR SITUATIONS.
C
  320 WRITE (IW,9999) I,J,K,DOTPJ
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
  340 WRITE (IW,9998) I,J,K
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
  360 WRITE (IW,9999) J,K,L,DOTPK
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
  380 WRITE (IW,9998) J,K,L
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
  400 WRITE (IW,9997) I,J,K,L
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
 9999 FORMAT(' COS ',I3,' - ',I3,' - ',I3,' = ',1PE12.5,'?',
     1       ' ROUNDOFF ERROR. RUN STOPPED.')
 9998 FORMAT(' ATOMS ',I3,1X,I3,1X,I3,' ARE COLLINEAR. RUN STOPPED.')
 9997 FORMAT(' ROUNDOFF ERROR IN CALCULATING TORSION ANGLE FOR',
     1       I3,' - ',I3,' - ',I3,' - ',I3,'. RUN STOPPED.')
      END
      SUBROUTINE HND_OPLA(EQVAL,NOINT,I,J,K,L,C,B,NCART,RTOD,PIO2,BOHR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     -----THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR AN
C          OUT OF PLANE WAGGING COORDINATE AS DEFINED BY WILSON.
C          SEE WDC P58. -----
C     -----I IS THE END ATOM, J IS THE APEX ATOM, AND
C          K AND L ARE THE ANCHOR ATOMS -----
C
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION B(NCART,1),C(3,1)
      DIMENSION RJI(3),RJK(3),RJL(3)
      DIMENSION EJI(3),EJK(3),EJL(3)
      DIMENSION C1(3),C2(3),C3(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','-OPLA  -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
C
      DJISQ = ZERO
      DJKSQ = ZERO
      DJLSQ = ZERO
      DO 100 M = 1,3
      RJI(M) = C(M,I)-C(M,J)
      DJISQ  = DJISQ+RJI(M)*RJI(M)
      RJK(M) = C(M,K)-C(M,J)
      DJKSQ  = DJKSQ+RJK(M)*RJK(M)
      RJL(M) = C(M,L)-C(M,J)
  100 DJLSQ  = DJLSQ+RJL(M)*RJL(M)
      DJI = SQRT(DJISQ)
      DJK = SQRT(DJKSQ)
      DJL = SQRT(DJLSQ)
      DO 120 M = 1,3
      EJI(M) = RJI(M)/DJI
      EJK(M) = RJK(M)/DJK
  120 EJL(M) = RJL(M)/DJL
      C1(1) = EJK(2)*EJL(3)-EJK(3)*EJL(2)
      C1(2) = EJK(3)*EJL(1)-EJK(1)*EJL(3)
      C1(3) = EJK(1)*EJL(2)-EJK(2)*EJL(1)
      C2(1) = EJL(2)*EJI(3)-EJL(3)*EJI(2)
      C2(2) = EJL(3)*EJI(1)-EJL(1)*EJI(3)
      C2(3) = EJL(1)*EJI(2)-EJL(2)*EJI(1)
      C3(1) = EJI(2)*EJK(3)-EJI(3)*EJK(2)
      C3(2) = EJI(3)*EJK(1)-EJI(1)*EJK(3)
      C3(3) = EJI(1)*EJK(2)-EJI(2)*EJK(1)
      DOTI = ZERO
      DO 140 M = 1,3
  140 DOTI = DOTI+EJK(M)*EJL(M)
      IF (ONE- ABS(DOTI)) 240,240,160
  160 SINI = SQRT(ONE-DOTI*DOTI)
      DET = EJI(1)*C1(1)+EJI(2)*C1(2)+EJI(3)*C1(3)
      SINT = DET/SINI
      IF (ONE- ABS(SINT)) 260,260,180
  180 COST = SQRT(ONE-SINT*SINT)
      TANT = SINT/COST
      DJI=DJI*BOHR
      DJK=DJK*BOHR
      DJL=DJL*BOHR
      DO 200 M = 1,3
      NOCOL1 = 3*(I-1)+M
      SMI = ((C1(M)/(COST*SINI))-(TANT*EJI(M)))/DJI
      B(NOCOL1,NOINT) = SMI
      NOCOL2 = 3*(K-1)+M
      SMK = ((C2(M)/(COST*SINI))-((TANT*(EJK(M)-DOTI*EJL(M)))/(SINI*
     1     SINI)))/DJK
      B(NOCOL2,NOINT) = SMK
      NOCOL3 = 3*(L-1)+M
      SML = ((C3(M)/(COST*SINI))-((TANT*(EJL(M)-DOTI*EJK(M)))/(SINI*
     1     SINI)))/DJL
      B(NOCOL3,NOINT) = SML
      NOCOL4 = 3*(J-1)+M
      B(NOCOL4,NOINT) = -(SMI+SMK+SML)
  200 CONTINUE
      DOT = ZERO
      DO 220 M = 1,3
  220 DOT = DOT+C1(M)*EJI(M)/SINI
      THETA = ACOS(DOT)
      THETA = PIO2-THETA
      EQVAL = THETA/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
      GO TO 260
  240 CONTINUE
      WRITE (IW,9008)
      CALL HND_HNDERR(3,ERRMSG)
  260 CONTINUE
      RETURN
 9008 FORMAT(' ROUNDOFF ERROR IN -OPLA- STOP ')
      END
      SUBROUTINE HND_LIBE(EQVAL,NOINT,I,J,K,L,C,B,NCART,RTOD,AA,NA,BOHR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- THIS ROUTINE COMPUTES THE B MATRIX ELEMENTS FOR A
C           PAIR OF PERPENDICULAR LINEAR BENDING COORDINATES.  SEE
C           FERIGLE AND MEISTER, J CHEM PHYS 19,982(1951) -----
C
C     ----- I AND K ARE THE END ATOMS AND J IS THE CENTRAL ATOM.
C           NO2 IS THE NUMBER OF THE SECOND INTERNAL COORDINATE,
C           I.E. THE ONE IN A PLANE PERPENDICULAR TO THE FIRST. A(I)
C           GIVES THE CARTESIAN COORDINATES OF A POINT IN SPACE SUCH
C           THAT THE VECTOR FROM ATOM J TO POINT A DEFINES
C           THE ORIENTATION OF THE COORDINATES IN SPACE.
C           THE FIRST OF THE TWO INTERNAL COORDINATES IS IN THE I-A-K
C           PLANE AND THE SECOND IN A PLANE PERPENDICULAR TO THE FIRST
C           THROUGH POINTS I,J, AND K -----
C
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION B(NCART,1),C(3,1)
      DIMENSION AA(3,1),A(3)
      DIMENSION RJI(3),RJK(3)
      DIMENSION UNIT(3),UP(3),UN(3)
      DIMENSION EJI(3),EJK(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','-LIBE  -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TOL,PIDEG /5.0D-05,180.00D+00/
      DATA PT0001    /0.0001D+00/
C
      IF(L.LE.NA) GO TO 120
      WRITE(IW,9048) L
      CALL HND_HNDERR(3,ERRMSG)
  120 CONTINUE
      A(1)=AA(1,L)
      A(2)=AA(2,L)
      A(3)=AA(3,L)
C
      DJISQ = ZERO
      DJKSQ = ZERO
      DAJSQ = ZERO
      DO 140 M = 1,3
      RJI(M) = C(M,I)-C(M,J)
      DJISQ  = DJISQ+RJI(M)*RJI(M)
      RJK(M) = C(M,K)-C(M,J)
      DJKSQ  = DJKSQ+RJK(M)*RJK(M)
      UN(M)  = A(M)-C(M,J)
  140 DAJSQ  = DAJSQ+UN(M)*UN(M)
      DJI = SQRT(DJISQ)
      DJK = SQRT(DJKSQ)
      DAJ = SQRT(DAJSQ)
      DOTJ = ZERO
      DOTP = ZERO
      DO 160 M = 1,3
      EJI(M) = RJI(M)/DJI
      EJK(M) = RJK(M)/DJK
      UNIT(M) = UN(M)/DAJ
      DOTJ = DOTJ+EJI(M)*EJK(M)
  160 DOTP = DOTP+EJI(M)*UNIT(M)
      TEST = ( ABS(DOTJ)-ONE)
      IF (PT0001- ABS(TEST)) 280,180,180
  180 CONTINUE
      TEST = ( ABS(DOTP)-ONE)
      IF (PT0001- ABS(TEST)) 200,300,300
  200 CONTINUE
      UP(1) = EJK(2)*UNIT(3)-EJK(3)*UNIT(2)
      UP(2) = EJK(3)*UNIT(1)-EJK(1)*UNIT(3)
      UP(3) = EJK(1)*UNIT(2)-EJK(2)*UNIT(1)
      DUM=ONE/ SQRT(UP(1)**2+UP(2)**2+UP(3)**2)
      UP(1)=UP(1)*DUM
      UP(2)=UP(2)*DUM
      UP(3)=UP(3)*DUM
      UNIT(1) = UP(2)*EJK(3)-UP(3)*EJK(2)
      UNIT(2) = UP(3)*EJK(1)-UP(1)*EJK(3)
      UNIT(3) = UP(1)*EJK(2)-UP(2)*EJK(1)
      DJI=DJI*BOHR
      DJK=DJK*BOHR
C
C     ----- FIRST COMPONENT -----
C
      DO 220 M = 1,3
      NOCOL1 = 3*(I-1)+M
      B(NOCOL1,NOINT) = -UNIT(M)/DJI
      NOCOL2 = 3*(K-1)+M
      B(NOCOL2,NOINT) = -UNIT(M)/DJK
      NOCOL3 = 3*(J-1)+M
  220 B(NOCOL3,NOINT) = (ONE/DJI+ONE/DJK)*UNIT(M)
      EQVAL = ACOS(DOTJ)/RTOD
      IF( ABS(EQVAL      ).LT.TOL) EQVAL=ZERO
      IF( ABS(EQVAL-PIDEG).LT.TOL) EQVAL=PIDEG
C
C     ----- SECOND COMPONENT -----
C
      NO2=NOINT+1
      DO 260 M = 1,3
      NOCOL1 = 3*(I-1)+M
      B(NOCOL1,NO2) = -UP(M)/DJI
      NOCOL2 = 3*(K-1)+M
      B(NOCOL2,NO2) = -UP(M)/DJK
      NOCOL3 = 3*(J-1)+M
  260 B(NOCOL3,NO2) = (ONE/DJI+ONE/DJK)*UP(M)
      GO TO 320
  280 CONTINUE
      WRITE (IW,9008) I,J,K,DOTJ
      CALL HND_HNDERR(3,ERRMSG)
  300 CONTINUE
      WRITE (IW,9028) I,J,K,L
      CALL HND_HNDERR(3,ERRMSG)
  320 CONTINUE
      RETURN
 9008 FORMAT(' IN -LIBE- I,J,K ARE NOT COLLINEAR. I,J,K,DOTJ = ',
     1 3I5,F15.9)
 9028 FORMAT(' IN -LIBE- DUMMY ATOM IS COLLINEAR WITH I,J,K.',
     1 ' I,J,K,L = ',4I5)
 9048 FORMAT(' IN -LIBE- DUMMY ATOM IS OUT OF RANGE. L = ',I5)
      END
      SUBROUTINE HND_DIHPLA(DIHANG,NOINT,IZMAT,CARTC,BMAT,NCART,
     1                      DTORAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION IZMAT(5), CARTC(3,1), BMAT(NCART,1),
     1          A(3), B(3), C(3), D(3), E1(3), E2(3), E3(3)
      DIMENSION ERRMSG(3)
      DATA ERRMSG   /'PROGRAM ','STOP IN ','-DIHPLA-'/
      DATA ZERO,ONE /0.0D+00,1.0D+00/
      DATA TOL      /1.0D-06/
C
C
C     COMPUTE THE B MATRIX AND DIHEDRAL ANGLE BETWEEN 5 ATOMS
C     FORMING TWO PLANES HAVING ONLY ONE ATOM IN COMMON.
C
C     WRITTEN BY STEPHEN T. ELBERT, AMES LABORATORY, ISU, 2/23/83
C
C     ON ENTRY,
C     IW     - PRINT MESSAGE FILE
C     NOINT  - NUMBER OF THE INTERNAL COORDINATE BEING EVALUATED
C     NCART  - NUMBER OF COORDINATES (3 * NATOMS)
C     IZMAT  - LIST OF THE FIVE ATOMS INVOLVED
C     CARTC  - CARTESIAN COORDINATES OF THE ATOMS
C     DTORAD - CONVERSION FACTOR (RADIANS TO DEGREES)
C
C     ON EXIT,
C     DIHANG - DIHEDRAL ANGLE, IN DEGREES
C     BMAT   - THE NOINTth COLUMN IS COMPUTED BY THIS ROUTINE
C
C     THE FIRST THREE ATOMS DEFINE THE FIRST PLANE,
C     THE LAST THREE ATOMS DEFINE THE SECOND PLANE,
C     THE MIDDLE ATOM IS IN BOTH PLANES.
C
      IATOM = IZMAT(1)
      JATOM = IZMAT(2)
      KATOM = IZMAT(3)
      LATOM = IZMAT(4)
      MATOM = IZMAT(5)
C
C     DEFINE VECTORS FROM THE PIVOT ATOM TO THE OTHER ATOMS
C     A= I-K, B=J-K, C=L-K, D=M-K
C
      DO 110 IXYZ = 1,3
         A(IXYZ) = CARTC(IXYZ,IATOM) - CARTC(IXYZ,KATOM)
         B(IXYZ) = CARTC(IXYZ,JATOM) - CARTC(IXYZ,KATOM)
         C(IXYZ) = CARTC(IXYZ,LATOM) - CARTC(IXYZ,KATOM)
         D(IXYZ) = CARTC(IXYZ,MATOM) - CARTC(IXYZ,KATOM)
  110 CONTINUE
C
C     E1 = A X B   IS PERPENDICULAR TO FIRST PLANE
C
      E1(1) = A(2)*B(3) - A(3)*B(2)
      E1(2) = A(3)*B(1) - A(1)*B(3)
      E1(3) = A(1)*B(2) - A(2)*B(1)
      E1MAG =  SQRT( E1(1)**2 + E1(2)**2 + E1(3)**2)
      IF(E1MAG.LT.TOL) GO TO 810
C
C     E2 = C X D   IS PERPENDICULAR TO SECOND PLANE
C
      E2(1) = C(2)*D(3) - C(3)*D(2)
      E2(2) = C(3)*D(1) - C(1)*D(3)
      E2(3) = C(1)*D(2) - C(2)*D(1)
      E2MAG =  SQRT( E2(1)**2 + E2(2)**2 + E2(3)**2)
      IF(E2MAG.LT.TOL) GO TO 820
C
C     NORMALIZE E1 AND E2
C
      E1MAGI = ONE/E1MAG
      E2MAGI = ONE/E2MAG
      DO 120 IXYZ=1,3
         E1(IXYZ) = E1(IXYZ) * E1MAGI
         E2(IXYZ) = E2(IXYZ) * E2MAGI
  120 CONTINUE
C
C     THE DIRECTION OF E2 IS DEFINED SUCH THAT WHEN THE PLANES ARE
C     ROTATED ALONG THE INTERSECTION E3 = E1 X E2 SO THEY COINCIDE
C     (DIHEDRAL ANGLE=0), THE ATOMS 1 AND 4 LIE ON THE SAME SIDE OF E3.
C     (IF ATOMS 1 AND/OR 4 LIE ON E3, THEN ATOMS 2 AND/OR 5 ARE USED)
C
      ADOTE2 = A(1)*E2(1) + A(2)*E2(2) + A(3)*E2(3)
      IF(ADOTE2.LT.TOL) ADOTE2 = B(1)*E2(1) + B(2)*E2(2) + B(3)*E2(3)
      CDOTE1 = C(1)*E1(1) + C(2)*E1(2) + C(3)*E1(3)
      IF(CDOTE1.LT.TOL) CDOTE1 = D(1)*E1(1) + D(2)*E1(2) + D(3)*E1(3)
C
C     POSSIBLY CHANGE SENSE OF E2 SO A (OR B) WILL BE ON THE SAME
C     SIDE OF E1 X E3 AS C (OR D)
C
      IF(ADOTE2*CDOTE1 .GT.ZERO) GO TO 140
         ADOTE2 = -ADOTE2
         DO 130 IXYZ=1,3
            E2(IXYZ) = -E2(IXYZ)
  130    CONTINUE
  140 CONTINUE
C
C      EVALUATE DIHEDRAL ANGLE    0 < ABS(DIHANG) < 180
C      ANGLE IS POSITIVE FOR SECOND PLANE 'ABOVE' ATOM I
C      ANGLE IS NEGATIVE FOR SECOND PLANE 'BELOW' ATOM I
C
      E1DE2 = E1(1)*E2(1) + E1(2)*E2(2) + E1(3)*E2(3)
      IF((      ABS(E1DE2)-ONE)  .GT.ZERO) GO TO 840
      IF(( ABS( ABS(E1DE2)-ONE)) .LT.TOL ) GO TO 830
      DIHANG = ACOS(E1DE2)
      IF(ADOTE2.GT.ZERO) DIHANG = -DIHANG
      SINDI = ONE / SIN(DIHANG)
      DIHANG = DIHANG/DTORAD
C
C     E3 = E1 X E2 IS ALONG INTERSECTION LINE OF THE PLANES
C
      E3(1) = SINDI * (E1(2)*E2(3) - E1(3)*E2(2))
      E3(2) = SINDI * (E1(3)*E2(1) - E1(1)*E2(3))
      E3(3) = SINDI * (E1(1)*E2(2) - E1(2)*E2(1))
C
      F1 =  E1MAGI * (B(1)*E3(1) + B(2)*E3(2) + B(3)*E3(3))
      F2 = -E1MAGI * (A(1)*E3(1) + A(2)*E3(2) + A(3)*E3(3))
      F4 = -E2MAGI * (D(1)*E3(1) + D(2)*E3(2) + D(3)*E3(3))
      F5 =  E2MAGI * (C(1)*E3(1) + C(2)*E3(2) + C(3)*E3(3))
C
      I = 3*IATOM - 3
      J = 3*JATOM - 3
      K = 3*KATOM - 3
      L = 3*LATOM - 3
      M = 3*MATOM - 3
C
      DO 160 N=1,3
         B1 = E1(N) * F1
         B2 = E1(N) * F2
         B4 = E2(N) * F4
         B5 = E2(N) * F5
         B3 = -(B1+B2+B4+B5)
         BMAT(NOINT,I+N) = B1
         BMAT(NOINT,J+N) = B2
         BMAT(NOINT,K+N) = B3
         BMAT(NOINT,L+N) = B4
         BMAT(NOINT,M+N) = B5
  160 CONTINUE
      RETURN
C                                        ERROR MESSAGES
  810 WRITE(IW,910)
      GO TO 880
  820 WRITE(IW,920)
      GO TO 880
  830 WRITE(IW,930)
      GO TO 880
  840 WRITE(IW,940) E1DE2
  880 CONTINUE
      WRITE(IW,950) NOINT,IATOM,JATOM,KATOM,LATOM,MATOM
      CALL HND_HNDERR(3,ERRMSG)
C
  910 FORMAT(1X,'---- ERROR, FIRST THREE ATOMS ARE COLLINEAR')
  920 FORMAT(1X,'---- ERROR, LAST THREE ATOMS ARE COLLINEAR ')
  930 FORMAT(1X,'---- ERROR, ALL FIVE ATOMS ARE COPLANAR ',
     1          '     THE -B- MATRIX COLUMN CAN NOT BE DEFINED. STOP')
  940 FORMAT(1X,'---- COSIN OF DIHEDRAL ANGLE=',F16.10,
     1          ' IS GREATER THAN ONE')
  950 FORMAT(1X,'THIS MESSAGE COMES TO YOU FROM -DIHPLA-, WHERE ',
     1          'INTERNAL COORDINATE',I4,' USES ATOMS',5I4)
      END
      SUBROUTINE HND_ZMTMAK(ZDONE,RTDB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "geom.fh"
      integer  rtdb
      logical  status
      logical  zcoord_input
      parameter (max_zcoord=16)
c
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (MXATOM=500)
      PARAMETER (MXCART=3*MXATOM)
      PARAMETER (MXZMAT=1500)    
      PARAMETER (MXCOOR=1500)    
      PARAMETER (MXIZMT=6000)
      PARAMETER (MXBOND= 64)
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXANGS= 8*MXATOM)
      PARAMETER (MXTORS= 8*MXATOM)
      PARAMETER (MXOOPA= 8*MXATOM)
      PARAMETER (MXLINB= 8*MXATOM)
      PARAMETER (MXLNBA=10       )
      PARAMETER (MXSEG=64)
      LOGICAL SOME,OUT,DBUG
      LOGICAL DBUG1
      LOGICAL BADZ
      LOGICAL REMOVE
      LOGICAL NUBOND
      LOGICAL MODBND
      LOGICAL MODANG
      LOGICAL MODTOR
      LOGICAL MODOOP
      LOGICAL MODLNB
      LOGICAL ENDATM
      LOGICAL ENDMOD
      LOGICAL  DONE
      LOGICAL IDONE
      LOGICAL ZDONE
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_MOLNUC/NUC(MXATOM)  
      COMMON/HND_MOLXYZ/C(3,MXATOM),ZAN(MXATOM),NAT
      COMMON/HND_ZMTPAR/NZMAT,NZVAR,NVAR
      COMMON/HND_ZMTDAT/ZMAT(MXCOOR),IZMAT(MXIZMT),NIZMAT
      DIMENSION XYZLNB(3,MXLNBA)
      DIMENSION IJBOND(2,MXBOND)
      DIMENSION IJKANG(3,MXANGS)
      DIMENSION IJKLTO(4,MXTORS)
      DIMENSION IJKLOP(4,MXOOPA)
      DIMENSION IJKLNB(4,MXLINB)
      DIMENSION IJBNDS(2,MXBNDS)
      DIMENSION IJMODS(2,MXBNDS)
      DIMENSION NIBOND(MXATOM),IIBOND(MXATOM)
      DIMENSION NIMODS(MXATOM),IIMODS(MXATOM)
      DIMENSION ENDATM(MXATOM)
      DIMENSION ENDMOD(MXATOM)
      DIMENSION  IDONE(MXATOM)
      DIMENSION IATSEG(MXATOM,MXSEG)
      DIMENSION NAMSEG(       MXSEG)
      DIMENSION LENSEG(       MXSEG)
      DIMENSION RCOV(103)
      DIMENSION ERRMSG(3)
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTMAK-'/
      DATA ZERO   /0.0D+00/
C
C      COVALENT RADII FROM J.EMSLEY, "THE ELEMENTS", 2ND EDITION, 1991
C      EXCEPT VAN DER WAALS RADII FOR HE,NE,AR,KR (SAME SOURCE),
C      AND GUESSES FOR NA,V,CR,RB,TC,PM,EU,YB,AT,RN
C
      DATA (RCOV(NUCZ),NUCZ=1,2)
     1  /0.30D+00,1.22D+00/
      DATA (RCOV(NUCZ),NUCZ=3,10)
     1  /1.23D+00,0.89D+00,0.88D+00,0.77D+00,
     2   0.70D+00,0.66D+00,0.58D+00,1.60D+00/
      DATA (RCOV(NUCZ),NUCZ=11,18)
     1  /1.66D+00,1.36D+00,1.25D+00,1.17D+00,
     2   1.10D+00,1.04D+00,0.99D+00,1.91D+00/
      DATA (RCOV(NUCZ),NUCZ=19,36)
     1  /2.03D+00,1.74D+00,
     2   1.44D+00,1.32D+00,1.22D+00,1.19D+00,1.17D+00,
     3  1.165D+00,1.16D+00,1.15D+00,1.17D+00,1.25D+00,
     4   1.25D+00,1.22D+00,1.21D+00,1.17D+00,1.14D+00,1.98D+00/
      DATA (RCOV(NUCZ),NUCZ=37,54)
     1  /2.22D+00,1.92D+00,
     2   1.62D+00,1.45D+00,1.34D+00,1.29D+00,1.27D+00,
     3   1.24D+00,1.25D+00,1.28D+00,1.34D+00,1.41D+00,
     4   1.50D+00,1.40D+00,1.41D+00,1.37D+00,1.33D+00,2.09D+00/
      DATA (RCOV(NUCZ),NUCZ=55,86)
     1  /2.35D+00,1.98D+00,
     2   1.69D+00,1.65D+00,1.65D+00,1.64D+00,1.65D+00,1.66D+00,1.65D+00,
     3   1.61D+00,1.59D+00,1.59D+00,1.58D+00,1.57D+00,1.56D+00,1.56D+00,
     4   1.56D+00,1.44D+00,1.34D+00,1.30D+00,1.28D+00,
     5   1.26D+00,1.26D+00,1.29D+00,1.34D+00,1.44D+00,
     6   1.55D+00,1.54D+00,1.52D+00,1.53D+00,1.50D+00,2.20D+00/
      DATA (RCOV(I),i=87,103) 
     1  /3.24d00,2.68d00,2.25d00,2.16d00,1.93d00,
     2   1.66d00,1.57d00,1.81d00,2.21d00,1.43d00,
     3   1.42d00,1.40d00,1.39d00,1.38d00,1.37d00,
     4   1.36d00,1.34d00/
C
      DIST(IAT,JAT)=SQRT((C(1,IAT)-C(1,JAT))**2+
     1                   (C(2,IAT)-C(2,JAT))**2+
     2                   (C(3,IAT)-C(3,JAT))**2) * TOANGS
      BADZ(NUCZ) = NUCZ.LE.0.OR.NUCZ.GT.103
C
      DBUG1=.FALSE.
      DBUG =.FALSE.
      OUT  =.FALSE.
      OUT  =OUT.OR.DBUG
      SOME =.TRUE.
      SOME =SOME.OR.OUT
      IF(SOME) THEN
         WRITE(IW,8888)
      ENDIF
      IF(DBUG) THEN
         DO IAT=1,NAT
            WRITE(IW,9976) IAT,NUC(IAT),C(1,IAT),C(2,IAT),C(3,IAT)
         ENDDO
      ENDIF
C
      NZVAR =0
      NIZMAT=0
      DO I=1,MXIZMT
         IZMAT(I)=0
      ENDDO
      DO I=1,MXCOOR
         ZMAT(I)=ZERO
      ENDDO
      DO MBNDS=1,MXBNDS
         IJBNDS(1,MBNDS)=0
         IJBNDS(2,MBNDS)=0
      ENDDO
      DO MBOND=1,MXBOND
         IJBOND(1,MBOND)=0
         IJBOND(2,MBOND)=0
      ENDDO
      DO MANGS=1,MXANGS
         DO I=1,3
            IJKANG(I,MANGS)=0
         ENDDO
      ENDDO
      DO MTORS=1,MXTORS
         DO I=1,4
            IJKLTO(I,MTORS)=0
         ENDDO
      ENDDO
      DO MOOPA=1,MXOOPA
         DO I=1,4
            IJKLOP(I,MOOPA)=0
         ENDDO
      ENDDO
      DO MLINB=1,MXLINB
         DO I=1,4
            IJKLNB(I,MLINB)=0
         ENDDO
      ENDDO
      DO MLNBA=1,MXLNBA
         DO I=1,3
            XYZLNB(I,MLNBA)=ZERO
         ENDDO
      ENDDO
      CVR_SCALING=1.20D+00
C
C     ----- READ IMPOSED BONDS, ANGLES, TORSIONS -----    
C     ----- O-O-PLANE, LINEAR BENDS, ... IF ANY  -----
C
c           read in zcoord data from -rtdb-
c    CVR_SCALING,IJBOND,IJKANG,IJKLTO,IJKLOP,IJKLNB
c
      status=rtdb_get(rtdb,'autoz:zcoord',mt_log,1,zcoord_input)
      status=status.or.        
     1       rtdb_get(rtdb,'autoz:cvr_scaling',mt_dbl,1,  
     2                                              cvr_scaling)
      status=status.or.        
     1       rtdb_get(rtdb,'autoz:ijbond',mt_int,2*max_zcoord,  
     2                                                 ijbond)
      status=status.or.
     1       rtdb_get(rtdb,'autoz:ijkang',mt_int,3*max_zcoord,  
     2                                                 ijkang)
      status=status.or.
     1       rtdb_get(rtdb,'autoz:ijklto',mt_int,4*max_zcoord,  
     2                                                 ijklto)
      status=status.or.
     1       rtdb_get(rtdb,'autoz:ijklop',mt_int,4*max_zcoord,  
     2                                                 ijklop)
      status=status.or.
     1       rtdb_get(rtdb,'autoz:ijklnb',mt_int,4*max_zcoord,  
     2                                                 ijklnb)
      if(status) then
         if(dbug) then
            write(iw,*) 'zcoord data read from rtdb'
         endif
      else
         if(dbug) then
            write(iw,*) 'no zcoord data found on rtdb'
         endif
      endif
      if(dbug) then
         write(iw,*) 'in hnd-zmtmak, cvr_scaling value = ',
     1                cvr_scaling
      endif
c
      CVFAC=CVR_SCALING
      IF(DBUG) THEN
         WRITE(IW,*) 'CVFAC = ',CVFAC
      ENDIF
C
      NUMBND=0
      DO MBOND=1,MXBOND
         ICON=IJBOND(1,MBOND)
         JCON=IJBOND(2,MBOND)
         IF(ICON.GT.0.AND.JCON.GT.0) THEN
            NUMBND=NUMBND+1
         ELSEIF(ICON.EQ.0.AND.JCON.EQ.0) THEN
         ELSE
            WRITE(IW,9977)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      MODBND=NUMBND.NE.0
      IF(DBUG) THEN
         IF(MODBND) THEN
            WRITE(IW,8887)
         ELSE
            WRITE(IW,8886)
         ENDIF
      ENDIF
C
      NUMANG=0
      DO MANGS=1,MXANGS
         ICON=IJKANG(1,MANGS)
         JCON=IJKANG(2,MANGS)
         KCON=IJKANG(3,MANGS)
         IF(ICON.GT.0.AND.JCON.GT.0.AND.KCON.GT.0) THEN
            NUMANG=NUMANG+1
         ELSEIF(ICON.EQ.0.AND.JCON.EQ.0.AND.KCON.EQ.0) THEN
         ELSE
            WRITE(IW,9977)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      MODANG=NUMANG.NE.0
      IF(DBUG) THEN
         IF(MODANG) THEN
            WRITE(IW,8885)
         ELSE
            WRITE(IW,8884)
         ENDIF
      ENDIF
C
      NUMTOR=0
      DO MTORS=1,MXTORS
         ICON=IJKLTO(1,MTORS)
         JCON=IJKLTO(2,MTORS)
         KCON=IJKLTO(3,MTORS)
         LCON=IJKLTO(4,MTORS)
         IF(ICON.GT.0.AND.JCON.GT.0.AND.
     1      KCON.GT.0.AND.LCON.GT.0     ) THEN
            NUMTOR=NUMTOR+1
         ELSEIF(ICON.EQ.0.AND.JCON.EQ.0.AND. 
     1          KCON.EQ.0.AND.LCON.EQ.0     ) THEN
         ELSE
            WRITE(IW,9977)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      MODTOR=NUMTOR.NE.0
      IF(DBUG) THEN
         IF(MODTOR) THEN
            WRITE(IW,8883)
         ELSE
            WRITE(IW,8882)
         ENDIF
      ENDIF
C
      NUMOOP=0
      DO MOOPA=1,MXOOPA
         ICON=IJKLOP(1,MOOPA)
         JCON=IJKLOP(2,MOOPA)
         KCON=IJKLOP(3,MOOPA)
         LCON=IJKLOP(4,MOOPA)
         IF(ICON.GT.0.AND.JCON.GT.0.AND.
     1      KCON.GT.0.AND.LCON.GT.0     ) THEN
            NUMOOP=NUMOOP+1
         ELSEIF(ICON.EQ.0.AND.JCON.EQ.0.AND.
     1          KCON.EQ.0.AND.LCON.EQ.0     ) THEN
         ELSE
            WRITE(IW,9977)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      MODOOP=NUMOOP.NE.0
      IF(DBUG) THEN
         IF(MODOOP) THEN
            WRITE(IW,8881)
         ELSE
            WRITE(IW,8880)
         ENDIF
      ENDIF
C
      NUMLNB=0
      DO MLINB=1,MXLINB
         ICON=IJKLNB(1,MLINB)
         JCON=IJKLNB(2,MLINB)
         KCON=IJKLNB(3,MLINB)
         LCON=IJKLNB(4,MLINB)
         IF(ICON.GT.0.AND.JCON.GT.0.AND.
     1      KCON.GT.0.AND.LCON.GT.0     ) THEN
            NUMLNB=NUMLNB+1
         ELSEIF(ICON.EQ.0.AND.JCON.EQ.0.AND.
     1          KCON.EQ.0.AND.LCON.EQ.0     ) THEN
         ELSE
            WRITE(IW,9977)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDDO
      MODLNB=NUMLNB.NE.0
      IF(DBUG) THEN
         IF(MODLNB) THEN
            WRITE(IW,8879)
         ELSE
            WRITE(IW,8878)
         ENDIF
      ENDIF
C
C     ----- SET UP CONNECTIVITY TABLE, INCLUDING FORCED BONDS -----
C
      IF(DBUG) THEN
         DO IAT=1,NAT
            WRITE(IW,9976) IAT,NUC(IAT),C(1,IAT),C(2,IAT),C(3,IAT)
         ENDDO
      ENDIF
C
      MBNDS=0
      DO IAT=1,NAT
      IIBOND(IAT)=MBNDS
C
C     ----- FIRST CHECK FORCED BONDS -----
C
         IF(MODBND) THEN
            DO MBOND=1,MXBOND
               NUBOND=.FALSE.
               IF(IJBOND(1,MBOND).EQ.IAT) THEN
                  MBNDS=MBNDS+1
                  IJBNDS(1,MBNDS)=IAT
                  IJBNDS(2,MBNDS)=IJBOND(2,MBOND)
                  NUBOND=.TRUE.
               ELSEIF(IJBOND(2,MBOND).EQ.IAT) THEN
                  MBNDS=MBNDS+1
                  IJBNDS(1,MBNDS)=IAT
                  IJBNDS(2,MBNDS)=IJBOND(1,MBOND)
                  NUBOND=.TRUE.
               ENDIF
               IF(NUBOND.AND.DBUG) THEN
                  WRITE(IW,9975) MBNDS,IJBNDS(1,MBNDS),IJBNDS(2,MBNDS)
               ENDIF
C
C     ----- CHECK FOR DUPLICATES AND ELIMINATE -----
C
               IF(NUBOND.AND.MBNDS.GT.1) THEN
                  M=MBNDS-1
                  DO I=1,M
                     IF(IJBNDS(1,I).EQ.IJBNDS(1,MBNDS).AND.
     1                  IJBNDS(2,I).EQ.IJBNDS(2,MBNDS)     ) THEN
                        MBNDS=MBNDS-1
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
C
C     ----- DISCOVER CONNECTIVITY TO OTHER ATOMS -----
C
         DO JAT=1,NAT
            IF(JAT.NE.IAT) THEN
               IZ=NUC(IAT)
               JZ=NUC(JAT)
               IF(BADZ(IZ).OR.BADZ(JZ)) THEN
                  WRITE(IW,9993) IAT,JAT
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
               NUBOND=.FALSE.
               RIJ=DIST(IAT,JAT)
               RCV=(RCOV(IZ)+RCOV(JZ))
               IF(RIJ.LE.RCV*CVFAC) THEN
                  IF(MBNDS.EQ.MXBNDS) THEN
                     WRITE(IW,*) 'TOO MANY BONDS FOUND',
     1                           ' MBNDS, MXBNDS = ',MBNDS, MXBNDS
                     CALL HND_HNDERR(3,ERRMSG)
                  ENDIF
                  MBNDS=MBNDS+1
                  IJBNDS(1,MBNDS)=IAT
                  IJBNDS(2,MBNDS)=JAT
                  NUBOND=.TRUE.
                  IF(DBUG1) THEN
                     WRITE(IW,9992) IAT,JAT,IZ,JZ,RIJ,RCV,MBNDS,
     1                              IJBNDS(1,MBNDS),IJBNDS(2,MBNDS),
     2                              NUBOND
                  ENDIF
C
C     ----- CHECK FOR DUPLICATES AND ELIMINATE -----
C
                  IF(NUBOND.AND.MBNDS.GT.1) THEN
                     REMOVE=.FALSE.
                     M=MBNDS-1
                     DO I=1,M
                        IF(IJBNDS(1,I).EQ.IJBNDS(1,MBNDS).AND.
     1                     IJBNDS(2,I).EQ.IJBNDS(2,MBNDS)     ) THEN
                           IF(DBUG) THEN
                              WRITE(IW,9974) MBNDS,I
                           ENDIF
                           REMOVE=.TRUE.
                        ENDIF
                     ENDDO
                     IF(REMOVE) THEN
                        IJBNDS(1,MBNDS)=0
                        IJBNDS(2,MBNDS)=0
                        MBNDS=MBNDS-1
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
C
         NIBOND(IAT)=MBNDS-IIBOND(IAT)
         IF(NIBOND(IAT).GT.1) THEN
            ENDATM(IAT)=.FALSE.
         ELSE
            ENDATM(IAT)=.TRUE.
         ENDIF
      ENDDO
      NBNDS=MBNDS
C
      IF(NBNDS.LE.0) THEN
         WRITE(IW,9999)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- DEBUG PRINTING OF CONNECTIVITY -----
C
      IF(DBUG1) THEN
         WRITE(IW,9998) NBNDS
         WRITE(IW,9997)
         DO IAT=1,NAT
            WRITE(IW,9996) IAT,NIBOND(IAT),IIBOND(IAT),ENDATM(IAT)
         ENDDO
         DO IAT=1,NAT
            WRITE(IW,9995) IAT,NIBOND(IAT)
            IF(NIBOND(IAT).GT.0) THEN
               IIBND1=IIBOND(IAT)+1
               IIBND2=IIBOND(IAT)+NIBOND(IAT)
               WRITE(IW,9994) (IJBNDS(1,I),IJBNDS(2,I),I=IIBND1,IIBND2)
            ENDIF
         ENDDO
      ENDIF
C
C
C     -----    TERMINAL ATOMS HAVE BEEN DETECTED   -----
C           TRY FOR SEGMENTS ( RINGS AND BRIDGES )
C
      DO I=1,NBNDS
         IJMODS(1,I)=IJBNDS(1,I)
         IJMODS(2,I)=IJBNDS(2,I)
      ENDDO
      MXCONN=0
      DO IAT=1,NAT
         ENDMOD(IAT)=ENDATM(IAT)
         IIMODS(IAT)=IIBOND(IAT)
         NIMODS(IAT)=NIBOND(IAT)
         IF(NIMODS(IAT).GT.MXCONN) THEN
            MXCONN=NIMODS(IAT)
         ENDIF
         DO MSEG=1,MXSEG
            IATSEG(IAT,MSEG)=0
         ENDDO
      ENDDO
      DO MSEG=1,MXSEG
         NAMSEG(MSEG)=0
         LENSEG(MSEG)=0
      ENDDO
C
C     -----     TAKE AWAY -END ATOM- CONNECTIVITY       -----
C           RING AND BRIDGE CONNECTIVITY SHOULD BE LEFT
C
      NPASS=MXCONN
      DO IPASS=1,NPASS
         IF(DBUG) THEN
            WRITE(IW,9991) IPASS,NPASS
         ENDIF
         DO IAT=1,NAT
            IF(NIMODS(IAT).EQ.1) THEN
               DO I=IIBOND(IAT)+1,IIBOND(IAT)+NIBOND(IAT)
                  KAT=IJMODS(1,I)
                  JAT=IJMODS(2,I)
                  IJMODS(1,I)=0
                  IJMODS(2,I)=0
                  IF(KAT.NE.0) THEN
                     IF(NIMODS(JAT).GE.1) THEN
                        DO J=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
                           IF(IJMODS(2,J).EQ.IAT) THEN
                              IJMODS(1,J)=0
                              IJMODS(2,J)=0
                              NIMODS(JAT)=NIMODS(JAT)-1
                           ENDIF
                        ENDDO
                     ELSE
                        WRITE(IW,9978) JAT,IAT
                        CALL HND_HNDERR(3,ERRMSG)
                     ENDIF
                     IF(NIMODS(JAT).EQ.1) THEN
                        ENDMOD(JAT)=.TRUE.
                     ENDIF
                  ENDIF
               ENDDO
               NIMODS(IAT)=NIMODS(IAT)-1
            ENDIF
         ENDDO
C
         IF(DBUG) THEN
            DO IAT=1,NAT
               WRITE(IW,9996) IAT,NIMODS(IAT),IIMODS(IAT),ENDMOD(IAT)
            ENDDO
            DO IAT=1,NAT
               WRITE(IW,9995) IAT,NIMODS(IAT)
               IF(NIBOND(IAT).GT.0) THEN
                  IIMOD1=IIBOND(IAT)+1
                  IIMOD2=IIBOND(IAT)+NIBOND(IAT)
                  WRITE(IW,9994) (IJMODS(1,I),IJMODS(2,I),
     1                                                 I=IIMOD1,IIMOD2)
               ENDIF
            ENDDO
         ENDIF
C
      ENDDO
C
C     ----- CHECK IF THERE WILL BE MORE TO DO ... -----
C
      DONE=.TRUE.
      DO IAT=1,NAT
         IF(NIMODS(IAT).GT.0) THEN
            IDONE(IAT)=.FALSE.
         ELSE
            IDONE(IAT)=.TRUE. 
         ENDIF
         DONE=DONE.AND.IDONE(IAT)
      ENDDO
      DONE=.TRUE.
      IF(DONE) THEN
         IF(DBUG) THEN
            WRITE(IW,9980)
         ENDIF
C
C     ----- CREATE -IZMAT- -----
C
         NZVAR =0
         NIZMAT=0
         DO I=1,MXIZMT
            IZMAT(I)=0
         ENDDO
         DO I=1,NBNDS
            IJMODS(1,I)=IJBNDS(1,I)
            IJMODS(2,I)=IJBNDS(2,I)
         ENDDO
         DO IAT=1,NAT
            ENDMOD(IAT)=ENDATM(IAT)
            IIMODS(IAT)=IIBOND(IAT)
            NIMODS(IAT)=NIBOND(IAT)
         ENDDO
         CALL HND_ZMTBLD(ZDONE,NZVAR,NIZMAT,IZMAT,C,
     1                   NAT,NBNDS,IJMODS,IIMODS,NIMODS,NUMBND,
     2                   IJKANG,NUMANG,IJKLTO,NUMTOR,
     3                   IJKLOP,NUMOOP,IJKLNB,NUMLNB,XYZLNB)
C
C     ----- SAVE -IZMAT(1:NZMAT)-   -----
C           SAVE -XYZLNB(1:MXLNBA)-
C
         NZMAT=NIZMAT
         NLNBA=3*MXLNBA
C
         IF(SOME) THEN
            WRITE(IW,8877) ZDONE
         ENDIF
         RETURN
      ENDIF
C
C     ----- NOW COMPRESS MODIFIED CONNECTIVITY TABLE -----
C
      DO IAT=1,NAT
         IF(NIMODS(IAT).GT.0) THEN
            NUMODS=0
            DO I=IIBOND(IAT)+1,IIBOND(IAT)+NIBOND(IAT)
               IF(IJMODS(1,I).NE.0.AND.IJMODS(2,I).NE.0) THEN
                  NUMODS=NUMODS+1
                  IJMODS(1,NUMODS)=IJBNDS(1,I)
                  IJMODS(2,NUMODS)=IJBNDS(2,I)
               ENDIF
            ENDDO
            IF(NIMODS(IAT).NE.NUMODS) THEN
               WRITE(IW,9985) IAT,NIMODS(IAT),NUMODS
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            IIBOND(IAT)=IIMODS(IAT)
            NIBOND(IAT)=IIMODS(IAT)
            ENDATM(IAT)=ENDMOD(IAT)
         ENDIF
      ENDDO
C
      IF(DBUG) THEN
         DO IAT=1,NAT
            WRITE(IW,9996) IAT,NIBOND(IAT),IIBOND(IAT),ENDATM(IAT)
         ENDDO
         DO IAT=1,NAT
            WRITE(IW,9995) IAT,NIMODS(IAT)
            IF(NIMODS(IAT).GT.0) THEN
               IIMOD1=IIMODS(IAT)+1
               IIMOD2=IIMODS(IAT)+NIMODS(IAT)
               WRITE(IW,9994) (IJMODS(1,I),IJMODS(2,I),I=IIMOD1,IIMOD2)
            ENDIF
         ENDDO
      ENDIF
C
C     ----- NOW WE HAVE ONLY ATOMS WITH REMAINING CONNECTIVITY -----
C           OF 2 OR MORE. THEY BELONG TO RINGS AND BRIDGES.
C
      NUCONN=0
      DO IAT=1,NAT
         IF(NIMODS(IAT).GT.NUCONN) THEN
            NUCONN=NIMODS(IAT)
         ENDIF
      ENDDO
C
      IF(NUCONN.LT.2) THEN
C
C     ----- WE SHOULD HAVE ONLY ( > 2 ) CONNECTIVITIES -----
C
         WRITE(IW,9989) NUCONN
         CALL HND_HNDERR(3,ERRMSG)
      ELSEIF(NUCONN.GT.3) THEN
C
C     ----- -NUCONN- > 3 .... NOT TREATED YET -----
C
         WRITE(IW,9989) NUCONN
         CALL HND_HNDERR(3,ERRMSG)
C
      ELSEIF(NUCONN.EQ.3) THEN
C
C     ----- -NUCONN- = 3 .... SEVERAL RINGS , FUSED OR NOT -----
C                             BRIDGES AS WELL.
C           FIRST , LIST OUT ALL THE CONNECTIVITY SEGMENTS
C
         IF(DBUG) THEN
            WRITE(IW,9984)
         ENDIF
C
         NSEG  =0
         IANCHR=0
  100    IANCHR=IANCHR+1
         IF(IANCHR.GT.NAT) THEN
            GO TO 120
         ENDIF
         IF(NIMODS(IANCHR).NE.3) THEN
            GO TO 100
         ENDIF
            DO I=1,3
               NSEG=NSEG+1
               IF(NSEG.GT.MXSEG) THEN
                  WRITE(IW,9983) NSEG
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
C
               JATOM=1
               IATSEG(JATOM,NSEG)=IANCHR
C
               IAT=IANCHR
               JAT=IJMODS(2,IIMODS(IAT)+I)
  110          JATOM=JATOM+1
               IATSEG(JATOM,NSEG)=JAT
               IF(NIMODS(JAT).EQ.2) THEN
                  IF(IJMODS(2,IIMODS(JAT)+1).NE.IAT) THEN
                     KAT=IJMODS(2,IIMODS(JAT)+1)
                  ELSE
                     KAT=IJMODS(2,IIMODS(JAT)+2)
                  ENDIF
                  IAT=JAT
                  JAT=KAT
                  GO TO 110
               ELSEIF(NIMODS(JAT).EQ.3) THEN
                  LENSEG(NSEG)=JATOM
                  IF(JAT.EQ.IANCHR) THEN
                     NAMSEG(NSEG)=1
                  ELSE
                     NAMSEG(NSEG)=2
                  ENDIF
               ELSE
                  WRITE(IW,9982) JAT
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
            ENDDO
C
  120    CONTINUE
         IF(DBUG) THEN
            DO MSEG=1,NSEG
               WRITE(IW,9981) MSEG,LENSEG(MSEG)
               WRITE(IW,9979) (IATSEG(IATOM,MSEG),IATOM=1,LENSEG(MSEG))
            ENDDO
         ENDIF
C
      ELSEIF(NUCONN.EQ.2) THEN
C
C     ----- -NUCONN- = 2 .... DISCONNECTED RINGS -----
C
         IF(DBUG) THEN
            WRITE(IW,9990)
         ENDIF
C
         NSEG  =0
         IANCHR=0
  200    IANCHR=IANCHR+1
         IF(IANCHR.GT.NAT) THEN
            GO TO 220
         ENDIF
         IF(NIMODS(IANCHR).NE.2) THEN
            GO TO 200
         ENDIF
            NSEG =NSEG+1
            IATOM=1
            IATSEG(IATOM,NSEG)=IANCHR
C
            IAT=IANCHR
            JAT=IJMODS(2,IIMODS(IAT)+1)
  210       IATOM=IATOM+1
            IATSEG(IATOM,NSEG)=JAT
            IF(NIMODS(JAT).EQ.2) THEN
               IF(IJMODS(2,IIMODS(JAT)+1).NE.IAT) THEN
                  KAT=IJMODS(2,IIMODS(JAT)+1)
               ELSE
                  KAT=IJMODS(2,IIMODS(JAT)+2)
               ENDIF
               IF(KAT.NE.IANCHR) THEN
                  IAT=JAT
                  JAT=KAT
                  GO TO 210
               ENDIF
            ELSE
               WRITE(IW,9986) JAT
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            LENSEG(NSEG)=IATOM
C
         IF(DBUG) THEN
            DO MSEG=1,NSEG
               WRITE(IW,9987) LENSEG(MSEG)
               WRITE(IW,9979) (IATSEG(IATOM,MSEG),IATOM=1,LENSEG(MSEG))
            ENDDO
         ENDIF
  220    CONTINUE
C
      ENDIF
C
      RETURN
 9999 FORMAT(' SOMETHING IS WRONG ... NO BONDS COULD BE DETECTED. STOP')
 9998 FORMAT(' TOTAL NUMBER OF BONDS = ',I5)
 9997 FORMAT(1X,'ATOM   NUM. OF BONDS   FIRST BOND  ENDATM',/,
     1       1X,41(1H-))
 9996 FORMAT(I5,6X,I5,5X,6X,I5,2X,4X,L1,3X)
 9995 FORMAT(' CONNECTIVITY FOR ATOM -IAT- = ',I4,' NUM.BONDS = ',I5)
 9994 FORMAT(6(1H(,I4,1H,,I4,1H)))
 9993 FORMAT(' NO COVALENT RADIUS SPECIFIED FOR ONE OF THESE ATOMS = ',
     1       2I5)
 9992 FORMAT(' IAT,JAT,IZ,JZ,RIJ,RCV,MBNDS= ',2I4,2I4,2F10.3,3I5,L4)
 9991 FORMAT(/,' PASS NUM. = ',I3,' OUT OF ',I3,' PASS(ES) ')
 9990 FORMAT(/,' SEARCHING FOR A SINGLE RING ...')
 9989 FORMAT(/,' ILLEGAL -MXCONN- HERE, -MXCONN- = ',I3,' . STOP.')
 9988 FORMAT(/,' THIS ROUTINE DOES NOT HANDLE -MXCONN- .GT. 3. STOP.')
 9987 FORMAT(/,' A ',I2,'-MEMBER RING WAS FOUND. ATOMS ARE = ',/,20I4)
 9986 FORMAT(/,' ATOM ',I3,' SHOULD HAVE TWO CONNECTIONS,',
     1         ' AND ONLY TWO. STOP.')
 9985 FORMAT(/,' ERROR DURING COMPRESSION OF CONNECTIVITY TABLE',
     1         ' FOR ATOM = ',I4,/,' NIMODS(IAT), NUMODS = ',2I5)
 9984 FORMAT(/,' SEARCHING FOR MULTIPLE RINGS AND BRIDGES ...')
 9983 FORMAT(/,' TOO MANY SEGMENTS FROM 3-CONNECTED ATOMS.',
     1         ' NSEG, MXSEG = ',2I4)
 9982 FORMAT(/,' ATOM ',I3,' SHOULD HAVE AT LEAST TWO CONNECTIONS.',
     1         ' STOP.')
 9981 FORMAT(/,' SEARCHING SEGMENT NO. = ',I4)
 9980 FORMAT(/,' ALL 1-CONNECTIVITY COMPLETED.')
 9979 FORMAT(16I4)
 9978 FORMAT(/,' ATOM ',I4,' HAS NO CONNECTIVITY LEFT,',/,
     1         ' BUT THE CODE IS TRYING TO REMOVE CONNECTION',
     2         ' TO ATOM ',I4,' . STOP.')
 9977 FORMAT(' INCORRECT DATA IN NAMELIST -COORDI- . STOP .')
 9976 FORMAT(2I5,3F10.6)
 9975 FORMAT(' FORCED BOND INSERTED = ',3I5)
 9974 FORMAT(' BOND REMOVED = ',I5,' EQUIVALENT TO = ',I5)
 8888 FORMAT(/,10X,6(1H-),/,10X,'auto-z',/,10X,6(1H-))
 8887 FORMAT(' SOME IMPOSED BONDS           IN -ZMTMAK- ')
 8886 FORMAT('   NO IMPOSED BONDS           IN -ZMTMAK- ')
 8885 FORMAT(' SOME IMPOSED ANGLES          IN -ZMTMAK- ')
 8884 FORMAT('   NO IMPOSED ANGLES          IN -ZMTMAK- ')
 8883 FORMAT(' SOME IMPOSED TORSIONS        IN -ZMTMAK- ')
 8882 FORMAT('   NO IMPOSED TORSIONS        IN -ZMTMAK- ')
 8881 FORMAT(' SOME IMPOSED O-O-PLANE-BENDS IN -ZMTMAK- ')
 8880 FORMAT('   NO IMPOSED O-O-PLANE-BENDS IN -ZMTMAK- ')
 8879 FORMAT(' SOME IMPOSED LINEAR-BENDS    IN -ZMTMAK- ')
 8878 FORMAT('   NO IMPOSED LINEAR-BENDS    IN -ZMTMAK- ')
 8877 FORMAT(' zdone = ',L4)
      END
      SUBROUTINE HND_ZMTBLD(ZDONE,NZVAR,NIZMAT,IZMAT,C,
     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,NUMBND,
     2                      IJKANG,NUMANG,IJKLTO,NUMTOR,
     3                      IJKLOP,NUMOOP,
     4                      IJKLNB,NUMLNB,XYZLNB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXANGS= 8*MXATOM)
      PARAMETER (MXTORS= 8*MXATOM)
      PARAMETER (MXOOPA= 8*MXATOM)
      PARAMETER (MXLINB= 8*MXATOM)
      LOGICAL     DBUG
      LOGICAL     OUT
      LOGICAL     ZDONE
      LOGICAL     DOOOPA
      LOGICAL     DOLINB
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,*)
      DIMENSION IZMAT(*)   
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION IJKANG(3,*),IJKLTO(4,*),IJKLOP(4,*)
      DIMENSION IJKLNB(4,*),XYZLNB(3,*)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTBLD-'/
C
      DBUG=.FALSE.
      OUT =.FALSE.
      OUT =OUT.OR.DBUG
      IF(DBUG) THEN
         WRITE(IW,9999)
      ENDIF
C
      IF(NAT.LE.1) THEN
         ZDONE=.TRUE.
         RETURN
      ENDIF
C
C     ----- BONDS -----
C
      CALL HND_ZMTYP1(NZVAR,NIZMAT,IZMAT,
     1                NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,NUMBND)
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,9995) NZVAR
      ENDIF
C
C     ----- SOME CHECKS FOR DIATOMICS AND TRIATOMICS -----
C
      IF(NAT.EQ.2) THEN
         IF(NZVAR.EQ.1) THEN
            ZDONE=.TRUE. 
            RETURN
         ELSE
            NZVAR=0
            ZDONE=.FALSE.
            RETURN
         ENDIF
      ENDIF
      IF(NAT.EQ.3) THEN
         IF(NZVAR.NE.2.AND.NZVAR.NE.3) THEN
            NZVAR=0
            ZDONE=.FALSE.
            RETURN
         ENDIF
      ENDIF
C
C     ----- ANGLES -----
C
      DOLINB=.FALSE.
      CALL HND_ZMTYP2(NZVAR,NIZMAT,IZMAT,C,
     1                NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,IJKANG,NUMANG,
     2                DOLINB)
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,9994) NZVAR
      ENDIF
C
C     ----- MORE CHECKS FOR TRIATOMICS -----
C
      IF(NAT.EQ.3) THEN
         IF(DOLINB) THEN
            NZVAR=0
            ZDONE=.FALSE.
            RETURN
         ENDIF
      ENDIF
C
C     ----- TORSIONS -----
C
      CALL HND_ZMTYP3(NZVAR,NIZMAT,IZMAT,C,
     1                NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,IJKLTO,NUMTOR)
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,9993) NZVAR
      ENDIF
C
      DOOOPA=.FALSE.
      IF(DOOOPA) THEN
C
C     ----- OUT-OF-PLANE BENDS -----
C
         CALL HND_ZMTYP4(NZVAR,NIZMAT,IZMAT,
     1                   NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,IJKLOP,NUMOOP)
         IF(DBUG) THEN
            WRITE(IW,9998)
            WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
            WRITE(IW,9996)
            WRITE(IW,9992) NZVAR
         ENDIF
      ENDIF
C
C     ----- LINEAR BENDS -----
C
      IF(DOLINB) THEN
         CALL HND_ZMTYP5(NZVAR,NIZMAT,IZMAT,C,
     1                   NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,
     2                   IJKLNB,NUMLNB,XYZLNB,IJKANG,NUMANG)
         IF(DBUG) THEN
            WRITE(IW,9998)
            WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
            WRITE(IW,9996)
            WRITE(IW,9991) NZVAR
         ENDIF
         ZDONE=.FALSE.
         RETURN
      ENDIF
C
C     ----- DONE -----
C
      IF(NZVAR.LT.3*NAT-6) THEN
         ZDONE=.FALSE.
      ELSE
         ZDONE=.TRUE.
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTBLD- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
 9995 FORMAT(' AFTER -ZMTYP1- , NZVAR = ',I5)
 9994 FORMAT(' AFTER -ZMTYP2- , NZVAR = ',I5)
 9993 FORMAT(' AFTER -ZMTYP3- , NZVAR = ',I5)
 9992 FORMAT(' AFTER -ZMTYP4- , NZVAR = ',I5)
 9991 FORMAT(' AFTER -ZMTYP5- , NZVAR = ',I5)
      END
      SUBROUTINE HND_ZMTYP1(NZVAR,NIZMAT,IZMAT,
     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,NUMBND)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IONE=1)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      LOGICAL     DBUG
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION IZMAT(*)   
      DIMENSION IJBOND(2,MXBNDS)
      DIMENSION IIBOND(  MXATOM)
      DIMENSION NIBOND(  MXATOM)
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP1-'/
C
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999)
      ENDIF
C
      NBOND=NBNDS
      DO IJBND=1,NBOND
         IJBOND(1,IJBND)=IJBNDS(1,IJBND)
         IJBOND(2,IJBND)=IJBNDS(2,IJBND)
      ENDDO
      DO IAT=1,NAT
         IIBOND(IAT)=IIBNDS(IAT)
         NIBOND(IAT)=NIBNDS(IAT)
      ENDDO
C
      IF(NAT.LE.1) RETURN
C
C     ----- BONDS -----
C
      NBND=NBOND
      DO IBND=1,NBND
         II=IJBOND(1,IBND)
         IJ=IJBOND(2,IBND)
         IF(II.NE.0) THEN
            IF(IBND.LT.NBND) THEN
               DO JBND=IBND+1,NBND
                  JI=IJBOND(1,JBND)
                  JJ=IJBOND(2,JBND)
                  IF(JI.NE.0) THEN
                     IF((JI.EQ.II.AND.JJ.EQ.IJ).OR.        
     1                  (JI.EQ.IJ.AND.JJ.EQ.II)    ) THEN
                        IJBOND(1,JBND)=0
                        IJBOND(2,JBND)=0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
             IZMAT(NIZMAT+1)=IONE
             IZMAT(NIZMAT+2)=II
             IZMAT(NIZMAT+3)=IJ   
            NIZMAT=NIZMAT+3
            NZVAR =NZVAR +1
            NUMBND=NUMBND+1
         ENDIF
      ENDDO
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,*)    'NZVAR = ',NZVAR
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTYP1- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
      END
      SUBROUTINE HND_ZMTYP2(NZVAR,NIZMAT,IZMAT,C,
     1           NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,IJKANG,NUMANG,
     2           DOLINB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ITWO=2)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXANGS= 8*MXATOM)
      LOGICAL     DBUG
      LOGICAL     DOLINB
      LOGICAL     LINBND
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,*)
      DIMENSION IZMAT(*)   
      DIMENSION IJKANG(3,*)     
      DIMENSION IJBOND(2,MXBNDS)
      DIMENSION IIBOND(  MXATOM)
      DIMENSION NIBOND(  MXATOM)
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP2-'/
      DATA TOL    /1.0D-08/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
C
      DISTSQ(IAT,JAT)=(C(1,IAT)-C(1,JAT))**2+
     1                (C(2,IAT)-C(2,JAT))**2+
     2                (C(3,IAT)-C(3,JAT))**2
      DOTPRD(IAT,JAT,KAT)=(C(1,IAT)-C(1,JAT))*(C(1,KAT)-C(1,JAT))+
     1                    (C(2,IAT)-C(2,JAT))*(C(2,KAT)-C(2,JAT))+
     2                    (C(3,IAT)-C(3,JAT))*(C(3,KAT)-C(3,JAT))
C
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999)
      ENDIF
C
C     ----- ANGLES -----
C
      NBOND=NBNDS
      DO IJBND=1,NBOND
         IJBOND(1,IJBND)=IJBNDS(1,IJBND)
         IJBOND(2,IJBND)=IJBNDS(2,IJBND)
      ENDDO
      DO IAT=1,NAT
         IIBOND(IAT)=IIBNDS(IAT)
         NIBOND(IAT)=NIBNDS(IAT)
      ENDDO
C
      IF(NAT.LE.2) RETURN
C
      NANG=0
      DO IAT=1,NAT
         DO IJBND=IIBOND(IAT)+1,IIBOND(IAT)+NIBOND(IAT)
            ICON=IJBOND(1,IJBND)
            JCON=IJBOND(2,IJBND)
            IF(ICON.NE.IAT) THEN
               WRITE(IW,9995) IJBND,IAT,ICON,JCON
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            JAT=JCON
            IF(NIBOND(JAT).GT.1) THEN
               DO JKBND=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
                  JCON=IJBOND(1,JKBND)
                  KCON=IJBOND(2,JKBND)
                  IF(JCON.NE.JAT) THEN
                     WRITE(IW,9995) JKBND,JAT,JCON,KCON
                     CALL HND_HNDERR(3,ERRMSG)
                  ENDIF
                  KAT=KCON
                  IF(KAT.NE.IAT) THEN
                     IF(NANG.EQ.MXANGS) THEN
                        WRITE(IW,*) 'TOO MANY BOND ANGLES',
     1                              ' NANG,MXANG = ',NANG,MXANGS
                        CALL HND_HNDERR(3,ERRMSG)
                     ENDIF
                     NANG=NANG+1
                     IJKANG(1,NANG)=IAT
                     IJKANG(2,NANG)=JAT
                     IJKANG(3,NANG)=KAT
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      IF(DBUG) THEN
         IF(NANG.GT.1) THEN
            WRITE(IW,9994)
            DO MANG=1,NANG
               WRITE(IW,9993) MANG,IJKANG(1,MANG),
     1                             IJKANG(2,MANG),IJKANG(3,MANG)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
C     ----- ELIMINATE DUPLICATES AND CLEAN UP -----
C
      IF(NANG.GE.2) THEN
         DO IANG=1,NANG-1
            II=IJKANG(1,IANG)
            IJ=IJKANG(2,IANG)
            IK=IJKANG(3,IANG)
            IF(IJKANG(1,IANG).NE.0) THEN
               DO JANG=IANG+1,NANG
                  JI=IJKANG(1,JANG)
                  JJ=IJKANG(2,JANG)
                  JK=IJKANG(3,JANG)
                  IF(JJ.EQ.IJ) THEN
                     IF((II.EQ.JI.AND.IK.EQ.JK).OR.
     1                  (II.EQ.JK.AND.IK.EQ.JI)    ) THEN
                        IJKANG(1,JANG)=0
                        IJKANG(2,JANG)=0
                        IJKANG(3,JANG)=0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
      IANG=0
  100 IANG=IANG+1
      IF(IANG.GT.NANG) GO TO 110
         IF(IJKANG(1,IANG).EQ.0) THEN
            IF(IANG.LT.NANG) THEN
               DO JANG=IANG+1,NANG
                  IJKANG(1,JANG-1)=IJKANG(1,JANG)
                  IJKANG(2,JANG-1)=IJKANG(2,JANG)
                  IJKANG(3,JANG-1)=IJKANG(3,JANG)
               ENDDO
               IANG=IANG-1
            ENDIF
            NANG=NANG-1
         ELSE
            RIJSQ=DISTSQ(IJKANG(1,IANG),IJKANG(2,IANG))
            RIJ  =SQRT(RIJSQ)                          
            RJKSQ=DISTSQ(IJKANG(2,IANG),IJKANG(3,IANG))
            RJK  =SQRT(RJKSQ)                           
            RIKSQ=DISTSQ(IJKANG(1,IANG),IJKANG(3,IANG))
            COSB =(RIJSQ+RJKSQ-RIKSQ)/(TWO*RIJ*RJK)
            DUM  =ABS(COSB)-ONE
            IF(DUM.GT.ZERO) THEN
               IF(DUM.LE.TOL) THEN
                  LINBND=.TRUE.
                  DOLINB=DOLINB.OR.LINBND
               ELSE
                  WRITE(IW,*) 'SOMETHING IS WRONG IN -ZMTYP2-'
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
            ELSE
               DUM=-DUM
               IF(DUM.LE.TOL) THEN
                   LINBND=.TRUE.
                   DOLINB=DOLINB.OR.LINBND
               ELSE
                  LINBND=.FALSE.
               ENDIF
            ENDIF
            IF(.NOT.LINBND) THEN
                IZMAT(NIZMAT+1)=ITWO
                IZMAT(NIZMAT+2)=IJKANG(1,IANG)
                IZMAT(NIZMAT+3)=IJKANG(2,IANG)
                IZMAT(NIZMAT+4)=IJKANG(3,IANG)
               NIZMAT=NIZMAT+4
               NZVAR =NZVAR +1
               NUMANG=NUMANG+1
            ENDIF
         ENDIF
         GO TO 100
  110 CONTINUE
C
      IF(DBUG) THEN
         IF(NANG.GT.1) THEN
            WRITE(IW,9994)
            DO MANG=1,NANG
               WRITE(IW,9993) MANG,IJKANG(1,MANG),
     1                             IJKANG(2,MANG),IJKANG(3,MANG)
            ENDDO
            IF(DOLINB) THEN
               WRITE(IW,9991)
            ENDIF
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,*)    'NZVAR = ',NZVAR
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTYP2- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
 9995 FORMAT(' SOMETHING WRONG WITH BOND ',I2,' OF ATOM ',I4,
     1       ' WITH -ICON-, -JCON- = ',2I5)
 9994 FORMAT(' TABLE OF ANGLES = ')
 9993 FORMAT(4I5)
 9992 FORMAT(' NO BOND ANGLES FOUND ')
 9991 FORMAT(' COLLINEAR ATOMS FOUND. DO -LINEAR BENDS- ... ')
      END
      SUBROUTINE HND_ZMTYP3(NZVAR,NIZMAT,IZMAT,C,
     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,
     2                                       IJKLTO,NUMTOR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ITHREE=3)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXTORS= 8*MXATOM)
      LOGICAL     DBUG
      LOGICAL     LINIJK,LINJKL
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,*)
      DIMENSION IZMAT(*)   
      DIMENSION IJKLTO(4,*)      
      DIMENSION IJBOND(2,MXBNDS)
      DIMENSION IIBOND(  MXATOM)
      DIMENSION NIBOND(  MXATOM)
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP3-'/
      DATA TOL    /1.0D-08/
      DATA ZERO   /0.0D+00/
      DATA ONE    /1.0D+00/
      DATA TWO    /2.0D+00/
C
      DISTSQ(IAT,JAT)=(C(1,IAT)-C(1,JAT))**2+
     1                (C(2,IAT)-C(2,JAT))**2+
     2                (C(3,IAT)-C(3,JAT))**2
      DOTPRD(IAT,JAT,KAT)=(C(1,IAT)-C(1,JAT))*(C(1,KAT)-C(1,JAT))+
     1                    (C(2,IAT)-C(2,JAT))*(C(2,KAT)-C(2,JAT))+
     2                    (C(3,IAT)-C(3,JAT))*(C(3,KAT)-C(3,JAT))
C
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999)
      ENDIF
C
C     ----- TORSIONS -----
C
      NBOND=NBNDS
      DO IJBND=1,NBOND
         IJBOND(1,IJBND)=IJBNDS(1,IJBND)
         IJBOND(2,IJBND)=IJBNDS(2,IJBND)
      ENDDO
      DO IAT=1,NAT
         IIBOND(IAT)=IIBNDS(IAT)
         NIBOND(IAT)=NIBNDS(IAT)
      ENDDO
C
      IF(NAT.LE.3) RETURN
C
      NTOR=0
      DO IAT=1,NAT
         DO IJBND=IIBOND(IAT)+1,IIBOND(IAT)+NIBOND(IAT)
            ICON=IJBOND(1,IJBND)
            JCON=IJBOND(2,IJBND)
            IF(ICON.NE.IAT) THEN
               WRITE(IW,9995) IJBND,IAT,ICON,JCON
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            JAT=JCON
            IF(NIBOND(JAT).GT.1) THEN
               DO JKBND=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
                  JCON=IJBOND(1,JKBND)
                  KCON=IJBOND(2,JKBND)
                  IF(JCON.NE.JAT) THEN
                     WRITE(IW,9995) JKBND,JAT,JCON,KCON
                     CALL HND_HNDERR(3,ERRMSG)
                  ENDIF
                  KAT=KCON
                  IF(KAT.NE.IAT) THEN
                     IF(NIBOND(KAT).GT.1) THEN
                        DO KLBND=IIBOND(KAT)+1,IIBOND(KAT)+NIBOND(KAT)
                           KCON=IJBOND(1,KLBND)
                           LCON=IJBOND(2,KLBND)
                           IF(KCON.NE.KAT) THEN
                              WRITE(IW,9995) KLBND,KAT,KCON,LCON
                              CALL HND_HNDERR(3,ERRMSG)
                           ENDIF
                           LAT=LCON
                           IF(LAT.NE.JAT) THEN
                              IF(NTOR.EQ.MXTORS) THEN
                                 WRITE(IW,*) 'TOO MANY TORSIONS',
     1                                       ' NTOR, MXTORS = ',NTOR,
     2                                       MXTORS
                                 CALL HND_HNDERR(3,ERRMSG)
                              ENDIF
                              NTOR=NTOR+1
                              IJKLTO(1,NTOR)=IAT
                              IJKLTO(2,NTOR)=JAT
                              IJKLTO(3,NTOR)=KAT
                              IJKLTO(4,NTOR)=LAT
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      IF(DBUG) THEN
         IF(NTOR.GT.1) THEN
            WRITE(IW,9994)
            DO MTOR=1,NTOR
               WRITE(IW,9993) MTOR,IJKLTO(1,MTOR),IJKLTO(2,MTOR),
     1                             IJKLTO(3,MTOR),IJKLTO(4,MTOR)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
C     ----- ELIMINATE DUPLICATES AND CLEAN UP -----
C
      IF(NTOR.GE.2) THEN
         DO ITOR=1,NTOR-1
            II=IJKLTO(1,ITOR)
            IJ=IJKLTO(2,ITOR)
            IK=IJKLTO(3,ITOR)
            IL=IJKLTO(4,ITOR) 
            IF(IJKLTO(1,ITOR).NE.0) THEN
               DO JTOR=ITOR+1,NTOR
                  JI=IJKLTO(1,JTOR)
                  JJ=IJKLTO(2,JTOR)
                  JK=IJKLTO(3,JTOR)
                  JL=IJKLTO(4,JTOR)
                  IF(JJ.EQ.IJ.AND.JK.EQ.IK) THEN
                     IF(II.EQ.JI.AND.IL.EQ.JL) THEN
                        IJKLTO(1,JTOR)=0
                        IJKLTO(2,JTOR)=0
                        IJKLTO(3,JTOR)=0
                        IJKLTO(4,JTOR)=0
                     ENDIF
                  ELSEIF(JJ.EQ.IK.AND.JK.EQ.IJ) THEN
                     IF(II.EQ.JL.AND.IL.EQ.JI) THEN
                        IJKLTO(1,JTOR)=0
                        IJKLTO(2,JTOR)=0
                        IJKLTO(3,JTOR)=0
                        IJKLTO(4,JTOR)=0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
      ITOR=0
  100 ITOR=ITOR+1
      IF(ITOR.GT.NTOR) GO TO 110
         IF(IJKLTO(1,ITOR).EQ.0) THEN
            IF(ITOR.LT.NTOR) THEN
               DO JTOR=ITOR+1,NTOR
                  IJKLTO(1,JTOR-1)=IJKLTO(1,JTOR)
                  IJKLTO(2,JTOR-1)=IJKLTO(2,JTOR)
                  IJKLTO(3,JTOR-1)=IJKLTO(3,JTOR)
                  IJKLTO(4,JTOR-1)=IJKLTO(4,JTOR)
               ENDDO
               ITOR=ITOR-1
            ENDIF
            NTOR=NTOR-1
         ELSE
            RIJSQ=DISTSQ(IJKLTO(1,ITOR),IJKLTO(2,ITOR))
            RIJ  =SQRT(RIJSQ)
            RJKSQ=DISTSQ(IJKLTO(2,ITOR),IJKLTO(3,ITOR))
            RJK  =SQRT(RJKSQ)
            RKLSQ=DISTSQ(IJKLTO(3,ITOR),IJKLTO(4,ITOR))
            RKL  =SQRT(RKLSQ)
            RIKSQ=DISTSQ(IJKLTO(1,ITOR),IJKLTO(3,ITOR))
            RJLSQ=DISTSQ(IJKLTO(2,ITOR),IJKLTO(4,ITOR))
            COSBJ=(RIJSQ+RJKSQ-RIKSQ)/(TWO*RIJ*RJK)
            COSBK=(RJKSQ+RKLSQ-RJLSQ)/(TWO*RJK*RKL)
            DUM  =ABS(COSBJ)-ONE
            IF(DUM.GT.ZERO) THEN
               IF(DUM.LE.TOL) THEN
                  LINIJK=.TRUE.
               ELSE
                  WRITE(IW,*) 'SOMETHING IS WRONG IN -ZMTYP2-'
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
            ELSE
               DUM=-DUM
               IF(DUM.LE.TOL) THEN
                   LINIJK=.TRUE.
               ELSE
                  LINIJK=.FALSE.
               ENDIF
            ENDIF
            DUM  =ABS(COSBK)-ONE
            IF(DUM.GT.ZERO) THEN
               IF(DUM.LE.TOL) THEN
                  LINJKL=.TRUE.
               ELSE
                  WRITE(IW,*) 'SOMETHING IS WRONG IN -ZMTYP2-'
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
            ELSE
               DUM=-DUM
               IF(DUM.LE.TOL) THEN
                   LINJKL=.TRUE.
               ELSE
                  LINJKL=.FALSE.
               ENDIF
            ENDIF
            IF(.NOT.LINIJK.AND..NOT.LINJKL) THEN
                IZMAT(NIZMAT+1)=ITHREE
                IZMAT(NIZMAT+2)=IJKLTO(1,ITOR)
                IZMAT(NIZMAT+3)=IJKLTO(2,ITOR)
                IZMAT(NIZMAT+4)=IJKLTO(3,ITOR)
                IZMAT(NIZMAT+5)=IJKLTO(4,ITOR)
               NIZMAT=NIZMAT+5
               NZVAR =NZVAR +1
            ENDIF
         ENDIF
         GO TO 100
  110 CONTINUE
C
      IF(DBUG) THEN
         IF(NTOR.GT.1) THEN
            WRITE(IW,9994)
            DO MTOR=1,NTOR
               WRITE(IW,9993) MTOR,IJKLTO(1,MTOR),IJKLTO(2,MTOR),
     1                             IJKLTO(3,MTOR),IJKLTO(4,MTOR)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,*)    'NZVAR = ',NZVAR
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTYP3- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
 9995 FORMAT(' SOMETHING WRONG WITH BOND ',I2,' OF ATOM ',I4,
     1       ' WITH -ICON-, -JCON- = ',2I5)
 9994 FORMAT(' TABLE OF TORSIONS = ')
 9993 FORMAT(5I5)
 9992 FORMAT(' NO TORSIONS FOUND ')
      END
      SUBROUTINE HND_ZMTYP4(NZVAR,NIZMAT,IZMAT,
     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,
     2                                       IJKLOP,NUMOOP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IFOUR=4)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXOOPA= 8*MXATOM)
      LOGICAL     DBUG
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION IZMAT(*)   
      DIMENSION IJKLOP(4,*)      
      DIMENSION IJBOND(2,MXBNDS)
      DIMENSION IIBOND(  MXATOM)
      DIMENSION NIBOND(  MXATOM)
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP4-'/
C
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999)
      ENDIF
C
C     ----- OUT-OF-PLANE-BENDS -----
C
      NBOND=NBNDS
      DO IJBND=1,NBOND
         IJBOND(1,IJBND)=IJBNDS(1,IJBND)
         IJBOND(2,IJBND)=IJBNDS(2,IJBND)
      ENDDO
      DO IAT=1,NAT
         IIBOND(IAT)=IIBNDS(IAT)
         NIBOND(IAT)=NIBNDS(IAT)
      ENDDO
C
      IF(NAT.LE.3) RETURN
C
      NOOP=0
      DO IAT=1,NAT
         DO IJBND=IIBOND(IAT)+1,IIBOND(IAT)+NIBOND(IAT)
            ICON=IJBOND(1,IJBND)
            JCON=IJBOND(2,IJBND)
            IF(ICON.NE.IAT) THEN
               WRITE(IW,9995) IJBND,IAT,ICON,JCON
               CALL HND_HNDERR(3,ERRMSG)
            ENDIF
            JAT=JCON
            IF(NIBOND(JAT).GT.2) THEN
               DO JKBND=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
                  JCON=IJBOND(1,JKBND)
                  KCON=IJBOND(2,JKBND)
                  IF(JCON.NE.JAT) THEN
                     WRITE(IW,9995) JKBND,JAT,JCON,KCON
                     CALL HND_HNDERR(3,ERRMSG)
                  ENDIF
                  KAT=KCON
                  IF(KAT.NE.IAT) THEN
                     DO JLBND=IIBOND(JAT)+1,IIBOND(JAT)+NIBOND(JAT)
                        JCON=IJBOND(1,JLBND)
                        LCON=IJBOND(2,JLBND)
                        IF(JCON.NE.JAT) THEN
                           WRITE(IW,9995) JLBND,JAT,JCON,LCON
                           CALL HND_HNDERR(3,ERRMSG)
                        ENDIF
                        LAT=LCON
                        IF(LAT.NE.IAT.AND.LAT.NE.KAT) THEN
                           NKCON=NIBOND(KAT)
                           NLCON=NIBOND(LAT)
                           IF((NKCON.EQ.1.AND.NLCON.EQ.1).OR.
     1                        (NKCON.GT.1.AND.NLCON.GT.1)    ) THEN
                              IF(NOOP.EQ.MXOOPA) THEN
                              WRITE(IW,*) 'TOO MANY OUT-OF-PLANE BENDS',
     1                                    ' NOOP, MXOOPA = ',NOOP,MXOOPA
                                 CALL HND_HNDERR(3,ERRMSG)
                              ENDIF
                              NOOP=NOOP+1
                              IJKLOP(1,NOOP)=IAT
                              IJKLOP(2,NOOP)=JAT
                              IJKLOP(3,NOOP)=KAT
                              IJKLOP(4,NOOP)=LAT
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      IF(DBUG) THEN
         IF(NOOP.GT.1) THEN
            WRITE(IW,9994)
            DO MOOP=1,NOOP
               WRITE(IW,9993) MOOP,IJKLOP(1,MOOP),IJKLOP(2,MOOP),
     1                             IJKLOP(3,MOOP),IJKLOP(4,MOOP)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
C     ----- ELIMINATE DUPLICATES AND CLEAN UP -----
C
      IF(NOOP.GE.2) THEN
         DO IOOP=1,NOOP-1
            II=IJKLOP(1,IOOP)
            IJ=IJKLOP(2,IOOP)
            IK=IJKLOP(3,IOOP)
            IL=IJKLOP(4,IOOP) 
            IF(IJKLOP(1,IOOP).NE.0) THEN
               DO JOOP=IOOP+1,NOOP
                  JI=IJKLOP(1,JOOP)
                  JJ=IJKLOP(2,JOOP)
                  JK=IJKLOP(3,JOOP)
                  JL=IJKLOP(4,JOOP)
                  IF(JI.EQ.II.AND.JJ.EQ.IJ) THEN
                     IF((JK.EQ.IK.AND.JL.EQ.IL).OR.
     1                  (JK.EQ.IL.AND.JL.EQ.IK)    ) THEN
                        IJKLOP(1,JOOP)=0
                        IJKLOP(2,JOOP)=0
                        IJKLOP(3,JOOP)=0
                        IJKLOP(4,JOOP)=0
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
C
      IOOP=0
  100 IOOP=IOOP+1
      IF(IOOP.GT.NOOP) GO TO 110
         IF(IJKLOP(1,IOOP).EQ.0) THEN
            IF(IOOP.LT.NOOP) THEN
               DO JOOP=IOOP+1,NOOP
                  IJKLOP(1,JOOP-1)=IJKLOP(1,JOOP)
                  IJKLOP(2,JOOP-1)=IJKLOP(2,JOOP)
                  IJKLOP(3,JOOP-1)=IJKLOP(3,JOOP)
                  IJKLOP(4,JOOP-1)=IJKLOP(4,JOOP)
               ENDDO
               IOOP=IOOP-1
            ENDIF
            NOOP=NOOP-1
         ELSE
             IZMAT(NIZMAT+1)=IFOUR 
             IZMAT(NIZMAT+2)=IJKLOP(1,IOOP)
             IZMAT(NIZMAT+3)=IJKLOP(2,IOOP)
             IZMAT(NIZMAT+4)=IJKLOP(3,IOOP)
             IZMAT(NIZMAT+5)=IJKLOP(4,IOOP)
            NIZMAT=NIZMAT+5
            NZVAR =NZVAR +1
         ENDIF
         GO TO 100
  110 CONTINUE
C
      IF(DBUG) THEN
         IF(NOOP.GT.1) THEN
            WRITE(IW,9994)
            DO MOOP=1,NOOP
               WRITE(IW,9993) MOOP,IJKLOP(1,MOOP),IJKLOP(2,MOOP),
     1                             IJKLOP(3,MOOP),IJKLOP(4,MOOP)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
         WRITE(IW,*)    'NZVAR = ',NZVAR
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTYP4- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
 9995 FORMAT(' SOMETHING WRONG WITH BOND ',I2,' OF ATOM ',I4,
     1       ' WITH -ICON-, -JCON- = ',2I5)
 9994 FORMAT(' TABLE OF OUT-OF-PLANE-BENDS = ')
 9993 FORMAT(5I5)
 9992 FORMAT(' NO OUT-OF-PLANE-ANGLE FOUND ')
      END
      SUBROUTINE HND_ZMTYP5(NZVAR,NIZMAT,IZMAT,C,
     1                      NAT,NBNDS,IJBNDS,IIBNDS,NIBNDS,
     2                      IJKLNB,NUMLNB,XYZLNB,IJKANG,NUMANG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (IFOUR=4)
      PARAMETER (MXATOM=500)        
      PARAMETER (MXBNDS= 8*MXATOM)
      PARAMETER (MXOOPA= 8*MXATOM)
      LOGICAL     DBUG
      CHARACTER*8 ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION C(3,*)
      DIMENSION IZMAT(*)   
      DIMENSION IJKLNB(4,*)      
      DIMENSION XYZLNB(3,*)
      DIMENSION IJKANG(3,*)           
      DIMENSION IJBNDS(2,*)
      DIMENSION IIBNDS(  *)
      DIMENSION NIBNDS(  *)
      DIMENSION ERRMSG(3)        
      DATA ERRMSG /'PROGRAM ','STOP IN ','-ZMTYP5-'/
C
      DBUG=.FALSE.
      IF(DBUG) THEN
         WRITE(IW,9999)
         WRITE(IW,9991)
      ENDIF
C
C     ----- LINEAR-BENDS -----
C
      NLBND=0
C
      IF(NAT.LE.3) RETURN
C
      IF(DBUG) THEN
         IF(NLBND.GT.1) THEN
            WRITE(IW,9994)
            DO MLBND=1,NLBND
               WRITE(IW,9993) MLBND,IJKLNB(1,MLBND),IJKLNB(2,MLBND),
     1                              IJKLNB(3,MLBND),IJKLNB(4,MLBND)
            ENDDO
         ELSE
            WRITE(IW,9992)
         ENDIF
      ENDIF
C
      IF(DBUG) THEN
         WRITE(IW,9998)
         WRITE(IW,9997) (IZMAT(I),I=1,NIZMAT)
         WRITE(IW,9996)
      ENDIF
C
      RETURN
 9999 FORMAT(' IN -ZMTYP5- ',/,' ----------- ')
 9998 FORMAT(' $ZMAT  ',/,' IZMAT= ')
 9997 FORMAT(12(I4,1H,))
 9996 FORMAT(' $END   ')
 9994 FORMAT(' TABLE OF LINEAR-BENDS = ')
 9993 FORMAT(5I5)
 9992 FORMAT(' NO LINEAR-BENDS FOUND ')
 9991 FORMAT(' LINEAR-BENDS NOT SET UP CURRENTLY . STOP . ')
      END
      SUBROUTINE GEOM_ZMT_GEO(COORDS,CHARGE,TAGS,NCENTER,
     1                        IZ,IZMAT,NZMOD,ICFRZ,NCFRZ,
     2                        UNITS)
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
     3             IZMAT,IZ,IZFRZ,MXIZMT,NZMOD,DBUG)
      ENDDO
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
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
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
#include "tcgmsg.fh"
#include "msgtypesf.h"
#include "global.fh"
#include "rtdb.fh"
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
      subroutine hnd_dparsc(a,la,c,lc)
      implicit double precision (a-h,o-z)
      character*(*) a
      character*(*) c
      character*1   blk
      data  blk /' '/
      do ic=1,lc
         c(ic:ic)=blk
      enddo
      if(lc.le.(la-2)) then
         c(1:lc)=a(2:lc+1)
      else
         c(1:lc)=a(2:la-1)
      endif
      return
      end
      subroutine hnd_dparsi(a,la,n)
      implicit double precision (a-h,o-z)
      character*(*) a
      character*1   char(12)
      data char     /'0','1','2','3','4','5','6','7','8','9',
     1               '+','-'/
c
      i1=1
      i2=la
c
c     ----- sign -----
c
      isign=1
      if(a(i1:i1).eq.char(12)) then
         isign=-1
         i1=i1+1
      elseif(a(i1:i1).eq.char(11)) then
         isign= 1
         i1=i1+1
      else
         isign= 1
      endif
      na=i2-i1+1
c
      n=0
      do ia=i1,i2
         ib=i2-ia
         do i=1,10
            if(a(ia:ia).eq.char(i)) then
               n=n+(i-1)*10**(ib)
            endif
         enddo
      enddo
      n=n*isign
c
      return
      end
      subroutine hnd_dparsr(a,la,x)
      implicit double precision (a-h,o-z)
      logical rep
      character*(*) a
      character*1   char(17)
      data char     /'0','1','2','3','4','5','6','7','8','9',
     1               '+','-','.','e','d','E','D'/
      data zero     /0.0d+00/
      data ten      /1.0d+01/
c     
      i1 =1            
      i2 =la
c
c     ----- sign -----
c
      isign=1
      if(a(i1:i1).eq.char(12)) then
         isign=-1
         i1=i1+1
      elseif(a(i1:i1).eq.char(11)) then
         isign= 1
         i1=i1+1
      else
         isign= 1
      endif
c
c     ----- exponent -----
c
      ie2=i2
      do ie=i1,i2
         if(a(ie:ie).eq.char(14).or.a(ie:ie).eq.char(15).or. 
     1      a(ie:ie).eq.char(16).or.a(ie:ie).eq.char(17)) go to 10
      enddo
      iexp=0
      go to 50
c
 10   ie2=i2
      ie1=ie+1
      i2 =ie-1
      iexp=1
      if(a(ie1:ie1).eq.char(12)                           ) then
         iexp=-1
      endif
      if(a(ie1:ie1).eq.char(12).or.a(ie1:ie1).eq.char(11)) then
         ie1=ie1+1
      endif
      itmp=0
      do i=ie1,ie2
         do j=1,10
            if(a(i:i).eq.char(j)) go to 30
         enddo
         go to 100
 30      itmp=itmp*10+j-1
      enddo
      iexp=iexp*itmp
c
c     ----- the number itself -----
c
 50   continue
      rep=.false.
      dum=zero
      do i=i1,i2
         if(a(i:i).ne.char(13)) then
            do j=1,10
               if(a(i:i).eq.char(j)) go to 70
            enddo
            go to 100
 70         dum=dum*ten+dble(j-1)
         else
            if(rep) go to 100
            iexp=iexp+i-i2
            rep=.true.
         endif
      enddo
      dum=dum*dble(isign)*ten**iexp
      x  =dum
      return
c
 100  continue       
      return
      end
      SUBROUTINE HND_HNDERR(LERR,ERRMSG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8  ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION    ERRMSG(LERR)
C
      WRITE(IW,9999)
      WRITE(IW,9998) ERRMSG
      WRITE(IW,9997)
C
      CALL ERRQUIT(' CALLS IT QUIT FROM HND_HNDERR ',0)
 9999 FORMAT(/,1X,72(1H-),/,31X,11HJOB STOPPED)
 9998 FORMAT(/,1X,9A8)
 9997 FORMAT(/,1X,72(1H-))
      END
      SUBROUTINE HND_GEOCLS(NFT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*8  ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION    ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-GEOCLS-'/
      CLOSE(UNIT=NFT,STATUS='KEEP',ERR=10)
      RETURN
   10 WRITE(IW,9999) NFT
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
 9999 FORMAT(/,' ----- ERROR CLOSING UNIT ',I3,' IN -GEOCLS- . STOP .')
      END
      SUBROUTINE HND_GEOOPN(NFT,GEOFIL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*80 GEOFIL
      CHARACTER*8  ERRMSG
      COMMON/HND_IOFILE/IR,IW
      DIMENSION    ERRMSG(3)
      DATA ERRMSG  /'PROGRAM ','STOP IN ','-GEOOPN-'/
      OPEN(UNIT=NFT,FILE=GEOFIL,STATUS='OLD',
     1     ACCESS='SEQUENTIAL',FORM=  'FORMATTED',ERR=10)
      RETURN
   10 WRITE(IW,9999) NFT
      CALL HND_HNDERR(3,ERRMSG)
      RETURN
 9999 FORMAT(/,' ----- ERROR OPENING FILE -',A80,'-',/,
     1         '       AS UNIT ',I3,' IN -GEOOPN- . STOP.')
      END
      SUBROUTINE HND_ZDAT(IAT,NAT,NVAR,ZVAL,ZLST,
     1                ZMT,NUMZMT,PRSZMT,FLGZMT,ZMTCHR,
     2                NUMVAR,PRSVAR,VARCHR,FRZVAR,FRZVAL,LST,
     3                IZMAT,IZ,IZFRZ,MXIZMT,NZMOD,DBUG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     THIS ROUTINE TAKES THE CHARACTER STRING VALUES FOR THE
C     -Z- MATRIX INPUT FOR THE INDICES, AND THE VALUES AND TRANSFORMS
C     THEM INTO INTEGER AND REAL VALUES. THE CHARACTER STRING
C     VALUES CAN BE EITHER NUMERIC (INTEGER AND REAL), OR A VARIABLE
C     THAT GETS REPLACED BY A VALUE FOR THE VARIABLE.
C
C     INPUT PARAMETERS:
C        IAT     - INTEGER    -  ATOM NUMBER BEING CONSIDERED
C        NAT     - INTEGER    -  TOTAL NUMBER OF ATOMS
C        NVAR    - INTEGER    -  TOTAL NUMBER OF VARIABLES
C        NUMZMT  - INTEGER    -  ARRAY OF NUMBERS OF WORDS PER LINE
C        PRSZMT  - CHARACTER*80 -  ARRAY OF PARSED -Z- DATA
C        FLGZMT  - INTEGER    -  ARRAY OF TYPES FOR -Z- DATA
C                                  0 ... BLANK
C                                  1 ... FLOATING
C                                  2 ... INTEGER
C                                  3 ... ALPHANUMERIC
C        ZMTCHR  - INTEGER    -  ARRAY OF LENGTHS FOR -Z- DATA
C        NUMVAR  - INTEGER    -  ARRAY OF NUM. OF PIECES OF INFO
C                                  FOR EACH VARIABLE FOR SUBSTITUTION
C        PRSVAR  - CHARACTER*80 -  ARRAY OF VARIABLES FOR SUBSTITUTION
C        VARCHR  - INTEGER    -  ARRAY OF LENGTHS FOR VARIABLES
C
C     OUTPUT PARAMETERS:
C        ZVAL    - REAL*8       -  ARRAY OF -Z- MATRIX VALUES
C        ZMT     - INTEGER    -  ARRAY OF -Z- MATRIX INDICES
C
C     S. CHIN: 11/08/90 - IBM KINGSTON, NY
C
      PARAMETER     (MXATOM=500)
      PARAMETER     (MAXGEO=MXATOM+1,MAXZMT=40,MAXVAR=256)
      LOGICAL       DBUG
      LOGICAL       LST
      LOGICAL       READY
      INTEGER     FLGZMT
      CHARACTER*80  PRSZMT
      INTEGER     ZMTCHR
      CHARACTER*80  PRSVAR
      INTEGER     VARCHR
      INTEGER     ZMT
      LOGICAL       FRZVAL
      LOGICAL       FRZVAR
      LOGICAL       CART
      LOGICAL       IZFRZ
      CHARACTER*8   ATNAME
      CHARACTER*1   PLUS
      CHARACTER*1   MINUS
      CHARACTER*8   ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_XYZGEO/XX(MAXGEO),YY(MAXGEO),ZZ(MAXGEO),
     1                  ATNAME(MAXGEO),ATNUM(MAXGEO),CART(MAXGEO)
      DIMENSION     ERRMSG(3)
      DIMENSION     NUMZMT(       MAXGEO)
      DIMENSION     PRSZMT(MAXZMT,MAXGEO)
      DIMENSION     FLGZMT(MAXZMT,MAXGEO)
      DIMENSION     ZMTCHR(MAXZMT,MAXGEO)
      DIMENSION     PRSVAR(MAXZMT,MAXVAR)
      DIMENSION     VARCHR(MAXZMT,MAXVAR)
      DIMENSION     FRZVAR(       MAXVAR)
      DIMENSION     NUMVAR(       MAXVAR)
      DIMENSION     FRZVAL(     3,MAXGEO)
      DIMENSION       ZVAL(     3,MAXGEO)
      DIMENSION        ZMT(     5,MAXGEO)
      DIMENSION       ZLST(     3,MAXGEO)
      DIMENSION      IZMAT(MXIZMT)
      DIMENSION      IZFRZ(MXIZMT)
      DATA ERRMSG /'PROGRAM ','STOP IN ','- ZDAT -'/
      DATA ZERO   /0.0D+00/
      DATA TWO    /2.0D+00/
      DATA PIDEG  /180.0D+00/
      DATA PLUS   /'+'/
      DATA MINUS  /'-'/
C
C     ---- THIS ATOM IS BEING INPUTED WITH CARTESIAN COORDINATES ----
C
      IF(FLGZMT(3,IAT).EQ.1.AND.FLGZMT(4,IAT).EQ.1) THEN
             IF(NUMZMT(  IAT).EQ.4.AND.FLGZMT(1,IAT).EQ.3.AND.
     1          FLGZMT(2,IAT).EQ.1                            ) THEN
c           READ(PRSZMT(1,IAT),*) ATNAME(IAT)
c           READ(PRSZMT(2,IAT),*) XX(IAT)
c           READ(PRSZMT(3,IAT),*) YY(IAT)
c           READ(PRSZMT(4,IAT),*) ZZ(IAT)
            call hnd_dparsc(prszmt(1,iat),zmtchr(1,iat),atname(iat),8)
            call hnd_dparsr(prszmt(2,iat),zmtchr(2,iat),xx(iat)      )
            call hnd_dparsr(prszmt(3,iat),zmtchr(3,iat),yy(iat)      )
            call hnd_dparsr(prszmt(4,iat),zmtchr(4,iat),zz(iat)      )
                 ZMTCHR(1,IAT)= ZMTCHR(1,IAT) - 2
         ELSEIF(NUMZMT(  IAT).EQ.5.AND.FLGZMT(1,IAT).EQ.3.AND.
     1          FLGZMT(2,IAT).EQ.2.AND.FLGZMT(5,IAT).EQ.1     ) THEN
c           READ(PRSZMT(1,IAT),*) ATNAME(IAT)
c           READ(PRSZMT(3,IAT),*) XX(IAT)
c           READ(PRSZMT(4,IAT),*) YY(IAT)
c           READ(PRSZMT(5,IAT),*) ZZ(IAT)
            call hnd_dparsc(prszmt(1,iat),zmtchr(1,iat),atname(iat),8)
            call hnd_dparsr(prszmt(3,iat),zmtchr(3,iat),xx(iat)      )
            call hnd_dparsr(prszmt(4,iat),zmtchr(4,iat),yy(iat)      )
            call hnd_dparsr(prszmt(5,iat),zmtchr(5,iat),zz(iat)      )
                 ZMTCHR(1,IAT)= ZMTCHR(1,IAT) - 2
         ELSE
            WRITE(IW,9994) IAT
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         CART(IAT) = .TRUE.
         RETURN
      ENDIF
C
C     ----- CONSTRUCT THE -Z- MATRIX -----
C
      IF(IAT.EQ.1.AND.(NUMZMT(IAT).NE.1.AND.NUMZMT(IAT).NE.4)) THEN
         WRITE(IW,9989)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(IAT.EQ.2.AND.(NUMZMT(IAT).NE.3.AND.NUMZMT(IAT).NE.4)) THEN
         WRITE(IW,9988)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(IAT.EQ.3.AND.(NUMZMT(IAT).NE.5.AND.NUMZMT(IAT).NE.4)) THEN
         WRITE(IW,9987)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- ANCHOR ATOM              -----
C
      ZMT(1,IAT) = IAT
      IF(FLGZMT(1,IAT) .EQ. 3) THEN
c        READ (PRSZMT(1,IAT),*) ATNAME(IAT)
         call hnd_dparsc(prszmt(1,iat),zmtchr(1,iat),atname(iat),8)
               ZMTCHR(1,IAT)= ZMTCHR(1,IAT) - 2
      ELSE
         WRITE(IW,9999)
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF(IAT .EQ. 1) GO TO 100
C
C     ----- BOND LENGTH -----
C
      IF(FLGZMT(2,IAT) .EQ. 2) THEN
c        READ (PRSZMT(2,IAT),*) ZMT(2,IAT)
         call hnd_dparsi(prszmt(2,iat),zmtchr(2,iat),zmt(2,iat)   )
      ELSE
         READY=.FALSE.
         DO JAT=1,IAT-1
            IF(PRSZMT(2,IAT) (1:8) .EQ. PRSZMT(1,JAT) (1:8)) THEN
               ZMT(2,IAT)=JAT
               READY=.TRUE.
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9992) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
      IF(FLGZMT(3,IAT) .EQ. 1) THEN
c        READ (PRSZMT(3,IAT),*) ZVAL(1,IAT)
         call hnd_dparsr(prszmt(3,iat),zmtchr(3,iat),zval(1,iat)  )
                                ZLST(1,IAT)=ZVAL(1,IAT)
      ELSEIF(FLGZMT(3,IAT) .EQ. 3) THEN
         READY=.FALSE.
         DO JVAR=1,NVAR
            IF((ZMTCHR(3, IAT).EQ.VARCHR(1,JVAR)).AND.
     1         (PRSZMT(3, IAT) (1:ZMTCHR(3, IAT)).EQ.
     2          PRSVAR(1,JVAR) (1:VARCHR(1,JVAR)))   ) THEN
c              READ (PRSVAR(2,JVAR),*) ZVAL(1,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(1,iat)  )
                                       ZLST(1,IAT)=ZVAL(1,IAT)
               FRZVAL(1,IAT)=FRZVAR(JVAR)
               READY=.TRUE.
               IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                 READ (PRSVAR(3,JVAR),*) ZLST(1,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(1,iat)  )
               ENDIF
            ELSEIF((ZMTCHR(3, IAT).EQ.VARCHR(1,JVAR)+1).AND.
     1             (PRSZMT(3, IAT) (3:ZMTCHR(3, IAT)).EQ.
     2              PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)))   ) THEN
               IF(PRSZMT(3, IAT) (2:2).EQ.PLUS) THEN
c                 READ (PRSVAR(2,JVAR),*) ZVAL(1,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(1,iat)  )
                                          ZLST(1,IAT)=ZVAL(1,IAT)
                  FRZVAL(1,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(1,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(1,iat)  )
                  ENDIF
               ELSEIF(PRSZMT(3, IAT) (2:2).EQ.MINUS) THEN
c                 READ (PRSVAR(2,JVAR),*) ZVAL(1,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(1,iat)  )
                                          ZLST(1,IAT)=ZVAL(1,IAT)
                  FRZVAL(1,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(1,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(1,iat)  )
                  ENDIF
                  ZVAL(1,IAT)=-ZVAL(1,IAT)
                  ZLST(1,IAT)=-ZLST(1,IAT)
               ENDIF
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9998) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ELSE
         WRITE(IW,9998) IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF(IAT .EQ. 2) GO TO 100
C
C     ----- BOND ANGLE -----
C
      IF(FLGZMT(4,IAT) .EQ. 2) THEN
c        READ (PRSZMT(4,IAT),*) ZMT(3,IAT)
         call hnd_dparsi(prszmt(4,iat),zmtchr(4,iat),zmt(3,iat)  )
      ELSE
         READY=.FALSE.
         DO JAT=1,IAT-1
            IF(PRSZMT(4,IAT) (1:8).EQ.PRSZMT(1,JAT) (1:8)) THEN
               ZMT(3,IAT) = JAT
               READY=.TRUE.
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9992) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
      IF(FLGZMT(5,IAT) .EQ. 1) THEN
c        READ (PRSZMT(5,IAT),*) ZVAL(2,IAT)
         call hnd_dparsr(prszmt(5,iat),zmtchr(5,iat),zval(2,iat)  )
                                ZLST(2,IAT)=ZVAL(2,IAT)
      ELSEIF(FLGZMT(5,IAT) .EQ. 3) THEN
         READY=.FALSE.
         DO JVAR=1,NVAR
            IF((ZMTCHR(5, IAT).EQ.VARCHR(1,JVAR)).AND.
     1         (PRSZMT(5, IAT) (1:ZMTCHR(5, IAT)).EQ.
     2          PRSVAR(1,JVAR) (1:VARCHR(1,JVAR)))    ) THEN
c              READ (PRSVAR(2,JVAR),*) ZVAL(2,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(2,iat)  )
                                       ZLST(2,IAT)=ZVAL(2,IAT)
               FRZVAL(2,IAT)=FRZVAR(JVAR)
               READY=.TRUE.
               IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                 READ (PRSVAR(3,JVAR),*) ZLST(2,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(2,iat)  )
               ENDIF
            ELSEIF((ZMTCHR(5, IAT).EQ.VARCHR(1,JVAR)+1).AND.
     1             (PRSZMT(5, IAT) (3:ZMTCHR(5, IAT)).EQ.
     2              PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)))    ) THEN
               IF(PRSZMT(5, IAT) (2:2).EQ.PLUS) THEN
c                 READ (PRSVAR(2,JVAR),*) ZVAL(2,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(2,iat)  )
                                          ZLST(2,IAT)=ZVAL(2,IAT)
                  FRZVAL(2,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(2,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(2,iat)  )
                  ENDIF
               ELSEIF(PRSZMT(5, IAT) (2:2).EQ.MINUS) THEN
c                 READ (PRSVAR(2,JVAR),*) ZVAL(2,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(2,iat)  )
                                          ZLST(2,IAT)=ZVAL(2,IAT)
                  FRZVAL(2,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(2,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(2,iat)  )
                  ENDIF
                  ZVAL(2,IAT)=-ZVAL(2,IAT)
                  ZLST(2,IAT)=-ZLST(2,IAT)
               ENDIF
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9997) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ELSE
         WRITE(IW,9997) IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
      IF(ZVAL(2,IAT).LT.ZERO.OR.ZVAL(2,IAT).GT.PIDEG) THEN
         WRITE(IW,9986) ZVAL(2,IAT),IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      IF(IAT .EQ. 3) GO TO 100
C
C     ----- TORSION -----
C
      IF(FLGZMT(6,IAT) .EQ. 2) THEN
c        READ (PRSZMT(6,IAT),*) ZMT(4,IAT)
         call hnd_dparsi(prszmt(6,iat),zmtchr(6,iat),zmt(4,iat)  )
      ELSE
         READY=.FALSE.
         DO JAT=1,IAT-1
            IF(PRSZMT(6,IAT) (1:8).EQ.PRSZMT(1,JAT) (1:8)) THEN
               ZMT(4,IAT)=JAT
               READY=.TRUE.
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9992) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
      IF(FLGZMT(7,IAT) .EQ. 1) THEN
c        READ (PRSZMT(7,IAT),*) ZVAL(3,IAT)
         call hnd_dparsr(prszmt(7,iat),zmtchr(7,iat),zval(3,iat)  )
                                ZLST(3,IAT)=ZVAL(3,IAT)
      ELSEIF(FLGZMT(7,IAT) .EQ. 3) THEN
         READY=.FALSE.
         DO JVAR=1,NVAR
            IF((ZMTCHR(7, IAT).EQ.VARCHR(1,JVAR)).AND.
     1         (PRSZMT(7, IAT) (1:ZMTCHR(7, IAT)).EQ.
     2          PRSVAR(1,JVAR) (1:VARCHR(1,JVAR)))    ) THEN
c              READ (PRSVAR(2,JVAR),*) ZVAL(3,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(3,iat)  )
                                       ZLST(3,IAT)=ZVAL(3,IAT)
               FRZVAL(3,IAT)=FRZVAR(JVAR)
               READY=.TRUE.
               IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                 READ (PRSVAR(3,JVAR),*) ZLST(3,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(3,iat)  )
               ENDIF
            ELSEIF((ZMTCHR(7, IAT).EQ.VARCHR(1,JVAR)+1).AND.
     1             (PRSZMT(7, IAT) (3:ZMTCHR(7, IAT)).EQ.
     2              PRSVAR(1,JVAR) (2:VARCHR(1,JVAR)))    ) THEN
               IF(PRSZMT(7, IAT) (2:2).EQ.PLUS) THEN
c                 READ (PRSVAR(2,JVAR),*) ZVAL(3,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(3,iat)  )
                                          ZLST(3,IAT)=ZVAL(3,IAT)
                  FRZVAL(3,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(3,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(3,iat)  )
                  ENDIF
               ELSEIF(PRSZMT(7, IAT) (2:2).EQ.MINUS) THEN
c                 READ (PRSVAR(2,JVAR),*) ZVAL(3,IAT)
         call hnd_dparsr(prsvar(2,jvar),varchr(2,jvar),zval(3,iat)  )
                                          ZLST(3,IAT)=ZVAL(3,IAT)
                  FRZVAL(3,IAT)=FRZVAR(JVAR)
                  READY=.TRUE.
                  IF(LST.AND.NUMVAR(JVAR).EQ.3) THEN
c                    READ (PRSVAR(3,JVAR),*) ZLST(3,IAT)
         call hnd_dparsr(prsvar(3,jvar),varchr(3,jvar),zlst(3,iat)  )
                  ENDIF
                  ZVAL(3,IAT)=-ZVAL(3,IAT)
                  ZLST(3,IAT)=-ZLST(3,IAT)
               ENDIF
            ENDIF
         ENDDO
         IF(.NOT.READY) THEN
            WRITE(IW,9996) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ELSE
         WRITE(IW,9996) IAT
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
C     ----- LAST INTEGER PARAMETER -----
C
      IF(NUMZMT(IAT).EQ.8) THEN
         IF(FLGZMT(8,IAT).EQ.2) THEN
c           READ(PRSZMT(8,IAT),*) ZMT(5,IAT)
         call hnd_dparsi(prszmt(8,iat),zmtchr(8,iat),zmt(5,iat)  )
         ELSE
            WRITE(IW,9993) IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
C     ----- POSSIBLE ERRORS -----
C
      IF(ZMT(5,IAT).EQ.0) THEN
         IF(ZVAL(3,IAT).LT.-TWO*PIDEG.OR.
     1      ZVAL(3,IAT).GT. TWO*PIDEG) THEN
            WRITE(IW,9985) ZVAL(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(ZVAL(3,IAT).LT.-PIDEG) ZVAL(3,IAT)=ZVAL(3,IAT)+TWO*PIDEG
         IF(ZVAL(3,IAT).GT. PIDEG) ZVAL(3,IAT)=ZVAL(3,IAT)-TWO*PIDEG
C
         IF(ZLST(3,IAT).LT.-TWO*PIDEG.OR.
     1      ZLST(3,IAT).GT. TWO*PIDEG) THEN
            WRITE(IW,9985) ZLST(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZLST(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
         IF(ZLST(3,IAT).LT.-PIDEG) ZLST(3,IAT)=ZLST(3,IAT)+TWO*PIDEG
         IF(ZLST(3,IAT).GT. PIDEG) ZLST(3,IAT)=ZLST(3,IAT)-TWO*PIDEG
      ELSE
         IF(ZVAL(3,IAT).LT. ZERO .OR.
     1      ZVAL(3,IAT).GT.PIDEG) THEN
            WRITE(IW,9985) ZVAL(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZVAL(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
C
         IF(ZLST(3,IAT).LT. ZERO .OR.
     1      ZLST(3,IAT).GT.PIDEG) THEN
            WRITE(IW,9985) ZLST(3,IAT),IAT
            WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(ZLST(I,IAT),I=1,3)
            CALL HND_HNDERR(3,ERRMSG)
         ENDIF
      ENDIF
C
  100 CONTINUE
      IF(DBUG) THEN
         WRITE(IW,9984) IAT,ATNAME(IAT)
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(  ZVAL(I,IAT),I=1,3),
     1                                         (FRZVAL(I,IAT),I=1,3)
         WRITE(IW,9995) IAT,(ZMT(I,IAT),I=1,5),(  ZLST(I,IAT),I=1,3),
     1                                         (FRZVAL(I,IAT),I=1,3)
      ENDIF
C
C     ----- CREATE -IZMAT- FOR -HONDO- INPUT -----
C
      IF(IAT.GT.1) THEN
         NZMOD=NZMOD+1
         IZ=IZ+1
         IZMAT(IZ)=1
         IZ=IZ+1
         IZMAT(IZ)=ZMT(1,IAT)
         IZ=IZ+1
         IZMAT(IZ)=ZMT(2,IAT)
         IZFRZ(NZMOD)=FRZVAL(1,IAT)
C
         IF(IAT.GT.2) THEN
            NZMOD=NZMOD+1
            IZ=IZ+1
            IZMAT(IZ)=2
            IZ=IZ+1
            IZMAT(IZ)=ZMT(1,IAT)
            IZ=IZ+1
            IZMAT(IZ)=ZMT(2,IAT)
            IZ=IZ+1
            IZMAT(IZ)=ZMT(3,IAT)
            IZFRZ(NZMOD)=FRZVAL(2,IAT)
C
            IF(IAT.GT.3) THEN
               IF(ZMT(5,IAT).EQ.0) THEN
                  NZMOD=NZMOD+1
                  IZ=IZ+1
                  IZMAT(IZ)=3
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(1,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(2,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(3,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(4,IAT)
                  IZFRZ(NZMOD)=FRZVAL(3,IAT)
               ELSE
                  NZMOD=NZMOD+1
                  IZ=IZ+1
                  IZMAT(IZ)=2
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(1,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(2,IAT)
                  IZ=IZ+1
                  IZMAT(IZ)=ZMT(4,IAT)
                  IZFRZ(NZMOD)=FRZVAL(3,IAT)
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF(IZ.GT.MXIZMT) THEN
         WRITE(IW,9991) IZ,MXIZMT,IAT
         CALL HND_HNDERR(3,ERRMSG)
      ENDIF
C
      RETURN
 9999 FORMAT(' FIRST PARSED -Z- DATA IS NOT A CHARACTER STRING. STOP')
 9998 FORMAT(' THE 3-D  PIECE OF -Z- DATA FOR ATOM = ',I4,
     1       ' IS NEITHER FLOATING POINT',/,' NOR ALPHANUMERIC',
     2       ' OR COULD NOT BE MATCHED WITH A VARIABLE. STOP')
 9997 FORMAT(' THE 5-D  PIECE OF -Z- DATA FOR ATOM = ',I4,
     1       ' IS NEITHER FLOATING POINT',/,' NOR ALPHANUMERIC',
     2       ' OR COULD NOT BE MATCHED WITH A VARIABLE. STOP')
 9996 FORMAT(' THE 7-D  PIECE OF -Z- DATA FOR ATOM = ',I4,
     1       ' IS NEITHER FLOATING POINT',/,' NOR ALPHANUMERIC',
     2       ' OR COULD NOT BE MATCHED WITH A VARIABLE. STOP')
 9995 FORMAT(' IAT=',I5,' ZMAT=',5I3,3F12.5,3(2X,L4))
 9994 FORMAT(' SOMETHING IS WRONG WITH THE CARTESIAN COORDINATES INPUT',
     1       ' DATA FOR THIS ATOM, IAT = ',I3,' . STOP')
 9993 FORMAT(' THE 8-TH (INTEGER) PARAMETER FOR ATOM = ',I4,
     1       ' IS OF THE WRONG TYPE. STOP')
 9992 FORMAT(' ATOM REFERED TO FOR IAT = ',I3,' IS INVALID. STOP')
 9991 FORMAT(' TOO MANY -IZMAT- ENTRIES, IZ,MXIZMT = ',2I5,
     1       ' FOR -IAT- = ',I5,' STOP.')
 9989 FORMAT(' THE -ZMATRIX- DATA FOR ATOM -1- IS INCORRECT. STOP')
 9988 FORMAT(' THE -ZMATRIX- DATA FOR ATOM -2- IS INCORRECT. STOP')
 9987 FORMAT(' THE -ZMATRIX- DATA FOR ATOM -3- IS INCORRECT. STOP')
 9986 FORMAT(' BOND ANGLE = ',F10.2,' FOR ATOM = ',I4,' IS OUT OF',
     1       ' ALLOWED RANGE ... 0.0 TO 180.0 . STOP')
 9985 FORMAT(' TORSION ANGLE OR BOND ANGLE = ',F10.2,' FOR ATOM = ',I4,
     1       ' IS OUT OF ALLOWED RANGE ... -360.0 TO 360.0 . STOP')
 9984 FORMAT(' IAT=',I5,' ATNAME=',1X,A8)
      END
      SUBROUTINE HND_ZXYZ(NAT,NZMT,ZVAL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ROUTINE COMPUTES :
C
C     THE CARTESIAN COORDINATES OF A MOLECULE FROM BONDING COORDINATES
C     (BOND LENGTHS, BOND AND DIHEDRAL ANGLES).
C
C     ATOMS NA, NB, NC, HAVE KNOWN COORDINATES AND ARE NOT COLLINEAR.
C     RCD   =ZVAL(1,I) : C-D BOND LENGTH (IN ANGSTROM)
C     THBCD =ZVAL(2,I) : B-C-D BOND ANGLE (IN DEGREES)
C     PHABCD=ZVAL(3,I) : DIHEDRAL ANGLE (IN DEGREES) OF C-D
C                        RELATIVE TO A-B MEASURED CLOCKWISE
C                        ( + ) ALONG THE DIRECTION B TO C
C
C     MEANING OF PARAMETER NZMT(5,I):
C          =0     ANGLES FROM DATA
C          =+/- 1 ZVAL(3,I) IS BOND ANGLE A-C-D
C          =2     BOND ANGLE B-C-D IS TETRAHEDRAL
C          =3     ATOMS B-C-D ARE COLLINEAR
C
      PARAMETER     (MXATOM=500)
      PARAMETER     (MAXGEO=MXATOM+1,MAXWRD=40,MAXVAR=256)
      LOGICAL       DBUG
      LOGICAL       CART
      CHARACTER*8   ATNAME
      CHARACTER*8   ERRMSG
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_XYZGEO/XX(MAXGEO),YY(MAXGEO),ZZ(MAXGEO),
     1                  ATNAME(MAXGEO),ATNUM(MAXGEO),CART(MAXGEO)
      DIMENSION   NZMT(5,MAXGEO)
      DIMENSION   ZVAL(3,MAXGEO)
      DIMENSION   ERRMSG(3)
      DATA ERRMSG    /'PROGRAM ','STOP IN ','- ZXYZ -'/
      DATA ZERO,ONE  /0.0D+00,1.0D+00/
      DATA TWO,THREE /2.0D+00,3.0D+00/
      DATA PT1,MAX   /0.1D+00,6/
C
      DBUG=.FALSE.
C
      PIFAC=3.1415926536D+00/180.0D+00
      COSTD=-ONE/THREE
      SINTD= TWO/THREE* SQRT(TWO)
      NA   =0
      NB   =0
      NC   =0
      ND   =0
      RCD  =ZERO
C
C     ----- ATOM -1- AND -2- -----
C
      IAT=1
      IF(CART(IAT)) THEN
         IF(DBUG) THEN
            WRITE(IW,9997) IAT
         ENDIF
      ELSE
         ND=NZMT(1,IAT)
         XX(ND)=ZERO
         YY(ND)=ZERO
         ZZ(ND)=ZERO
         CART(IAT)=.TRUE.
         IF(DBUG) THEN
            WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
         ENDIF
      ENDIF
C
      IF(NAT.GT.1) THEN
         IAT=2
         IF(CART(IAT)) THEN
            IF(DBUG) THEN
               WRITE(IW,9997) IAT
            ENDIF
         ELSE
            ND=NZMT(1,IAT)
            NC=NZMT(2,IAT)
            XX(ND)=ZERO
            YY(ND)=ZERO
            ZZ(ND)=ZVAL(1,IAT)
            RCD   =ZVAL(1,IAT)
            CART(IAT)=.TRUE.
            IF(DBUG) THEN
               WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
            ENDIF
         ENDIF
C
         IF(NAT.GT.2) THEN
C
C     ----- ATOM -3- -----
C
            IAT=3
            IF(CART(IAT)) THEN
               IF(DBUG) THEN
                  WRITE(IW,9997) IAT
               ENDIF
            ELSE
                   IF(NZMT(5,IAT).EQ.3) THEN
                  CCOS=-ONE
                  SSIN= ZERO
               ELSEIF(NZMT(5,IAT).EQ.2) THEN
                  CCOS=COSTD
                  SSIN=SINTD
               ELSEIF(NZMT(5,IAT).EQ.0) THEN
                  THBCD=ZVAL(2,IAT)*PIFAC
                  CCOS= COS(THBCD)
                  SSIN= SIN(THBCD)
               ELSE
                  WRITE(IW,9999) IAT
                  CALL HND_HNDERR(3,ERRMSG)
               ENDIF
               ND =NZMT(1,IAT)
               NC =NZMT(2,IAT)
               NB =NZMT(3,IAT)
               RCD=ZVAL(1,IAT)
               IF(DBUG) THEN
                  WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
               ENDIF
               XX(ND)=         RCD*SSIN
               YY(ND)=ZERO
               ZZ(ND)=ZZ(NC) - RCD*CCOS
               IF(ZZ(NC).LT.ZZ(NB)) THEN
                  ZZ(ND)=ZZ(NC) + RCD*CCOS
               ENDIF
               CART(IAT) = .TRUE.
            ENDIF
C
C     ----- ATOM -4- AND HIGHER -----
C
            IF(NAT.GT.3) THEN
               DO 20 IAT=4,NAT
                  IF(CART(IAT)) THEN
                     IF(DBUG) THEN
                        WRITE(IW,9997) IAT
                     ENDIF
                     GO TO 20
                  ENDIF
                  NA =NZMT(4,IAT)
                  NB =NZMT(3,IAT)
                  NC =NZMT(2,IAT)
                  ND =NZMT(1,IAT)
                  RCD=ZVAL(1,IAT)
                  IF(DBUG) THEN
                     WRITE(IW,9998) IAT,NA,NB,NC,ND,RCD
                  ENDIF
C
C     ----- LINEAR CONNECTION FOR THIS ATOM : D...C...B -----
C
                  IF(NZMT(5,IAT).EQ.3) THEN
                     RCB= SQRT((XX(NB)-XX(NC))**2+
     1                         (YY(NB)-YY(NC))**2+(ZZ(NB)-ZZ(NC))**2)
                     XX(ND)=XX(NC)+RCD*(XX(NC)-XX(NB))/RCB
                     YY(ND)=YY(NC)+RCD*(YY(NC)-YY(NB))/RCB
                     ZZ(ND)=ZZ(NC)+RCD*(ZZ(NC)-ZZ(NB))/RCB
                  ELSE
C
C     ----- OTHER CONNECTIONS FOR THIS ATOM : DEFINE LOCAL FRAME -----
C
                     XB=XX(NB)-XX(NC)
                     YB=YY(NB)-YY(NC)
                     ZB=ZZ(NB)-ZZ(NC)
                     RCB= SQRT(XB*XB+YB*YB+ZB*ZB)
                     XB=XB/RCB
                     YB=YB/RCB
                     ZB=ZB/RCB
                     T11=XB
                     T21=YB
                     T31=ZB
                     XA=XX(NA)-XX(NC)
                     YA=YY(NA)-YY(NC)
                     ZA=ZZ(NA)-ZZ(NC)
                     RCA= SQRT(XA*XA+YA*YA+ZA*ZA)
                     XA=XA/RCA
                     YA=YA/RCA
                     ZA=ZA/RCA
                     DOT=XA*XB+YA*YB+ZA*ZB
                     XA=XA-DOT*XB
                     YA=YA-DOT*YB
                     ZA=ZA-DOT*ZB
                     DUM= SQRT(XA*XA+YA*YA+ZA*ZA)
                     XA=XA/DUM
                     YA=YA/DUM
                     ZA=ZA/DUM
                     T12=XA
                     T22=YA
                     T32=ZA
                     T13=T21*T32-T31*T22
                     T23=T31*T12-T11*T32
                     T33=T11*T22-T21*T12
C
C     ----- GET CARTESIAN COORDINATES IN LOCAL  FRAME -----
C
                         IF(NZMT(5,IAT).EQ.0) THEN
                        PHI  = ZVAL(2,IAT)*PIFAC
                        THETA=-ZVAL(3,IAT)*PIFAC
                        CCOSP= COS(PHI)
                        SSINP= SIN(PHI)
                        CCOST= COS(THETA)
                        SSINT= SIN(THETA)
                     ELSEIF(NZMT(5,IAT).EQ.-1.OR.
     1                      NZMT(5,IAT).EQ. 1    ) THEN
                        RAB  =SQRT((XX(NB)-XX(NA))**2+(YY(NB)-YY(NA))**2
     1                            +(ZZ(NB)-ZZ(NA))**2)
                        PHI=ZVAL(2,IAT)*PIFAC
                        BET=ZVAL(3,IAT)*PIFAC
                        ALP=PHI
                        GAM  =ACOS((RCB**2+RCA**2-RAB**2)/(TWO*RCB*RCA))
                        THETA=ACOS(( COS(BET)- COS(ALP)* COS(GAM))/
     1                             (           SIN(ALP)* SIN(GAM))  )
                        IF(NZMT(5,IAT).EQ.-1) THEN
                           THETA=-THETA
                        ENDIF
                        CCOSP= COS(PHI)
                        SSINP= SIN(PHI)
                        CCOST= COS(THETA)
                        SSINT= SIN(THETA)
                     ELSEIF(NZMT(5,IAT).EQ.2) THEN
                        CCOSP= COSTD
                        SSINP= SINTD
                        THETA=ZVAL(3,IAT)*PIFAC
                        CCOST= COS(THETA)
                        SSINT= SIN(THETA)
                     ELSE
                        WRITE(IW,9999) IAT
                        CALL HND_HNDERR(3,ERRMSG)
                     ENDIF
                     XXD=RCD*CCOSP
                     YYD=RCD*SSINP*CCOST
                     ZZD=RCD*SSINP*SSINT
C
C     ----- GET CARTESIAN COORDINATES IN MASTER FRAME -----
C
                     XX(ND)=T11*XXD+T12*YYD+T13*ZZD + XX(NC)
                     YY(ND)=T21*XXD+T22*YYD+T23*ZZD + YY(NC)
                     ZZ(ND)=T31*XXD+T32*YYD+T33*ZZD + ZZ(NC)
                  ENDIF
                  CART(IAT)=.TRUE.
   20             CONTINUE
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
 9999 FORMAT(' ILLEGAL 5-TH -ZMT- PARAMETER FOR -IAT- = ',I3,' STOP.')
 9998 FORMAT(' IAT,NA,NB,NC,ND,RCD = ',5I4,F12.5)
 9997 FORMAT(' IAT = ',I4,' IS ALREADY SPECIFIED IN CARTESIAN SPACE.')
      END
      SUBROUTINE HND_PRSQ(V,M,N,NDIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ----- PRINT OUT A SQUARE MATRIX -----
C
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_LISTNG/LIST
      DIMENSION V(NDIM,1)
      DIMENSION IC(5),C(5)
      DATA VTOL  /1.5D-01/
      DATA ICMAX /5/
C
      LIST=1
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IF(LIST.EQ.0.AND.N.GT.0) GO TO 200
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 120 J = 1,N
      IF(LIST.EQ.0) WRITE (IW,9048) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9148) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9248) J,(V(J,I),I = IMIN,IMAX)
  120 CONTINUE
      IF (IMAX .LT. M) GO TO 100
      RETURN
C
  200 CONTINUE
C
      DO 240 J=1,M
      WRITE(IW,9328) J
C
      II=0
      DO 220 I=1,N
      IF( ABS(V(I,J)).LT.VTOL) GO TO 220
      II=II+1
      IC (II)=I
      C  (II)=V(I,J)
      IF(II.LT.ICMAX) GO TO 220
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
      II=0
  220 CONTINUE
      IF(II.EQ.0) GO TO 240
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
  240 CONTINUE
C
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(6X,10(4X,I3,4X))
 9048 FORMAT(I5,1X,10F11.5)
 9128 FORMAT(6X,7(6X,I3,6X))
 9148 FORMAT(I5,1X,7F15.10)
 9228 FORMAT(6X,7(6X,I3,6X))
 9248 FORMAT(I5,1X,7E15.8)
 9328 FORMAT(/,' VEC =',I5,/,1X,10(1H-))
 9348 FORMAT(5(I5,F11.5))
      END
      SUBROUTINE HND_PREV(V,E,M,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- PRINT OUT E AND V-MATRICES
C
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_LISTNG/LIST
      DIMENSION V(NDIM,1),E(1)
      DIMENSION IC(5),C(5)
      DATA VTOL  /1.5D-01/
      DATA ICMAX /5/
C
      LIST=1
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IF(LIST.EQ.0.AND.N.GT.0) GO TO 200
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. M) IMAX = M
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9068) (E(I),I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9168) (E(I),I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9268) (E(I),I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 120 J = 1,N
      IF(LIST.EQ.0) WRITE (IW,9048) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9148) J,(V(J,I),I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9248) J,(V(J,I),I = IMIN,IMAX)
  120 CONTINUE
      IF (IMAX .LT. M) GO TO 100
      RETURN
C
  200 CONTINUE
C
      DO 240 J=1,M
      WRITE(IW,9328) J,E(J)
C
      II=0
      DO 220 I=1,N
      IF( ABS(V(I,J)).LT.VTOL) GO TO 220
      II=II+1
      IC (II)=I
      C  (II)=V(I,J)
      IF(II.LT.ICMAX) GO TO 220
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
      II=0
  220 CONTINUE
      IF(II.EQ.0) GO TO 240
      WRITE(IW,9348) (IC(IDUM),C(IDUM),IDUM=1,II)
  240 CONTINUE
C
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(15X,10(4X,I3,4X))
 9048 FORMAT(I5,1X,10F11.5)
 9068 FORMAT(15X,10F11.5)
 9128 FORMAT(15X,7(6X,I3,6X))
 9148 FORMAT(10X,I5,7F15.10)
 9168 FORMAT(15X,7F15.10)
 9228 FORMAT(15X,7(6X,I3,6X))
 9248 FORMAT(I5,1X,7E15.8)
 9268 FORMAT(15X,7E15.8)
 9328 FORMAT(/,' VEC =',I5,' EIG/OCC =',F11.5,/,1X,31(1H-))
 9348 FORMAT(5(I5,F11.5))
      END
      SUBROUTINE HND_PRTR(D,N)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- PRINT OUT A TRIANGULAR MATRIX -----
C
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_LISTNG/LIST
      DIMENSION D(1),DD(10)
C
      LIST=1
      IF(LIST.EQ.0) MAX=10
      IF(LIST.EQ.1) MAX=7
      IF(LIST.EQ.2) MAX=7
C
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. N) IMAX = N
      WRITE (IW,9008)
      IF(LIST.EQ.0) WRITE (IW,9028) (I,I = IMIN,IMAX)
      IF(LIST.EQ.1) WRITE (IW,9128) (I,I = IMIN,IMAX)
      IF(LIST.EQ.2) WRITE (IW,9228) (I,I = IMIN,IMAX)
      WRITE (IW,9008)
      DO 160 J = 1,N
      K = 0
      DO 140 I = IMIN,IMAX
      K = K+1
      II = MAX0( I, J)
      JJ = MIN0( I, J)
      IJ = (II*(II-1))/2 + JJ
  140 DD(K) = D(IJ)
      IF(LIST.EQ.0) WRITE (IW,9048) J,(DD(I),I = 1,K)
      IF(LIST.EQ.1) WRITE (IW,9148) J,(DD(I),I = 1,K)
      IF(LIST.EQ.2) WRITE (IW,9248) J,(DD(I),I = 1,K)
  160 CONTINUE
      IF (IMAX .LT. N) GO TO 100
      RETURN
 9008 FORMAT(/)
 9028 FORMAT(6X,10(4X,I3,4X))
 9048 FORMAT(I5,1X,10F11.5)
 9128 FORMAT(6X,7(6X,I3,6X))
 9148 FORMAT(I5,1X,7F15.10)
 9228 FORMAT(6X,7(6X,I3,6X))
 9248 FORMAT(I5,1X,7E15.8)
      END
      SUBROUTINE HND_TFTR(H,F,Q,T,IA,M,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- H(M,M) = Q(DAGGER)(N,M) * F(N,N) * Q(N,M) -----
C
      DIMENSION H(1),F(1),Q(NDIM,1),T(1),IA(1)
      DATA SMALL /1.0D-11/
      DATA ZERO  /0.0D+00/
      IJ = 0
      DO 180 J = 1,M
      IK = 0
      DO 140 I = 1,N
      DUM = ZERO
      QIJ = Q(I,J)
      MAX = I-1
      IF (MAX .EQ. 0) GO TO 120
      DO 100 K = 1,MAX
      IK = IK+1
      T(K) = T(K)+ F(IK)*QIJ
      DUM  = DUM + F(IK)*Q(K,J)
  100 CONTINUE
  120 IK = IK+1
      T(I) = DUM + F(IK)*QIJ
  140 CONTINUE
      DO 160 I = 1,J
      IJ = IJ+1
      HIJ = ddot(N,Q(1,I),1,T,1)
      IF( ABS(HIJ).LT.SMALL) HIJ=ZERO
      H(IJ)=HIJ
  160 CONTINUE
  180 CONTINUE
      RETURN
      END
      SUBROUTINE HND_DIAGIV(A,VEC,EIG,IA,NVEC,N,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ----- GENERAL CALLING  ROUTINE FOR -----
C               -GIVEIS- OR -GIVDIA-
C
#include "mafdecls.fh"
#include "global.fh"
      logical status
      COMMON/HND_IOFILE/IR,IW
      COMMON/HND_DIAOPT/MODDIA
      DIMENSION A(1),VEC(NDIM,1),EIG(1),IA(1)
      dimension x(1),ix(1)
      equivalence (dbl_mb(1),x(1)),(int_mb(1),ix(1))
c
c     ----- memory -----
c
      status=ma_push_get(MT_DBL,4*ndim,'mem_dbl dia',i_idbl,idbl)
      status=ma_push_get(MT_INT,3*ndim,'mem_int dia',i_iint,iint)
c
c     ----- call -----
c
      CALL HND_GIVDIA(A,VEC,EIG,IA,N,NDIM,
     1     x(idbl), x(idbl+  ndim), x(idbl+2*ndim), x(idbl+3*ndim),
     2                              x(idbl+2*ndim), x(idbl+2*ndim),
     3    ix(iint),ix(iint+  ndim),ix(iint+2*ndim)                )
c
c     ----- memory -----
c
      status=ma_pop_stack(i_iint)
      if(.not.status) then
         call errquit('hnd_diagiv,ma_pop_stack of -iint- failed',911)
      endif
      status=ma_pop_stack(i_idbl)
      if(.not.status) then
         call errquit('hnd_diagiv,ma_pop_stack of -idbl- failed',911)
      endif
c
      RETURN
      END
      SUBROUTINE HND_GIVDIA(A,VEC,EIG,IA,N,NDIM,
     1           w,gamma,beta,betasq,p,q,iposv,ivpos,iord)
      IMPLICIT REAL*8 (A-H,O-Z)
C        -----  A GIVENS HOUSHOLDER MATRIX DIAGONALIZATION  -----
C        -----  ROUTINE SAME AS EIGEN BUT WORKS WITH A      -----
C        -----  LINEAR ARRAY.                               -----
      DIMENSION A(1),VEC(NDIM,1),EIG(1),IA(1)
      DIMENSION W(*),GAMMA(*),BETA(*),BETASQ(*)
      DIMENSION P(*),Q(*),IPOSV(*),IVPOS(*),IORD(*)
      DATA ZERO,PT5,ONE,TWO /0.0D+00,0.5D+00,1.0D+00,2.0D+00/
      DATA RHOSQ /1.0D-22/
C
      IF(N.EQ.0) GO TO 560
      N1=N-1
      N2=N-2
      GAMMA(1)=A(1)
      IF(N2) 280,270,120
  120 DO 260 NR=1,N2
      IK=IA(NR+1)+NR
      B=A(IK)
      S=ZERO
      DO 130 I=NR,N2
      IJ=IA(I+2)+NR
  130 S=S+A(IJ)**2
C
C        -----  PREPARE FOR POSSIBLE BYPASS OF TRANSFORMATION ----
C
      A(IK)=ZERO
      IF(S.LE.ZERO) GO TO 250
      S=S+B*B
      SGN=+ONE
      IF(B.GE.ZERO) GO TO 160
      SGN=-ONE
  160 SQRTS= SQRT(S)
      D=SGN/(SQRTS+SQRTS)
      TEMP= SQRT(PT5+B*D)
      W(NR)=TEMP
      A(IK)=TEMP
      D=D/TEMP
      B=-SGN*SQRTS
C
C        -----  -D- IS FACTOR OF PROPORTIONALITY. NOW       -----
C        -----  COMPUTE AND SAVE -W- VECTOR. EXTRA SINGLY   -----
C        -----  SUBSCRIPTED -W- VECTOR FOR SPEED.           -----
C
      DO 170 I=NR,N2
      IJ=IA(I+2)+NR
      TEMP=D*A(IJ)
      W(I+1)=TEMP
  170 A(IJ)=TEMP
C
C        -----  PREMULTIPLY VECTOR -W- BY MATRIX -A- TO     -----
C        -----  OBTAIN -P- VECTOR. SIMULTANEOUSLY ACCUMULATE ----
C        -----  DOT PRODUCT -WP- -- SCALR -K-.              -----
C
      WTAW=ZERO
      DO 220 I=NR,N1
      SUM=ZERO
      II=IA(I+1)
      DO 180 J=NR,I
      IJ=II+J+1
  180 SUM=SUM+A(IJ)*W(J)
      I1=I+1
      IF(N1.LT.I1) GO TO 210
      DO 200 J=I1,N1
      IJ=IA(J+1)+I+1
  200 SUM=SUM+A(IJ)*W(J)
  210 P(I)=SUM
  220 WTAW=WTAW+SUM*W(I)
      DO 230 I=NR,N1
  230 Q(I)=P(I)-WTAW*W(I)
C
C        -----  NOW FORM -PAP- MATRIX, REQUIRED PART        -----
C
      DO 240 J=NR,N1
      QJ=Q(J)
      WJ=W(J)
      JJ=J+1
      DO 240 I=J,N1
      IJ=IA(I+1)+JJ
  240 A(IJ)=A(IJ)-TWO*(W(I)*QJ+WJ*Q(I))
  250 BETA(NR)=B
      BETASQ(NR)=B*B
      IL=IK+1
  260 GAMMA(NR+1)=A(IL)
  270 IJ=IA(N)+N-1
      B=A(IJ)
      BETA(N-1)=B
      BETASQ(N-1)=B*B
      IJ=IJ+1
      GAMMA(N)=A(IJ)
  280 BETASQ(N)=ZERO
C
C        -----  ADJOIN AN IDENTYTY MATRIX TO BE POST-       -----
C        -----  MULTIPLIED BY ROTATIONS                     -----
C
      DO 300 I=1,N
      DO 299 J=1,N
  299 VEC(I,J)=ZERO
  300 VEC(I,I)=ONE
      M=N
      SUM=ZERO
      NPAS=1
      GO TO 400
  310 SUM=SUM+SHIFT
      COSA=ONE
      G=GAMMA(1)-SHIFT
      PP=G
      PPBS=PP*PP+BETASQ(1)
      PPBR= SQRT(PPBS)
      DO 370 J=1,M
      COSAP=COSA
      IF(PPBS.NE.ZERO) GO TO 320
      SINA=ZERO
      SINA2=ZERO
      COSA=ONE
      GO TO 350
  320 SINA=BETA(J)/PPBR
      SINA2=BETASQ(J)/PPBS
      COSA=PP/PPBR
C
C        -----  POSTMULTIPLY IDENTITY BY -P- TRANSPOSE .    -----
C
      NT=J+NPAS
      IF(NT.LT.N) GO TO 330
      NT=N
  330 CONTINUE
      DO 340 I=1,NT
      TEMP=COSA*VEC(I,J)+SINA*VEC(I,J+1)
      VEC(I,J+1)=-SINA*VEC(I,J)+COSA*VEC(I,J+1)
  340 VEC(I,J)=TEMP
  350 DIA=GAMMA(J+1)-SHIFT
      U=SINA2*(G+DIA)
      GAMMA(J)=G+U
      G=DIA-U
      PP=DIA*COSA-SINA*COSAP*BETA(J)
      IF(J.NE.M) GO TO 360
      BETA(J)=SINA*PP
      BETASQ(J)=SINA2*PP*PP
      GO TO 380
  360 PPBS=PP*PP+BETASQ(J+1)
      PPBR= SQRT(PPBS)
      BETA(J)=SINA*PPBR
  370 BETASQ(J)=SINA2*PPBS
  380 GAMMA(M+1)=G
C
C        -----  TEST FOR CONVERGENCE OF LAST DIAGONAL ELEMENT ----
C
      NPAS=NPAS+1
      IF(BETASQ(M).GT.RHOSQ) GO TO 410
  390 EIG(M+1)=GAMMA(M+1)+SUM
  400 BETA(M)=ZERO
      BETASQ(M)=ZERO
      M=M-1
      IF(M.EQ.0) GO TO 430
      IF(BETASQ(M).LE.RHOSQ) GO TO 390
C
C        -----  TAKE ROOT OF CORMER 2 BY 2 NEAREST TO       -----
C        -----  LOWER DIAGONAL IN VALUE AS ESTIMATE OF      -----
C        -----  EIGENVALUE TO USE FOR SHIFT                 -----
C
  410 A2=GAMMA(M+1)
      R2=PT5*A2
      R1=PT5*GAMMA(M)
      R12=R1+R2
      DIF=R1-R2
      TEMP= SQRT(DIF*DIF+BETASQ(M))
      R1=R12+TEMP
      R2=R12-TEMP
      DIF=  ABS(A2-R1)-  ABS(A2-R2)
      IF(DIF.LT.ZERO) GO TO 420
      SHIFT=R2
      GO TO 310
  420 SHIFT=R1
      GO TO 310
  430 EIG(1)=GAMMA(1)+SUM
      DO 440 J=1,N
      IPOSV(J)=J
      IVPOS(J)=J
  440  IORD(J)=J
      M=N
      GO TO 470
  450 DO 460 J=1,M
      IF(EIG(J).LE.EIG(J+1)) GO TO 460
      TEMP=EIG(J)
      EIG(J)=EIG(J+1)
      EIG(J+1)=TEMP
      ITEMP=IORD(J)
      IORD(J)=IORD(J+1)
      IORD(J+1)=ITEMP
  460 CONTINUE
  470 M=M-1
      IF(M.NE.0) GO TO 450
      IF(N1.EQ.0) GO TO 500
      DO 490 L=1,N1
      NV=IORD(L)
      NP=IPOSV(NV)
      IF(NP.EQ.L) GO TO 490
      LV=IVPOS(L)
      IVPOS(NP)=LV
      IPOSV(LV)=NP
      DO 480 I=1,N
      TEMP=VEC(I,L)
      VEC(I,L)=VEC(I,NP)
  480 VEC(I,NP)=TEMP
  490 CONTINUE
C
C        -----  BACK TRANSFORM THE VECTORS OF THE TRIPLE    -----
C        -----  DIAGONAL MATRIX.                            -----
C
  500 DO 550 NRR=1,N
      K=N1
  510 K=K-1
      IF(K.LE.0) GO TO 550
      SUM=ZERO
      DO 520 I=K,N1
      IJ=IA(I+1)+K
  520 SUM=SUM+VEC(I+1,NRR)*A(IJ)
      SUM=SUM+SUM
      DO 530 I=K,N1
      IJ=IA(I+1)+K
  530 VEC(I+1,NRR)=VEC(I+1,NRR)-SUM*A(IJ)
      GO TO 510
  550 CONTINUE
  560 CONTINUE
      RETURN
      END
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

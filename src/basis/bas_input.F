      subroutine bas_input(rtdb)
c $Id: bas_input.F,v 1.37 1997-12-10 01:00:53 d3e129 Exp $
      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "inp.fh"
#include "stdio.fh"
c
c   basis [<name>] [library [<standard set>] [file <filename>] \
c         [spherical|cartesian] [segment||nosegment] [print|noprint]\
c         [ecp ecp_name]
c
c     tag library [tag in library] <standard set> [file <filename>]
c     tag <shell type>
c       <exponent> <contraction coefficients>
c       ...
c     end basis
c
c also used for ecp input
c
c   ecp [<name>] [library [<standard set>] [file <filename>] \
c       [spherical|cartesian] [segment||nosegment] [print] [noprint]
c
c     tag nelec <integer>
c     tag library [tag in library] <standard set> [file <filename>]
c     tag <shell type>
c       <r-exponent> <exponent> <contraction coefficients>
c       ...
c     end ecp
c
c     parse the main directive
c
c::functions
      logical  bas_do_destroy
      external bas_do_destroy
c::passed
      integer rtdb              ! [input] handle to database
c::local
      integer nopt
      parameter (nopt = 10)
      character*10 opts(nopt)
      character*255 test, name, filename, standard, ecpname
      character*1000 errmsg
      logical status, ospherical, osegment, oprint, o_add_ecpname
      logical oIs_ecp
      integer ind, basis
      logical bas_add_ucnt_init, bas_set_spherical
      external bas_add_ucnt_init, bas_set_spherical
      data opts /'spherical', 'cartesian', 'segment', 'nosegment',
     $     'library', 'file', 'print', 'noprint', 'ecpset', 'version'/
c
      ecpname = '                                                '
c
c     Check is a basis/ecp directive and read in name of the basis
c
      call inp_set_field(0)
      status = inp_a(test)
      oIs_ecp = inp_compare(.false.,test,'ecp')
      status = status.and.
     &    (inp_compare(.false.,test,'basis') .or.
     &     oIs_ecp)
      if (.not.status) goto 10000
*:debug-s
*debug:      write(luout,*)' debug output - start'
*debug:      write(luout,*)'      test?        ',test
*debug:      write(luout,*)' oIs_ecp is        :',oIs_ecp
*debug:      write(luout,*)' debug output - end  '
*:debug-e
c
c     Parse rest of basis/ecp directive line
c
      name = ' '
      filename = ' '
      test = ' '
      standard = ' '
      ospherical = .false. ! Default is cartesian
      osegment   = .true.  ! Default is to force segmentation
      oprint     = .true.  ! Default is to print the basis on input
      o_add_ecpname = .false. ! Default is to not associate a name
c
 10   if (inp_a(test)) then
c
         if (.not. inp_match(nopt, .false., test, opts, ind)) then
c
c     Not a recognized option ... the name of the basis or an error
c
            if ((name.ne.' ') .or. (inp_cur_field() .ne. 2)) then
               write(6,*)
     &            ' bas_input: basis/ecp name must be first option'
               goto 10000
            endif
            name = test
            goto 10
         endif
c
         goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000) ind
         goto 10000
c
c     spherical
c
 100     ospherical = .true.
         goto 10
c
c     cartesian
c     
 200     ospherical = .false.
         goto 10
c
c     segment
c
 300     osegment = .true.
         goto 10
c
c     nosegment
c
 400     osegment = .false.
         goto 10
c
c     library
c
 500     if (.not. inp_a(standard)) goto 10000
         goto 10
c
c     file
c
 600     if (.not. inp_a(filename)) goto 10000
         goto 10
c
c     print
c
 700     oprint = .true.
         goto 10
c
c     noprint
c
 800     oprint = .false.
         goto 10
c
c     ecpset "name"
c
00900    o_add_ecpname = .true.
         if (.not. inp_a(ecpname)) goto 10000
         if (oIs_ecp) then 
           errmsg = 'bas_input: you are attaching an ecp to '
           errmsg = errmsg//'an ecp basis.  Not possible'
           call errquit(errmsg, 911)
         endif
         goto 10
c
c print basis version info
c
01000    continue
         if (.not.bas_version()) call errquit
     &       ('bas_input: bas_version failed ',911)
         goto 10
c
      endif
c
c     Now check reality against input
c
c     Open a new basis set to receive the new data
c
      if (oIs_ecp) then
        if (name .eq. ' ') name = 'ecp basis'
      else
        if (name .eq. ' ') name = 'ao basis'
      endif
      if (.not. bas_create(basis, name))
     $     call errquit('bas_input: failed to create basis', 0)
      if (oIs_ecp) then
        if (.not.bas_set_ecp_basis(basis)) call errquit
     &      ('bas_input: ecp_set_basis failed',911)
      endif
      if (o_add_ecpname) then
        write(6,*)' ecp basis set associated is ',ecpname
        if (.not.bas_set_ecp_name(basis,ecpname)) call errquit
     &      ('bas_input: bas_set_ecp_name failed',911)
      endif
c
c     Process standard basis sets directive ... not yet done
c
      if (standard .ne. ' ') call errquit
     $     ('bas_input: specify standard basis sets per tag', 0)
c
c    Now left with reading in from the input additional specifications
c    for basis functions or standard sets on specific tags
c
      if (.not.bas_add_ucnt_init(basis)) call errquit
     &    ('bas_input: failed to init add_ucnt ',911)
      call bas_input_body(basis, osegment, oIs_ecp)
c
      if (.not. bas_set_spherical(basis, ospherical))
     &      call errquit
     &      (' bas_set_spherical failed ',911)
      
c
c     Now have processed the entire basis directive.  Print out
c     info if desired, write it to the data base, tidy up and go home
c
      if (oprint) then
         if (.not. bas_print(basis))
     $        call errquit('bas_input: print failed', 0)
         if (.not.bas_summary_print(basis))
     &        call errquit('scf: basis summary print failed',911)
      endif
c
      if (.not. bas_rtdb_store(rtdb, name, basis)) call errquit
     $     ('bas_input: failed to store basis', 0)
c
      if (.not. bas_do_destroy(basis)) call errquit
     $     ('bas_input: bas_destroy failed', 0)
c
      return
c
10000 write(6,10001)
10001 format(' basis|ecp [<name>] [library <standard set>] \\',/,
     $    '      [file <filename>] [spherical|cartesian] [segment]\\',
     &    '      [print|noprint] [ecp name_of_ecp]')
      call errquit('bas_input: invalid format for basis directive', 0)
c
      end
*.....................................................................
      subroutine bas_input_body(basis, osegment, oIs_ecp)
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "inp.fh"
c::passed
      integer basis             ! [input]
      logical osegment          ! [input]
      logical oIs_ecp            ! [input]
c
c     Read the body of a basis directive that describes the
c     tags/exponents/contraction coefficients 
c
c
c     tag library [<tag in library>] <standard set> [file <filename>]
c     tag <contraction type>
c       <exponent> <contraction coefficients>
c       ...
c     end basis
c
c for ecp
c     tag nelec <integer>      
c     tag library [<tag in library>] <standard set> [file <filename>]
c     tag <contraction type>
c       <r-exponent> <exponent> <contraction coefficients>
c       ...
c     end ecp
c
c
c     
      character*16 tag, tag_in_lib
      character*16 cont_type
      character*255 field, standard, filename
      integer nltypes           ! No. of known angular momentum types
      integer nsptypes          ! No. of known sp type shells
      integer nopts             ! No. of options
      parameter (nltypes = 7, nsptypes = 3, nopts = 3)
      character*1 ltypes(nltypes)
      character*2 sptypes(nsptypes)
      character*8 opts(nopts)
      integer spvalues(nsptypes)
      integer num_elec
      integer l_value, ind
      data ltypes /'s','p','d','f','g','h','i'/
      data sptypes / 'sp', 'l ', 'ul'/
      data spvalues/  -1 , -1 ,   -1 /
      data opts / 'library', 'file', 'nelec' /
c
c     Input a new line
c
 10   if (.not. inp_read()) call errquit
     $     ('bas_input_body: premature EOF', 0)
c
c     Start parsing current line
c
 20   call inp_set_field(0)
      standard = ' '            ! Must reset these for every tag
      filename =
     &BASIS_LIBRARY
      if (.not. inp_a(tag)) goto 10000
c
*. . . . . . . . . . . . . . . . . . . . . . ! End of basis directive
      if (inp_compare(.false., 'end', tag)) goto 9000 
c
      if (.not. inp_a(cont_type)) goto 10000
c     
      if (inp_match(nltypes, .false., cont_type, ltypes, ind)) then
c     
c  The contraction is a simple shell or an ecp projector shell of the 
c  given type
c     
        l_value = ind - 1
      else if (
     &      inp_match(nsptypes, .false., cont_type,sptypes,ind)) then
c     
c     The contraction is an sp-type shell or an ecp Ulocal shell
c     
         l_value = spvalues(ind)
      else
c     
c     It might be an option
c     
         call inp_prev_field
 22      if (inp_a(field)) then
            if (inp_match(nopts, .false., field, opts, ind)) then
               goto (24, 26, 28) ind
c     
c     library
c     
 24            if (.not. inp_a(standard)) goto 10000
	       tag_in_lib = tag
	       if (inp_a(field)) then
                 if (inp_match(nopts, .false., field, opts, ind)) then
		   call inp_prev_field
                 else
		   tag_in_lib = standard
		   standard   = field
		   write(6,*) ' tag_in_lib, standard ', 
     $                tag_in_lib, standard
                 endif
               endif
               goto 22
c     
c     file
c     
 26            if (.not. inp_a(filename)) goto 10000
               goto 22
c
c...  nelec  ...  =>  number of electrons for ecp
c
 28            if (.not.inp_i(num_elec)) goto 10000
*debug:               write(6,*)' basis    = ',basis
*debug:               write(6,*)' tag      = ',tag
*debug:               write(6,*)' num_elec = ',num_elec
               if (.not.ecp_set_num_elec(basis,tag,num_elec,'unknown'))
     &             call errquit
     &             ('bas_input_body: ecp_set_num_elec failed',911)
*. . . . no more input allowed on line so no goto 22 structure needed?
               goto 10
            else
               goto 10000       ! Unknown option
            endif
         endif
      endif
c
c     Line with tag on has been parsed ... either a standard set
c     or explicit input
c
      if (standard .ne. ' ') then
         call bas_tag_lib(basis, osegment, tag, tag_in_lib, standard,
     $        filename, oIs_ecp)
         goto 10
      endif
c
c     Fall thru to here to read in a set of contraction coefficients
c
      call bas_input_cont(basis, osegment, tag, l_value,
     &    'user specified', oIs_ecp)
c
c     Have already read in the next line ... parse it
c
      goto 20
c
c     Have read in all of the basis set/ecp info. 
c     
 9000 return
c
10000 write(6,1)
 1    format(' basis directive body format is:',/,
     $       '       tag library <standard set> [file <filename>]',/,
     $       '       tag <contraction type>',/,
     $       '           <exponent> <contraction coefficients>',/,
     $       '           ... ',/,
     $       ' end basis',/,/,/,'or',/,
     &       ' ecp directive body format is:',/,
     $       '       tag library <standard set> [file <filename>]',/,
     $       '       tag nelec <number of electrons ',
     &            'replaced on tag>',/,
     $       '       tag <contraction type>',/,
     $       '           <exponent> <contraction coefficients>',/,
     $       '           ... ',/,
     $       ' end ecp')
      call errquit('bas_input_body: format error in the input', 0)
c
      end
*.....................................................................
      subroutine bas_input_cont(basis, osegment, tag, l_value, stdtag,
     &    oIs_ecp)
      implicit none
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "inp.fh"
c::functions
      logical bas_add_ucnt
      external bas_add_ucnt
c
      integer basis             ! [input] basis set handle
      logical osegment          ! [input] segment the basis set?
      character*16 tag          ! [input] tag string
      integer l_value           ! [input] l value of contraction
      character*(*) stdtag      ! [input] standard name
      logical oIs_ecp            ! [input] is this an ecp basis
c
      integer ngen, iprim, nprim, i, j, nptmp
      integer cont_max          ! Max. no. of contractions
      integer prim_max          ! Max. no. of primitives
      logical status
      character*1000 errmsg
      logical rex_okay
      logical ostat_f, ostat_field, ostat_split
      parameter (cont_max = 30, prim_max = 30)
      double precision expnt(prim_max), coeff(prim_max,cont_max),
     $     etmp(prim_max), ctmp(prim_max),
     &     rex(prim_max), rex_tmp(prim_max)
c
c
c    The SUBSEQUENT lines contain coefficients and exponents
c    for a contraction ... read them in and add the contraction to the
c    basis set, segmenting if required.  The input file is left
c    positioned having read the end directive
c
      if (.not. inp_read()) goto 10000
c
      if (oIs_ecp) then
        ngen = inp_n_field() - 2   ! rex, expn, cont coeff(s)
      else
        ngen = inp_n_field() - 1   ! expn, cont coeff(s)
      endif
c
      if (ngen .lt. 1) goto 10000
      if (ngen .gt. cont_max) then
        errmsg = 'bas_input_cont: '//
     &        'too many contractions - increase cont_max'
        call errquit(errmsg,cont_max)
      endif
      if ((.not.oIs_ecp).and.(l_value.eq.-1 .and. ngen.ne.2))
     &    call errquit
     &    ('bas_input_cont: sp shell requires exactly 2 coeffs',ngen)
c
      do iprim = 1, prim_max
c
        if (oIs_ecp) then
          ostat_f = inp_f(rex(iprim))
          ostat_f = ostat_f .and. inp_f(expnt(iprim))
        else
          ostat_f = inp_f(expnt(iprim))
        endif
        if (.not. ostat_f) then
c     
c     If cannot read the first field as an exponent then 
c     it is the end of this contraction
c     
          goto 30
        else if (expnt(iprim) .le. 0.0d0) then
          call errquit('bas_input_cont: invalid exponent', 0)
        else
          if (oIs_ecp) then
            ostat_field = (inp_n_field() - 2) .ne. ngen
          else
            ostat_field = (inp_n_field() - 1) .ne. ngen
          endif
          if (Ostat_field) then
            write(luout,*) ' bas_input_cont: no. of coefficients?'
            goto 10000
          endif
          do i= 1, ngen
            if (.not. inp_f(coeff(iprim,i))) then
              write(luout,*) ' bas_input_cont: failed reading ',
     $            'coefficient'
              goto 10000
            else if (coeff(iprim,i) .eq. 0.0d0) then
              write(luout,10100)tag(1:inp_strlen(tag)),
     &            l_value,stdtag(1:inp_strlen(stdtag))
              call inp_outrec()
              write(luout,*)' '
            endif
          enddo
          if (.not. inp_read()) goto 10000
        endif
      enddo
      call errquit
     &    ('bas_input_cont: too many primitives in contraction',
     $     prim_max)
 30   nprim = iprim-1
      if (nprim .le. 0) call errquit('bas_input_cont: no primitives?',
     $     nprim)
c check to make sure any contraction is not all zero's
c
c
      call bas_check_contractions(nprim,ngen,prim_max,coeff,expnt)
c
c   Now have tag, contraction type, no. of contractions, no. of prims,
c   exponents, coeffs.  Shove this lot into the basis set.
c
c
c   bas_add_ucnt -> adds a new general contraction on the specified 
c                   tag.
c
c   If the tag is not present it will also add that.
c
      status = .true.
      if (osegment) then
c     
c     Add contractions one-at-a-time to force segmentation ... remove
c     functions with zero coefficients
c     
        ostat_split = (l_value.ge.0).or.oIs_ecp
        if (ostat_split) then
c     simple shell or all ecp shells
          do i = 1, ngen
            nptmp = 0
            do j = 1, nprim
              if (coeff(j,i).ne.0.0d0) then
                nptmp = nptmp + 1
                ctmp(nptmp) = coeff(j,i)
                etmp(nptmp) = expnt(j)
                if (oIs_ecp) then
                  rex_okay =
     &                ((abs(rex(j)-0.0d00)).lt.1.0d-9).or.
     &                ((abs(rex(j)-1.0d00)).lt.1.0d-9).or.
     &                ((abs(rex(j)-2.0d00)).lt.1.0d-9)
                  if (rex_okay) then
                    rex_tmp(nptmp) = rex(j)
                  else
                    write(luout,*)'component',j,
     &                  ' of the current contraction has ',
     &                  'an r-exponent of',rex(j)
                    call errquit
     &                 ('*ERR* bas_input_cont: Fatal r-exponent error',
     &                 911)
                  endif
                endif
              endif
            enddo
            status = status .and.
     $          bas_add_ucnt(basis, tag, l_value, 1, nptmp, 
     $          rex_tmp, etmp, ctmp, prim_max,
     &          stdtag, oIs_ecp)
          enddo
        else
c     sp shell
          status = status .and.
     $        bas_add_ucnt(basis, tag, 0, 1, nprim, 
     $        rex, expnt, coeff(1,1), prim_max,
     &        stdtag, oIs_ecp)
          status = status .and.
     $        bas_add_ucnt(basis, tag, 1, 1, nprim, 
     $        rex, expnt, coeff(1,2), prim_max,
     &        stdtag, oIs_ecp)
        endif
      else
c     
c     Add as a single general contraction or sp shell
c     
        status = status .and.
     $      bas_add_ucnt(basis, tag, l_value, ngen, nprim, 
     $      rex, expnt, coeff, prim_max,
     &      stdtag, oIs_ecp)
      endif
      if (.not. status) call errquit
     $        ('bas_input_cont: bas_add_ucnt failed!!', 0)
c
      return
c
10000 write(6,1)
 1    format(' basis/ecp contraction format is:',/,
     $       '       tag <contraction type>',/,
     $       '           <exponent> <contraction coefficients>',/,
     $       '           ... ',/,
     $       ' end basis/ecp')
      call errquit('bas_input_body: format error in the input', 0)
10100 format(/,' **** WARNING Zero Coefficient **** on atom "', a16,
     &    '"',/,' angular momentum value:',i2,
     &    '  standard basis set name: "',a,'"',/,
     &    1x,'input line that generated warning:')
c
      end
*.....................................................................
      subroutine bas_tag_lib(basis, osegment, tag, tag_in_lib, 
     $     standard, filename, oIs_ecp)
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "inp.fh"
      integer basis             ! [input] basis handle
      logical osegment          ! [input] segment if true
      character*16 tag          ! [input] tag
      character*16 tag_in_lib   ! [input] tag to lookup in library
      character*(*) standard    ! [input] name of standard basis set
      character*(*) filename    ! [input] name of library file
      logical oIs_ecp           ! [input] ecp library function read?
c
      character*5 id_string
      character*2 symbol
      character*16 element
      character*255 string
      character*4 cont_type
      integer atn, lenn
      integer nltypes           ! No. of known angular momentum types
      integer nsptypes          ! No. of known sp type shells
      parameter (nltypes = 7, nsptypes = 3)
      character*1 ltypes(nltypes)
      character*2 sptypes(nsptypes)
      character*255 field
      logical status
      logical inp_search_fast
      external inp_search_fast
      integer spvalues(nsptypes)
      integer l_value, ind
      integer num_elec
      data ltypes /'s','p','d','f','g','h','i'/
      data sptypes / 'sp', 'l ', 'ul'/
      data spvalues/  -1 , -1  ,  -1 /
c
c     Try to read in a standard basis set for the atom type
c     associated with this tag
c
c     Note that this routine temporarily resets the inp package
c     to read from a different file.
c
c     Translate tag to symbol of atom and assemble a search string of
c     the form 'basis <symbol>_<basis name>'
c      or form 'ecp <symbol>_<basis name>'
c
      if (.not. geom_tag_to_element(tag_in_lib, symbol, element, atn))
     $     call errquit
     &    ('bas_tag_lib: tag does not refer to an atom', 0)

      string = symbol
      lenn = inp_strlen(string)
      string(lenn+1:lenn+1) = '_'
      string(lenn+2:) = standard
      call inp_lcase(string)
c
c     replace all non-trailing spaces with underscore
c
      lenn = inp_strlen(string)
      do ind = 1, lenn
         if (string(ind:ind).eq.' ') string(ind:ind) = '_'
      enddo
c
c     Open the library file and try to locate the basis
c
      if (filename .eq. ' ') filename = 'library'
      open(33, file=filename, status='old', err=1000)

      call inp_save_state()     ! Save state for unit 5.

      call inp_init(33,6)
      if (oIs_ecp) then
        id_string = 'ecp'
      else
        id_string = 'basis'
      endif
 33   if (inp_search_fast(id_string)) then
         status = inp_a(field)
         if (inp_a(field)) then
            if (inp_compare(.false., string, field)) goto 34
         else
            goto 35
         endif
         goto 33
      endif
 35   call inp_restore_state    ! Failure ... restore unit 5 and abort
      call inp_clear_err
      write(6,*) ' bas_tag_lib: failed to locate ', string
      call errquit('bas_tag_lib: no such basis available', 0)
c
 34   continue                  ! success
c
c     Now read the basis for this magical atom type in using
c     the given tag
c
 10   if (.not. inp_read()) call errquit
     $     ('bas_tag_lib: premature EOF', 0)
c     
c     Start parsing current line
c     
 20   call inp_set_field(0)
      if (.not. inp_a(field)) goto 10000 ! tag
      if (inp_compare(.false., 'end', field)) goto 9000 ! End
      if (.not. inp_a(field)) goto 10000 ! cont_type or 'nelec'
      if (inp_compare(.false.,'nelec', field)) then
         if (.not.inp_i(num_elec)) goto 10000
         if (.not.ecp_set_num_elec(basis,tag,num_elec,'unknown'))
     &       call errquit
     &       ('bas_tag_lib: ecp_set_num_elec failed',911)
         if (.not.inp_read())
     &       call errquit
     &       ('bas_tag_lib: premature end of input in library',911)
         goto 20
      else
         call inp_prev_field()
         if (.not. inp_a(cont_type)) goto 10000
      endif
c
      if (inp_match(nltypes, .false., cont_type, ltypes, ind)) then
c     
c     The contraction is a simple shell
c     
         l_value = ind - 1
      else if (
     &       inp_match(nsptypes, .false., cont_type,sptypes,ind)) then
c     
c     The contraction is an sp-type shell
c     
         l_value = spvalues(ind)
      else
         call errquit('bas_tag_lib: invalide shell type?', 0)
      endif
c
c     Have tag and l_value ... read in the contraction coeffcients
c
      call bas_input_cont(basis, osegment, tag, l_value,
     &    standard, oIs_ecp)
      goto 20
c
 9000 close(33)
      call inp_restore_state    ! Restore state for unit 5
c
      return
c
 1000 write(6,*) ' bas_tag_lib: looking for ', tag_in_lib, standard,
     $     filename
      call errquit('bas_tag_lib: failed to open basis library', 0)
c
10000 call errquit('bas_tag_lib: problem with library format?', 0)
c
      end
      subroutine bas_check_contractions(np,ng,ldc,cc,ee)
      implicit none
#include "stdio.fh"
      integer np  ! [input] number of primatives
      integer ng  ! [input] number of general contractions
      integer ldc ! [input] leading dimension of cc
      double precision ee(np)    ! [input] exponents
      double precision cc(ldc,ng) ! [input] coefficients
c
      double precision thresh
      parameter (thresh = 1.0d-15)
      logical okay
      integer ig, ip, nzero
c
*debug:      write(luout,*)'number of primitives          :',np
*debug:      write(luout,*)'number of general contractions:',ng
*debug:      write(luout,*)' input exponents, coeffs'
*debug:      do ip = 1,np
*debug:        write(luout,10000)ee(ip),(cc(ip,ig),ig=1,ng)
*debug:      enddo
*debug:      
      okay = .true.
      do ig = 1,ng
        nzero = 0
        do ip = 1,np
          if (abs(cc(ip,ig)).lt.thresh) nzero = nzero + 1
        enddo
        okay = okay.and.nzero.lt.np
      enddo
      if (okay) return
      write(luout,*) ' bad basis set input :'
      write(luout,*) ' one contraction has all zero coefficients'
      do ip = 1,np
        write(luout,10000)ee(ip),(cc(ip,ig),ig=1,ng)
      enddo
      call errquit('bas_check_contractions: fatal error',911)
10000 format(10(1pd14.6))
      end

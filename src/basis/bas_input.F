      subroutine bas_input(rtdb)
c $Id: bas_input.F,v 1.18 1996-01-19 19:16:06 d3g681 Exp $
      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "inp.fh"
c
c     basis [<name>] [library [<standard set>] [file <filename>] \
c           [spherical|cartesian] [segment||nosegment] [print] [noprint]
c
c     tag library [tag in library] <standard set> [file <filename>]
c     tag <shell type>
c       <exponent> <contraction coefficients>
c       ...
c     end basis
c
c     parse the main directive
c
      integer rtdb              ! [input] handle to database
      integer nopt
      parameter (nopt = 8)
      character*10 opts(nopt)
      character*255 test, name, filename, standard
      logical status, ospherical, osegment, oprint
      integer ind, basis
      data opts /'spherical', 'cartesian', 'segment', 'nosegment',
     $     'library', 'file', 'print', 'noprint'/
c
c     Check is a basis directive and read in name of the basis
c
      call inp_set_field(0)
      status = inp_a(test)
      if ((.not. status) .or.
     $     (.not. inp_compare(.false., test, 'basis'))) goto 10000
c
c     Parse rest of basis directive line
c
      name = ' '
      filename = ' '
      test = ' '
      standard = ' '
      ospherical = .false.      ! Default is cartesian
      osegment   = .true.       ! Default is to force segmentation
      oprint     = .true.       ! Default is to print the basis on input
c
 10   if (inp_a(test)) then
c
         if (.not. inp_match(nopt, .false., test, opts, ind)) then
c
c     Not a recognized option ... the name of the basis or an error
c
            if ((name.ne.' ') .or. (inp_cur_field() .ne. 2)) then
               write(6,*) ' bas_input: basis name must be first option'
               goto 10000
            endif
            name = test
            goto 10
         endif
c
         goto (100, 200, 300, 400, 500, 600, 700, 800) ind
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
      endif
c
c     Now check reality against input
c
      if (ospherical) call errquit
     $     ('bas_input: spherical harmonics not yet', 0)
c
c     Open a new basis set to receive the new data
c
      if (name .eq. ' ') name = 'ao basis'
      if (.not. bas_create(basis, name))
     $     call errquit('bas_input: failed to create basis', 0)
c
c     Process standard basis sets directive ... not yet done
c
      if (standard .ne. ' ') call errquit
     $     ('bas_input: specify standard basis sets per tag', 0)
c
c     Now left with reading in from the input additional specifications
c     for basis functions or standard sets on specific tags
c
      call bas_input_body(basis, osegment)
c
c     Now have processed the entire basis directive.  Print out
c     info if desired, write it to the data base, tidy up and go home
c
      if (oprint) then
         if (.not. bas_print(basis))
     $        call errquit('bas_input: print failed', 0)
      endif
c
      if (.not. bas_rtdb_store(rtdb, name, basis)) call errquit
     $     ('bas_input: failed to store basis', 0)
c
      if (.not. bas_destroy(basis)) call errquit
     $     ('bas_input: bas_destroy failed', 0)
c
      return
c
10000 write(6,10001)
10001 format(' basis [<name>] [library <standard set>] \\'/
     $       '       [file <filename>] [spherical|cartesian] [segment]')
      call errquit('bas_input: invalid format for basis directive', 0)
c
      end
      subroutine bas_input_body(basis, osegment)
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "inp.fh"
      integer basis             ! [input]
      logical osegment          ! [input]
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
      character*16 tag, tag_in_lib
      character*16 cont_type
      character*255 field, standard, filename
      integer nltypes           ! No. of known angular momentum types
      integer nsptypes          ! No. of known sp type shells
      integer nopts             ! No. of options
      parameter (nltypes = 7, nsptypes = 2, nopts = 2)
      character*1 ltypes(nltypes)
      character*2 sptypes(nsptypes)
      character*8 opts(nopts)
      integer spvalues(nsptypes)
      integer l_value, ind
      data ltypes /'s','p','d','f','g','h','i'/
      data sptypes / 'sp', 'l '/
      data spvalues/  -1 , -1 /
      data opts / 'library', 'file' /
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
      if (inp_compare(.false., 'end', tag)) goto 9000 ! End of basis directive
c
      if (.not. inp_a(cont_type)) goto 10000
c     
      if (inp_match(nltypes, .false., cont_type, ltypes, ind)) then
c     
c     The contraction is a simple shell
c     
         l_value = ind - 1
      else if (inp_match(nsptypes, .false., cont_type,sptypes,ind)) then
c     
c     The contraction is an sp-type shell
c     
         l_value = spvalues(ind)
      else
c     
c     It might be an option
c     
         call inp_prev_field
 22      if (inp_a(field)) then
            if (inp_match(nopts, .false., field, opts, ind)) then
               goto (24, 26) ind
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
         call bas_tag_lib(basis,osegment,tag_in_lib,standard,filename)
         goto 10
      endif
c
c     Fall thru to here to read in a set of contraction coefficients
c
      call bas_input_cont(basis, osegment, tag, l_value)
c
c     Have already read in the next line ... parse it
c
      goto 20
c
c     Have read in all of the basis set info. 
c     
 9000 return
c
10000 write(6,1)
 1    format(' basis directive body format is:'/
     $       '       tag library <standard set> [file <filename>]'/
     $       '       tag <contraction type>'/
     $       '           <exponent> <contraction coefficients>'/
     $       '           ... '/
     $       ' end basis')
      call errquit('bas_input_body: format error in the input', 0)
c
      end
      subroutine bas_input_cont(basis, osegment, tag, l_value)
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "inp.fh"
c::functions
      logical bas_add_ucnt
      external bas_add_ucnt
c
      integer basis             ! [input]
      logical osegment          ! [input]
      character*16 tag          ! [input]
      integer l_value           ! [input]
c
      integer ngen, iprim, nprim, i, j, nptmp
      integer cont_max          ! Max. no. of contractions
      integer prim_max          ! Max. no. of primitives
      logical status
      parameter (cont_max = 30, prim_max = 30)
      double precision expnt(prim_max), coeff(prim_max,cont_max),
     $     etmp(prim_max), ctmp(prim_max)
c
c     The SUBSEQUENT lines contain coefficients and exponents
c     for a contraction ... read them in and add the contraction to the
c     basis set, segmenting if required.  The input file is left
c     positioned having read the end directive
c
      if (.not. inp_read()) goto 10000
c
      ngen = inp_n_field() - 1
      if (ngen .lt. 1) goto 10000
      if (ngen .gt. cont_max) call errquit
     $     ('bas_input_cont: too many contractions - increase cont_max',
     $     cont_max)
      if (l_value.eq.-1 .and. ngen.ne.2) call errquit
     $     ('bas_input_cont: sp shell requires exactly 2 coeffs',ngen)
c
      do iprim = 1, prim_max
c
         if (.not. inp_f(expnt(iprim))) then
c     
c     If cannot read the first field as an exponent then 
c     it is the end of this contraction
c     
            goto 30
         else if (expnt(iprim) .le. 0.0d0) then
            call errquit('bas_input_cont: invalid exponent', 0)
         else
            if ((inp_n_field()-1) .ne. ngen) then
               write(6,*) ' bas_input_cont: no. of coefficients?'
               goto 10000
            endif
            do i= 1, ngen
               if (.not. inp_f(coeff(iprim,i))) then
                  write(6,*) ' bas_input_cont: failed reading ',
     $                 'coefficient'
                  goto 10000
               else if (coeff(iprim,i) .eq. 0.0d0) then
                 write(6,*)' **** WARNING Zero Coefficient **** '
*                call errquit('bas_input_cont: invalid coefficient', 0)
               endif
            enddo
            if (.not. inp_read()) goto 10000
         endif
      enddo
      call errquit('bas_input_cont: too many primitives in contraction',
     $     prim_max)
 30   nprim = iprim-1
      if (nprim .le. 0) call errquit('bas_input_cont: no primitives?',
     $     nprim)
c
c     Now have tag, contraction type, no. of contractions, no. of prims,
c     exponents, coeffs.  Shove this lot into the basis set.
c
c     bas_add_ucnt -> adds a new general contraction on the specified tag.
c     If the tag is not present it will also add that.
c
      status = .true.
      if (osegment) then
c     
c     Add contractions one-at-a-time to force segmentation ... remove
c     functions with zero coefficients
c     
         if (l_value .ne. -1) then
c     simple shell
            do i = 1, ngen
               nptmp = 0
               do j = 1, nprim
                  if (coeff(j,i).ne.0.0d0) then
                     nptmp = nptmp + 1
                     ctmp(nptmp) = coeff(j,i)
                     etmp(nptmp) = expnt(j)
                  endif
               enddo
               status = status .and.
     $              bas_add_ucnt(basis, tag, l_value, 1, nptmp, 
     $              etmp, ctmp, prim_max)
            enddo
         else
c     sp shell
            status = status .and.
     $           bas_add_ucnt(basis, tag, 0, 1, nprim, 
     $           expnt, coeff(1,1), prim_max)
            status = status .and.
     $           bas_add_ucnt(basis, tag, 1, 1, nprim, 
     $           expnt, coeff(1,2), prim_max)
         endif
      else
c     
c     Add as a single general contraction or sp shell
c     
         status = status .and.
     $        bas_add_ucnt(basis, tag, l_value, ngen, nprim, 
     $        expnt, coeff, prim_max)
      endif
      if (.not. status) call errquit
     $        ('bas_input_cont: bas_add_ucnt failed!!', 0)
c
      return
c
10000 write(6,1)
 1    format(' basis contraction format is:'/
     $       '       tag <contraction type>'/
     $       '           <exponent> <contraction coefficients>'/
     $       '           ... '/
     $       ' end basis')
      call errquit('bas_input_body: format error in the input', 0)
c
      end
      subroutine bas_tag_lib(basis, osegment, tag, standard, filename)
      implicit none
#include "bas.fh"
#include "geom.fh"
#include "inp.fh"
      integer basis             ! [input] basis handle
      logical osegment          ! [input] segment if true
      character*16 tag          ! [input] tag
      character*(*) standard    ! [input] name of standard basis set
      character*(*) filename    ! [input] name of library file
c
      character*2 symbol
      character*16 element
      character*30 string
      character*4 cont_type
      integer atn, lenn
      integer nltypes           ! No. of known angular momentum types
      integer nsptypes          ! No. of known sp type shells
      parameter (nltypes = 7, nsptypes = 2)
      character*1 ltypes(nltypes)
      character*2 sptypes(nsptypes)
      character*255 field
      logical status
      integer spvalues(nsptypes)
      integer l_value, ind
      data ltypes /'s','p','d','f','g','h','i'/
      data sptypes / 'sp', 'l '/
      data spvalues/  -1 , -1 /
c
c     Try to read in a standard basis set for the atom type
c     associated with this tag
c
c     Note that this routine temporarily resets the inp package
c     to read from a different file.  This has the effect of messing
c     up line counts and discarding the remainder of the current
c     input line.
c
c     Translate tag to symbol of atom and assemble a search string of
c     the form 'basis <symbol>_<basis name>'
c
      if (.not. geom_tag_to_element(tag, symbol, element, atn))
     $     call errquit('bas_tag_lib: tag does not refer to an atom', 0)
      string = symbol
      lenn = inp_strlen(string)
      string(lenn+1:lenn+1) = '_'
      string(lenn+2:) = standard
      call inp_lcase(string)
c
c     Open the library file and try to locate the basis
c
      if (filename .eq. ' ') filename = 'library'
      open(33, file=filename, status='old', err=1000)

      call inp_save_state()     ! Save state for unit 5.

      call inp_init(33,6)
 33   if (inp_search(.false., 'basis')) then
         status = inp_a(field)
         if (inp_a(field)) then
            if (inp_compare(.false., string, field)) goto 34
         else
            goto 35
         endif
         goto 33
      endif
 35   write(6,*) ' bas_tag_lib: failed to locate ', string
      call errquit('bas_tag_lib: no such basis available', 0)
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
      if (.not. inp_a(field)) goto 10000
c     
      if (inp_compare(.false., 'end', field)) goto 9000 ! End
      if (.not. inp_a(cont_type)) goto 10000
      if (inp_match(nltypes, .false., cont_type, ltypes, ind)) then
c     
c     The contraction is a simple shell
c     
         l_value = ind - 1
      else if (inp_match(nsptypes, .false., cont_type,sptypes,ind)) then
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
      call bas_input_cont(basis, osegment, tag, l_value)
      goto 20
c
 9000 close(33)
      call inp_restore_state    ! Restore state for unit 5
      return
c
 1000 write(6,*) ' bas_tag_lib: looking for ', tag, standard, filename
      call errquit('bas_tag_lib: failed to open basis library', 0)
c
10000 call errquit('bas_tag_lib: problem with library format?', 0)
c
      end

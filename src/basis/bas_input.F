      subroutine bas_input(rtdb)
      implicit none
#include "rtdb.fh"
#include "context.fh"
#include "geom.fh"
#include "mafdecls.h"
#include "bas.fh"
#include "inp.fh"
c
c     basis [<name>] [library <standard set>] [file <filename>] \
c           [spherical|cartesian] [segment] [print]
c
c     tag library <standard set> [file <filename>]
c     tag <shell type>
c       <exponent> <contraction coefficients>
c       ...
c     end basis
c
c     parse the main directive
c
      integer rtdb              ! [input] handle to database
      integer nopt
      parameter (nopt = 6)
      character*10 opts(nopt)
      character*255 test, name, filename, standard
      logical status, ospherical, osegment, oprint
      integer ind, basis
      data opts /'spherical', 'cartesian', 'segment', 'library', 'file',
     $     'print'/
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
      osegment   = .false.      ! Default is to NOT force segmentation
      oprint     = .false.      ! Default is NOT to print basis on input
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
         goto (100, 200, 300, 400, 500, 600) ind
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
c     library
c
 400     if (.not. inp_a(standard)) goto 10000
         goto 10
c
c     file
c
 500     if (.not. inp_a(filename)) goto 10000
         goto 10
c
c     print
c
 600     oprint = .true.
         goto 10
c
      endif
c
c     Now check reality against input
c
*      if (.not. osegment) then
*         write(6,*) ' Integral code requires segmentation'
*         osegment = .true.
*      endif
      if (standard .ne. ' ') then
         write(6,*) ' Standard basis ', standard(1:inp_strlen(standard))
         call errquit('bas_input: standard basis set not yet' ,0)
      endif
      if (ospherical) call errquit
     $     ('bas_input: spherical harmonics not yet', 0)
c
c     Open a new basis set to receive the new data
c
      if (name .eq. ' ') name = 'mo basis'
      if (.not. bas_create(basis, name))
     $     call errquit('bas_input: failed to create basis', 0)
c
c     Here will soon process reading standard basis sets
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
#include "inp.fh"
#include "bas.fh"
#include "geom.fh"
      integer basis             ! [input]
      logical osegment          ! [input]
c
c     Read the body of a basis directive that describes the
c     tags/exponents/contraction coefficients 
c
c
c     tag library <standard set> [file <filename>]
c     tag <contraction type>
c       <exponent> <contraction coefficients>
c       ...
c     end basis
c     
      character*16 tag
      character*16 cont_type
      integer nltypes           ! No. of known angular momentum types
      integer nsptypes          ! No. of known sp type shells
      integer nopts             ! No. of options
      integer cont_max          ! Max. no. of contractions
      integer prim_max          ! Max. no. of primitives
      parameter (nltypes = 7, nsptypes = 2, nopts = 2)
      parameter (cont_max = 30, prim_max = 30)
      double precision expnt(prim_max), coeff(prim_max,cont_max)
      character*1 ltypes(nltypes)
      character*2 sptypes(nsptypes)
      character*8 opts(nopts)
      integer spvalues(nsptypes)
      integer l_value, ngen, iprim, nprim, i, ind
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
      else if (inp_match(nopts, .false., cont_type, opts, ind)) then
c     
c     It is actually an option to input a standard basis
c     
c     Don't bother parsing this yet
c     
         call errquit('bas_input_body: no standard basis sets yet',0)
c
         goto 10                ! Process the next input line
c     
      else
c     
c     Only god and the user knows what was intended
c     
         goto 10000
      endif
c
c     Fall thru to here to read in a set of contraction coefficients
c
      if (.not. inp_read()) goto 10000
c
      ngen = inp_n_field() - 1
      if (ngen .lt. 1) goto 10000
      if (ngen .gt. cont_max) call errquit
     $     ('bas_input_body: too many contractions - increase cont_max',
     $     cont_max)
c
      do iprim = 1, prim_max
c
         if (.not. inp_f(expnt(iprim))) then
c     
c     If cannot read the first field as an exponent then 
c     it is the end of this contraction
c     
            goto 30
         else
            if ((inp_n_field()-1) .ne. ngen) then
               write(6,*) ' bas_input_body: no. of coefficients?'
               goto 10000
            endif
            do i= 1, ngen
               if (.not. inp_f(coeff(iprim,i))) then
                  write(6,*) ' bas_input_body: failed reading ',
     $                 'coefficient'
                  goto 10000
               endif
            enddo
            if (.not. inp_read()) goto 10000
         endif
      enddo
      call errquit('bas_input_body: too many primitives in contraction',
     $     prim_max)
 30   nprim = iprim-1
      if (nprim .le. 0) call errquit('bas_input_body: no primitives?',
     $     nprim)
c
c     Now have tag, contraction type, no. of contractions, no. of prims,
c     exponents, coeffs.  Shove this lot into the basis set.
c
c     bas_add_ucnt -> adds a new general contraction on the specified tag.
c     If the tag is not present it will also add that.
c
      if (osegment) then
c
c     Add contractions one-at-a-time to force segmentation
c
         do i = 1, ngen
            if (.not. bas_add_ucnt(basis, tag, l_value, 1, nprim, 
     $           expnt, coeff(1,i), prim_max)) call errquit
     $           ('bas_input_body: bas_add_ucnt failed!!', 0)
         enddo
      else
c
c     Add as a single general contraction
c
         if (.not. bas_add_ucnt(basis, tag, l_value, ngen, nprim, 
     $        expnt, coeff, prim_max)) call errquit
     $        ('bas_input_body: bas_add_ucnt failed!!', 0)
      endif
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



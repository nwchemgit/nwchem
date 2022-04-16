      subroutine bas_input(rtdb)
c $Id$
      implicit none
#include "errquit.fh"
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
c         [ecpset ecp_name] [soset so_name] [rel] [bse]
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
c     tag <shell type> <channel>
c       <r-exponent> <exponent> <contraction coefficients>
c       ...
c     end ecp
c
c also used for so potential input
c
c   so [<name>] [library [<standard set>] [file <filename>] \
c       [spherical|cartesian] [segment||nosegment] [print] [noprint]
c
c     tag nelec <integer>
c     tag library [tag in library] <standard set> [file <filename>]
c     tag <shell type>
c       <r-exponent> <exponent> <contraction coefficients>
c       ...
c     end so
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
      parameter (nopt = 13)
      character*10 opts(nopt)
      character*255 test, name, filename, standard
      character*255 ecpname, soname
      character*256 mydir
      character*1000 errmsg
      logical status, ospherical, osegment, oprint
      logical o_add_ecpname, o_add_soname
      logical oHas_Star
      logical oIs_ecp, oIs_so, oIs_rel
      integer ind, basis
      logical bas_add_ucnt_init, bas_set_spherical, bas_add_ucnt_tidy
      external bas_add_ucnt_init, bas_set_spherical
      data opts /
     &    'spherical', 'cartesian', 'segment', 'nosegment', 'library',
     &    'file',      'print',     'noprint', 'ecpset',    'soset',
     &    'version', 'rel','bse'/
c
      ecpname = '                                                '
      soname = '                                                '
      mydir = ' '
c
c     Check is a basis/ecp directive and read in name of the basis
c
      call inp_set_field(0)
      status = inp_a(test)
      oIs_ecp = inp_compare(.false.,test,'ecp')
      oIs_so = inp_compare(.false.,test,'so')
      status = status.and.
     &    (inp_compare(.false.,test,'basis') .or.
     &     oIs_ecp .or. oIs_so)
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
      o_add_ecpname = .false. ! Default is to not associate an ecp name
      o_add_soname = .false.  ! Default is to not associate a so name
      oIs_rel = .false.    ! Default is nonrelativistic
      oHas_Star = .false.  ! Default, not using star-tags 
c
 10   if (inp_a(test)) then
c
         if (.not. inp_match(nopt, .false., test, opts, ind)) then
c
c     Not a recognized option ... the name of the basis or an error
c
            if ((name.ne.' ') .or. (inp_cur_field() .ne. 2)) then
               write(LuOut,*)
     &            ' bas_input: basis/ecp/so name must be first option'
               goto 10000
            endif
            name = test
            goto 10
         endif
c
         goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
     &         1100,1200,1300) ind
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
         call errquit('bas_input: nosegment option disabled',
     S        0,BASIS_ERR)
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
         if (oIs_ecp.or.oIs_so) then 
           errmsg = 'bas_input: you are attaching an ecp to '
     $           // 'an so or ecp basis.  Not possible'
           call errquit(errmsg, 911, INPUT_ERR)
         endif
         goto 10
c
c     soset "name"
c
01000    o_add_soname = .true.
         if (.not. inp_a(soname)) goto 10000
         if (oIs_so.or.oIs_ecp) then 
           errmsg = 'bas_input: you are attaching an so to '
     $           // 'an so or ecp basis.  Not possible'
           call errquit(errmsg, 911, INPUT_ERR)
         endif
         goto 10
c
c print basis version info
c
01100    continue
         if (.not.bas_version()) call errquit
     &       ('bas_input: bas_version failed ',911, INPUT_ERR)
         goto 10
c
c set entire basis set to be relativistic
c
01200    continue
         oIs_rel = .true.
         if (oIs_so.or.oIs_ecp) then 
           errmsg = 'bas_input: rel keyword incompatible with '
     $           // 'so or ecp basis.'
           call errquit(errmsg, 911, INPUT_ERR)
         endif
         goto 10
c
c use libraries.bse
c
01300    continue
         mydir='libraries.bse'
         goto 10
c
      endif
c
c initialize the default library name based on rules in bas_library_file
c the file parameter can overide this setting for a given standard basis set
c 
      call bas_set_library_name(mydir)
c
      
c
c     Now check reality against input
c
c     Open a new basis set to receive the new data
c
      if (oIs_ecp) then
        if (name .eq. ' ') name = 'ecp basis'
      else if (oIs_so) then
        if (name .eq. ' ') name = 'so potential'
      else
        if (name .eq. ' ') name = 'ao basis'
      endif
      if (.not. bas_create(basis, name))
     $     call errquit('bas_input: failed to create basis', 0,
     &       BASIS_ERR)
      if (oIs_ecp) then
        if (.not.bas_set_ecp_basis(basis)) call errquit
     &      ('bas_input: ecp_set_basis failed',911, BASIS_ERR)
      endif
      if (oIs_so) then
        if (.not.bas_set_so_basis(basis)) call errquit
     &      ('bas_input: so_set_basis failed',911, BASIS_ERR)
      endif
      if (o_add_ecpname) then
        write(LuOut,*)' ecp basis set associated is ',ecpname
        if (.not.bas_set_ecp_name(basis,ecpname)) call errquit
     &      ('bas_input: bas_set_ecp_name failed',911,
     &       BASIS_ERR)
      endif
      if (o_add_soname) then
        write(LuOut,*)' so basis set associated is ',soname
        if (.not.bas_set_so_name(basis,soname)) call errquit
     &      ('bas_input: bas_set_so_name failed',911, BASIS_ERR)
      endif
c
c     Process standard basis sets directive ... not yet done
c
      if (standard .ne. ' ') call errquit
     $     ('bas_input: specify standard basis sets per tag', 0,
     &       INPUT_ERR)
c
c    Now left with reading in from the input additional specifications
c    for basis functions or standard sets on specific tags
c
      if (.not.bas_add_ucnt_init(basis)) call errquit
     &    ('bas_input: failed to init add_ucnt ',911,
     &       BASIS_ERR)
      if (.not. bas_set_spherical(basis, ospherical))
     &      call errquit
     &      (' bas_set_spherical failed ',911, BASIS_ERR)
      
      call bas_input_body(basis, osegment, oIs_rel, oHas_Star)
c
c
c     Now have processed the entire basis directive.  Print out
c     info if desired, write it to the data base, tidy up and go home
c
c     Always print summary. Only print full basis when it does not
c     contain any star tags. Do store the info that the basis should be 
c     printed on rtdb, so taht it can be printed at the bas_rtdb_load 
c     stage
c
      if (oprint) then
         test = 'basisprint:'//name(1:inp_strlen(name))
         if (.not. rtdb_put(rtdb,test,mt_log,1,oHas_Star))
     $          call errquit('bas_input: rtdb_put * failed', 0,
     &       RTDB_ERR)
         if (.not. oHas_Star) then
            if (.not. bas_print(basis))
     $          call errquit('bas_input: print failed', 0, BASIS_ERR)
         endif
         if (.not.bas_summary_print(basis))
     &      call errquit('scf: basis summary print failed',911,
     &       BASIS_ERR)
      endif
c
      if (.not. bas_rtdb_store(rtdb, name, basis)) call errquit
     $     ('bas_input: failed to store basis', 0, BASIS_ERR)
c
      if (.not. bas_add_ucnt_tidy(basis))
     $     call errquit('bas_input: tidy failed',0, BASIS_ERR)
c
      if (.not. bas_do_destroy(basis)) call errquit
     $     ('bas_input: bas_destroy failed', 0, BASIS_ERR)
c
      return
c
10000 write(LuOut,10001)
10001 format(' basis|ecp|so [<name>] [library <standard set>] \\',/,
     $    '      [file <filename>] [segment]\\',
     $    '      [spherical|cartesian]\\',
     &    '      [print|noprint] [rel] [ecpset name_of_ecp]',
     &    ' [soset name_of_so_potential] [bse]')
      call errquit('bas_input: invalid format for basis directive', 0,
     &       INPUT_ERR)
c
      end
*.....................................................................
      subroutine bas_input_body(basis, osegment, oIs_rel, oHas_Star)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "inp.fh"
#include "baslibraryP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "ecpso_decP.fh"
#include "stdio.fh"
#include "bas_starP.fh"
#include "util.fh"
c::passed
      integer basis             ! [input]
      logical osegment          ! [input]
      logical oHas_Star         ! [output] true if basis input contains star-tags
c
c     Read the body of a basis directive that describes the
c     tags/exponents/contraction coefficients 
c
c
c     tag library [<tag in library>] <standard set> [file <filename>]
c         [rel]
c     tag <contraction type> [rel]
c       <exponent> <contraction coefficients>
c       ...
c     end basis
c
c for ecp or so 
c     tag nelec <integer>      ! only for ecp
c     tag library [<tag in library>] <standard set> [file <filename>]
c     tag <contraction type> <channel>
c       <r-exponent> <exponent> <contraction coefficients>
c       ...
c     end ecp
c
c     Returning oIs_star to indicate that the basis had a star tag 
c     for printing purposes
c     
      character*16 tag, tag_in_lib
      character*16 type_of_cont, rel_shell
      character*255 field, standard, filename
      integer nltypes           ! No. of known angular momentum types
      integer nsptypes          ! No. of known sp type shells
      integer nopts             ! No. of options
      parameter (nltypes = 8, nsptypes = 3, nopts = 6)
      character*1 ltypes(nltypes)
      character*2 sptypes(nsptypes)
      character*8 opts(nopts)
      logical oIs_ecp            
      logical oIs_so
      logical oIs_rel
      logical oshell_is_rel
      logical oIs_star
      integer spvalues(nsptypes)
      integer num_elec
      integer l_value, ind, icount
cc AJL/Begin/SPIN-POLARISED ECPs
      integer channel           ! Both (Default)=0; Alpha=1; Beta=2
cc AJL/End
      data ltypes /'s','p','d','f','g','h','i','k'/
      data sptypes / 'sp', 'l ', 'ul'/
      data spvalues/  -1 , -1 ,   -1 /
      data opts / 'except', 'library', 'file', 'rel', 'nelec' ,
     A     'details' /
c
#include "ecpso_sfnP.fh"
c
      star_details=.false.
c
c     if user_library_name ends in /libraries.bse/ star_details=t
c
      if(user_library_name(
     /     (inp_strlen(user_library_name)-14):
     /     inp_strlen(user_library_name)).eq.
     /     '/libraries.bse/') then
         star_details=.true.
      endif
      oIs_ecp = Is_ECP_in(basis)
      oIs_so  = Is_SO_in(basis)
c
c     Reset data structures for reading '*' tags
c
      star_nr_tags = 0
      star_tot_excpt = 0
      do icount = 1, max_star_tag
         star_nr_excpt(icount) = 0
      enddo
c
c     Input a new line
c
 10   if (.not. inp_read()) call errquit
     $     ('bas_input_body: premature EOF', 0, INPUT_ERR)
c
c     Start parsing current line
c
 20   call inp_set_field(0)
      standard = ' '            ! Must reset these for every tag
      filename = user_library_name
      oshell_is_rel = oIs_rel
      oIs_star = .false.
      if (.not. inp_a(tag)) goto 10000
c
c  Check if we have a tag that contains a '*'. If so, we have to store
c  the data "as is" in the RTDB until we know the full geometry without 
c  processing
c
      ind = inp_strlen(tag)
      if (tag(ind:ind) .eq. '*') then
        oIs_star = .true.
        star_nr_tags = star_nr_tags + 1
        star_tag(star_nr_tags) = tag
        if (star_nr_tags .gt. max_star_tag) 
     &     call errquit('bas_input_body: max_star_tag limit',211,
     &       INPUT_ERR)
      endif
c
*. . . . . . . . . . . . . . . . . . . . . . ! End of basis directive
      if (inp_compare(.false., 'end', tag)) goto 9000 
c
      if (.not. inp_a(type_of_cont)) goto 10000
c     
      if (inp_match(nltypes, .false., type_of_cont, ltypes, ind)) then
c     
c  The contraction is a simple shell or an ecp projector shell of the 
c  given type
c  Check if we do a * basis, if so bail out as we only do this for 
c  basis sets read from the library
c     
        if (oIs_star) goto 10000
c
        l_value = ind - 1

CC AJL/Begin/SPIN-POLARISED ECPs
c
c     Check if it has a channel specified for spin-polarised ECPs
c     This is only implemented for user-defined values, not Library
c
        if (oIs_ecp) then
          if (.not.inp_i(channel)) then
            channel = 0 
          else if (channel.lt.0.or.channel.gt.2) then
            call errquit('bas_input_body: channel choice for ECP 
     &        must be either 0: Both (default); 1: Alpha; 2: Beta',
     &        99, INPUT_ERR)
          end if
        else
          channel = 0
        end if
c
CC AJL/End

c
c     check if it is a relativistic shell
c
        if (inp_a(rel_shell)) then
          if (.not. inp_compare(.false., 'rel', rel_shell)) call errquit
     &        ('bas_input_body: option on simple shell must be "rel"',
     &        99, INPUT_ERR)
          oshell_is_rel = .true.
        end if
      else if (
     &      inp_match(nsptypes, .false., type_of_cont,sptypes,ind)) then
c     
c  The contraction is an sp-type shell or an ecp Ulocal shell
c     
        l_value = spvalues(ind)

CC AJL/Begin
c
c     check if it has a channel specified for spin-polarised ECPs
c
        if (oIs_ecp) then
          if (.not.inp_i(channel)) then
            channel = 0 
          else if (channel.lt.0.or.channel.gt.2) then 
            call errquit('bas_input_body: channel choice for ECP 
     &        must be either 0: Both (default); 1: Alpha; 2: Beta',
     &        99, INPUT_ERR)
          end if
        else
          channel = 0
        end if
c
CC AJL/End

c
c     check for relativistic shell - this is illegal.
c
        if (oIs_rel) call errquit
     &      ('bas_input_body: sp shell illegal with relativistic basis',
     &      99, INPUT_ERR)
        if (inp_a(rel_shell)) then
          if (inp_compare(.false., 'rel', rel_shell)) call errquit
     &        ('bas_input_body: sp shell cannot be "rel"', 99,
     &       INPUT_ERR)
          call inp_prev_field
        end if
      else
c     
c  It might be an option
c     
         call inp_prev_field
 22      if (inp_a(field)) then
            if (inp_match(nopts, .false., field, opts, ind)) then
               goto (23, 25, 26, 27, 28, 101) ind
c     
c     except
c     
 23            if (.not. oIs_star) goto 10000
 24            if (inp_a(field)) then
                 if (inp_match(nopts, .false., field, opts, ind)) then
                   call inp_prev_field
                   star_nr_excpt(star_nr_tags) = star_tot_excpt
                   goto 22
                 else
                   star_tot_excpt = star_tot_excpt + 1
                   if (star_tot_excpt .gt. max_star_excpt) call 
     &                errquit('bas_input_body: max_star_excpt',211,
     &       INPUT_ERR)
                   star_excpt(star_tot_excpt) = field
                   goto 24
                 endif
               endif
               star_nr_excpt(star_nr_tags) = star_tot_excpt
               goto 22
c     
c     library
c     
 25            if (.not. inp_a(standard)) goto 10000
               tag_in_lib = tag
               if (inp_a(field)) then
                 if (inp_match(nopts, .false., field, opts, ind)) then
                   call inp_prev_field
                 else
                   tag_in_lib = standard
                   standard   = field
*                   write(LuOut,*) ' tag_in_lib, standard ', 
*     $                tag_in_lib, standard
                 endif
               endif
               goto 22
c     
c     file
c     
 26            if (.not. inp_a(filename)) goto 10000
               goto 22
c     
c     relativistic shell
c    
 27            oshell_is_rel = .true.
               goto 22
c
c     details, ie read spherical from library
c
 101           star_details = .true.
               goto 22
c
c...  nelec  ...  =>  number of electrons for ecp
c
 28            if (.not.inp_i(num_elec)) goto 10000
*debug:               write(LuOut,*)' basis    = ',basis
*debug:               write(LuOut,*)' tag      = ',tag
*debug:               write(LuOut,*)' num_elec = ',num_elec
               if (.not.ecp_set_num_elec(basis,tag,num_elec,'unknown'))
     &             call errquit
     &             ('bas_input_body: ecp_set_num_elec failed',911,
     &       INPUT_ERR)
               if (oshell_is_rel) call errquit
     &             ('bas_input_body: rel keyword incompatible with ECP',
     &             99, INPUT_ERR)
*. . . . no more input allowed on line so no goto 22 structure needed?
               goto 10
            else
               goto 10000       ! Unknown option
            endif
         endif
      endif
c
c     Line with tag on has been parsed ... either a standard set
c     or explicit input or a *-tag
c
      if (oIs_star) then
         if (tag_in_lib .ne. ' ' .and. standard .ne. ' ') then
            star_in_lib(star_nr_tags) = tag_in_lib
            star_bas_typ(star_nr_tags) = standard
         else
            goto 10000
         endif
         star_file(star_nr_tags) = filename
         star_rel(star_nr_tags) = oshell_is_rel
         star_segment = osegment
         oHas_Star = .true.
         goto 10
      endif
c
      if (standard .ne. ' ') then
         call bas_tag_lib(basis, osegment, tag, tag_in_lib, standard,
     $        filename, oshell_is_rel, star_details)
         goto 10
      endif
c
c     Fall thru to here to read in a set of contraction coefficients
c
      call bas_input_cont(basis, osegment, tag, l_value,
     &    'user specified', oshell_is_rel, channel)
c
c     Have already read in the next line ... parse it
c
      goto 20
c
c     Have read in all of the basis set/ecp info. 
c     
 9000 return
c
10000 write(LuOut,1)
 1    format(' basis directive body format is:',/,
     $'       tag library [details] <standard set> [file <filename>]',/,
     $       '                   [except <exceptions list>] [rel]',/,
     $       '       tag <contraction type> [rel]',/,
     $       '           <exponent> <contraction coefficients>',/,
     $       '           ... ',/,
     $       ' end basis',/,/,/,'or',/,
     &       ' ecp directive body format is:',/,
     $       '       tag library <standard set> [file <filename>]',/,
     $       '                   [except <exceptions list>] [rel]',/,
     $       '       tag nelec <number of electrons ',
     &            'replaced on tag>',/,
     $       '       tag <contraction type> <channel>',/,
     $       '           <exponent> <contraction coefficients>',/,
     $       '           ... ',/,
     $       ' end ecp',/,/,/,'or',/,
     &       ' so directive body format is:',/,
     $       '       tag library <standard set> [file <filename>]',/,
     $       '                   [except <exceptions list>] [rel]',/,
     $       '       tag <contraction type>',/,
     $       '           <exponent> <contraction coefficients>',/,
     $       '           ... ',/,
     $       ' end so')
      call errquit('bas_input_body: format error in the input', 0,
     &       INPUT_ERR)
c
      end
*.....................................................................
      subroutine bas_input_cont(basis, osegment, tag, l_value, stdtag,
     &    oshell_is_rel, channel)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basdeclsP.fh"
#include "basP.fh"
#include "geom.fh"
#include "inp.fh"
#include "ecpso_decP.fh"
#include "util.fh"
c::functions
      logical bas_add_ucnt0
      external bas_add_ucnt0
c
      integer basis             ! [input] basis set handle
      logical osegment          ! [input] segment the basis set?
      character*16 tag          ! [input] tag string
      integer l_value           ! [input] l value of contraction
      character*(*) stdtag      ! [input] standard name
      logical oshell_is_rel     ! [input] flag for relativistic shell
cc AJL/Begin/SPIN-POLARISED-ECPs
      integer channel           ! [input] ecp channel [both/alpha/beta]
cc AJL/End
c
      integer ngen, iprim, nprim, i, j, nptmp
      integer cont_max          ! Max. no. of contractions
      integer prim_max          ! Max. no. of primitives
      logical status
      character*1000 errmsg
      logical oIs_ecp            ! is this an ecp basis
      logical oIs_so             ! is this an so basis
      logical o_basis
      logical rex_okay
      logical ostat_f, ostat_field, ostat_split
      parameter (cont_max = nw_max_gen_con, prim_max = nw_max_prim)
      double precision expnt(prim_max), coeff(prim_max,cont_max),
     $     etmp(prim_max), ctmp(prim_max),
     &     rex(prim_max), rex_tmp(prim_max)
c
#include "ecpso_sfnP.fh"
c
      oIs_ecp = Is_ECP_in(basis)
      oIs_so  = Is_SO_in(basis)
      o_basis = .not.(oIs_ecp.or.oIs_so)
c
c    The SUBSEQUENT lines contain coefficients and exponents
c    for a contraction ... read them in and add the contraction to the
c    basis set, segmenting if required.  The input file is left
c    positioned having read the end directive
c
      if (.not. inp_read()) goto 10000
c
      if (oIs_ecp.or.oIs_so) then
        ngen = inp_n_field() - 2   ! rex, expn, cont coeff(s)
      else
        ngen = inp_n_field() - 1   ! expn, cont coeff(s)
      endif
c
      if (ngen .lt. 1) goto 10000
      if (ngen .gt. cont_max) then
        errmsg = 'bas_input_cont: '//
     &        'too many contractions - increase cont_max='
        call errquit(errmsg,cont_max, INPUT_ERR)
      endif
      if (o_basis.and.(l_value.eq.-1 .and. ngen.ne.2))
     &    call errquit
     &    ('bas_input_cont: sp shell requires exactly 2 coeffs',ngen,
     &       INPUT_ERR)
c
      do iprim = 1, (prim_max+1)
c
        if (oIs_ecp.or.oIs_so) then
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
        else if (iprim.gt.prim_max) then
          continue  ! drop out of loop without any assignment
        else if (expnt(iprim) .le. 0.0d0) then
          call errquit('bas_input_cont: invalid exponent', 0,
     &       INPUT_ERR)
        else
          if (oIs_ecp.or.oIs_so) then
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
     $     prim_max, BASIS_ERR)
 30   nprim = iprim-1
      if (nprim .le. 0) call errquit('bas_input_cont: no primitives?',
     $     nprim, BASIS_ERR)
c check to make sure any contraction is not all zero's
c
c
      if (oIs_ecp.and.nprim.eq.1.and.ngen.eq.1.and.coeff(1,1).eq.0) then
         write(luout,*) 'Local ECP potential is zero and ignored'
         write(luout,*) ' ' 
         return
      endif
      call bas_check_contractions(nprim,ngen,prim_max,coeff,expnt)
c
c   Now have tag, contraction type, no. of contractions, no. of prims,
c   exponents, coeffs.  Shove this lot into the basis set.
c
c
c   bas_add_ucnt -> adds a new general contraction on the specified 
c                   tag.
c AJL/SPIN-POLARISED ECPs
c   bas_add_ucnt0 -> adds the contraction to a specified channel
c
c   If the tag is not present it will also add that.
c
      status = .true.
      if (osegment) then
c     
c     Add contractions one-at-a-time to force segmentation ... remove
c     functions with zero coefficients
c     
        ostat_split = (l_value.ge.0).or.oIs_ecp.or.oIs_so
        if (ostat_split) then
c     simple shell or all ecp/so shells
          do i = 1, ngen
            nptmp = 0
            do j = 1, nprim
              if (coeff(j,i).ne.0.0d0) then
                nptmp = nptmp + 1
                ctmp(nptmp) = coeff(j,i)
                etmp(nptmp) = expnt(j)
                if (oIs_ecp.or.oIs_so) then
                  rex_okay =
     &                ((abs(rex(j)-0.0d00)).lt.1.0d-9).or.
     &                ((abs(rex(j)-1.0d00)).lt.1.0d-9).or.
     &                ((abs(rex(j)-2.0d00)).lt.1.0d-9).or.
     &                ((abs(rex(j)-3.0d00)).lt.1.0d-9).or.
     &                ((abs(rex(j)-4.0d00)).lt.1.0d-9)
                  if (rex_okay) then
                    rex_tmp(nptmp) = rex(j)
                  else
                    write(luout,*)'component',j,
     &                  ' of the current contraction has ',
     &                  'an r-exponent of',rex(j)
                    call errquit
     &                 ('*ERR* bas_input_cont: Fatal r-exponent error',
     &                 911, BASIS_ERR)
                  endif
                endif
              endif
            enddo
cc AJL/Begin/SPIN-POLARISED ECPs
c            status = status .and.
c     $          bas_add_ucnt(basis, tag, l_value, 1, nptmp, 
c     $          rex_tmp, etmp, ctmp, prim_max,
c     $          stdtag, oshell_is_rel)
            status = status .and.
     $          bas_add_ucnt0(basis, tag, l_value, 1, nptmp, 
     $          rex_tmp, etmp, ctmp, prim_max,
     $          stdtag, oshell_is_rel, channel)
cc AJL/End
          enddo
        else
c     sp shell
cc AJL/Begin
c          status = status .and.
c     $        bas_add_ucnt(basis, tag, 0, 1, nprim, 
c     $        rex, expnt, coeff(1,1), prim_max,
c     &        stdtag, oshell_is_rel)
c          status = status .and.
c     $        bas_add_ucnt(basis, tag, 1, 1, nprim, 
c     $        rex, expnt, coeff(1,2), prim_max,
c     &        stdtag, oshell_is_rel)
          status = status .and.
     $        bas_add_ucnt0(basis, tag, 0, 1, nprim, 
     $        rex, expnt, coeff(1,1), prim_max,
     &        stdtag, oshell_is_rel,channel) ! Added Channel
          status = status .and.
     $        bas_add_ucnt0(basis, tag, 1, 1, nprim, 
     $        rex, expnt, coeff(1,2), prim_max,
     &        stdtag, oshell_is_rel,channel) ! Added Channel
        endif
      else
c     
c     Add as a single general contraction or sp shell
c     
c        status = status .and.
c     $      bas_add_ucnt(basis, tag, l_value, ngen, nprim, 
c     $      rex, expnt, coeff, prim_max,
c     &      stdtag, oshell_is_rel)
        status = status .and.
     $      bas_add_ucnt0(basis, tag, l_value, ngen, nprim, 
     $      rex, expnt, coeff, prim_max,
     &      stdtag, oshell_is_rel, channel)
      endif
cc AJL/End
      if (.not. status) call errquit
     $        ('bas_input_cont: bas_add_ucnt0 failed!!', 0, BASIS_ERR)
c
      return
c
10000 write(LuOut,1)
 1    format(' basis/ecp contraction format is:',/,
     $       '       tag <contraction type> [rel] [channel]',/,
     $       '           <exponent> <contraction coefficients>',/,
     $       '           ... ',/,
     $       ' end basis/ecp')
      call errquit('bas_input_cont: format error in the input', 0,
     &       INPUT_ERR)
10100 format(/,' **** WARNING Zero Coefficient **** on atom "', a16,
     &    '"',/,' angular momentum value:',i2,
     &    '  standard basis set name: "',a,'"',/,
     &    1x,'input line that generated warning:')
c
      end
*.....................................................................
      subroutine bas_tag_lib(basis, osegment, tag, tag_in_lib, 
     $     standard, filename, oshell_is_rel,read_sphere)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "geom.fh"
#include "inp.fh"
#include "ecpso_decP.fh"
#include "stdio.fh"
#include "util.fh"
      integer basis             !< [Input] basis handle
      logical osegment          !< [Input] segment if true
      character*16 tag          !< [Input] tag
      character*16 tag_in_lib   !< [Input] tag to lookup in library
      character*(*) standard    !< [Input] name of standard basis set
      character*(*) filename    !< [Input] name of library file
      logical oshell_is_rel     !< [Input] flag for relativistic shell(s)
      logical read_sphere       !< [Input] flag for reading spherical from library
c
      character*5 id_string
      character*14 id_str_ecp
      character*2 symbol
      character*16 element
      character*255 string
      character*255 filestring
      character*4 type_of_cont
      integer atn, lenn, lena
      integer nltypes           ! No. of known angular momentum types
      integer nsptypes          ! No. of known sp type shells
      parameter (nltypes = 7, nsptypes = 3)
      character*1 ltypes(nltypes)
      character*2 sptypes(nsptypes)
      character*255 field
      logical status
      logical oIs_ecp           ! ecp library function read?
      logical oIs_so            ! so library function read?
      logical file_exists
      logical inp_search_fast
      external inp_search_fast
      integer spvalues(nsptypes)
      integer l_value, ind
      integer num_elec
      logical bas_set_spherical
      external bas_set_spherical
      data ltypes /'s','p','d','f','g','h','i'/
      data sptypes / 'sp', 'l ', 'ul'/
      data spvalues/  -1 , -1  ,  -1 /
c
#include "ecpso_sfnP.fh"
c
      id_str_ecp = ""
      oIs_ecp = Is_ECP_in(basis)
      oIs_so  = Is_SO_in(basis)
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
     &    ('bas_tag_lib: tag does not refer to an atom', 0, INPUT_ERR)
c
c     replace all underscores with spaces in basis name in standard
c
      lenn = inp_strlen(standard)
      do ind = 1, lenn
         if (standard(ind:ind).eq.'_') standard(ind:ind) = ' '
      enddo
c
c     check for user typing multiple spaces in basis set by accident
c
  21  lenn = inp_strlen(standard)
      do ind = 1, lenn-1
         if (standard(ind:ind+1).eq.'  ') then
            standard(ind+1:lenn-1) = standard(ind+2:lenn)
            standard(lenn:lenn) = ' '
            goto 21
         endif
      enddo
c
      string = symbol
      lenn = inp_strlen(string)
      string(lenn+1:lenn+1) = '_'
      string(lenn+2:) = standard
      call inp_lcase(string)
c
c two possibilities for filename:
c 1. user defined file -> do not do anything, just read file
c 2. default path or user defined path -> append a file, check if it exists
c
      filestring = filename
      lenn = inp_strlen(filename)
#ifdef WIN32
      if (filestring(lenn:lenn) .eq. '\') then
#else
      if (filestring(lenn:lenn) .eq. '/') then
#endif
c
c we have to append the name of a basis set file to the path
c append, but remove some problem characters from basis name in standard:
c ( and ) -> remove
c blank , and / -> replace with underscore
c * -> replace with s
c
          lena = lenn 
          do 31 ind = 1, inp_strlen(standard)
             if (standard(ind:ind) .eq. '(' .or. 
     &           standard(ind:ind) .eq. ')' .or.     
     &           standard(ind:ind) .eq. '[' .or. 
     &           standard(ind:ind) .eq. ']' ) goto 31
             lena = lena + 1
             if (standard(ind:ind) .eq. ',' .or.
     &           standard(ind:ind) .eq. '/' .or.
     &           standard(ind:ind) .eq. ' ' ) then
                 filestring(lena:lena) = '_'
             elseif (standard(ind:ind) .eq. '*') then
                 filestring(lena:lena) = 's'
             else
                 filestring(lena:lena) = standard(ind:ind)
             endif
 31       continue
          call inp_lcase(filestring(lenn+1:))
          file_exists = .false.
          inquire(file=filestring,exist=file_exists)
          if (.not. file_exists) then
             write(LuOut,32) filestring
 32          format(' Attempting to read a basis set from a ',
     &              ' non-existing file:'/' ',A)
             call errquit('bas_tag_lib: failed opening basis file',0,
     &       BASIS_ERR)
          endif
      endif
c
c     Open the library file and try to locate the basis
c
      if (filename .eq. ' ') filename = 'library'
      open(33, file=filestring, status='old', err=1000)

      call inp_save_state()     ! Save state for unit 5.

      call inp_init(33,6)
      if (oIs_ecp) then
        id_string = 'ecp'
      else if (oIs_so) then
        id_string = 'so'
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
c
c     We did not find what we were looking for in the current file!
c
c     If we were looking for an ECP or SO the basis file might provide
c     an alternative location with the ASSOCIATED_ECP keyword, so let us
c     see if we can find that.
c
      if ((oIs_ecp.or.oIs_so).and.
     +    inp_compare(.false.,"",id_str_ecp)) then
c
c       Restart the I/O on the current basis set file to get ready to
c       look for the "ASSOCIATED_ECP" keyword.
c
        rewind(33)
        call inp_init(33,6)
        id_str_ecp = "ASSOCIATED_ECP"
 333    if (inp_search_fast(id_str_ecp)) then
c
c         Found the "ASSOCIATED_ECP" keyword, now read the basis set
c         name it is pointing to.
c
          status = inp_a(field)
          if (inp_a(field)) then
c
c           Got the basis set name. Now replace the original basis set
c           name with the ecp basis set name and re-issue the I/O
c
            standard = field
            lenn = inp_strlen(standard)
            do ind = 1, lenn
               if (standard(ind:ind).eq.'_') standard(ind:ind) = ' '
            enddo
c
c           Restore the I/O as if we are leaving the current basis set
c           file.
c
            call inp_restore_state
            call inp_clear_err
            goto 21 ! Re-issue to I/O on the ecp basis set file.
          else
            goto 35
          endif
          goto 333
        endif
      endif
 35   call inp_restore_state    ! Failure ... restore unit 5 and abort
      call inp_clear_err
      write(LuOut,*) ' bas_tag_lib: failed to locate ', 
     C     id_string(1:inp_strlen(id_string)),' ',
     C     string(1:inp_strlen(string)), 
     C     ' in file ',filename(1:inp_strlen(filename))
      call errquit('bas_tag_lib: no such basis available', 0,
     &       INPUT_ERR)
c
 34   continue                  ! success
c
c     Now read the basis for this magical atom type in using
c     the given tag
c
 10   if (.not. inp_a(field)) field=' '
c
c     ignore spherical vs cartesian from bse for He/H
c     since they do assign spherical when the other elements are cartesian
c
      if (.not.inp_compare(.false., 'h', string(1:1))) then
      if(read_sphere) then
         if (inp_compare(.false., 'spherical', field)) then
            if (.not. bas_set_spherical(basis, .true.))
     &           call errquit
     &           (' bas_set_spherical failed ',911, BASIS_ERR)
         endif
      endif
      endif
      if (.not. inp_read()) call errquit
     $     ('bas_tag_lib: premature EOF', 0, INPUT_ERR)
c     
c     Start parsing current line
c     
 20   call inp_set_field(0)
      id_str_ecp = ""
      if (.not. inp_a(field)) goto 10000 ! tag
      if (inp_compare(.false., 'end', field)) goto 9000 ! End
      if (.not. inp_a(field)) goto 10000 ! type_of_cont or 'nelec'
      if (inp_compare(.false.,'nelec', field)) then
         if (.not.inp_i(num_elec)) goto 10000
         if (.not.ecp_set_num_elec(basis,tag,num_elec,'unknown'))
     &       call errquit
     &       ('bas_tag_lib: ecp_set_num_elec failed',911, BASIS_ERR)
         if (.not.inp_read())
     &       call errquit
     &       ('bas_tag_lib: premature end of input in library',911,
     &       BASIS_ERR)
         goto 20
      else
         call inp_prev_field()
         if (.not. inp_a(type_of_cont)) goto 10000
      endif
c
      if (inp_match(nltypes, .false., type_of_cont, ltypes, ind)) then
c     
c     The contraction is a simple shell
c     
         l_value = ind - 1
      else if (
     &       inp_match(nsptypes,.false.,type_of_cont,sptypes,ind)) then
c     
c     The contraction is an sp-type shell
c     
         l_value = spvalues(ind)
         if (oshell_is_rel) call errquit(
     &       'bas_tag_lib: sp shell illegal with relativistic basis',
     &       99, INPUT_ERR)
      else
         call errquit('bas_tag_lib: invalid shell type?', 0, INPUT_ERR)
      endif
c
c     Have tag and l_value ... read in the contraction coeffcients
c
      call bas_input_cont(basis, osegment, tag, l_value,
cc AJL/Begin/SPIN ECPs
c     &    standard, oshell_is_rel)
     &    standard, oshell_is_rel, 0)
c                                  ^ Indicates that this applies to both
c                                    channels i.e. no spin-polarised ecps
c                                    when reading from libraries
cc AJL/End
      goto 20
c
 9000 close(33)
      call inp_restore_state    ! Restore state for unit 5
c
      return
c
 1000 write(LuOut,*) ' bas_tag_lib: looking for ', tag_in_lib, standard,
     $     filename
      call errquit('bas_tag_lib: failed to open basis library', 0,
     &       BASIS_ERR)
c
10000 call errquit('bas_tag_lib: problem with library format?', 0,
     &       BASIS_ERR)
c
      end
      subroutine bas_check_contractions(np,ng,ldc,cc,ee)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "util.fh"
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
      call errquit('bas_check_contractions: fatal error',911,
     &       BASIS_ERR)
10000 format(10(1pd14.6))
      end
      subroutine bas_library_file(libname)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
*passed
      character*(*) libname ! [output] filename for library file
*local
      character*(nw_max_path_len) compiled_name
      character*128 test
      character*24  basis_library
      integer length
      integer unitrc
      logical from_environment
      logical from_nwchemrc
      logical from_compile
      logical does_it_exist
      logical is_it_open
      logical debug
      logical foundit
      logical nwrcopen
      logical util_find_dir
      logical noslash
      external util_find_dir
      integer calls
      data calls/0/
      save calls
*
      nwrcopen=.false.
      from_environment = .false.
      from_nwchemrc    = .false.
      from_compile     = .false.
      debug = .false.
      noslash=.false.
      calls = calls + 1
      if (debug) write(luout,*)' calls:calls: ',calls
      call util_nwchem_srcdir(compiled_name)
c
c     add here basis/libraries bit
c
      if(libname.ne.' ') then
      compiled_name=compiled_name(1:inp_strlen(compiled_name))
     $     //"/basis/"//libname(1:inp_strlen(libname))//"/"
      else
      compiled_name=compiled_name(1:inp_strlen(compiled_name))
     $     //"/basis/libraries/"
      endif
*
* order of precedence for choosing name
* 1) value of NWCHEM_BASIS_LIBRARY environment variable
* 2) value of NWCHEM_BASIS_LIBRARY set in $HOME/.nwchemrc file
* 3) value of the compiled in library name
*      
*1: check for NWCHEM_BASIS_LIBRARY environment variable
      call util_getenv('NWCHEM_BASIS_LIBRARY',libname)
      if (debug) then
        write(luout,*)
     &      'env return value of NWCHEM_BASIS_LIBRARY <',
     &      libname(1:inp_strlen(libname)),'>'
      endif
      length = inp_strlen(libname)
      if (length.gt.0) then
        does_it_exist = .false.
        if (libname(length:length).eq.'/') then
           does_it_exist=util_find_dir(libname)
        else
           inquire(file=libname,exist=does_it_exist)
           noslash=does_it_exist
        endif
        if (does_it_exist)then
          from_environment = .true.
          goto 99999
        else
          write(luout,*)' warning:::::::::::::: from_environment'
          write(luout,*)' NWCHEM_BASIS_LIBRARY set to: <',
     &        libname(1:inp_strlen(libname)),'>'
          write(luout,*)' but file does not exist !'
          write(luout,*)' using .nwchemrc or compiled library'
        endif
      endif
*2: check for NWCHEM_BASIS_LIBRARY defined in nwchemrc config file
*
      basis_library='nwchem_basis_library'
      if(.not.util_nwchemrc_get(basis_library,libname)) then
        if (debug) then
          write(luout,*)'util_nwchemrc_get failed'
        endif
      else
        length = inp_strlen(libname)
        does_it_exist = .false.
        if (libname(length:length).eq.'/') then
           does_it_exist=util_find_dir(libname)
        else
           inquire(file=libname,exist=does_it_exist)
           noslash=does_it_exist
        endif
        if (does_it_exist)then
          from_nwchemrc = .true.
          goto 99999
        else
          write(luout,*)' warning:::::::::::::: from_nwchemrc'
          write(luout,*)' NWCHEM_BASIS_LIBRARY set to: <',
     &        libname(1:inp_strlen(libname)),'>'
          write(luout,*)' but file does not exist !'
          write(luout,*)' using compiled in library'
        endif
      endif
      libname = compiled_name
      does_it_exist = .false.
#ifndef FUJITSU
      inquire(file=libname,exist=does_it_exist)
#endif
c
c     check if it is a directory
c
      if(.not.does_it_exist) then
         does_it_exist=util_find_dir(libname)
      endif
      if (does_it_exist)then
        from_compile = .true.
      else
        write(luout,*)' warning:::::::::::::: from_compile'
        write(luout,*)' NWCHEM_BASIS_LIBRARY is: <',
     &      libname(1:inp_strlen(libname)),'>'
        write(luout,*)' but file does ','not exist or you ',
     &      'do not have ','access to it !'
      endif
99999 continue
      if (debug.or.(from_environment.or.from_nwchemrc)) then
        if (from_environment)
     &      write(luout,*)' library name resolved from: environment'
        if (from_nwchemrc)
     &      write(luout,*)' library name resolved from: .nwchemrc'
        if (from_compile)
     &      write(luout,*)
     &      ' library name resolved from: compiled reference'
        write(luout,*)' library file name is: <',
     &      libname(1:inp_strlen(libname)),'>'
        write(luout,*) ' '
      endif
      if(noslash) then
         write(luout,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
         write(luout,*) "! DANGEROUS: nwchem_basis_library set equal !"
         write(luout,*) "! to a filename and not to a directory      !"
         write(luout,*) "! since you did not use the ending slash.   !"
         write(luout,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      endif
c     if (nwrcopen)
c    .     close(unit=unitrc,status='keep') ! close file after library name found
      end
      subroutine bas_set_library_name(mydir)
      implicit none
      character*(*) mydir
#include "baslibraryP.fh"
#include "inp.fh"
      user_library_name=mydir
      call bas_library_file(user_library_name)
      end
*.....................................................................
      subroutine bas_ecce_print_basis(basisin,module_id)
      implicit none
#include "errquit.fh"
#include "basdeclsP.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "ecpso_decP.fh"
#include "bas_exndcf_dec.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
      logical ecp_get_num_elec
      external ecp_get_num_elec
      logical is_ecce_print_on
      external is_ecce_print_on
*
* print basis set details to ecce output file
*
      integer basisin  ! [input] basis set handle
      character*(*) module_id
*
      integer basis  ! lexical basis set index
      character*128 myline  ! line for printing
      integer ntags, itag
      integer ncont, icont, fcont, lcont
      integer nelec
      logical oIs_ecp, oIs_so, oIs_basis
      integer j, mylen
      integer type, nprim, iprim, ngen, igen, iexp, icoef, irexp, bptr
      character*2 stype
      character*16 mytag
      double precision bufl(32)
#include "bas_exndcf_sfn.fh"
#include "ecpso_sfnP.fh"
*
      if (is_ecce_print_on()) then
*        write(LuOut,*)' ecce print is on '
        continue
      else
*        write(LuOut,*)' ecce print is off '
        return
      endif
      basis = basisin + Basis_Handle_Offset 
      oIs_ecp = Is_ECP_in(basisin)
      oIs_so  = Is_SO_in(basisin)
      oIs_basis = (.not.(oIs_ecp.or.oIs_so))
      call ecce_print_module_entry(module_id)
      myline = 'b-a-d'
      if (oIs_ecp)   myline = 'ecp   '
      if (oIs_so)    myline = 'so    '
      if (oIs_basis) myline = 'basis '
*        
*      write(LuOut,*)'len basis',len_bs_name(basis)
*      write(LuOut,*)'len trans',len_bs_trans(basis)
*
      if (len_bs_trans(basis).eq.0) then
        myline = myline(1:6)//'"'//
     &      bs_name(basis)(1:len_bs_name(basis))//'"'
      else
        myline = myline(1:6)//'"'//
     &      bs_trans(basis)(1:len_bs_trans(basis))//'"'
      endif
*
* add spherical or cartesean tag
*
      mylen = inp_strlen(myline)
      if (bas_spherical(basis)) then
         myline = myline(1:mylen)//' spherical'
      else
         myline = myline(1:mylen)//' cartesian'
      endif
      call ecce_print_echo_string(myline)
      ntags = infbs_head(Head_Ntags,basis)
      do itag = 1,ntags
        ncont = infbs_tags(Tag_Ncont,itag,basis)
        fcont = infbs_tags(Tag_Fcont,itag,basis)
        lcont = infbs_tags(Tag_Lcont,itag,basis)
        if ((lcont-fcont+1).ne.ncont) then
          write(luout,*)' Number of tag contractions strange '
          write(luout,*)' First contraction: ',fcont
          write(luout,*)' Last  contraction: ',lcont
          write(luout,*)' computed number  : ',(lcont-fcont+1)
          write(luout,*)' stored number    : ',ncont
          call errquit(
     &        'bas_ecce_print_basis: fatal error:ncont',911,
     &       BASIS_ERR)
          
        endif
        mytag = bs_tags(itag,basis)
        if (oIs_ecp) then
          if (.not.ecp_get_num_elec(basisin,mytag,nelec))
     &        call errquit(
     &        'bas_ecce_print_basis: fatal error:nelec',911,
     &       BASIS_ERR)
          if (nelec.gt.0) then
             myline = ' '
             write(myline,10000)mytag,nelec
             call ecce_print_echo_string(myline)
           endif
        endif
        do icont = fcont,lcont
          type = infbs_cont(Cont_Type,icont,basis)
          if (type.eq.-1.and.oIs_ecp) then
            stype = 'ul'
          elseif (type.eq.-1.and.oIs_basis) then
            stype = 'sp'
          elseif (type.eq.0) then
            stype = 's'
          elseif (type.eq.1) then
            stype = 'p'
          elseif (type.eq.2) then
            stype = 'd'
          elseif (type.eq.3) then
            stype = 'f'
          elseif (type.eq.4) then
            stype = 'g'
          elseif (type.eq.5) then
            stype = 'h'
          elseif (type.eq.6) then
            stype = 'i'
          else
            call errquit(
     &        'bas_ecce_print_basis: fatal error:stype',911,
     &       BASIS_ERR)
          endif
          myline = ' '
          write(myline,10001)mytag,stype
          call ecce_print_echo_string(myline)
          nprim = infbs_cont(Cont_Nprim,icont,basis)
          ngen  = infbs_cont(Cont_Ngen,icont,basis)
          iexp  = infbs_cont(Cont_Iexp,icont,basis)
          icoef = infbs_cont(Cont_Icfp,icont,basis)
          if (.not.oIs_basis) then
            irexp = infbs_cont(Cont_Irexp,icont,basis)
          endif
          do iprim = 0,(nprim-1)
            bptr = 1
            if (.not.oIs_basis) then
              bufl(bptr) = sf_exndcf((irexp+iprim),basis)
              bptr = bptr + 1
            endif
            bufl(bptr) = sf_exndcf((iexp+iprim),basis)
            bptr = bptr + 1
            do igen = 0,(ngen-1)
              bufl(bptr) = sf_exndcf((icoef+iprim+(igen*nprim)),basis)
              bptr = bptr + 1
            enddo
            myline = ' '
            if (.not.oIs_basis) then
              write(myline,10002)(bufl(j),j=1,(ngen+2))
            else
              write(myline,10003)(bufl(j),j=1,(ngen+1))
            endif
            call ecce_print_echo_string(myline)
          enddo
        enddo
      enddo
      myline = 'end'
      call ecce_print_echo_string(myline)
      call ecce_print_module_exit(module_id, 'ok')
10000 format(1x,a16,1x,'nelec',1x,i3)
10001 format(1x,a16,1x,a2)
10002 format(1x,f5.0,31(2x,1pe15.6))
10003 format(31(2x,1pe15.6))
      end

      

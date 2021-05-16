c $Id: md_xs_input.F 25011 2013-12-19 17:45:25Z d3y133 $

*     *********************************************
*     *                                           *
*     *              md_xs_input                  *
*     *                                           *
*     *********************************************
      subroutine md_xs_input(rtdb)
      implicit none
      integer rtdb

#include "inp.fh"
#include "bafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"

*     **** local variables ****
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 15)
      character*18 dirs(num_dirs)
      data dirs / 'spectroscopy',
     >            'absorber',
     >            'edge',
     >            'center',
     >            'rmax',
     >            'stride',
     >            'rcut',
     >            'Fourier',
     >            'DW_Harmonic', 
     >            'first',
     >            'last',
     >            'hydrogen',
     >            'e0',
     >            's0',
     >            'end'/

c
      logical       value,md_xs,hydrogen,nohydrogen
      character*80  rtdb_name
      character*255 test,test2
      real*8 rmax,rcut,ft(3),e0,s0
      integer icenter,istride,iw,isteps(2)
      integer ind               ! Index of matched directive

*     **** external functions ****
      logical  nwpw_parse_boolean
      external nwpw_parse_boolean

      md_xs = .true.
      if (inp_a(test)) then
          md_xs = nwpw_parse_boolean(test,.true.)
          go to 9999
      end if

 10   if (.not. inp_read()) 
     >  call errquit(
     >  'md_xs_input: inp_read failed',0,INPUT_ERR)
      if (.not. inp_a(test))
     >  call errquit(
     >  'md_xs_input: failed to read keyword',0,INPUT_ERR)

      if (.not. inp_match(num_dirs,.false.,test,dirs,ind))
     >  call errquit(
     >  'md_xs_input: unknown_directive',ind,INPUT_ERR)

      goto ( 100, 200,300, 400, 500, 600, 700,800,900,
     >      1000, 1100,1200,1300,1400,
     >      9999) ind
      call errquit(
     >     'md_xs_input: unimplemented directive',ind,INPUT_ERR)

c
c  spectroscopy
c
 100  if (inp_a(test2)) then
         rtdb_name = 'md_xs:spectroscopy'
         value = rtdb_cput(rtdb,rtdb_name,1,test2)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing spectroscopy',100,RTDB_ERR)
      end if

      goto 10

c
c  absorbers
c
 200  if (inp_a(test2)) then
         rtdb_name = 'md_xs:absorber'
         value = rtdb_cput(rtdb,rtdb_name,1,test2)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',200,RTDB_ERR)
      end if

      goto 10

c
c  edge
c
 300  if (inp_a(test2)) then
         rtdb_name = 'md_xs:edge'
         value = rtdb_cput(rtdb,rtdb_name,1,test2)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',300, RTDB_ERR)
      end if

      goto 10

c
c  center
c
 400  if (inp_i(icenter)) then
         rtdb_name = 'md_xs:center'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,icenter)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',400, RTDB_ERR)
      end if

      goto 10

c
c  rmax
c
 500  if (inp_f(rmax)) then
         rtdb_name = 'md_xs:rmax'
         value = rtdb_put(rtdb,rtdb_name,mt_dbl,1,rmax)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',500, RTDB_ERR)
      end if

      goto 10

c
c  stride
c
 600  if (inp_i(istride)) then
         rtdb_name = 'md_xs:stride'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,istride)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',600, RTDB_ERR)
      end if

      goto 10

c
c  rcut
c
 700  if (inp_f(rcut)) then
         rtdb_name = 'md_xs:rcut'
         value = rtdb_put(rtdb,rtdb_name,mt_dbl,1,rcut)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',700, RTDB_ERR)
      end if

      goto 10

c
c  Fourier Kmin Kmax dK iw
c
 800  if (.not.inp_f(ft(1))) ft(1) = 2.3
      if (.not.inp_f(ft(2))) ft(2) = 17.0
      if (.not.inp_f(ft(3))) ft(3) = 1.0
      if (.not.inp_i(iw))    iw = 2
         rtdb_name = 'md_xs:Fourier'
         value = rtdb_put(rtdb,rtdb_name,mt_dbl,3,ft)
         rtdb_name = 'md_xs:Fourier_iw'
         value = value.and.
     >           rtdb_put(rtdb,rtdb_name,mt_int,1,iw)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',800, RTDB_ERR)

      goto 10
c 
c DW_Harmonic temp equil_steps steps 
c
 900  if (.not.inp_f(ft(1)))     ft(1)     = 298.15d0
      if (.not.inp_i(isteps(1))) isteps(1) = 500
      if (.not.inp_i(isteps(2))) isteps(2) = 1000

         rtdb_name = 'md_xs:DW_Harmonic_Temperature'
         value = rtdb_put(rtdb,rtdb_name,mt_dbl,1,ft)
         rtdb_name = 'md_xs:DW_Harmonic_steps'
         value = value.and.
     >           rtdb_put(rtdb,rtdb_name,mt_int,2,isteps)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',900, RTDB_ERR)

      goto 10

c
c  first
c
 1000 if (inp_i(istride)) then
         rtdb_name = 'md_xs:first'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,istride)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',1000, RTDB_ERR)
      end if

      goto 10

c
c  last
c
 1100 if (inp_i(istride)) then
         rtdb_name = 'md_xs:last'
         value = rtdb_put(rtdb,rtdb_name,mt_int,1,istride)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',1100, RTDB_ERR)
      end if

      goto 10

c
c  hydrogen
c
 1200 if (inp_a(test2)) then
         hydrogen = nwpw_parse_boolean(test2,.true.)
         rtdb_name = 'md_xs:nohydrogen'
         nohydrogen = .not.hydrogen
         value = rtdb_put(rtdb,rtdb_name,mt_log,1,nohydrogen)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',1200, RTDB_ERR)
      end if

      goto 10

c
c  e0
c
 1300 if (inp_f(e0)) then
         rtdb_name = 'md_xs:e0'
         value = rtdb_put(rtdb,rtdb_name,mt_dbl,1,e0)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',1300, RTDB_ERR)
      end if

      goto 10

c
c  s0
c
 1400 if (inp_f(s0)) then
         rtdb_name = 'md_xs:s0'
         value = rtdb_put(rtdb,rtdb_name,mt_dbl,1,s0)
         if (.not. value) call errquit(
     >   'md_xs_input:failed writing update',1400, RTDB_ERR)
      end if

      goto 10



 9999 continue

         if (.not.rtdb_put(rtdb,'nwpw:md_xs',mt_log,1,md_xs))
     >    call errquit(
     >   'md_xs_input:failed writing md_xs:',0,
     >    RTDB_ERR)

      
      return
      end

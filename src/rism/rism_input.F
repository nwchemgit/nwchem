      subroutine rism_input(rtdb)
c
      implicit none
#include "errquit.fh"
c
#include "stdio.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      character*32 tag
      character*72 keyword
      character*32 pname
      character*255 token
      character*32 mtoken(10)
      integer  itoken(10)
      double precision  eatoms
      double precision bqzone
      logical bq_dynamical
      integer linkatm,bq_exclude
      logical status
      integer ip,np
      logical ltoken
      logical ignore
      double precision f
c
      pname = "rism_input: "
c
      write(luout,*) "in ",pname
c
      call inp_set_field(0)
c
c     start parsing input
c     ------------------
      if (.not.inp_a(token))
     + call errquit(pname//': no input available',0, INPUT_ERR)
      if (.not.inp_compare(.false.,token,'rism'))
     + call errquit(pname//': no input available',0, INPUT_ERR)
      goto 2
    1 continue
      if (.not.inp_read()) call errquit(pname//': premature EOF',0,
     &       INPUT_ERR)
    2 continue
      if(.not.inp_a(token)) goto 1
c
c     permittivity
c     -----------
      keyword="permittivity"
      if (inp_compare(.false.,keyword, token)) then
        if (.not.inp_f(f))
     >    call errquit(pname//'failed to read '//keyword,0, INPUT_ERR)
        tag = "rism:"//keyword(1:inp_strlen(keyword))
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_dbl,1,f)) 
     >      goto 999
        goto 2
      endif
c
c     tau
c     tau separation between short and long range functions in 1/Angstroms
c     -----------
      keyword="tau"
      if (inp_compare(.false.,keyword, token)) then
        if (.not.inp_f(f))
     >    call errquit(pname//'failed to read '//keyword,0, INPUT_ERR)
        tag = "rism:"//keyword(1:inp_strlen(keyword))
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_dbl,1,f)) 
     >      goto 999
        goto 2
      endif
c
c     Temperature
c     -----------
      keyword="temp"
      if (inp_compare(.false.,keyword, token)) then
        if (.not.inp_f(f))
     >    call errquit(pname//'failed to read '//keyword,0, INPUT_ERR)
        tag = "rism:"//keyword(1:inp_strlen(keyword))
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_dbl,1,f)) 
     >      goto 999
        goto 2
      endif
c
c     tolerance
c     -----------
      keyword="tol"
      if (inp_compare(.false.,keyword, token)) then
        if (.not.inp_f(f))
     >    call errquit(pname//'failed to read '//keyword,0, INPUT_ERR)
        tag = "rism:"//keyword(1:inp_strlen(keyword))
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_dbl,1,f)) 
     >      goto 999
        goto 2
      endif
c     
c
c     lambda
c     -----------
      keyword="lambda"
      if (inp_compare(.false.,keyword, token)) then
        if (.not.inp_f(f))
     >    call errquit(pname//'failed to read '//keyword,0, INPUT_ERR)
        tag = "rism:"//keyword(1:inp_strlen(keyword))
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_dbl,1,f)) 
     >      goto 999
        goto 2
      endif
c
c     diis
c     -----------
      keyword="diis"
      if (inp_compare(.false.,keyword, token)) then
        if (.not.inp_i(itoken(1)))
     >    call errquit(pname//'failed to read '//keyword,0, INPUT_ERR)
        tag = "rism:"//keyword(1:inp_strlen(keyword))
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_int,
     +      1,itoken(1))) 
     >      goto 999
        goto 2
      endif
c     
c
c     vdw mixing rule
c     -----------
      keyword="vdw"
      if (inp_compare(.false.,keyword, token)) then
        if (.not.inp_a(token))
     >    call errquit(pname//'failed to read '//keyword,0, INPUT_ERR)
        tag = "rism:"//keyword(1:inp_strlen(keyword))
        if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),
     +      1,token)) 
     >      goto 999
        goto 2
      endif
c     
c
c     solvent
c     -----------
      keyword="solvent"
      if (inp_compare(.false.,keyword, token)) then
        if (.not.inp_a(token))
     >    call errquit(pname//'failed to read '//keyword,0, INPUT_ERR)
        tag = "rism:"//keyword(1:inp_strlen(keyword))
        if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),
     +      1,token)) 
     >      goto 999
        goto 2
      endif
cc     
c
c     closure
c     -----------
      keyword="closure"
      if (inp_compare(.false.,keyword, token)) then
        if (.not.inp_a(token))
     >    call errquit(pname//'failed to read '//keyword,0, INPUT_ERR)
        tag = "rism:"//keyword(1:inp_strlen(keyword))
        if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),
     +      1,token)) 
     >      goto 999
        goto 2
      endif
c     
c
c     number of grid points
c     -----------
      keyword="ngrid"
      if (inp_compare(.false.,keyword, token)) then
        if (.not.inp_i(itoken(1)))
     >    call errquit(pname//'failed to read '//keyword,0, INPUT_ERR)
        tag = "rism:"//keyword(1:inp_strlen(keyword))
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_int,
     +      1,itoken(1))) 
     >      goto 999
        goto 2
      endif
c     
cc     
c     print level 
c     -------------------
      if (inp_compare(.false.,'print', token)) then
         call util_print_input(rtdb, "rism")
         go to 2
      end if
c
      if (token.eq.'end') then
        write(luout,*) "out of ",pname
        return
      endif
c
      write(luout,*)pname//' unrecognized token:',
     + token(1:inp_strlen(token))
       call errquit(pname//'failed ',0,
     >        RTDB_ERR)
      return
998   call errquit(pname//'no token found '//token,0,
     >        RTDB_ERR)
999   call errquit(pname//'failed to store '//tag,0,
     >        RTDB_ERR)

      end


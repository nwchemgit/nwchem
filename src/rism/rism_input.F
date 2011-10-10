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
#ifdef RISM_DEBUG      
      write(luout,*) "in ",pname
#endif      
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
c     solute
c     -----------
      keyword="solute"
      if (inp_compare(.false.,keyword, token)) then
        call rism_solute_input(rtdb)
c        if (.not.inp_a(token))
c     >    call errquit(pname//'failed to read '//keyword,0, INPUT_ERR)
c        tag = "rism:"//keyword(1:inp_strlen(keyword))
c        if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),
c     +      1,token)) 
c     >      goto 999
        goto 2
      endif
c
c     vdw
c     -----------
      keyword="vdw"
      if (inp_compare(.false.,keyword, token)) then
32       continue
         if (.not.inp_a(token)) goto 2
         if(inp_compare(.false.,"rule", token)) then
           if (.not.inp_a(token))
     >       call errquit(pname//'vdw rule',0,INPUT_ERR)
           tag="rism:vdw:rule"
           if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),1,token)) 
     >         goto 999
           goto 32
         else if(inp_compare(.false.,"parameters", token)) then
           if (.not.inp_a(token))
     >       call errquit(pname//'vdw parameters',0,INPUT_ERR)
           tag="rism:vdw:parameters"
           if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),1,token)) 
     >         goto 999
           goto 32
         else
            call errquit(pname//'unknown tag '//token,0,INPUT_ERR)
         end if
        goto 2
      endif
c     
cc
c     solvent
c     -----------
      keyword="solvent"
      if (inp_compare(.false.,keyword, token)) then
              call rism_solvent_input(rtdb)
c31       continue
c         if (.not.inp_a(token)) goto 2
c         if(inp_compare(.false.,"model", token)) then
c           if (.not.inp_a(token))
c     >       call errquit(pname//'solvent model',0,INPUT_ERR)
c           tag="rism:solvent:model"
c           if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),1,token)) 
c     >         goto 999
c           goto 31
c         else if(inp_compare(.false.,"rdf", token)) then
c           if (.not.inp_a(token))
c     >       call errquit(pname//'solvent rdf',0,INPUT_ERR)
c           tag="rism:solvent:rdf"
c           if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),1,token)) 
c     >         goto 999
c           goto 31
c         end if
         goto 2
      endif
c     
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
#ifdef RISM_DEBUG              
        write(luout,*) "out of ",pname
#endif        
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

      subroutine rism_solute_input(rtdb)
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
      character*255 message
      character*32 mtoken(10)
      integer  itoken(10)
      logical status
      integer ip,np
      logical ltoken
      logical ignore
      double precision f
c
      pname = "rism_solute_input: "
c
#ifdef RISM_DEBUG      
      write(luout,*) "in ",pname
#endif      
c
c     read initial token
c     ------------------
      if (.not.inp_a(token)) goto 1000
#ifdef RISM_DEBUG      
      write(*,*) "found token",token
#endif      
c      
c     solute model
c     -------------
      if(inp_compare(.false.,"configuration", token)) then
        message = "processing solute configuration"
        if (.not.inp_a(token)) goto 998
        tag="rism:solute:configuration"
        if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),1,token)) 
     >      goto 999
        if (.not.inp_a(token)) goto 1000
      end if
      message = "unknown tokens"
      call rism_message(pname(1:inp_strlen(pname))//
     > ": ignoring the rest of input starting from "//
     > token)
c
c     entry point for sucsessfull return
c     -----------------------------------
1000  continue
      return
c     -----------------------------------
997   write(luout,*)pname//' unrecognized token:',
     + token(1:inp_strlen(token))
       call errquit(pname//'failed ',0,
     >        RTDB_ERR)
      return
998   call errquit(pname//message(1:inp_strlen(message)),
     >             0,0)
999   call errquit(pname//'failed to store '//tag,0,
     >        RTDB_ERR)

      end

      subroutine rism_solvent_input(rtdb)
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
      character*255 message
      character*32 mtoken(10)
      integer  itoken(10)
      logical status
      integer ip,np
      logical ltoken
      logical ignore
      double precision f
c
      pname = "rism_solvent_input: "
c
#ifdef RISM_DEBUG      
      write(luout,*) "in ",pname
#endif      
c
c     read initial token
c     ------------------
      if (.not.inp_a(token)) goto 1000
#ifdef RISM_DEBUG      
      write(*,*) "found token",token
#endif      
c      
c     solvent model
c     -------------
      if(inp_compare(.false.,"model", token)) then
        message = "processing for solvent model"
        if (.not.inp_a(token)) goto 998
        tag="rism:solvent:model"
        if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),1,token)) 
     >      goto 999
        if (.not.inp_a(token)) goto 1000
      end if
c      
c     solvent permittivity
c     -------------------------------------
      if(inp_compare(.false.,"permittivity", token)) then
        message = "processing solvent permittivity"
        if (.not.inp_f(f)) goto 998
        tag="rism:solvent:permittivity"
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),
     >  mt_dbl,1,f)) 
     >  goto 999
        if (.not.inp_a(token)) goto 1000
      end if
c      
c     solvent rdf with optional kspace flag
c     -------------------------------------
      if(inp_compare(.false.,"rdf", token)) then
        message = "processing solvent rdf"
        if (.not.inp_a(token)) goto 998
        if(inp_compare(.false.,"kspace", token)) then
          tag="rism:solvent:kspace"
          if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),
     >    mt_log,1,token)) 
     >    goto 999
          if (.not.inp_a(token)) goto 998
        end if
        tag="rism:solvent:rdf"
        if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),1,token)) 
     >      goto 999
        if (.not.inp_a(token)) goto 1000
      end if
      message = "unknown tokens"
      call rism_message(pname(1:inp_strlen(pname))//
     > ": ignoring the rest of input starting from "//
     > token)
c
c     entry point for sucsessfull return
c     -----------------------------------
1000  continue
      return
c     -----------------------------------
997   write(luout,*)pname//' unrecognized token:',
     + token(1:inp_strlen(token))
       call errquit(pname//'failed ',0,
     >        RTDB_ERR)
      return
998   call errquit(pname//message(1:inp_strlen(message)),
     >             0,0)
999   call errquit(pname//'failed to store '//tag,0,
     >        RTDB_ERR)

      end
c $Id$

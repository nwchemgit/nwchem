      subroutine qmmm_input(rtdb)
c
      implicit none
#include "errquit.fh"
c
#include "stdio.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "qmmm_params.fh"
c
      integer rtdb
c
      character*32 tag
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
c
      pname = "qmmm_input: "
c
      call inp_set_field(0)
c
c     start parsing input
c     ------------------
      if (.not.inp_a(token))
     + call errquit('qmmm_input: no input available',0, INPUT_ERR)
      if (.not.inp_compare(.false.,token,'qmmm'))
     + call errquit('qmmm_input: no input available',0, INPUT_ERR)
      goto 2
    1 continue
      if (.not.inp_read()) call errquit('qmmm_input: premature EOF',0,
     &       INPUT_ERR)
    2 continue
      if(.not.inp_a(token)) goto 1
c
c     zero point energy
c     -----------------
      if (inp_compare(.false.,"eatoms", token) .or.
     >    inp_compare(.false.,"eref", token)) then
        if (.not.inp_f(eatoms))
     >    call errquit('qmmm_input: failed to read eatoms',0, INPUT_ERR)
        tag = "qmmm:uqmatm"
        if (.not.rtdb_put(rtdb,"qmmm:uqmatm",mt_dbl,1,eatoms)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
c     number of cycles
c     ----------------
      if (inp_compare(.false.,"ncycles", token) ) then
        if (.not.inp_i(itoken(1)))
     >    call errquit('failed to read ncycles',0, INPUT_ERR)
        tag = "qmmm:ncycles"
        if (.not.rtdb_put(rtdb,"qmmm:ncycles",mt_int,1,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
c
c     link atom type
c     --------------
       if(inp_compare(.false.,'link_atoms',token) .or.
     >   inp_compare(.false.,'link_atom',token)) then
        if(.not.inp_a(token)) call errquit('no link_atoms option',0,
     &         INPUT_ERR)
        if(inp_compare(.false.,'hydrogen',token)) then
          linkatm=qmmm_hydrogen
        else if (inp_compare(.false.,'halogen',token)) then
          linkatm=qmmm_halogen
        else
         call errquit('unknown option for link_atoms',0,INPUT_ERR)
        end if
        tag="qmmm:linkatm"
        if (.not.rtdb_put(rtdb,tag,mt_int,1,linkatm)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
      endif
c
c     link atom basis set
c     -------------------
      if(inp_compare(.false.,'link_ecps',token) .or.
     >   inp_compare(.false.,'link_ecp',token)) then
        if(.not.inp_a(token)) then
         write(*,*) "please provide an option for link_ecps"
         call errquit('qmmm_input: link-ecp',0,
     &         INPUT_ERR)
        end if
        status = inp_compare(.false.,'auto',token)
        status = status .or. inp_compare(.false.,'user',token)
        if(.not.status) then
         call errquit('unknown option for link_ecps',0,INPUT_ERR)
        end if
        tag="qmmm:link_ecp"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
      if(inp_compare(.false.,'geom',token)) then
        if(.not.inp_a(token)) then
         write(*,*) "please provide an option for geom"
         call errquit('qmmm_input: geom',0,
     &         INPUT_ERR)
        end if
        tag="qmmm:geom"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
      if(inp_compare(.false.,'maxiter',token)) then
        np = inp_n_field()
        do ip = 1,np-1
        itoken(ip) = 0
        if(.not.inp_i(itoken(ip))) then
         call errquit('qmmm_input: optimization',0,
     &         INPUT_ERR)
        end if
        end do
        tag="qmmm:maxiter_n"
        if (.not.rtdb_put(rtdb,tag,mt_int,1,np-1)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        tag="qmmm:maxiter"
        if (.not.rtdb_put(rtdb,tag,mt_int,np-1,itoken)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
      if(inp_compare(.false.,'method',token)) then
        np = inp_n_field()
        do ip = 1,np-1
        mtoken(ip) = " "
        if(.not.inp_a(mtoken(ip))) then
         call errquit('qmmm_input: optimization',0,
     &         INPUT_ERR)
        end if
        end do
        tag="qmmm:method_n"
        if (.not.rtdb_put(rtdb,tag,mt_int,1,np-1)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        tag="qmmm:method"
        if (.not.rtdb_cput(rtdb,tag,np-1,mtoken)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
      if(inp_compare(.false.,'region',token)) then
        np = inp_n_field()
        do ip = 1,np-1
        mtoken(ip) = " "
        if(.not.inp_a(mtoken(ip))) then
         call errquit('qmmm_input: optimization',0,
     &         INPUT_ERR)
        end if
        end do
        tag="qmmm:region_n"
        if (.not.rtdb_put(rtdb,tag,mt_int,1,np-1)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        tag="qmmm:region"
        if (.not.rtdb_cput(rtdb,tag,np-1,mtoken)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c    
       if(inp_compare(.false.,'bq_exclude',token)) then
        tag="qmmm:bq_exclude"
        if(.not.inp_a(token)) call errquit('qmmm_input:  bq_exclude ',0,
     &       INPUT_ERR)
        if(token.eq.'linkbond_H') then
          bq_exclude=qmmm_linkbond_H
        else if(token.eq.'linkbond') then
          bq_exclude=qmmm_linkbond
        else if(token.eq.'linkres') then
          bq_exclude=qmmm_linkres
        else if(token.eq.'all') then
          bq_exclude=qmmm_all_atoms
        else if(token.eq.'none') then
          bq_exclude=qmmm_no_atoms
        else
         call errquit('unknown option for bq_exclude',0,INPUT_ERR)
        end if
        if (.not.rtdb_put(rtdb,tag,mt_int,1,bq_exclude)) 
     >      call errquit(pname//'failed to store bq_exclude',0,
     >       RTDB_ERR)
        goto 2
      endif
c
      if(inp_compare(.false.,'load',token)) then
        if (.not.inp_a(token))
     +    call errquit("please provide an option for load",
     +                0, INPUT_ERR)
        if(inp_compare(.false.,'bq',token)) then
          if(.not.rtdb_put(rtdb,"qmmm:readbq",mt_log,1,.true.))
     >     call errquit(pname//'failed to store '//tag,0,
     >        RTDB_ERR)
        else
          call errquit('unknown option for load',0,INPUT_ERR)
        end if
        goto 2
      end if
c
      if(inp_compare(.false.,'bq_zone',token) .or.
     >   inp_compare(.false.,'bqzone',token)) then
        tag="qmmm:bqzone"
        if (.not.inp_f(bqzone))
     +    call errquit('qmmm_input: failed to read bqzone',0, INPUT_ERR)
        if(.not.rtdb_put(rtdb,tag,mt_dbl,1,bqzone))
     >  call errquit(pname//'failed to store '//tag,0,
     >        RTDB_ERR)

        goto 2
      end if

      if(inp_compare(.false.,'bq_update',token)) then
        if(.not.inp_a(token)) then
         write(*,*) "please provide an option for "//token
         call errquit('qmmm_input: '//token,0,
     &         INPUT_ERR)
        end if
        if(inp_compare(.false.,'static',token)) then
          bq_dynamical=.false.
        else if(inp_compare(.false.,'dynamic',token)) then
          bq_dynamical=.true.
        else
          call errquit('unknown option for bq_update',0,INPUT_ERR)
        end if
        tag="qmmm:bq_dynamical"
        if (.not.rtdb_put(rtdb,tag,mt_log,1,bq_dynamical)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      end if

c     
c     print level 
c     -------------------
      if (inp_compare(.false.,'print', token)) then
         call util_print_input(rtdb, "qmmm")
         go to 2
      end if
c
      if (token.eq.'end') then
        return
      endif
c
      write(luout,*)' unrecognized token in qmmm input:',
     + token(1:inp_strlen(token))
       call errquit(pname//'failed ',0,
     >        RTDB_ERR)
      return
      end

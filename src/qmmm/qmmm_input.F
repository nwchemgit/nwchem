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
      character*32 interface
      real*8 eatoms
      logical fixlinks,fixqm
      integer linkatm,bq_exclude
      logical status
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
      if (inp_compare(.false.,"eatoms", token)) then
        if (.not.inp_f(eatoms))
     +    call errquit('qmmm_input: failed to read eatoms',0, INPUT_ERR)
        tag = "qmmm:uqmatm"
        if (.not.rtdb_put(rtdb,"qmmm:uqmatm",mt_dbl,1,eatoms)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
      if (inp_compare(.false.,"link_bonds", token) ) then
        if(.not.inp_a(token)) call errquit(' no link_bonds option',0,
     &         INPUT_ERR)
        if(inp_compare(.false.,"free", token)) then
          fixlinks=.false.
        else if(inp_compare(.false.,"fixed", token)) then
          fixlinks=.true.
        else
          call errquit('unknown option for linkbonds',0,INPUT_ERR)
        end if
        tag="qmmm:fixlinks"
        if (.not.rtdb_put(rtdb,tag,mt_log,1,fixlinks)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
       end if
c
      if(inp_compare(.false.,'link_atoms',token)) then
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
      if(inp_compare(.false.,'link_ecps',token)) then
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
      if(inp_compare(.false.,'optimization',token)) then
        if(.not.inp_a(token)) then
         write(*,*) "please provide an option for optimization"
         call errquit('qmmm_input: optimization',0,
     &         INPUT_ERR)
        end if
        status = inp_compare(.false.,'all',token)
        status = status .or. inp_compare(.false.,'qm',token)
        status = status .or. inp_compare(.false.,'mm',token)
        if(.not.status) then
         call errquit('unknown option for optimization',0,INPUT_ERR)
        end if
        tag="qmmm:optimization"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
      if(inp_compare(.false.,'bq_exclude',token)) then
        tag="qmmm:bq_exclude"
        if(.not.inp_a(token)) call errquit('qmmm_input:  bq_exclude ',0,
     &       INPUT_ERR)
        if(token.eq.'hydrogen') then
          bq_exclude=qmmm_hydrogen
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
      if (token.eq.'end') then
        return
      endif
c
      write(luout,*)' unrecognized token in qmmm input:',
     + token(1:inp_strlen(token))
      return
      end

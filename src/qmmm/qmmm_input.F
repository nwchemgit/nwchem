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
      character*32 optimize
      character*32 interface
      real*8 eatoms
      logical xyz,xyzq,fixlinks,fixqm
      integer linkatm,bq_exclude
c
      pname = "qmmm_input: "
c
      call inp_set_field(0)
c
c**** default values for parameters
      xyz  = .false.
      xyzq = .false.
      fixlinks = .false.
      fixqm = .false.
      linkatm=qmmm_halogen
      bq_exclude=qmmm_hydrogen
      optimize='qm'
      interface='qm'
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
      if (token.eq.'eatoms') then
        if (.not.inp_f(eatoms))
     +    call errquit('qmmm_input: failed to read eatoms',0, INPUT_ERR)
        tag = "qmmm:uqmatm"
        if (.not.rtdb_put(rtdb,"qmmm:uqmatm",mt_dbl,1,eatoms)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
 
      endif
c
      if (token.eq.'xyz') then
        xyz = .true.
        tag="qmmm:xyz"
        if (.not.rtdb_put(rtdb,tag,mt_log,1,xyz)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)

        goto 2 

      end if
      if (token.eq.'fix') then
        if(.not.inp_a(token)) call errquit(' no fix option',0,
     &         INPUT_ERR)
        if(token.eq.'linkbonds') then
          fixlinks=.true.
          tag="qmmm:fixlinks"
          if (.not.rtdb_put(rtdb,tag,mt_log,1,fixlinks)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

          goto 2
        elseif(token.eq.'qmregion') then
          fixqm=.true.
          tag="qmmm:fixqm"
          if (.not.rtdb_put(rtdb,tag,mt_log,1,fixqm)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

          goto 2
        else
          call errquit('qmmm_input: invalid fix option',0, INPUT_ERR)
        endif
       goto 2 
      end if
      if (token.eq.'xyzq') then
        xyzq = .true.
        tag="qmmm:xyz"
        if (.not.rtdb_put(rtdb,tag,mt_log,1,xyzq)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)

        goto 2 
      end if
c
      if(token.eq.'link') then
        if(.not.inp_a(token)) call errquit('qmmm_input: no link atom',0,
     &         INPUT_ERR)
        if(token.eq.'hydrogen') then
          linkatm=qmmm_hydrogen
          tag="qmmm:linkatm"
          if (.not.rtdb_put(rtdb,tag,mt_int,1,linkatm)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

          goto 2
        elseif(token.eq.'halogen') then
          linkatm=qmmm_halogen
          tag="qmmm:linkatm"
          if (.not.rtdb_put(rtdb,tag,mt_int,1,linkatm)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

          goto 2
        else
          call errquit('qmmm_input: invalid link atom',0, INPUT_ERR)
        endif
      endif
c
      if(token.eq.'link_ecp') then
        if(.not.inp_a(token)) then
         write(*,*) "please provide an option for link-ecp"
         call errquit('qmmm_input: optimize',0,
     &         INPUT_ERR)
        end if
        tag="qmmm:link_ecp"
          if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c
      if(token.eq.'optimize') then
        if(.not.inp_a(token)) then
         write(*,*) "please provide an option for optimize"
         call errquit('qmmm_input: optimize',0,
     &         INPUT_ERR)
        end if
        if(token.eq.'qm') then
          optimize='qm'
          tag="qmmm:optimize"
          if (.not.rtdb_cput(rtdb,tag,1,optimize)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
        elseif(token.eq.'mm') then
          optimize='mm'
          tag="qmmm:optimize"
          if (.not.rtdb_cput(rtdb,tag,1,optimize)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
        else
          call errquit('qmmm_input: invalid option',0, INPUT_ERR)
        endif
      endif
c
      if(token.eq.'interface') then
        if(.not.inp_a(token)) then
         write(*,*) "please provide an option for interface"
         call errquit('qmmm_input: optimize',0,
     &         INPUT_ERR)
        end if
        if(token.eq.'qm') then
          interface='qm'
          tag = "qmmm:interface"
          if (.not.rtdb_cput(rtdb,tag,1,interface)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)
          goto 2
        elseif(token.eq.'mm') then
          interface='mm'
          tag = "qmmm:interface"
          if (.not.rtdb_cput(rtdb,tag,1,interface)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)
           goto 2
        else
          call errquit('qmmm_input: invalid option',0, INPUT_ERR)
        endif
      endif
cc
      if(token.eq.'bq_exclude') then
      tag="qmmm:bq_exclude"
      if(.not.inp_a(token)) call errquit('qmmm_input:  bq_exclude ',0,
     &       INPUT_ERR)
      if(token.eq.'hydrogen') then
      bq_exclude=qmmm_hydrogen
        if (.not.rtdb_put(rtdb,tag,mt_int,1,bq_exclude)) 
     >      call errquit(pname//'failed to store bq_exclude',0,
     >       RTDB_ERR)

      goto 2
      endif
      if(token.eq.'all') then
      bq_exclude=qmmm_all_atoms
        if (.not.rtdb_put(rtdb,tag,mt_int,1,bq_exclude)) 
     >      call errquit(pname//'failed to store bq_exclude',0,
     >       RTDB_ERR)
      goto 2
      endif
      if(token.eq.'none') then
      bq_exclude=qmmm_no_atoms
        if (.not.rtdb_put(rtdb,tag,mt_int,1,bq_exclude)) 
     >      call errquit(pname//'failed to store bq_exclude',0,
     >       RTDB_ERR)
      goto 2
      endif
      call errquit('qmmm_input: invalid bq_exclude option',0, INPUT_ERR)
      endif
c
      if (token.eq.'end') then
      if(optimize.eq."mm") interface="mm"
       if (.not.rtdb_cput(rtdb,'qmmm:interface',1,interface))
     + call errquit('qmmm_input: failed to store qmmm:interface',
     &       0,RTDB_ERR)
c      if (.not.rtdb_put(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
c     + call errquit('qmmm_input: failed to store qmmm:uqmatm',0,
c     &       RTDB_ERR)
c      if (.not.rtdb_put(rtdb,'qmmm:eatoms',mt_dbl,1,eatoms))
c     + call errquit('qmmm_input: failed to store qmmm:eatoms',0,
c     &       RTDB_ERR)
c       if (.not.rtdb_put(rtdb,'qmmm:linkatm',mt_int,1,linkatm))
c     + call errquit('qmmm_input: failed to store qmmm:linkatm',linkatm,
c     &       RTDB_ERR)
c      if (.not.rtdb_put(rtdb,'qmmm:xyz',mt_log,1,xyz))
c     + call errquit('qmmm_input: failed to store qmmm:xyz',xyz,
c     &       RTDB_ERR)
c      if (.not.rtdb_put(rtdb,'qmmm:xyzq',mt_log,1,xyzq))
c     + call errquit('qmmm_input: failed to store qmmm:xyzq',xyzq,
c     &       RTDB_ERR)
c      if (.not.rtdb_put(rtdb,'qmmm:bq_exclude',mt_int,1,bq_exclude))
c     + call errquit('qmmm_input: failed to store qmmm:bq_exclude',
c     &       bq_exclude,RTDB_ERR)
c      if (.not.rtdb_put(rtdb,'qmmm:fixlinks',mt_log,1,fixlinks))
c     + call errquit('qmmm_input: failed to store qmmm:fixlinks',
c     &       0,RTDB_ERR)
c      if (.not.rtdb_put(rtdb,'qmmm:fixqm',mt_log,1,fixqm))
c     + call errquit('qmmm_input: failed to store qmmm:fixqm',
c     &       0,RTDB_ERR)
c       if (.not.rtdb_cput(rtdb,'qmmm:optimize',1,optimize))
c     + call errquit('qmmm_input: failed to store qmmm:optimize',
c     &       0,RTDB_ERR)
        return
      endif
c
      write(luout,*)' unrecognized token in qmmm input:',
     + token(1:inp_strlen(token))
      write(luout,1000)
 1000 format(/,/,'QM/MM input format:',/,'qmmm',/,'eatoms <real>',/,
     + 'end',///)
      call errquit ('qmmm_input: called with bad input',911, INPUT_ERR)
      return
      end

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
      character*255 token
      character*32 optimize
      character*32 interface
      real*8 eatoms,qatoms
      logical xyz,xyzq,fixlinks,fixqm
      integer linkatm,bq_exclude
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
     +  call errquit('qmmm_input: failed to read eatoms',0, INPUT_ERR)
      write(luout,*)' QM/MM separated atom energy:',eatoms
      call util_flush(luout)
      goto 2
      endif
c
      if (token.eq.'xyz') then
       xyz = .true.
       write(*,*) " found xyz",xyz
       goto 2 
      end if
      if (token.eq.'fix') then
        if(.not.inp_a(token)) call errquit(' no fix option',0,
     &         INPUT_ERR)
        if(token.eq.'linkbonds') then
          fixlinks=.true.
          goto 2
        elseif(token.eq.'qmregion') then
          fixqm=.true.
          goto 2
        else
          call errquit('qmmm_input: invalid link atom',0, INPUT_ERR)
        endif
       goto 2 
      end if
      if (token.eq.'xyzq') then
       xyzq = .true.
       write(*,*) " found xyzq",xyzq
       goto 2 
      end if
      if (token.eq.'qatoms') then
      if (.not.inp_f(qatoms))
     +  call errquit('qmmm_input: failed to read qatoms',0, INPUT_ERR)
      write(luout,*)' QM/MM atom charge:',qatoms
      call util_flush(luout)
      if (.not.rtdb_put(rtdb,'qmmm:qqmatm',mt_dbl,1,qatoms))
     + call errquit('qmmm_input: failed to store qmmm:uqmatm',0,
     &       INPUT_ERR)
      goto 2
      endif
c
      if(token.eq.'link') then
        if(.not.inp_a(token)) call errquit('qmmm_input: no link atom',0,
     &         INPUT_ERR)
        if(token.eq.'hydrogen') then
          linkatm=qmmm_hydrogen
          goto 2
        elseif(token.eq.'halogen') then
          linkatm=qmmm_halogen
          goto 2
        else
          call errquit('qmmm_input: invalid link atom',0, INPUT_ERR)
        endif
      endif
c
      if(token.eq.'optimize') then
        if(.not.inp_a(token)) then
         write(*,*) "please provide an option for optimize"
         call errquit('qmmm_input: optimize',0,
     &         INPUT_ERR)
        end if
        if(token.eq.'qm') then
          optimize='qm'
          goto 2
        elseif(token.eq.'mm') then
          optimize='mm'
          goto 2
        else
          call errquit('qmmm_input: invalid option',0, INPUT_ERR)
        endif
      endif
c
      if(token.eq.'interface') then
        if(.not.inp_a(token)) then
         write(*,*) "please provide an option for optimize"
         call errquit('qmmm_input: optimize',0,
     &         INPUT_ERR)
        end if
        if(token.eq.'qm') then
          interface='qm'
          goto 2
        elseif(token.eq.'mm') then
          interface='mm'
          goto 2
        else
          call errquit('qmmm_input: invalid option',0, INPUT_ERR)
        endif
      endif
cc
      if(token.eq.'bq_exclude') then
      if(.not.inp_a(token)) call errquit('qmmm_input:  bq_exclude ',0,
     &       INPUT_ERR)
      if(token.eq.'hydrogen') then
      bq_exclude=qmmm_hydrogen
      goto 2
      endif
      if(token.eq.'all') then
      bq_exclude=qmmm_all_atoms
      goto 2
      endif
      if(token.eq.'none') then
      bq_exclude=qmmm_no_atoms
      goto 2
      endif
      call errquit('qmmm_input: invalid bq_exclude option',0, INPUT_ERR)
      endif
c
      if (token.eq.'end') then
      if (.not.rtdb_put(rtdb,'qmmm:uqmatm',mt_dbl,1,eatoms))
     + call errquit('qmmm_input: failed to store qmmm:uqmatm',0,
     &       RTDB_ERR)
      if (.not.rtdb_put(rtdb,'qmmm:eatoms',mt_dbl,1,eatoms))
     + call errquit('qmmm_input: failed to store qmmm:eatoms',0,
     &       RTDB_ERR)
       if (.not.rtdb_put(rtdb,'qmmm:linkatm',mt_int,1,linkatm))
     + call errquit('qmmm_input: failed to store qmmm:linkatm',linkatm,
     &       RTDB_ERR)
      if (.not.rtdb_put(rtdb,'qmmm:xyz',mt_log,1,xyz))
     + call errquit('qmmm_input: failed to store qmmm:xyz',xyz,
     &       RTDB_ERR)
      if (.not.rtdb_put(rtdb,'qmmm:xyzq',mt_log,1,xyzq))
     + call errquit('qmmm_input: failed to store qmmm:xyzq',xyzq,
     &       RTDB_ERR)
      if (.not.rtdb_put(rtdb,'qmmm:bq_exclude',mt_int,1,bq_exclude))
     + call errquit('qmmm_input: failed to store qmmm:bq_exclude',
     &       bq_exclude,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'qmmm:fixlinks',mt_log,1,fixlinks))
     + call errquit('qmmm_input: failed to store qmmm:fixlinks',
     &       0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'qmmm:fixqm',mt_log,1,fixqm))
     + call errquit('qmmm_input: failed to store qmmm:fixqm',
     &       0,RTDB_ERR)
       if (.not.rtdb_cput(rtdb,'qmmm:optimize',1,optimize))
     + call errquit('qmmm_input: failed to store qmmm:optimize',
     &       0,RTDB_ERR)
       if (.not.rtdb_cput(rtdb,'qmmm:interface',1,interface))
     + call errquit('qmmm_input: failed to store qmmm:interface',
     &       0,RTDB_ERR)
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

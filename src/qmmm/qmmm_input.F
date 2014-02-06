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
      logical ltoken
      logical ignore
      double precision f
c
      pname = "qmmm_input: "
c
      write(luout,*) "in ",pname
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
c     convergence
c     -----------
      if (inp_compare(.false.,"convergence", token)) then
        if (.not.inp_f(f))
     >    call errquit(pname//'failed to read convergence',0, INPUT_ERR)
        tag = "qmmm:convergence"
        if (.not.rtdb_put(rtdb,tag(1:inp_strlen(tag)),mt_dbl,1,f)) 
     >      goto 999
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
c     nsamples
c     ----------
      if(inp_compare(.false.,'nsamples',token)) then
        if(.not.inp_i(itoken(1))) 
     +      call errquit('qmmm_input: nsamples',0,
     +         INPUT_ERR)
        tag = "qmmm:nsamples"
        if (.not.rtdb_put(rtdb,tag,mt_int,1,itoken)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
c     trajectory
c     ----------
      if(inp_compare(.false.,'trajectory',token)) then
        if(.not.inp_a(token)) then
         call errquit('qmmm_input: trajectory filename',0,
     &         INPUT_ERR)
        end if
        do ip=1,3
          itoken(ip) = ip
        end do
        np = inp_n_field()
        do ip=1,min(np-2,3)
          if(.not.inp_i(itoken(ip))) 
     +      call errquit('qmmm_input: trajectory',0,
     +         INPUT_ERR)
        end do
        tag = "qmmm:trajectory"
        if (.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        tag = "qmmm:trajectory_file"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        tag = "qmmm:trajectory_n"
        if (.not.rtdb_put(rtdb,tag,mt_int,3,itoken)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
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
      if(inp_compare(.false.,'density',token)) then
        if(.not.inp_a(token)) then
         call errquit('qmmm_input: density',0,
     &         INPUT_ERR)
        end if
        tag="qmmm:density_treatment"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        if(inp_compare(.false.,'static',token)) then
          tag="qmmm:frozen_density"
          if (.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        else if (inp_compare(.false.,'espfit',token)) then
          tag="qmmm:esp_density"
          if (.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        else if (inp_compare(.false.,'dynamical',token)) then
          tag="qmmm:frozen_density"
          ignore = rtdb_delete(rtdb,tag)
          tag="qmmm:esp_density"
          ignore = rtdb_delete(rtdb,tag)
        else
           call errquit(pname//'unknown token'//token,0,
     >       RTDB_ERR)
        end if
        goto 2
      end if  
c
c     option to load custom names
c     ---------------------------
      if(inp_compare(.false.,'rename',token)) then
        np = inp_n_field()
        if(np.ne.2) 
     >    call errquit(pname//'rename needs one argument',0,
     >       RTDB_ERR)
        do ip = 1,np-1
          if(.not.inp_a(mtoken(ip))) then
             call errquit(pname//" rename" ,0,
     &           INPUT_ERR)
          end if
        end do
        tag="qmmm:rename"
        if (.not.rtdb_cput(rtdb,tag,1,mtoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     option to load coordinates from rtdb geometry
c     ---------------------------------------------
      if(inp_compare(.false.,'geomrtdb',token)) then
        np = inp_n_field()
        if(np.ne.2) 
     >    call errquit(pname//'geomrtdb needs one argument',0,
     >       RTDB_ERR)
        do ip = 1,np-1
          if(.not.inp_a(mtoken(ip))) then
             call errquit(pname//" geom" ,0,
     &           INPUT_ERR)
          end if
        end do
        tag="qmmm:geomrtdb"
        if (.not.rtdb_cput(rtdb,tag,1,mtoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     option to load coordinates from xyzi file
c     -----------------------------------------
      if(inp_compare(.false.,'geomfile',token)) then
        np = inp_n_field()
        if(np.ne.2) 
     >    call errquit(pname//'geomfile needs one argument',0,
     >       RTDB_ERR)
        do ip = 1,np-1
          if(.not.inp_a(mtoken(ip))) then
             call errquit(pname//" geom" ,0,
     &           INPUT_ERR)
          end if
        end do
        tag="qmmm:geom_file"
        if (.not.rtdb_cput(rtdb,tag,1,mtoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif

c      if(inp_compare(.false.,'geom',token)) then
c        np = inp_n_field()
c        if(np.ne.2) 
c     >    call errquit(pname//'geom needs one argument',0,
c     >       RTDB_ERR)
c        do ip = 1,np-1
c          if(.not.inp_a(mtoken(ip))) then
c             call errquit(pname//" geom" ,0,
c     &           INPUT_ERR)
c          end if
c        end do
c        tag="qmmm:geom"
c        if (.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
c     >      call errquit(pname//'failed to store'//tag,0,
c     >       RTDB_ERR)
c        tag="qmmm:geom_file"
c        if (.not.rtdb_cput(rtdb,tag,1,mtoken(1))) 
c     >      call errquit(pname//'failed to store'//tag,0,
c     >       RTDB_ERR)
cc        tag="qmmm:geom_region"
cc        if (.not.rtdb_cput(rtdb,tag,1,mtoken(2))) 
cc     >      call errquit(pname//'failed to store'//tag,0,
cc     >       RTDB_ERR)
c          goto 2
c      endif
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
      if(inp_compare(.false.,'xyz',token)) then
        np = inp_n_field()
        do ip = 1,np-1
        mtoken(ip) = " "
        if(.not.inp_a(mtoken(ip))) then
         call errquit('qmmm_input',0,
     &         INPUT_ERR)
        end if
        if(mtoken(ip).eq."none") 
     $     mtoken(ip) = " "
        end do
        tag="qmmm:xyz_n"
        if (.not.rtdb_put(rtdb,tag,mt_int,1,np-1)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        tag="qmmm:xyz"
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
       if(inp_compare(.false.,'mm_charges',token)) then
600     continue
        if(.not.inp_a(token)) goto 1
        if(inp_compare(.false.,'exclude',token)) then
         tag="qmmm:bq_exclude"
         if(.not.inp_a(token)) goto 998
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
          call errquit('unknown option for charge',0,INPUT_ERR)
         end if
         if (.not.rtdb_put(rtdb,tag,mt_int,1,bq_exclude)) 
     >       call errquit(pname//'failed to store bq_exclude',0,
     >        RTDB_ERR)
        else if (inp_compare(.false.,'expand',token)) then
         if(.not.inp_a(token)) goto 998
         tag="qmmm:bq_expand"
         if(.not.rtdb_cput(rtdb,tag,1,token)) 
     >       goto 999
        else if (inp_compare(.false.,'save',token)) then
         if(.not.inp_a(token)) goto 998
         tag="qmmm:bq:save:file"
         if(.not.rtdb_cput(rtdb,tag,1,token)) 
     >       goto 999
         write(*,*) "found qmmm:bq:save:file",token
        else if (inp_compare(.false.,'load',token)) then
         if(.not.inp_a(token)) goto 998
         tag="qmmm:bq:load:file"
         if(.not.rtdb_cput(rtdb,tag,1,token)) 
     >       goto 999
        else if(inp_compare(.false.,'update',token)) then
          if(.not.inp_i(itoken(1))) goto 998
          tag="qmmm:bq_update"
          if (.not.rtdb_put(rtdb,tag,mt_int,1,itoken(1))) 
     >        goto 999
        end if
        goto 600
      end if
c
      if(inp_compare(.false.,'abc',token)) then
        call qmmm_input_abc(rtdb)     
        go to 1
      end if
c
      if(inp_compare(.false.,'pes',token)) then
        call qmmm_input_pes(rtdb)     
        go to 1
      end if
c
      if(inp_compare(.false.,'esp',token)) then
        call qmmm_input_esp(rtdb)     
        go to 1
      end if
c      if(inp_compare(.false.,'load',token)) then
c        if (.not.inp_a(token))
c     +    call errquit("please provide an option for load",
c     +                0, INPUT_ERR)
c        if(inp_compare(.false.,'bq',token)) then
c          tag = "qmmm:readbq"
c          if(.not.rtdb_put(rtdb,tag,mt_log,1,.true.))
c     >     goto 999
c          if(inp_a(token)) then
c            tag = "qmmm:bqfilename"
c            if(.not.rtdb_cput(rtdb,tag,1,token)) 
c     >       goto 999
c          end if
c        else if(inp_compare(.false.,'esp',token)) then
c          tag = "qmmm:readesp"
c          if(.not.rtdb_put(rtdb,tag,mt_log,1,.true.))
c     >     goto 999
c          if(inp_a(token)) then
c            tag = "qmmm:espfilename"
c            if(.not.rtdb_cput(rtdb,tag,1,token)) 
c     >       goto 999
c          end if
c        else
c          call errquit('unknown option for load',0,INPUT_ERR)
c        end if
c        goto 2
c      end if
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
c     
c     print level 
c     -------------------
      if (inp_compare(.false.,'print', token)) then
         call util_print_input(rtdb, "qmmm")
         go to 2
      end if
c
      if (token.eq.'end') then
        write(luout,*) "out of ",pname
        return
      endif
c
      write(luout,*)' unrecognized token in qmmm input:',
     + token(1:inp_strlen(token))
       call errquit(pname//'failed ',0,
     >        RTDB_ERR)
      return
998   call errquit(pname//'no token found '//token,0,
     >        RTDB_ERR)
999   call errquit(pname//'failed to store '//tag,0,
     >        RTDB_ERR)

      end

      subroutine qmmm_input_pes(rtdb)
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
      integer  itok
      integer  itoken(10)
      double precision   ftoken(10)
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
      pname = "qmmm_input_pes "
c
      write(luout,*) "in ",pname
c
      call inp_set_field(0)
c
c     start parsing input
c     ------------------
      if (.not.inp_a(token))
     + call errquit(pname,0, INPUT_ERR)
      if (.not.inp_compare(.false.,token,'pes'))
     + call errquit(pname,0, INPUT_ERR)
      goto 2
    1 continue
      if (.not.inp_read()) call errquit('qmmm_input: premature EOF',0,
     &       INPUT_ERR)
    2 continue
      if(.not.inp_a(token)) goto 1
c
c     number of points
c     -----------------
      if (inp_compare(.false.,"npoints",token)) then
        if (.not.inp_i(itok))
     >    call errquit(pname//'failed to read npoints',0, INPUT_ERR)
        tag = "qmmm:pes_npoints"
        if (.not.rtdb_put(rtdb,tag,mt_int,1,itok)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
c     starting point
c     --------------
      if (inp_compare(.false.,"istart",token)) then
        if (.not.inp_i(itok))
     >    call errquit(pname//'failed to read istart',0, INPUT_ERR)
        tag = "qmmm:pes_istart"
        if (.not.rtdb_put(rtdb,tag,mt_int,1,itok)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
c     end point
c     --------------
      if (inp_compare(.false.,"iend",token)) then
        if (.not.inp_i(itok))
     >    call errquit(pname//'failed to read iend',0, INPUT_ERR)
        tag = "qmmm:pes_iend"
        if (.not.rtdb_put(rtdb,tag,mt_int,1,itok)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
c     direction
c     ---------
      if (inp_compare(.false.,"reverse",token)) then
        tag = "qmmm:pes_backward"
        if (.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
      if (inp_compare(.false.,"forward",token) ) then
        tag = "qmmm:pes_backward"
        if (.not.rtdb_put(rtdb,tag,mt_log,1,.false.)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
c     rst files
c     --------
      if(inp_compare(.false.,'rst',token) ) then
        if(.not.inp_a(token)) 
     >   call errquit(pname//' rst ',0,INPUT_ERR)
        tag="qmmm:pes_rst_input"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
      endif
c
c     orbital files
c     ------------
      if(inp_compare(.false.,'vectors',token) ) then
        if(.not.inp_a(token)) 
     >   call errquit(pname//' vectors ',0,INPUT_ERR)
        tag="qmmm:pes_vectors_input"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
      endif
c
c
      if(inp_compare(.false.,'spring',token)) then
        np = inp_n_field()
        if(np.ne.3)
     >   call errquit(pname," spring",0,INPUT_ERR)
        do ip = 1,np-1
        ftoken(ip) = 0.0d0
        if(.not.inp_f(ftoken(ip))) then
         call errquit(pname//'spring',0,INPUT_ERR)
        end if
        end do
        tag="qmmm:pes_bounds"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,2,ftoken)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
      if (token.eq.'end') then
        write(luout,*) "out of ",pname
        return
      endif
c
      write(luout,*)' unrecognized token in qmmm input:',
     + token(1:inp_strlen(token))
       call errquit(pname//'failed ',0,
     >        RTDB_ERR)
      return
998   call errquit(pname//'no token found '//token,0,
     >        RTDB_ERR)
999   call errquit(pname//'failed to store '//tag,0,
     >        RTDB_ERR)

      end

      subroutine qmmm_input_abc(rtdb)
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
      integer  itok
      integer  itoken(10)
      double precision   ftoken(10)
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
      pname = "qmmm_input_abc "
c
      write(luout,*) "in ",pname
c
      call inp_set_field(0)
c
c     start parsing input
c     ------------------
      if (.not.inp_a(token))
     + call errquit(pname,0, INPUT_ERR)
      if (.not.inp_compare(.false.,token,'abc'))
     + call errquit(pname,0, INPUT_ERR)
      goto 2
    1 continue
      if (.not.inp_read()) call errquit('qmmm_input: premature EOF',0,
     &       INPUT_ERR)
    2 continue
      if(.not.inp_a(token)) goto 1
c
c     number of cycles
c     ----------------
      if (inp_compare(.false.,"ncycles",token)) then
        if (.not.inp_i(itok))
     >    call errquit(pname//'failed to read npoints',0, INPUT_ERR)
        tag = "qmmm:abc_ncycles"
        if (.not.rtdb_put(rtdb,tag,mt_int,1,itok)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
c     new or extend
c     -------------
      if (inp_compare(.false.,"new",token)) then
        tag = "qmmm:abc_extend"
        if (.not.rtdb_put(rtdb,tag,mt_log,1,.false.)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
      if (inp_compare(.false.,"extend",token)) then
        tag = "qmmm:abc_extend"
        if (.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
        goto 2
      endif
c
c     theory
c     ------
      if(inp_compare(.false.,'theory',token) ) then
        do ip=1,2
          if(.not.inp_a(mtoken(ip))) 
     >     call errquit(pname//' rst ',0,INPUT_ERR)
        end do
        tag="qmmm:abc_theory"
        if (.not.rtdb_cput(rtdb,tag,2,mtoken)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
      endif
c
c     topology
c     --------
      if(inp_compare(.false.,'topology',token) ) then
        do ip=1,2
          if(.not.inp_a(mtoken(ip))) 
     >     call errquit(pname//' rst ',0,INPUT_ERR)
        end do
        tag="qmmm:abc_top"
        if (.not.rtdb_cput(rtdb,tag,2,mtoken)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
      endif
c
      if (token.eq.'end') then
        write(luout,*) "out of ",pname
        return
      endif
c
      write(luout,*)' unrecognized token in qmmm input:',
     + token(1:inp_strlen(token))
       call errquit(pname//'failed ',0,
     >        RTDB_ERR)
      return
998   call errquit(pname//'no token found '//token,0,
     >        RTDB_ERR)
999   call errquit(pname//'failed to store '//tag,0,
     >        RTDB_ERR)

      end

      subroutine qmmm_input_esp(rtdb)
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
      character*255 token0
c
      pname = "qmmm_input_esp"
c
600   continue
      if(.not.inp_a(token)) goto 20
      if (inp_compare(.false.,'save',token)) then
       token0 = token
       tag="qmmm:esp:save"
       if(.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
     >     goto 999
      else if (inp_compare(.false.,'load',token)) then
       token0 = token
       tag="qmmm:readesp"
       if(.not.rtdb_put(rtdb,tag,mt_log,1,.true.)) 
     >     goto 999
      else
       if(inp_compare(.false.,'load',token0)) then
         tag="qmmm:espfilename"
         if(.not.rtdb_cput(rtdb,tag,1,token)) 
     >     goto 999
       else if (inp_compare(.false.,'save',token0)) then
         tag="qmmm:esp:save:file"
         if(.not.rtdb_cput(rtdb,tag,1,token)) 
     >     goto 999
       else
        goto 997
       end if
      end if
      goto 600
20    continue
      return
998   call errquit(pname//'no token found '//token,0,
     >        RTDB_ERR)
999   call errquit(pname//'failed to store '//tag,0,
     >        RTDB_ERR)
997   call errquit(pname//'invalid token '//token,0,
     >        RTDB_ERR)

      end
c $Id$

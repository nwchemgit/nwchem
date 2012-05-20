      subroutine smd_input(rtdb)
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
      character*32 pname
      character*255 token
      character*80 mtoken(10)
      integer  itoken(10)
      double precision  ftoken(10)
      integer ip,np
c
      pname = "smd_input: "
c
c      write(luout,*) "in ",pname
c
      call inp_set_field(0)
c
c     start parsing input
c     ------------------
      if (.not.inp_a(token))
     + call errquit(pname//'no input available',0, INPUT_ERR)
      if (.not.inp_compare(.false.,token,'smd'))
     + call errquit('smd_input: no input available',0, INPUT_ERR)
      goto 2
    1 continue
      if (.not.inp_read()) call errquit('smd_input: premature EOF',0,
     &       INPUT_ERR)
    2 continue
      if(.not.inp_a(token)) goto 1
c
c     coordinates
c     ----------
      if(inp_compare(.false.,'coord',token)) then
        np = 1
        do ip = 1,np
          if(.not.inp_a(mtoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:coordfile"
        if (.not.rtdb_cput(rtdb,tag,1,mtoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     parameters
c     ----------
      if(inp_compare(.false.,'param',token)) then
        np = 1
        do ip = 1,np
          if(.not.inp_a(mtoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:paramfile"
        if (.not.rtdb_cput(rtdb,tag,1,mtoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     kvec
c     -----
      if(inp_compare(.false.,'kvec',token)) then
        np = 3
        do ip = 1,np
          if(.not.inp_i(itoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:kvec"
        if (.not.rtdb_put(rtdb,tag,mt_int,3,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     lat_a
c     -----
      if(inp_compare(.false.,'lat_a',token)) then
        np = 3
        do ip = 1,np
          if(.not.inp_f(ftoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:lat_a"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,3,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     lat_b
c     -----
      if(inp_compare(.false.,'lat_b',token)) then
        np = 3
        do ip = 1,np
          if(.not.inp_f(ftoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:lat_b"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,3,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     lat_c
c     -----
      if(inp_compare(.false.,'lat_c',token)) then
        np = 3
        do ip = 1,np
          if(.not.inp_f(ftoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:lat_c"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,3,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     ndata
c     -----
      if(inp_compare(.false.,'ndata',token)) then
        np = 1
        do ip = 1,np
          if(.not.inp_i(itoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:ndata"
        if (.not.rtdb_put(rtdb,tag,mt_int,np,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     nequil
c     -----
      if(inp_compare(.false.,'nequil',token)) then
        np = 1
        do ip = 1,np
          if(.not.inp_i(itoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:nequil"
        if (.not.rtdb_put(rtdb,tag,mt_int,np,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     nprint
c     -----
      if(inp_compare(.false.,'nprint',token)) then
        np = 1
        do ip = 1,np
          if(.not.inp_i(itoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:nprint"
        if (.not.rtdb_put(rtdb,tag,mt_int,np,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     step
c     -----
      if(inp_compare(.false.,'step',token)) then
        np = 1
        do ip = 1,np
          if(.not.inp_f(ftoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:step"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     rcut
c     ----
      if(inp_compare(.false.,'rcut',token)) then
        np = 1
        do ip = 1,np
          if(.not.inp_f(ftoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:rcut"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     temp
c     -----
      if(inp_compare(.false.,'temp',token)) then
        np = 1
        do ip = 1,np
          if(.not.inp_f(ftoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:temp_target"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     ewald
c     -----
      if(inp_compare(.false.,'ewald',token)) then
        np = 1
        do ip = 1,np
          if(.not.inp_f(ftoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:ewald"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c
c     verlet
c     ------
      if(inp_compare(.false.,'verlet',token)) then
        np = 1
        do ip = 1,np
          if(.not.inp_f(ftoken(ip))) then
             call errquit(pname//token,0,
     &           INPUT_ERR)
          end if
        end do
        tag="smd:rcut_verlet"
        if (.not.rtdb_put(rtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          goto 2
      endif
c     
c     print level 
c     -------------------
      if (inp_compare(.false.,'print', token)) then
         call util_print_input(rtdb, "smd")
         go to 2
      end if
c
      if (token.eq.'end') then
c        write(luout,*) "out of ",pname
c        if(.not.rtdb_print(rtdb,.true.))
c     >      call errquit(pname//'failed to print rtdb',0,
c     >       RTDB_ERR)
        return
      endif
c
      write(luout,*)' unrecognized token in smd input:',
     + token(1:inp_strlen(token))
       call errquit(pname//'failed ',0,
     >        RTDB_ERR)
      return
998   call errquit(pname//'no token found '//token,0,
     >        RTDB_ERR)
999   call errquit(pname//'failed to store '//tag,0,
     >        RTDB_ERR)

      end


c $Id$

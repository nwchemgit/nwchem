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
c     charges
c     ----------
      if(inp_compare(.false.,'charge',token)) then
        call smd_input_charge(rtdb)
        goto 2
      endif
c
c     velocities
c     ----------
      if(inp_compare(.false.,'veloc',token)) then
        call smd_input_veloc(rtdb)
        goto 2
      endif
c
c     coordinates
c     ----------
      if(inp_compare(.false.,'coord',token)) then
        call smd_input_coord(rtdb)
        goto 2
      endif
cc      if(inp_compare(.false.,'coord',token)) then
c        np = 1
c        do ip = 1,np
c          if(.not.inp_a(mtoken(ip))) then
c             call errquit(pname//token,0,
c     &           INPUT_ERR)
c          end if
c        end do
c        tag="smd:coordfile"
c        if (.not.rtdb_cput(rtdb,tag,1,mtoken(1))) 
c     >      call errquit(pname//'failed to store'//tag,0,
c     >       RTDB_ERR)
c          goto 2
c      endif
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
        tag="smd:verlet"
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
c        write(*,*) "RTDB after smd"
c        if(.not.rtdb_print(rtdb,.true.))
c     >      call errquit(pname//'failed to print rtdb',0,
c     >       RTDB_ERR)
        return
      endif
c
        if(.not.rtdb_print(rtdb,.true.))
     >      call errquit(pname//'failed to print rtdb',0,
     >       RTDB_ERR)
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

      subroutine smd_input_coord(rtdb)
c
      implicit none
#include "errquit.fh"
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
      integer ip,np
c
      pname = "smd_input_coord"
c
      write(luout,*) "in ",pname
c
      np = inp_n_field()
c
      if(np.eq.1) goto 200
c
      call inp_set_field(1)
c
c     start parsing input
c     ------------------
2     continue
      if(inp_cur_field().eq.np) return
      if (.not.inp_a(token))
     + call errquit(pname,0, INPUT_ERR)
c
c     input files
c     -----------
      if(inp_compare(.false.,'input',token) ) then
        if(.not.inp_a(token)) 
     >   call errquit(pname//' input ',0,INPUT_ERR)
        tag="smd:coord:input"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
      endif
c
c     output files
c     -----------
      if(inp_compare(.false.,'output',token) ) then
        if(.not.inp_a(token)) 
     >   call errquit(pname//' input ',0,INPUT_ERR)
        tag="smd:coord:output"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
      endif


200   continue
      write(luout,*) "out ",pname
      return
      end

      subroutine smd_input_veloc(rtdb)
c
      implicit none
#include "errquit.fh"
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
      integer ip,np
c
      pname = "smd_input_vel"
c
      write(luout,*) "in ",pname
c
      np = inp_n_field()
c
      if(np.eq.1) goto 200
c
      call inp_set_field(1)
c
c     start parsing input
c     ------------------
2     continue
      if(inp_cur_field().eq.np) return
      if (.not.inp_a(token))
     + call errquit(pname,0, INPUT_ERR)
c
c     input files
c     -----------
      if(inp_compare(.false.,'input',token) ) then
        if(.not.inp_a(token)) 
     >   call errquit(pname//' input ',0,INPUT_ERR)
        tag="smd:veloc:input"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
      endif
c
c     output files
c     -----------
      if(inp_compare(.false.,'output',token) ) then
        if(.not.inp_a(token)) 
     >   call errquit(pname//' input ',0,INPUT_ERR)
        tag="smd:veloc:output"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
      endif


200   continue
      write(luout,*) "out ",pname
      return
      end

      subroutine smd_input_charge(rtdb)
c
      implicit none
#include "errquit.fh"
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
      integer ip,np
c
      pname = "smd_input_charge"
c
      write(luout,*) "in ",pname
c
      np = inp_n_field()
c
      if(np.eq.1) goto 200
c
      call inp_set_field(1)
c
c     start parsing input
c     ------------------
2     continue
      if(inp_cur_field().eq.np) return
      if (.not.inp_a(token))
     + call errquit(pname,0, INPUT_ERR)
c
c     input files
c     -----------
      if(inp_compare(.false.,'input',token) ) then
        if(.not.inp_a(token)) 
     >   call errquit(pname//' input ',0,INPUT_ERR)
        tag="smd:charge:input"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
      endif
c
c     output files
c     -----------
      if(inp_compare(.false.,'output',token) ) then
        if(.not.inp_a(token)) 
     >   call errquit(pname//' input ',0,INPUT_ERR)
        tag="smd:charge:output"
        if (.not.rtdb_cput(rtdb,tag,1,token)) 
     >        call errquit(pname//'failed to store'//tag,0,
     >         RTDB_ERR)

        goto 2
      endif


200   continue
      write(luout,*) "out ",pname
      return
      end
c $Id$

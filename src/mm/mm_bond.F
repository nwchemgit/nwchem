      block data initial_mm_bond_data
      implicit none
#include "mm_bond_data.fh"
 
      data nbon /0/
      data nang /0/
      data ndih /0/

      end

      subroutine mm_bndparm_init(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh" 
      integer rtdb 

      character*180 filename
      character*180 message
      character*30  pname
      integer       fn
      integer       nb, na, nd     

      pname = "mm_bndparm_load"
c     write(*,*) pname

c     deallocate all previous allocated arrays just in case
      call mm_bonded_deallocate()
          
      call mm_get_nbon(rtdb,nb)
      call mm_get_nang(rtdb,na)
      call mm_get_ndih(rtdb,nd)

c     write(*,*) 'nbon = ', nb, 'nang = ', na, 'ndih = ', nd

      if(nb.gt.0) call mm_bond_allocate(nb)
      if(na.gt.0) call mm_angl_allocate(na)
      if(nd.gt.0) call mm_dihe_allocate(nd)

      call mm_bond_parm_load(rtdb)
      call mm_angl_parm_load(rtdb)
      call mm_dihe_parm_load(rtdb)
 
c     call mm_bond_test()
c     call mm_angl_test()
c     call mm_dihe_test()

      return

911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_get_nbon(rtdb,nbon0)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb 
      integer nbon0

      character*180 filename
      character*255 field
      character*180 message
      character*30  pname
      integer       fn
      integer       nb, na, nd     
        
      pname = "mm_get_nbon"   
   
      if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,filename))
     & goto 911

      call mm_open_file(filename,fn)
      
      call inp_init(fn,LuOut)
      nbon0 = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'bond',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.4) then
                nbon0 = nbon0 + 1
                goto 10
              end if
            end if
            goto 10
          end if
        end if
        goto 5
      end if
20    continue
      close(fn)
      call inp_init(LuIn,LuOut)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_get_nang(rtdb,nang0)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb 
      integer nang0

      character*180 filename
      character*255 field
      character*180 message
      character*30  pname
      integer       fn
      integer       nb, na, nd     
        
      pname = "mm_get_nang"   
   
      if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,filename))
     & goto 911

      call mm_open_file(filename,fn)
      
      call inp_init(fn,LuOut)
      nang0 = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'angle',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.5) then
                nang0 = nang0 + 1
                goto 10
              end if
            end if
            goto 10
          end if
        end if
        goto 5
      end if
20    continue
      close(fn)
      call inp_init(LuIn,LuOut)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_get_ndih(rtdb,ndih0)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb 
      integer ndih0

      character*180 filename
      character*255 field
      character*180 message
      character*30  pname
      integer       fn
      integer       nb, na, nd     
        
      pname = "mm_get_ndih"   
   
      if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,filename))
     & goto 911

      call mm_open_file(filename,fn)
      
      call inp_init(fn,LuOut)
      ndih0 = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'dihedral',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.7) then
                ndih0 = ndih0 + 1
                goto 10
              end if
            end if
            goto 10
          end if
        end if
        goto 5
      end if
20    continue
      close(fn)
      call inp_init(LuIn,LuOut)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_bond_parm_load(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh" 
      integer rtdb 

      character*180 filename
      character*255 field
      character*180 message
      character*30  pname
      integer       fn
      integer       nb, na, nd     
        
      pname = "mm_bond_parm_load"   
   
      if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,filename))
     & goto 911

      call mm_open_file(filename,fn)
      
      call inp_init(fn,LuOut)
      nb = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'bond',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.4) then
                message = "checking for nbon bounds"
                if(nb.gt.nbon) goto 911
                nb = nb + 1
                if(.not.inp_i(int_mb(i_ibon+nb-1))) goto 911 
                if(.not.inp_i(int_mb(i_jbon+nb-1))) goto 911 
                if(.not.inp_f(dbl_mb(i_kb+nb-1)))   goto 911 
                if(.not.inp_f(dbl_mb(i_rb+nb-1)))   goto 911
                goto 10
              end if
            end if
            goto 10
          end if
        end if
        goto 5
      end if
20    continue
      close(fn)
      call inp_init(LuIn,LuOut)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_angl_parm_load(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh" 
      integer rtdb 

      character*180 filename
      character*255 field
      character*180 message
      character*30  pname
      integer       fn
      integer       nb, na, nd     
        
      pname = "mm_angle_parm_load"   
   
      if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,filename))
     & goto 911

      call mm_open_file(filename,fn)
      
      call inp_init(fn,LuOut)
      na = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'angle',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.5) then
                message = "checking for nang bounds"
                if(na.gt.nang) goto 911
                na = na + 1
                if(.not.inp_i(int_mb(i_iang+na-1))) goto 911 
                if(.not.inp_i(int_mb(i_jang+na-1))) goto 911 
                if(.not.inp_i(int_mb(i_kang+na-1))) goto 911 
                if(.not.inp_f(dbl_mb(i_kthe+na-1))) goto 911 
                if(.not.inp_f(dbl_mb(i_thet+na-1))) goto 911 
                goto 10
              end if
            end if
            goto 10
          end if
        end if
        goto 5
      end if
20    continue
      close(fn)
      call inp_init(LuIn,LuOut)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_dihe_parm_load(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh" 
      integer rtdb 

      character*180 filename
      character*255 field
      character*180 message
      character*30  pname
      integer       fn
      integer       nb, na, nd     
        
      pname = "mm_dihe_parm_load"   
   
      if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,filename))
     & goto 911

      call mm_open_file(filename,fn)
      
      call inp_init(fn,LuOut)
      nd = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'dihedral',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.7) then
                message = "checking for ndih bounds"
                if(nd.gt.ndih) goto 911
                nd = nd + 1
                if(.not.inp_i(int_mb(i_idih+nd-1))) goto 911 
                if(.not.inp_i(int_mb(i_jdih+nd-1))) goto 911 
                if(.not.inp_i(int_mb(i_kdih+nd-1))) goto 911 
                if(.not.inp_i(int_mb(i_ldih+nd-1))) goto 911 
                if(.not.inp_f(dbl_mb(i_kphi+nd-1))) goto 911 
                if(.not.inp_f(dbl_mb(i_peri+nd-1))) goto 911 
                if(.not.inp_f(dbl_mb(i_phas+nd-1))) goto 911 
                goto 10
              end if
            end if
            goto 10
          end if
        end if
        goto 5
      end if
20    continue
      close(fn)
      call inp_init(LuIn,LuOut)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_bond_allocate(n)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh" 
      integer n
      
      character*30 pname

      pname = "mm_bond_allocate"
 
      if(n.ne.nbon) then
        call mm_bond_end()
        nbon = n
 
        if (.not. ma_alloc_get(MT_INT, nbon,
     &       'bond scratch i',
     &       h_ibon, i_ibon) ) call errquit(
     &       'failed scratch space '//pname,
     &       nbon, MA_ERR)
 
        if (.not. ma_alloc_get(MT_INT, nbon,
     &       'bond scratch j',
     &       h_jbon, i_jbon) ) call errquit(
     &       'failed scratch space '//pname,
     &       nbon, MA_ERR)
 
        if (.not. ma_alloc_get(MT_DBL, nbon,
     &       'bond scratch kb',
     &       h_kb, i_kb) ) call errquit(
     &       'failed scratch space '//pname,
     &       nbon, MA_ERR)
 
        if (.not. ma_alloc_get(MT_DBL, nbon,
     &       'bond scratch rb',
     &       h_rb, i_rb) ) call errquit(
     &       'failed scratch space '//pname,
     &       nbon, MA_ERR)
      end if
      call ifill(n,0.0,int_mb(i_ibon),1)
      call ifill(n,0.0,int_mb(i_jbon),1)
      call dfill(n,0.0d0,dbl_mb(i_kb),1)
      call dfill(n,0.0d0,dbl_mb(i_rb),1)
   
      return

      end

      subroutine mm_angl_allocate(n)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh" 
      integer n
      
      character*30 pname

      pname = "mm_angl_allocate"
 
      if(n.ne.nang) then
        call mm_angl_end()
        nang = n
 
        if (.not. ma_alloc_get(MT_INT, nang,
     &       'angle scratch i',
     &       h_iang, i_iang) ) call errquit(
     &       'failed scratch space '//pname,
     &       nang, MA_ERR)
 
        if (.not. ma_alloc_get(MT_INT, nang,
     &       'angle scratch j',
     &       h_jang, i_jang) ) call errquit(
     &       'failed scratch space '//pname,
     &       nang, MA_ERR)
 
        if (.not. ma_alloc_get(MT_INT, nang,
     &       'angle scratch k',
     &       h_kang, i_kang) ) call errquit(
     &       'failed scratch space '//pname,
     &       nang, MA_ERR)
 
        if (.not. ma_alloc_get(MT_DBL, nang,
     &       'angle scratch kthe',
     &       h_kthe, i_kthe) ) call errquit(
     &       'failed scratch space '//pname,
     &       nang, MA_ERR)
 
        if (.not. ma_alloc_get(MT_DBL, nang,
     &       'angle scratch thet',
     &       h_thet, i_thet) ) call errquit(
     &       'failed scratch space '//pname,
     &       nang, MA_ERR)
      end if
      call ifill(n,0.0,int_mb(i_iang),1)
      call ifill(n,0.0,int_mb(i_jang),1)
      call ifill(n,0.0,int_mb(i_kang),1)
      call dfill(n,0.0d0,dbl_mb(i_kthe),1)
      call dfill(n,0.0d0,dbl_mb(i_thet),1)
   
      return
      end

      subroutine mm_dihe_allocate(n)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh" 
      integer n
      
      character*30 pname

      pname = "mm_dihe_allocate"
 
      if(n.ne.ndih) then
        call mm_dihe_end()
        ndih = n
 
        if (.not. ma_alloc_get(MT_INT, ndih,
     &       'dihedral scratch i',
     &       h_idih, i_idih) ) call errquit(
     &       'failed scratch space '//pname,
     &       ndih, MA_ERR)
 
        if (.not. ma_alloc_get(MT_INT, ndih,
     &       'dihedral scratch j',
     &       h_jdih, i_jdih) ) call errquit(
     &       'failed scratch space '//pname,
     &       ndih, MA_ERR)
 
        if (.not. ma_alloc_get(MT_INT, ndih,
     &       'dihedral scratch k',
     &       h_kdih, i_kdih) ) call errquit(
     &       'failed scratch space '//pname,
     &       ndih, MA_ERR)
 
        if (.not. ma_alloc_get(MT_INT, ndih,
     &       'dihedral scratch l',
     &       h_ldih, i_ldih) ) call errquit(
     &       'failed scratch space '//pname,
     &       ndih, MA_ERR)
 
        if (.not. ma_alloc_get(MT_DBL, ndih,
     &       'dihedral scratch kphi',
     &       h_kphi, i_kphi) ) call errquit(
     &       'failed scratch space '//pname,
     &       ndih, MA_ERR)
 
        if (.not. ma_alloc_get(MT_DBL, ndih,
     &       'dihedral scratch periodicity',
     &       h_peri, i_peri) ) call errquit(
     &       'failed scratch space '//pname,
     &       ndih, MA_ERR)
 
        if (.not. ma_alloc_get(MT_DBL, ndih,
     &       'dihedral scratch phase',
     &       h_phas, i_phas) ) call errquit(
     &       'failed scratch space '//pname,
     &       ndih, MA_ERR)
      end if
      call ifill(n,0.0,int_mb(i_idih),1)
      call ifill(n,0.0,int_mb(i_jdih),1)
      call ifill(n,0.0,int_mb(i_kdih),1)
      call ifill(n,0.0,int_mb(i_ldih),1)
      call dfill(n,0.0d0,dbl_mb(i_kphi),1)
      call dfill(n,0.0d0,dbl_mb(i_peri),1)
      call dfill(n,0.0d0,dbl_mb(i_phas),1)
   
      return
      end

      subroutine mm_bonded_deallocate()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh" 
        call mm_bond_end()
        call mm_angl_end()
        call mm_dihe_end()
      end

      subroutine mm_bond_end()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh"
      
      character*30 pname
      pname = "mm_bond_end"

      if (nbon .gt. 0) then
        if (.not.ma_free_heap(h_ibon)) goto 911
        if (.not.ma_free_heap(h_jbon)) goto 911
        if (.not.ma_free_heap(h_kb))   goto 911
        if (.not.ma_free_heap(h_rb))   goto 911
        nbon = 0
      end if
      return
911   call errquit("error "//trim(pname),0,-1)
      return
      end

      subroutine mm_angl_end()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh"
      
      character*30 pname
      pname = "mm_angl_end"

      if (nang .gt. 0) then
        if (.not.ma_free_heap(h_iang)) goto 911
        if (.not.ma_free_heap(h_jang)) goto 911
        if (.not.ma_free_heap(h_kang)) goto 911
        if (.not.ma_free_heap(h_kthe)) goto 911
        if (.not.ma_free_heap(h_thet)) goto 911
        nang = 0
      end if
      return
911   call errquit("error "//trim(pname),0,-1)
      return
      end

      subroutine mm_dihe_end()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh"
      
      character*30 pname
      pname = "mm_dihe_end"

      if (ndih .gt. 0) then
        if (.not.ma_free_heap(h_idih)) goto 911
        if (.not.ma_free_heap(h_jdih)) goto 911
        if (.not.ma_free_heap(h_kdih)) goto 911
        if (.not.ma_free_heap(h_ldih)) goto 911
        if (.not.ma_free_heap(h_kphi)) goto 911
        if (.not.ma_free_heap(h_peri)) goto 911
        if (.not.ma_free_heap(h_phas)) goto 911
        ndih = 0
      end if
      return
911   call errquit("error "//trim(pname),0,-1)
      return
      end

      subroutine mm_bond_test()
      implicit none
#include "mm_bond_data.fh"
#include "mafdecls.fh"

      integer i

      do i=1,nbon
         write(6,'(2(1X,I3),1X,2(1X,ES12.6))') 
     $        (int_mb(i_ibon+i-1)),
     $        (int_mb(i_jbon+i-1)),
     $        (dbl_mb(i_kb+i-1)),
     $        (dbl_mb(i_rb+i-1))
      end do
      end

      subroutine mm_angl_test()
      implicit none
#include "mm_bond_data.fh"
#include "mafdecls.fh"

      integer i

      do i=1,nang
         write(6,'(3(1X,I3),1X,2(1X,ES12.6))') 
     $        (int_mb(i_iang+i-1)),
     $        (int_mb(i_jang+i-1)),
     $        (int_mb(i_kang+i-1)),
     $        (dbl_mb(i_kthe+i-1)),
     $        (dbl_mb(i_thet+i-1))
      end do
      end

      subroutine mm_dihe_test()
      implicit none
#include "mm_bond_data.fh"
#include "mafdecls.fh"

      integer i

      do i=1,ndih
         write(6,'(4(1X,I3),1X,3(1X,ES12.6))') 
     $        (int_mb(i_idih+i-1)),
     $        (int_mb(i_jdih+i-1)),
     $        (int_mb(i_kdih+i-1)),
     $        (int_mb(i_ldih+i-1)),
     $        (dbl_mb(i_kphi+i-1)),
     $        (dbl_mb(i_peri+i-1)),
     $        (dbl_mb(i_phas+i-1))
      end do
      end

      subroutine mm_bond_add_energy(rtdb,e)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
      integer rtdb
      double precision e

      double precision scale_energy,scale_length
      double precision scale_force
      double precision ebon
      character*30 message
      character*30 pname
      character*16 thetag
      double precision crd(3),q
      integer i,j

      pname = "mm_bond_add_energy"
c     write(*,*) pname 

      call mm_bond_energy(nqm,nbtot,dbl_mb(i_br),log_mb(i_lqm),
     &                   int_mb(i_indx),nbon,int_mb(i_ibon),
     &                   int_mb(i_jbon),dbl_mb(i_kb),
     &                   dbl_mb(i_rb),ebon)

      call util_convert_units("kcal","au",scale_energy)
 
      ebon = ebon*scale_energy
      e = e + ebon

      return
911   call errquit("error "//trim(message),0,-1)

      end

      subroutine mm_bond_add_egrad(rtdb,e,n,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
      integer rtdb
      double precision e
      integer n
      double precision g(3,n)

      integer i_g, h_g
      double precision scale_energy,scale_length
      double precision scale_force
      double precision ebon
      character*30 message
      character*30 pname
      character*16 thetag
      double precision crd(3),q
      integer i,j

      pname = "mm_bond_add_egrad"
c     write(*,*) pname 

      if ( .not. ma_push_get( MT_DBL, 3*n, 
     &      'bond grad scratch ',
     &      h_g, i_g) ) call errquit(
     &      'failed scratch space ',
     &      nqm, MA_ERR)

      call dfill(3*n,0.0d0,dbl_mb(i_g),1)
            
      call mm_bond_egrad(n,nbtot,dbl_mb(i_br),log_mb(i_lqm),
     &                   int_mb(i_indx),nbon,int_mb(i_ibon),
     &                   int_mb(i_jbon),dbl_mb(i_kb),
     &                   dbl_mb(i_rb),ebon,dbl_mb(i_g))

      call util_convert_units("kcal","au",scale_energy)
      call util_convert_units("ang","au",scale_length)
 
      ebon = ebon*scale_energy
      e = e + ebon
      scale_force = scale_energy/scale_length      

      call daxpy(3*n,scale_force,dbl_mb(i_g),1,g,1)

      message = "memory deallocation"
      if(.not.ma_pop_stack(h_g)) goto 911
      return
911   call errquit("error "//trim(message),0,-1)

      end

      subroutine mm_bond_energy(ng,nbtot,br,lqm,indx,
     &                               nbon,ibon,jbon,kb,rb,
     &                               ebon)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
      integer ng
      integer nbtot
      double precision br(3,nbtot)
      logical lqm(nbtot)
      integer indx(nbtot)
      integer nbon
      integer ibon(nbon)
      integer jbon(nbon)
      double precision kb(nbon)
      double precision rb(nbon)
      double precision ebon

      character*16 tag
      character*30 pname
      character*30 message
      integer i, n, m
      integer iqm
      integer i_indx, j_indx
      integer i1, i2
      double precision xij,yij,zij
      double precision da,de,df
      double precision rr
      double precision e_bon
      logical is_qm

      pname = "mm_bond_energy"
c     write(*,*) pname
      
      ebon = 0.0d0
      i1 = 0 
      i2 = 0

      do n=1,nbon
        i_indx = ibon(n)
        j_indx = jbon(n)
c       first obtain indx i1, i2 in bond coord array
        do m=1,nbtot
          if(i_indx.eq.indx(m)) i1 = m
          if(j_indx.eq.indx(m)) i2 = m
        end do

c       Now doing energy and grad calculations
        message = "checking bounds for bond"
        if(i1.eq.0 .or. i2.eq.0) then
c         write(*,*) i1,i2
          goto 911
        end if

        xij = br(1,i2) - br(1,i1)
        yij = br(2,i2) - br(2,i1)
        zij = br(3,i2) - br(3,i1)
 
        rr = xij*xij + yij*yij + zij*zij
        rr = sqrt(rr)
          
        da = rr - rb(n)
        de = kb(n)*da
        e_bon = de*da
        ebon = ebon + e_bon
      end do
       
c     write(*,*) 'ebon = ', ebon, ' kcal/mol'
      return
911   call errquit("error "//trim(message),0,-1)

      end

      subroutine mm_bond_egrad(ng,nbtot,br,lqm,indx,
     &                               nbon,ibon,jbon,kb,rb,
     &                               ebon,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_geom_data.fh"
#include "mm_coords_data.fh"
      integer ng
      integer nbtot
      double precision br(3,nbtot)
      logical lqm(nbtot)
      integer indx(nbtot)
      integer nbon
      integer ibon(nbon)
      integer jbon(nbon)
      double precision kb(nbon)
      double precision rb(nbon)
      double precision ebon
      double precision g(3,ng)

      character*16 tag
      character*30 pname
      character*30 message
      integer i, n, m
      integer iqm
      integer i_indx, j_indx
      integer i1, i2
      double precision xij,yij,zij
      double precision da,de,df
      double precision rr
      double precision e_bon
      logical is_qm

      pname = "mm_bond_egrad"
c     write(*,*) pname
      
      ebon = 0.0d0

      do n=1,nbon
        i_indx = ibon(n)
        j_indx = jbon(n)
        i1 = 0 
        i2 = 0
c       first obtain indx i1, i2 in bond coord array
        do m=1,nbtot
          if(i_indx.eq.indx(m)) i1 = m
          if(j_indx.eq.indx(m)) i2 = m
        end do

c       Now doing energy and grad calculations
        message = "checking bounds for bond"
        if(i1.eq.0 .or. i2.eq.0) then
c         write(*,*) i1,i2
          goto 911
        end if

        xij = br(1,i2) - br(1,i1)
        yij = br(2,i2) - br(2,i1)
        zij = br(3,i2) - br(3,i1)
 
        rr = xij*xij + yij*yij + zij*zij
        rr = sqrt(rr)
          
        da = rr - rb(n)
        de = kb(n)*da
        e_bon = de*da
        ebon = ebon + e_bon
        
        df = (de+de)/rr

        do i=1,nact
          iqm = int_mb(i_iact+i-1)
          if(iqm.eq.i_indx) then
            g(1,i) = g(1,i) - df*xij     
            g(2,i) = g(2,i) - df*yij     
            g(3,i) = g(3,i) - df*zij
          end if
          if(iqm.eq.j_indx) then
            g(1,i) = g(1,i) + df*xij     
            g(2,i) = g(2,i) + df*yij     
            g(3,i) = g(3,i) + df*zij
          end if
        end do     
      end do
  
c     write(*,*) 'ebon = ', ebon , ' kcal/mol'
      return
911   call errquit("error "//trim(message),0,-1)

      end

      subroutine mm_angl_add_energy(rtdb,e)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
      integer rtdb
      double precision e

      integer i_g, h_g
      double precision scale_energy,scale_length
      double precision scale_force
      double precision eang
      character*30 message
      character*30 pname
      character*16 thetag
      double precision crd(3),q
      integer i,j

      pname = "mm_ang_add_energy"
c     write(*,*) pname
 
      call mm_angl_energy(nqm,nbtot,dbl_mb(i_br),log_mb(i_lqm), 
     &                   int_mb(i_indx),nang,int_mb(i_iang),
     &                   int_mb(i_jang),int_mb(i_kang),
     &                   dbl_mb(i_kthe),dbl_mb(i_thet),
     &                   eang)

      call util_convert_units("kcal","au",scale_energy)
 
      eang = eang*scale_energy
      e = e + eang

      return
911   call errquit("error "//trim(message),0,-1)

      end     
 
      subroutine mm_angl_add_egrad(rtdb,e,n,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
      integer rtdb
      double precision e
      integer n
      double precision g(3,n)

      integer i_g, h_g
      double precision scale_energy,scale_length
      double precision scale_force
      double precision eang
      character*30 message
      character*30 pname
      character*16 thetag
      double precision crd(3),q
      integer i,j

      pname = "mm_ang_add_egrad"
c     write(*,*) pname

      if ( .not. ma_push_get( MT_DBL, 3*n, 
     &      'angl grad scratch ',
     &      h_g, i_g) ) call errquit(
     &      'failed scratch space ',
     &      nqm, MA_ERR)

      call dfill(3*n,0.0d0,dbl_mb(i_g),1)
            
      call mm_angl_egrad(n,nbtot,dbl_mb(i_br),log_mb(i_lqm), 
     &                   int_mb(i_indx),nang,int_mb(i_iang),
     &                   int_mb(i_jang),int_mb(i_kang),
     &                   dbl_mb(i_kthe),dbl_mb(i_thet),
     &                   eang,dbl_mb(i_g))

      call util_convert_units("kcal","au",scale_energy)
      call util_convert_units("ang","au",scale_length)
 
      eang = eang*scale_energy
      e = e + eang
      scale_force = scale_energy/scale_length      

      call daxpy(3*n,scale_force,dbl_mb(i_g),1,g,1)

      message = "memory deallocation"
      if(.not.ma_pop_stack(h_g)) goto 911
      return
911   call errquit("error "//trim(message),0,-1)

      end

      subroutine mm_angl_energy(ng,nbtot,br,lqm,indx,
     &                          nang,iang,jang,kang,kthe,thet,eang)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
      integer ng
      integer nbtot
      double precision br(3,nbtot)
      logical lqm(nbtot)
      integer indx(nbtot)
      integer nang
      integer iang(nang)
      integer jang(nang)
      integer kang(nang)
      double precision kthe(nang)
      double precision thet(nang)
      double precision eang

      character*30 pname
      character*16 tag
      character*30 message
      integer i, n, m
      integer iqm
      integer i_indx, j_indx, k_indx
      integer i1, i2, i3
      double precision xij,yij,zij
      double precision xkj,ykj,zkj
      double precision da,de,df
      double precision rij, rkj, rik
      double precision e_ang
      double precision pt999
      double precision ct0, ct1, ct2, cst
      double precision ant
      double precision sth
      double precision cii, cik, ckk
      double precision dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9
    
      data pt999 /0.9990d0/

      pname = "mm_angl_energy"
c     write(*,*) pname
        
      eang = 0.0d0
      do n=1,nang
        i_indx = iang(n)
        j_indx = jang(n)
        k_indx = kang(n)

        do m=1,nbtot
          if(i_indx.eq.indx(m)) i1 = m
          if(j_indx.eq.indx(m)) i2 = m
          if(k_indx.eq.indx(m)) i3 = m
        end do

        xij = br(1,i1) - br(1,i2)
        yij = br(2,i1) - br(2,i2)
        zij = br(3,i1) - br(3,i2)
        xkj = br(1,i3) - br(1,i2)
        ykj = br(2,i3) - br(2,i2)
        zkj = br(3,i3) - br(3,i2)
        
        rij = xij**2 + yij**2 + zij**2      
        rkj = xkj**2 + ykj**2 + zkj**2      
        rik = sqrt(rij*rkj)
        ct0 = (xij*xkj+yij*ykj+zij*zkj)/rik
        ct1 = max(-pt999,ct0)
        ct2 = min(pt999,ct1)
        cst = ct2
        ant = acos(ct2)

        da = ant - thet(n)
        de = kthe(n)*da
        e_ang = de*da
        eang = eang + e_ang
      end do
 
c     write(*,*) 'eang = ', eang, ' kcal/mol'
         
      end

      subroutine mm_angl_egrad(ng,nbtot,br,lqm,indx,
     &                          nang,iang,jang,kang,kthe,thet,eang,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_geom_data.fh"
#include "mm_coords_data.fh"
      integer ng
      integer nbtot
      double precision br(3,nbtot)
      logical lqm(nbtot)
      integer indx(nbtot)
      integer nang
      integer iang(nang)
      integer jang(nang)
      integer kang(nang)
      double precision kthe(nang)
      double precision thet(nang)
      double precision g(3,ng)
      double precision eang

      character*30 pname
      character*16 tag
      character*30 message
      integer i, n, m
      integer iqm
      integer i_indx, j_indx, k_indx
      integer i1, i2, i3
      double precision xij,yij,zij
      double precision xkj,ykj,zkj
      double precision da,de,df
      double precision rij, rkj, rik
      double precision e_ang
      double precision pt999
      double precision ct0, ct1, ct2, cst
      double precision ant
      double precision sth
      double precision cii, cik, ckk
      double precision dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8,dt9
      double precision sphi     
    
      data pt999 /0.9990d0/

      pname = "mm_angl_egrad"
c     write(*,*) pname
      
      eang = 0.0d0
      
      do n=1,nang
        i_indx = iang(n)
        j_indx = jang(n)
        k_indx = kang(n)

        do m=1,nbtot
          if(i_indx.eq.indx(m)) i1 = m
          if(j_indx.eq.indx(m)) i2 = m
          if(k_indx.eq.indx(m)) i3 = m
        end do

        xij = br(1,i1) - br(1,i2)
        yij = br(2,i1) - br(2,i2)
        zij = br(3,i1) - br(3,i2)
        xkj = br(1,i3) - br(1,i2)
        ykj = br(2,i3) - br(2,i2)
        zkj = br(3,i3) - br(3,i2)
        
        rij = xij**2 + yij**2 + zij**2      
        rkj = xkj**2 + ykj**2 + zkj**2      
        rik = sqrt(rij*rkj)
        ct0 = (xij*xkj+yij*ykj+zij*zkj)/rik
        ct1 = max(-pt999,ct0)
        ct2 = min(pt999,ct1)
        cst = ct2
        ant = acos(ct2)
     
        da = ant - thet(n)
        de = kthe(n)*da
        e_ang = de*da
        eang = eang + e_ang

C       doing gradients and add to QM atoms
C       
        sphi = sin(ant)
        df = -(de+de)/sphi
        sth = df*cst
        cik = df/rik
        cii = sth/rij
        ckk = sth/rkj
        dt1 = cik*xkj-cii*xij
        dt2 = cik*ykj-cii*yij
        dt3 = cik*zkj-cii*zij
        dt7 = cik*xij-ckk*xkj
        dt8 = cik*yij-ckk*ykj
        dt9 = cik*zij-ckk*zkj
        dt4 = -dt1-dt7
        dt5 = -dt2-dt8
        dt6 = -dt3-dt9

        do i=1,nact
          iqm = int_mb(i_iact+i-1)
          if(i_indx.eq.iqm) then        
            g(1,i) = g(1,i) + dt1
            g(2,i) = g(2,i) + dt2
            g(3,i) = g(3,i) + dt3
          end if
          if(j_indx.eq.iqm) then        
            g(1,i) = g(1,i) + dt4
            g(2,i) = g(2,i) + dt5
            g(3,i) = g(3,i) + dt6
          end if
          if(k_indx.eq.iqm) then        
            g(1,i) = g(1,i) + dt7
            g(2,i) = g(2,i) + dt8
            g(3,i) = g(3,i) + dt9
          end if
        end do
      end do
 
c     write(*,*) 'eang = ', eang, ' kcal/mol'
         
      end

      subroutine mm_dihe_add_energy(rtdb,e)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
      integer rtdb
      double precision e

      integer i_g, h_g
      double precision scale_energy,scale_length
      double precision scale_force
      double precision edih
      character*30 message
      character*30 pname
      character*16 thetag
      double precision crd(3),q
      integer i,j

      pname = "mm_dihe_add_energy"
c     write(*,*) pname
 
      call mm_dihe_energy(nqm,nbtot,dbl_mb(i_br),log_mb(i_lqm), 
     &                   int_mb(i_indx),ndih,int_mb(i_idih),
     &                   int_mb(i_jdih),int_mb(i_kdih),
     &                   int_mb(i_ldih),dbl_mb(i_kphi),
     &                   dbl_mb(i_peri),dbl_mb(i_phas),
     &                   edih)

      call util_convert_units("kcal","au",scale_energy)
 
      edih = edih*scale_energy
      e = e + edih

      return
911   call errquit("error "//trim(message),0,-1)

      end

      subroutine mm_dihe_add_egrad(rtdb,e,n,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_bond_data.fh"
#include "mm_coords_data.fh"
#include "mm_bond_coords_data.fh"
      integer rtdb
      double precision e
      integer n
      double precision g(3,n)

      integer i_g, h_g
      double precision scale_energy,scale_length
      double precision scale_force
      double precision edih
      character*30 message
      character*30 pname
      character*16 thetag
      double precision crd(3),q
      integer i,j

      pname = "mm_dihe_add_egrad"
c     write(*,*) pname
 
      if ( .not. ma_push_get( MT_DBL, 3*n, 
     &      'dih grad scratch ',
     &      h_g, i_g) ) call errquit(
     &      'failed scratch space ',
     &      nqm, MA_ERR)

      call dfill(3*n,0.0d0,dbl_mb(i_g),1)
            
      call mm_dihe_egrad(n,nbtot,dbl_mb(i_br),log_mb(i_lqm), 
     &                   int_mb(i_indx),ndih,int_mb(i_idih),
     &                   int_mb(i_jdih),int_mb(i_kdih),
     &                   int_mb(i_ldih),dbl_mb(i_kphi),
     &                   dbl_mb(i_peri),dbl_mb(i_phas),
     &                   edih,dbl_mb(i_g))

      call util_convert_units("kcal","au",scale_energy)
      call util_convert_units("ang","au",scale_length)
 
      edih = edih*scale_energy
      e = e + edih
      scale_force = scale_energy/scale_length      

      call daxpy(3*n,scale_force,dbl_mb(i_g),1,g,1)

      message = "memory deallocation"
      if(.not.ma_pop_stack(h_g)) goto 911
      return
911   call errquit("error "//trim(message),0,-1)

      end

      subroutine mm_dihe_energy(ng,nbtot,br,lqm,indx,
     &                          ndih,idih,jdih,kdih,ldih,
     &                          kphi,peri,phas,edih)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
      integer nbtot
      double precision br(3,nbtot)
      logical lqm(nbtot)
      integer indx(nbtot)
      integer ndih
      integer idih(ndih)
      integer jdih(ndih)
      integer kdih(ndih)
      integer ldih(ndih)
      double precision kphi(ndih)
      double precision peri(ndih)
      double precision phas(ndih)
      integer ng
      double precision edih

      character*30 pname
      character*16 tag
      character*30 message
      integer i, n, m
      integer iqm
      integer i_indx, j_indx, k_indx, l_indx
      integer i1, i2, i3, i4
      double precision xij,yij,zij
      double precision xkj,ykj,zkj
      double precision xkl,ykl,zkl
      double precision dx,dy,dz
      double precision gx,gy,gz
      double precision da,de,df
      double precision rij, rkj, rik
      double precision e_dih
      double precision z10,z20,z11,z12,z22,z1,z2 
      double precision ct,ct0,ct1
      double precision s
      double precision ftem
      double precision ap0,ap1
      double precision cphi,sphi
      double precision cosnp,sinnp
      double precision pi
      double precision zero,one,tm06,tenm3
      double precision df0,df1,dums,dflim
      double precision inc
      double precision gmul(10)
      double precision dc1,dc2,dc3,dc4,dc5,dc6
      double precision dr1,dr2,dr3,dr4,dr5,dr6
      double precision drx,dry,drz
      double precision fxi,fyi,fzi
      double precision fxj,fyj,fzj
      double precision fxk,fyk,fzk
      double precision fxl,fyl,fzl
      
       
      data gmul/0.0d+00,2.0d+00,0.0d+00,4.0d+00,0.0d+00,6.0d+00,
     &          0.0d+00,8.0d+00,0.0d+00,10.0d+00/

      data zero,one,tm06,tenm3 /0.0d0,1.0d0,1.0d-06,1.0d-03/

      parameter (pi = 3.1415926535897932384626433832795d0)

      pname = "mm_dihe_energy"
c     write(*,*) pname

      edih = 0.0d0
      do n=1,ndih
        i_indx = idih(n)
        j_indx = jdih(n)
        k_indx = kdih(n)
        l_indx = ldih(n)

        do m=1,nbtot
          if(i_indx.eq.indx(m)) i1 = m
          if(j_indx.eq.indx(m)) i2 = m
          if(k_indx.eq.indx(m)) i3 = m
          if(l_indx.eq.indx(m)) i4 = m
        end do
        
        xij = br(1,i1)-br(1,i2)
        yij = br(2,i1)-br(2,i2)
        zij = br(3,i1)-br(3,i2)
        xkj = br(1,i3)-br(1,i2)
        ykj = br(2,i3)-br(2,i2)
        zkj = br(3,i3)-br(3,i2)
        xkl = br(1,i3)-br(1,i4)
        ykl = br(2,i3)-br(2,i4)
        zkl = br(3,i3)-br(3,i4)

        dx = yij*zkj - zij*ykj
        dy = zij*xkj - xij*zkj
        dz = xij*ykj - yij*xkj
        gx = zkj*ykl - ykj*zkl
        gy = xkj*zkl - zkj*xkl
        gz = ykj*xkl - xkj*ykl

        fxi = sqrt(dx**2+dy**2+dz**2+1.0d-18)
        fyi = sqrt(gx**2+gy**2+gz**2+1.0d-18)
        ct  = dx*gx+dy*gy+dz*gz
        
        z10 = one/fxi
        z20 = one/fyi
        if (tenm3 > fxi) z10 = zero
        if (tenm3 > fyi) z20 = zero
        z12 = z10*z20
        z1 = z10
        z2 = z20
        ftem = zero
        if (z12 /= zero) ftem = one
        fzi = ftem
        ct0 = min(one,ct*z12)
        ct1 = max(-one,ct0)
        s = xkj*(dz*gy-dy*gz)+ 
     &      ykj*(dx*gz-dz*gx)+ 
     &      zkj*(dy*gx-dx*gy)
        ap0 = acos(ct1)
        ap1 = pi-sign(ap0,s)
        ct = ap1
        cphi = cos(ap1)
        sphi = sin(ap1)
        
        ct0 = peri(n)*ct
        cosnp = cos(ct0)
        sinnp = sin(ct0)

        e_dih = (kphi(n)+cosnp*kphi(n)*cos(phas(n)) +
     &                   sinnp*kphi(n)*sin(phas(n)))*fzi

        edih = edih + e_dih
      end do

c     write(*,*) 'edih = ', edih, ' kcal/mol'

      end

      subroutine mm_dihe_egrad(ng,nbtot,br,lqm,indx,
     &                          ndih,idih,jdih,kdih,ldih,
     &                          kphi,peri,phas,edih,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_geom_data.fh"
#include "mm_coords_data.fh"
      integer nbtot
      double precision br(3,nbtot)
      logical lqm(nbtot)
      integer indx(nbtot)
      integer ndih
      integer idih(ndih)
      integer jdih(ndih)
      integer kdih(ndih)
      integer ldih(ndih)
      double precision kphi(ndih)
      double precision peri(ndih)
      double precision phas(ndih)
      integer ng
      double precision edih
      double precision g(3,ng)

      character*30 pname
      character*16 tag
      character*30 message
      integer i, n, m
      integer iqm
      integer i_indx, j_indx, k_indx, l_indx
      integer i1, i2, i3, i4
      integer inc
      double precision xij,yij,zij
      double precision xkj,ykj,zkj
      double precision xkl,ykl,zkl
      double precision dx,dy,dz
      double precision gx,gy,gz
      double precision da,de,df
      double precision rij, rkj, rik
      double precision e_dih
      double precision z10,z20,z11,z12,z22,z1,z2 
      double precision ct,ct0,ct1
      double precision s
      double precision ftem
      double precision ap0,ap1
      double precision cphi,sphi
      double precision cosnp,sinnp
      double precision pi
      double precision zero,one,tm06,tenm3
      double precision df0,df1,dums,dflim
      double precision gmul(10)
      double precision dc1,dc2,dc3,dc4,dc5,dc6
      double precision dr1,dr2,dr3,dr4,dr5,dr6
      double precision drx,dry,drz
      double precision fxi,fyi,fzi
      double precision fxj,fyj,fzj
      double precision fxk,fyk,fzk
      double precision fxl,fyl,fzl
      
       
      data gmul/0.0d+00,2.0d+00,0.0d+00,4.0d+00,0.0d+00,6.0d+00,
     &          0.0d+00,8.0d+00,0.0d+00,10.0d+00/

      data zero,one,tm06,tenm3 /0.0d0,1.0d0,1.0d-06,1.0d-03/

      parameter (pi = 3.1415926535897932384626433832795d0)

      pname = "mm_dihe_egrad"
c     write(*,*) pname

      edih = 0.0d0
      do n=1,ndih
        i_indx = idih(n)
        j_indx = jdih(n)
        k_indx = kdih(n)
        l_indx = ldih(n)

        do m=1,nbtot
          if(i_indx.eq.indx(m)) i1 = m
          if(j_indx.eq.indx(m)) i2 = m
          if(k_indx.eq.indx(m)) i3 = m
          if(l_indx.eq.indx(m)) i4 = m
        end do
        
        xij = br(1,i1)-br(1,i2)
        yij = br(2,i1)-br(2,i2)
        zij = br(3,i1)-br(3,i2)
        xkj = br(1,i3)-br(1,i2)
        ykj = br(2,i3)-br(2,i2)
        zkj = br(3,i3)-br(3,i2)
        xkl = br(1,i3)-br(1,i4)
        ykl = br(2,i3)-br(2,i4)
        zkl = br(3,i3)-br(3,i4)

        dx = yij*zkj - zij*ykj
        dy = zij*xkj - xij*zkj
        dz = xij*ykj - yij*xkj
        gx = zkj*ykl - ykj*zkl
        gy = xkj*zkl - zkj*xkl
        gz = ykj*xkl - xkj*ykl

        fxi = sqrt(dx**2+dy**2+dz**2+1.0d-18)
        fyi = sqrt(gx**2+gy**2+gz**2+1.0d-18)
        ct  = dx*gx+dy*gy+dz*gz
        
        z10 = one/fxi
        z20 = one/fyi
        if (tenm3 > fxi) z10 = zero
        if (tenm3 > fyi) z20 = zero
        z12 = z10*z20
        z1 = z10
        z2 = z20
        ftem = zero
        if (z12 /= zero) ftem = one
        fzi = ftem
        ct0 = min(one,ct*z12)
        ct1 = max(-one,ct0)
        s = xkj*(dz*gy-dy*gz)+ 
     &      ykj*(dx*gz-dz*gx)+ 
     &      zkj*(dy*gx-dx*gy)
        ap0 = acos(ct1)
        ap1 = pi-sign(ap0,s)
        ct = ap1
        cphi = cos(ap1)
        sphi = sin(ap1)
        
        ct0 = peri(n)*ct
        cosnp = cos(ct0)
        sinnp = sin(ct0)

        e_dih = (kphi(n)+cosnp*kphi(n)*cos(phas(n)) +
     &                   sinnp*kphi(n)*sin(phas(n)))*fzi

        edih = edih + e_dih

C       doing gradients and add to qm atoms    
C
        inc = int(peri(n))
        df0 = peri(n)*(cos(phas(n))*kphi(n)*sinnp -
     &                 sin(phas(n))*kphi(n)*cosnp)
        dums = sphi+sign(1.0d-18,sphi)         
        dflim = cos(phas(n))*kphi(n)*(peri(n)-gmul(inc)+gmul(inc)*cphi)
        df1 = df0/dums
        if(tm06 > abs(dums)) df1 = dflim
   
        df = df1*fzi

        z11 = z1*z1
        z12 = z1*z2
        z22 = z2*z2
        dc1 = -gx*z12-cphi*dx*z11
        dc2 = -gy*z12-cphi*dy*z11
        dc3 = -gz*z12-cphi*dz*z11
        dc4 =  dx*z12+cphi*gx*z22
        dc5 =  dy*z12+cphi*gy*z22
        dc6 =  dz*z12+cphi*gz*z22

        dr1 = df*(dc3*ykj - dc2*zkj)
        dr2 = df*(dc1*zkj - dc3*xkj)
        dr3 = df*(dc2*xkj - dc1*ykj)
        dr4 = df*(dc6*ykj - dc5*zkj)
        dr5 = df*(dc4*zkj - dc6*xkj)
        dr6 = df*(dc5*xkj - dc4*ykj)
        drx = df*(-dc2*zij + dc3*yij + 
     &        dc5*zkl - dc6*ykl)
        dry = df*( dc1*zij - dc3*xij - 
     &        dc4*zkl + dc6*xkl)
        drz = df*(-dc1*yij + dc2*xij + 
     &        dc4*ykl - dc5*xkl)
        fxi = - dr1
        fyi = - dr2
        fzi = - dr3
        fxj = - drx + dr1
        fyj = - dry + dr2
        fzj = - drz + dr3
        fxk = + drx + dr4
        fyk = + dry + dr5
        fzk = + drz + dr6
        fxl = - dr4
        fyl = - dr5
        fzl = - dr6
      
        do i=1,nact
          iqm = int_mb(i_iact+i-1)
          if(i_indx.eq.iqm) then
            g(1,i) = g(1,i) - fxi
            g(2,i) = g(2,i) - fyi
            g(3,i) = g(3,i) - fzi
          end if
          if(j_indx.eq.iqm) then
            g(1,i) = g(1,i) - fxj
            g(2,i) = g(2,i) - fyj
            g(3,i) = g(3,i) - fzj
          end if 
          if(k_indx.eq.iqm) then
            g(1,i) = g(1,i) - fxk
            g(2,i) = g(2,i) - fyk
            g(3,i) = g(3,i) - fzk
          end if
          if(l_indx.eq.iqm) then
            g(1,i) = g(1,i) - fxl
            g(2,i) = g(2,i) - fyl
            g(3,i) = g(3,i) - fzl
          end if
        end do
         
      end do

c     write(*,*) 'edih = ', edih, ' kcal/mol'
      
      end 

      subroutine mm_bond_call_test()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_bond_coords_data.fh"
#include "mm_bond_data.fh"
#include "mm_coords_data.fh"

      double precision ebon, eang, edih
      
      call mm_bond_energy(nqm,nbtot,dbl_mb(i_br),log_mb(i_lqm),
     &                   int_mb(i_indx),nbon,int_mb(i_ibon),
     &                   int_mb(i_jbon),dbl_mb(i_kb),dbl_mb(i_rb),ebon)

      call mm_angl_energy(nqm,nbtot,dbl_mb(i_br),log_mb(i_lqm),
     &                   int_mb(i_indx),nang,int_mb(i_iang),
     &                   int_mb(i_jang),int_mb(i_kang),
     &                   dbl_mb(i_kthe),dbl_mb(i_thet),eang)
   
      call mm_dihe_energy(nqm,nbtot,dbl_mb(i_br),log_mb(i_lqm),
     &                   int_mb(i_indx),ndih,int_mb(i_idih),
     &                   int_mb(i_jdih),int_mb(i_kdih),int_mb(i_ldih),
     &                   dbl_mb(i_kphi),dbl_mb(i_peri),
     &                   dbl_mb(i_phas),edih)
 
c     write(*,*) 'ebon =', ebon, ' kcal/mol'
c     write(*,*) 'eang =', eang, ' kcal/mol'
c     write(*,*) 'edih =', edih, ' kcal/mol'

      end

      subroutine mm_write_bond(lunit)
      implicit none
#include "mafdecls.fh"
#include "mm_bond_data.fh"

      integer lunit
      character*30 pname 

      integer i
      integer indx_i, indx_j
      double precision kb
      double precision rb

      pname = "mm_write_bond"
c     write(*,*) pname

      write(lunit,'(a)') 'bond'
      write(lunit,'(a)') '#      i       j       k_ij             r0'
      write(lunit,'(a)') '#'

      do i=1,nbon
        indx_i = int_mb(i_ibon+i-1)
        indx_j = int_mb(i_jbon+i-1)
        kb = dbl_mb(i_kb+i-1)
        rb = dbl_mb(i_rb+i-1)
        write(lunit,'(2(x,i7),2ES16.8)') indx_i, indx_j, kb, rb
      end do

      write(lunit,'(a,/)') 'end'

      end

      subroutine mm_write_angl(lunit)
      implicit none
#include "mafdecls.fh"
#include "mm_bond_data.fh"

      integer lunit
      character*30 pname 

      integer i
      integer indx_i, indx_j, indx_k
      double precision kthe
      double precision thet

      pname = "mm_write_angl"
c     write(*,*) pname

      write(lunit,'(a)') 'angle'
      write(lunit,'(a)') '#      i       j       k      k_ijk'//
     >                   '           theta0'
      write(lunit,'(a)') '#'

      do i=1,nang
        indx_i = int_mb(i_iang+i-1)
        indx_j = int_mb(i_jang+i-1)
        indx_k = int_mb(i_kang+i-1)
        kthe = dbl_mb(i_kthe+i-1)
        thet = dbl_mb(i_thet+i-1)
        write(lunit,'(3(x,i7),2ES16.8)') indx_i, indx_j, indx_k, kthe,
     >                                   thet
      end do

      write(lunit,'(a,/)') 'end'

      end

      subroutine mm_write_dihe(lunit)
      implicit none
#include "mafdecls.fh"
#include "mm_bond_data.fh"

      integer lunit
      character*30 pname 

      integer i
      integer indx_i, indx_j, indx_k, indx_l
      double precision kphi
      double precision peri
      double precision phase

      pname = "mm_write_dihe"
c     write(*,*) pname

      write(lunit,'(a)') 'dihedral'
      write(lunit,'(a)') '#      i       j       k       l'//
     >            '      k_ijkl       periodicity        phase'
      write(lunit,'(a)') '#'
      do i=1,ndih
        indx_i = int_mb(i_idih+i-1)
        indx_j = int_mb(i_jdih+i-1)
        indx_k = int_mb(i_kdih+i-1)
        indx_l = int_mb(i_ldih+i-1)
        kphi = dbl_mb(i_kphi+i-1)
        peri = dbl_mb(i_peri+i-1)
        phase = dbl_mb(i_phas+i-1)
        write(lunit,'(4(x,i7),3ES16.8)') indx_i, indx_j, indx_k, indx_l,
     >                               kphi, peri, phase
      end do

      write(lunit,'(a,/)') 'end'

      end

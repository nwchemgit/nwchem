c
c     empirical dispersion: hessian
c
      subroutine nwxc_vdw_hess(delta,s6,s8,sr6,sr8,n,x,z)
c
      implicit none
c
c This function makes vdw empirical correction to the hessian
c must be called before thermochemical data and vibrational
c analysis are generated.
c
#include "inp.fh"
#include "util.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
c
      double precision s6,s8,sr6,sr8
      integer n                
      double precision x(3,n)
      integer z(n)
      double precision l_force(3,n)
      double precision r_force(3,n)
      double precision hessvdw(3,n,3,n)
      double precision dval0
      double precision dval1
      double precision delta
      integer i, j, A, B
      integer n3xyz 
      integer nat2
      integer nhesst
      integer lenhess 
      integer l_exybs,k_exybs 
      integer l_exybt,k_exybt
c
      character*(nw_max_path_len) filehess
c
c  dispersion contribution to hessian 
c
c in principle this task is very fast,
c so only master node works, and read-write to disk.
c 
      if (ga_nodeid().eq.0) then
c
      call util_file_name('hess', .false., .false.,filehess)

      lenhess=inp_strlen(filehess)
      n3xyz=3*n
      nhesst=n3xyz*(n3xyz+1)/2
      nat2=n3xyz*n3xyz
c
      if (.not.ma_push_get(mt_dbl,nat2,'xcvdwhess exybs ',
     *  l_exybs,k_exybs))
     &   call errquit('xcvdwhess: cannot allocate exybs',0, MA_ERR)
      call dfill(nat2,0.0d00,dbl_mb(k_exybs),1)
c
      if (.not.ma_push_get(mt_dbl,nhesst,'xcvdwhess exybt ',
     *  l_exybt,k_exybt))
     &   call errquit('xcvdwhess: cannot allocate exybt',0, MA_ERR)
      call dfill(nhesst,0.0d00,dbl_mb(k_exybt),1)
c
      write(LuOut,* ) 'Read old Hessian file : ', filehess 
c
      call dfill(nhesst,0.d0,dbl_mb(k_exybt),1)
      !write(6,* ) 'leee '
      open(unit=69,file=filehess,form='formatted',status='old',
     &    err=99900,access='sequential')
      do j = 0,(nhesst-1)
        read(69,*,err=99901,end=99902) dval0
        dbl_mb(k_exybt+j) = dval0
      !write(6,* ) 'dval ', dval0
      enddo
      close(unit=69,status='keep')
c
c     complete the square matrix from triangle values
c
      call trin2squa(dbl_mb(k_exybs),dbl_mb(k_exybt),n3xyz)
c
      write(LuOut,* ) 'vdW contribution to Hessian '
c
      call output(dbl_mb(k_exybs), 1, n3xyz, 1, n3xyz, n3xyz, n3xyz, 1)
c
       do A = 1, n
         do i = 1, 3
           do B=1, n
            do j=1, 3
               r_force(j,B)=0.0d0
               l_force(j,B)=0.0d0
            enddo
           enddo
c right displacement
            x(i,A) = x(i,A) + delta
            call nwxc_vdw_der(s6,s8,sr6,sr8,n,x,z,r_force)
c left displacement
            x(i,A) = x(i,A) - 2.0d0*delta
            call nwxc_vdw_der(s6,s8,sr6,sr8,n,x,z,l_force)
c back to original position
            x(i,A) = x(i,A) + delta
c
            do B=1,n
              do  j=1, 3
                dval1 = (r_force(j,B)-l_force(j,B)) / (2.0d00*delta)
                hessvdw(i,A,j,B)=dval1
              enddo
            enddo
         enddo
      enddo
c
      call output(hessvdw,1,n3xyz,1,n3xyz,
     &            n3xyz,n3xyz,1)
c
      call daxpy(nat2,1d0,hessvdw,1,dbl_mb(k_exybs),1)
c
c:write the final Hessian
c
      call output(dbl_mb(k_exybs),1,n3xyz,1,n3xyz,
     &            n3xyz,n3xyz,1)
      call  stpr_wrt_fd_from_sq(dbl_mb(k_exybs),n3xyz,filehess)
      write(LuOut,* ) 'New Hessian file vdw corrected has been 
     .                      written:', filehess 
      if (.not.ma_chop_stack(l_exybs)) 
     C   call errquit('xcvdwhess: cannot pop stack exybs',4, MA_ERR)

      endif
c
      call ga_sync
c
      return
99900 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('xc_vdw_hess 99900', 911, DISK_ERR)
99901 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('xc_vdw_hess 99901', 911, DISK_ERR)
99902 continue
      write(luout,*)'hess_file => ',filehess
      call errquit('xc_vdw_hess 99902', 911, DISK_ERR)
      end
c
      subroutine nwxc_vdw_to_hessian(rtdb)
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "geom.fh"
c
      integer rtdb
      integer  geom
      double precision dum
      double precision dum1(1)
      character*255 name
c
      if (.not. rtdb_cget(rtdb,'geometry', 1, name))
     $     name = 'geometry'
c
      if (.not. geom_create(geom, name))
     $ call errquit('nwxc_vdw_to_hessian: geom_create failed !', 
     $                                                0,GEOM_ERR)
c
      if (.not. geom_rtdb_load(rtdb, geom, name))
     $call errquit('nwxc_vdw_to_hessian: no geometry load form rtdb',0,
     $        GEOM_ERR)
c
      call nwxc_vdw(rtdb, geom, dum, dum1, 'hessian')
c
      if(.not. geom_destroy(geom))
     $ call errquit('xc_vdw_to_hessian: geom_create failed !', 
     $                                                0,GEOM_ERR)
c
      return
      end
c $Id$

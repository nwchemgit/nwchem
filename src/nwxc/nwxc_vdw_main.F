c
c     empirical dispersion: main wrapper
c
C> \ingroup nwxc_api
C> @{
C>
      subroutine nwxc_vdw(rtdb,geom,exc,force,what)
c
      implicit none
      character *(*) what
      integer geom,rtdb
      double precision exc,force(*),s6,s8,sr6,sr8,a1,a2
c
#include "geom.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "nwxc_vdw.fh"
c
      integer n
      integer itags,ltags,i_xyz,l_xyz,icharge,lcharge,
     I     l_fvdw,i_fvdw
      external nwxc_vdw_e
      double precision nwxc_vdw_e,evdw,scalea
      integer iz,lz,i
      logical nwxc_vdw_init
      external nwxc_vdw_init
      logical oprint,oprinth
      logical       stat
c
c      logical geom_cart_get2 
c      external geom_cart_get2 
c
      double precision delta,delta_default
c
      oprint = util_print('vdw', print_medium)
      oprinth = util_print('vdw high', print_high)
c
c     Allocate memory blocks
c
      if (.not. geom_ncent(geom, n)) 
     &   call errquit('nwxc_vdw: geom_ncent failed',geom, GEOM_ERR)
      if (.not.MA_push_get(MT_Dbl,n*3,'xyz',l_xyz,i_xyz))
     &   call errquit('nwxc_vdw: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_int,n,'atns',lz,iz))
     &   call errquit('nwxc_vdw: cannot allocate atns',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,n,'charge',lcharge,icharge))
     &   call errquit('nwxc_vdw: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,n*16,'tags',ltags,itags))
     &   call errquit('nwxc_vdw: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get2(geom, n, Byte_MB(itags),
     &        Dbl_MB(i_xyz), Dbl_MB(icharge), int_mb(iz)))
     &   call errquit('nwxc_vdw: geom_cart_get failed',74, GEOM_ERR)
      if (.not.ma_pop_stack(ltags))
     &   call errquit('nwxc_vdw: cannot pop stack',2, MA_ERR)
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('nwxc_vdw: cannot pop stack',3, MA_ERR)
c
c     Which Grimme dispersion version
c
      if (.not.rtdb_get(rtdb,'dft:ivdw',MT_INT,1,ivdw))
     &      ivdw = 2
c
c     conversion factor angs 2 au
c
       if(.not.geom_get_ang2au(geom, scalea)) call
     S     errquit('nwxc_vdw: gang2au failed',0,0)
c
c     Initialize are variables
c
      if(.not.nwxc_vdw_init(scalea))
     &   call errquit('nwxc_vdw: vwdinit failed',0, 0)
c
c     Read in some user defined parameters
c
      call nwxc_get_disp(ivdw,s6,s8,sr6,sr8,alpha,a1,a2)
cDEBUG
      write(*,*)'nwxc_vdw: after nwxc_get_disp: s6,s8,sr6,sr8,a1,a2=',
     &          s6,s8,sr6,sr8,a1,a2
cDEBUG
c     if (.not.rtdb_get(rtdb,'dft:vdwalpha',MT_DBL,1,alpha))
c    &      alpha = 20.0d0
c     if (ivdw.eq.3) alpha = 14.0d0
c
c     Get proper scaling factors depending on Grimme dispersion version
c
c     stat = rtdb_get(rtdb, 'dft:vdw', mt_dbl, 1, s6)
cc    if (.not.rtdb_get(rtdb, 'dft:vdw', mt_dbl, 1, s6))
cc     &    call get_scaling_fac(s6,s8,sr6,sr8)
cc   &      call nwxc_get_disp(ivdw,s6,s8,alpha)
c
      if(what.eq.'energy') then
c
c     Compute energy contribution
c
        if(oprinth.and.ga_nodeid().eq.0) then
           write(luout,*) ' s6     =',s6
           write(luout,*) ' s8     =',s8
           write(luout,*) ' sr6    =',sr6
           write(luout,*) ' sr8    =',sr8
           write(luout,*) ' alpha  =',alpha
           write(luout,*) ' a1     =',a1
           write(luout,*) ' a2     =',a2
           write(luout,*) ' ivdw   =',ivdw
           write(luout,*) ' k1     = ',k1
           write(luout,*) ' k2     = ',k2
           write(luout,*) ' k3     = ',k3
         endif
c
         evdw=nwxc_vdw_e(s6,s8,sr6,sr8,n,dbl_mb(i_xyz),int_mb(iz))
c
         if(oprint.and.ga_nodeid().eq.0) then
           if (ivdw.eq.1.or.ivdw.eq.2) then
              write(luout,222) s6, evdw
 222  format(//
     &     '          s6 scale factor =', f22.12/
     &     '              vdW contrib =', f22.12//)
           endif
           if (ivdw.eq.3) then
              write(luout,223) s6, s8, sr6, sr8, evdw
 223  format(//
     &     '          s6 scale factor  =', f22.12/
     &     '          s8 scale factor  =', f22.12/
     &     '          sr6 scale factor =', f22.12/
     &     '          sr8 scale factor =', f22.12/
     &     '              vdW contrib  =', f22.12//)
           endif
           if (ivdw.eq.4) then
              write(luout,224) s6, s8, a1, a2, evdw
 224  format(//
     &     '          s6 scale factor =', f22.12/
     &     '          s8 scale factor =', f22.12/
     &     '          a1 parameter    =', f22.12/
     &     '          a2 parameter    =', f22.12/
     &     '             vdW contrib  =', f22.12//)
           endif
         endif
c
c        Add contribution to Exc
c
         Exc=Exc+evdw
c
      elseif(what.eq.'forces') then
c
c     Gradient calculation
c
      if (.not.MA_push_get(MT_Dbl,n*3,'xyz',l_fvdw,i_fvdw))
     &   call errquit('nwxc_vdw: cannot allocate forcev',0, MA_ERR)
c
         call nwxc_vdw_der(s6,s8,sr6,sr8,n,dbl_mb(i_xyz),int_mb(iz),
     D     dbl_mb(i_fvdw))
c
         if(oprinth.and.ga_nodeid().eq.0) then
             write(luout,*) ' vdW contrib for S6=',s6
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               dbl_mb(i_fvdw+3*(i-1)),
     Y               dbl_mb(i_fvdw+3*(i-1)+1),
     Z        dbl_mb(i_fvdw+3*(i-1)+2)
             enddo
             write(luout,*) ' before vdw contr @@@@@'
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               force(1+3*(i-1)),
     Y               force(1+3*(i-1)+1),
     Z               force(1+3*(i-1)+2)
             enddo

          endif
c
c         Add to force matrix
c
          call daxpy(3*n,1d0,dbl_mb(i_fvdw),1,force,1)
c
          if(oprinth.and.ga_nodeid().eq.0) then
             write(luout,*) ' after vdw contr @@@@@'
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               force(1+3*(i-1)),
     Y               force(1+3*(i-1)+1),
     Z               force(1+3*(i-1)+2)
             enddo
          endif
c
      else if(what.eq.'hessian') then
c
c     Hessian calculation, numerical from analytical gradients
c
c     Get delta as used in a numerical hessian DFT calculation
c
      delta_default  =  0.01d0
      if (.not.rtdb_get(rtdb,'stpr_gen:delta',MT_DBL,1,delta))
     &      delta = delta_default 
c
      call nwxc_vdw_hess(delta,s6,s8,sr6,sr8,n,dbl_mb(i_xyz),int_mb(iz))
c
        if (oprint.and.ga_nodeid().eq.0)  then
            write(luout,*) ' s6 = ',s6
            write(luout,*) ' vdw to hessian contribution is done'
        endif
      endif ! what
c
c     Clean up
c
      if (.not.ma_chop_stack(l_xyz)) 
     C   call errquit('nwxc_vdw: cannot pop stack',4, MA_ERR)
c
      return
      end
C>
C> @}
c $Id$

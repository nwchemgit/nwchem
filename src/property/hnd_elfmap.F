      subroutine hnd_elfmap(rtdb,basis,geom)
c
c $Id: hnd_elfmap.F,v 1.1 2004-07-21 17:28:27 bert Exp $
c
c     This routine calculates the electric field
c     for a given density at the atomic positions.
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "geom.fh"
c
      integer rtdb      ! [Input] rtdb
      integer basis     ! [Input] Basis set
      integer geom      ! [Input] Geometry
c
      character*2  symbol
      character*16 element, at_tag
      integer iat, atn, nat, i, j, ij
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt, l_efld, k_efld
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      character*3 scftyp
      double precision xp, yp, zp, xn, yn, zn, zan
      double precision elf(3), elfld
      double precision rr, rr3
c
c     Initialize integrals
c
      call int_init(rtdb,1, basis)
      call schwarz_init(geom, basis)
c
c     Get density matrix
c
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
c
c     ----- calculate electric field -----
c
      if (ga_nodeid().eq.0) write(luout,9999)
      if (ga_nodeid().eq.0) write(luout,9994)
c
       call ecce_print_module_entry('EField')
c
c     ----- define points for calculation -----
c           1. nuclei
c
      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_elfmap: geom_ncent',911,GEOM_ERR)
c
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_elfmap: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,6*nat,'efld pnt',l_efld,k_efld))
     &    call errquit('hnd_elfmap: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_elfmap: ma failed',911,MA_ERR)
c
      do 30 iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_elfmap: geom_cent_get',911,GEOM_ERR)
   30 continue
c
c     ----- calculate electronic contribution at all points -----
c
      call hnd_elfcon(basis,geom,g_dens(ndens),dbl_mb(k_xyzpt),nat,
     &                dbl_mb(k_efld),1)
c
c     ----- collect and output results of all points -----
c
      if (ga_nodeid().gt.0) goto 300
c
      write(luout,9997)
      do 230  iat=1,nat
         xp = dbl_mb(k_xyzpt  +3*(iat-1))
         yp = dbl_mb(k_xyzpt+1+3*(iat-1))
         zp = dbl_mb(k_xyzpt+2+3*(iat-1))
c
c     ----- add nuclear contribution -----
c
         elf(1) = dbl_mb(k_efld  +3*(iat-1))
         elf(2) = dbl_mb(k_efld+1+3*(iat-1))
         elf(3) = dbl_mb(k_efld+2+3*(iat-1))
         do 210 i = 1,nat
            xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
            yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
            zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
            zan = dbl_mb(k_zanpt+i-1)
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            if(rr.lt.1.0d-3) go to 210
            rr3=rr*rr*rr
            elf(1) = elf(1) - zan*xn/rr3
            elf(2) = elf(2) - zan*yn/rr3
            elf(3) = elf(3) - zan*zn/rr3
  210    continue
c
c     ----- add -efc- contribution -----
c
!        if(iefc.eq.keyefc) then
!           do 220 i = 1,nefc
!              xn = efcc(1,i) - xp
!              yn = efcc(2,i) - yp
!              zn = efcc(3,i) - zp
!              rr =  sqrt(xn*xn + yn*yn + zn*zn)
!              if(rr.lt.tenm3) then
!                 write(luout,9993) xp,yp,zp,i
!                 go to 220
!              else
!                 rr3=rr*rr*rr
!                 elf(1) = elf(1) - efcz(i)*xn/rr3
!                 elf(2) = elf(2) - efcz(i)*yn/rr3
!                 elf(3) = elf(3) - efcz(i)*zn/rr3
!              endif
! 220       continue
!        endif
         if (.not. geom_cent_tag(geom,iat,at_tag)) call
     &      errquit('hnd_elfmap: geom_cent_tag failed',0,GEOM_ERR)
         if (.not. geom_tag_to_element(at_tag,symbol,element,atn)) call
     &      errquit('hnd_elfmap: geom_tag_to_element failed',0,GEOM_ERR)
         elfld   = sqrt(elf(1)*elf(1) + elf(2)*elf(2) + elf(3)*elf(3))
         write(luout,9995) iat,symbol,xp,yp,zp,elf(1),elf(2),elf(3),
     &                     elfld
c
c        ----- store ecce data -----
c
         call ecce_print1_char('atom name',symbol,1)
         call ecce_print1('EField components',MT_DBL,elf,3)
         call ecce_print1('EField',MT_DBL,elfld,1)
c
  230 continue ! Assembling and printing next atom
c
      call ecce_print_module_exit('EField','ok')
      call util_flush(luout)
c
c     ----- release memory block -----
c
  300 call ga_sync()
c
c     ------- Deallocate MA memory ------
c
      if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_elfmap, ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_efld)) call errquit
     &   ('hnd_elfmap, ma_pop_stack of l_efld failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_elfmap, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('elfmap: ga_destroy failed g_dens',0,GA_ERR)
      enddo
c
c     Terminate integrals
c
      call schwarz_tidy()
      call int_terminate()
c
      return
 9999 format(/,10x,14(1H-),/,10x,'Electric field',
     1       /,10x,14(1H-),/)
 9997 format(3x,'Atom ',6x,'X',9x,'Y',9x,'Z',8x,1x,
     1 15x,'Electric field (a.u.)',/,
     2 39x,7x,'X',14x,'Y',14x,'Z',11x,'Field',/,
     3 2x,96(1h-))
 9995 format(i5,1x,a2,3f10.5,4f15.6)
 9996 format(' --- Warning - electric field at ',
     1 3f10.5,' . contribution from nucleus ',i3,' ignored')
 9994 format(' 1 a.u. = 0.171524 10**(-8) dyn/esu ',/)
 9993 format(' --- warning - electric field at ',
     1 3f10.5,' . contribution from  -EFC-  ',i3,' ignored')
      end

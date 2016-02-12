      subroutine hnd_mtpole(rtdb,basis,geom,mpole)
c
c $Id$
c
c     This routine calculates the multipole moments around
c     choosen center of expansion
c
c     mpole = 1 : Dipole
c     mpole = 2 : Quadrupole
c     mpole = 3 : Octupole
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "cosmo.fh"
c
      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
      integer mpole   ! [input] order of multipole
c
c     Setting a max_pole. We have up to Octupole.
c
      integer max_pole, maxcomp
      parameter (max_pole = 3, maxcomp = (max_pole+1)*(max_pole+2)/2)
c
      double precision cm(3), mptval(maxcomp,3)
      double precision mptval2(maxcomp,3)
      double precision octo, debye, buck, cx, cy, cz
      double precision dipol, dipefc, diptot, rsquar, dum1, dum2
      logical status
      integer Nxyz(3)
      integer ncomp, nat, iat, icomp, i,j
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt
      integer nefc, l_efcc, k_efcc, l_efcz, k_efcz
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      character*3 scftyp
      character*16 at_tag
      data octo    /0.711688d+00/
      data buck    /1.344911d+00/
      data debye   /2.54176568d+00/
c
      if (ga_nodeid().eq.0.and.mpole.eq.1) write(luout,2001)
      if (ga_nodeid().eq.0.and.mpole.eq.2) write(luout,3001)
      if (ga_nodeid().eq.0.and.mpole.eq.3) write(luout,9001)
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
c     Determine number of components and zero mptval
c
      ncomp = (mpole+1)*(mpole+2)/2
      call dcopy(2*maxcomp,0.0d0,0,mptval,1)
c
c     Get center of expansion
c
      call prp_cent(rtdb,geom,cm)
c
c     Electronic contribution
c
      call hnd_mtpcon(basis,geom,g_dens(ndens),mptval,mpole,cm)
c
c     Only node 0 writing to rtdb
c
      status = rtdb_parallel(.false.)
c
      if (ga_nodeid().gt.0) goto 1000     
c
c     Allocate some local memory
c
      if (.not.geom_ncent(geom,nat)) call
     &   errquit('hnd_efgmap: geom_ncent failed',911,GEOM_ERR)
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_mtpole: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_mtpole: ma failed',911,MA_ERR)
c
      do iat=1,nat
         if (.not.geom_cent_get(geom,iat,at_tag,
     &       dbl_mb(k_xyzpt+3*(iat-1)),dbl_mb(k_zanpt+iat-1))) call
     &       errquit('hnd_mtpole: geom_cent_get failed',911,GEOM_ERR)
      enddo ! iat
c
c     For quadrupole case, need electronic contribution for 
c     diamagnetic susceptibility
c
      if (mpole.eq.2) rsquar=mptval(1,1)+mptval(4,1)+mptval(6,1)
c
c     Nuclear contribution 
c
      do icomp = 1, ncomp
        mptval(icomp,1) = -mptval(icomp,1)
        call getNxyz(mpole,icomp,Nxyz)
        do iat = 1, nat
           cx = 1d0
           cy = 1d0
           cz = 1d0
           if(Nxyz(1).ne.0)
     &         cx = (dbl_mb(k_xyzpt  +3*(iat-1))-cm(1))**Nxyz(1)
           if(Nxyz(2).ne.0)
     &         cy = (dbl_mb(k_xyzpt+1+3*(iat-1))-cm(2))**Nxyz(2)
           if(Nxyz(3).ne.0)
     &         cz = (dbl_mb(k_xyzpt+2+3*(iat-1))-cm(3))**Nxyz(3)
           mptval(icomp,1) = mptval(icomp,1) + 
     &                       dbl_mb(k_zanpt+iat-1) * cx * cy * cz
        enddo ! iat
      enddo   ! icomp
c
c     ----- form -efc- contribution -----
c           from cosmo point charges !!!!
c
      if (cosmo_last) then
         if (.not.rtdb_get(rtdb,'cosmo:nefc',mt_int,1,nefc))
     &      call errquit('hnd_mtpole: rtdb get failed for nefc ',911,
     &      rtdb_err)
         if (.not.ma_push_get(mt_dbl,nefc*3,'efcc',l_efcc,k_efcc))
     &      call errquit('hnd_mtpole: malloc k_efcc fail',911,ma_err)
         if (.not.ma_push_get(mt_dbl,nefc,'efcz',l_efcz,k_efcz))
     &      call errquit('hnd_mtpole: malloc k_efcz fail',911,ma_err)
         if (.not.rtdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,
     &      dbl_mb(k_efcc))) call
     &      errquit('hnd_mtpole: rtdb get failed efcc',912,rtdb_err)
         if (.not.rtdb_get(rtdb,'cosmo:efcz',mt_dbl,nefc,
     &      dbl_mb(k_efcz))) call
     &      errquit('hnd_mtpole: rtdb get failed efcz',913,rtdb_err)
         do icomp = 1, ncomp
            call getNxyz(2,icomp,Nxyz)
            do iat=1,nefc
               cx=1d0
               cy=1d0
               cz=1d0
               if(Nxyz(1).ne.0)
     &             cx = (dbl_mb(k_efcc+3*(icomp-1)  )-cm(1))**Nxyz(1)
               if(Nxyz(2).ne.0)
     &             cy = (dbl_mb(k_efcc+3*(icomp-1)+1)-cm(2))**Nxyz(2)
               if(Nxyz(3).ne.0)
     &             cz = (dbl_mb(k_efcc+3*(icomp-1)+2)-cm(3))**Nxyz(3)
               mptval(icomp,2) = mptval(icomp,2)
     &                         + dbl_mb(k_efcz+iat-1)*cx*cy*cz
            enddo ! iat
         enddo ! icomp
         if (.not.ma_chop_stack(l_efcc)) call
     &      errquit('hnd_mtpole: chop stack l_efcc',913,ma_err)
      endif
c
c     Print output for
c     mpole = 1 : Dipole
c     mpole = 2 : Quadrupole
c     mpole = 3 : Octupole
c
      goto (101,102,103) mpole
c
c     Dipole output
c
  101 dipol  = sqrt(mptval(1,1)*mptval(1,1) +
     &              mptval(2,1)*mptval(2,1) +
     &              mptval(3,1)*mptval(3,1))
      dipefc = sqrt(mptval(1,2)*mptval(1,2) +
     &              mptval(2,2)*mptval(2,2) +
     &              mptval(3,2)*mptval(3,2))
      diptot = sqrt((mptval(1,1)+mptval(1,2))**2 +
     &              (mptval(2,1)+mptval(2,2))**2 +
     &              (mptval(3,1)+mptval(3,2))**2)
      write(luout,2995) dipol,mptval(1,1),mptval(1,2),mptval(2,1),
     &                  mptval(2,2),mptval(3,1),mptval(3,2),dipefc,
     &                  diptot
c
      do icomp = 1, ncomp
         mptval(icomp,1) = mptval(icomp,1) * debye
         mptval(icomp,2) = mptval(icomp,2) * debye
         mptval(icomp,3) = mptval(icomp,1) + mptval(icomp,2)
      enddo
      dipol  =  dipol*debye
      dipefc = dipefc*debye
      diptot = diptot*debye
      write(luout,2996) dipol,mptval(1,1),mptval(1,2),mptval(2,1),
     &                  mptval(2,2),mptval(3,1),mptval(3,2),dipefc,
     &                  diptot
      write(luout,2994)
      if (.not.rtdb_put(rtdb,'prop:dipval',mt_dbl,3,mptval(1,3)))
     &     call errquit('prop: rtdb_put of dipval failed',555,
     &       RTDB_ERR)
c
c     Done, goto cleanup and sync with other nodes
c
      goto 999 
c
c     Quadrupole output
c
c     Ordering in mptval (from defNxyz):
c     xx, xy, xz, yy, yz, zz
c
c     First < r**2 > = diamagnetic susceptibility
c
  102 write(luout,3992) rsquar
c
c     Second moments
c
      write(luout,3999)
      write(luout,3996) mptval(1,1),mptval(1,2),mptval(1,1)+mptval(1,2),
     1                  mptval(4,1),mptval(4,2),mptval(4,1)+mptval(4,2),
     2                  mptval(6,1),mptval(6,2),mptval(6,1)+mptval(6,2),
     3                  mptval(2,1),mptval(2,2),mptval(2,1)+mptval(2,2),
     4                  mptval(3,1),mptval(3,2),mptval(3,1)+mptval(3,2),
     5                  mptval(5,1),mptval(5,2),mptval(5,1)+mptval(5,2)
c
      do j = 1, 2
        do i = 1, 6
          mptval2(i,j) = mptval(i,j)*buck
        enddo
      enddo
c
      write(luout,3995)
      write(luout,3996) 
     1    mptval2(1,1),mptval2(1,2),mptval2(1,1)+mptval2(1,2),
     2    mptval2(4,1),mptval2(4,2),mptval2(4,1)+mptval2(4,2),
     3    mptval2(6,1),mptval2(6,2),mptval2(6,1)+mptval2(6,2),
     4    mptval2(2,1),mptval2(2,2),mptval2(2,1)+mptval2(2,2),
     5    mptval2(3,1),mptval2(3,2),mptval2(3,1)+mptval2(3,2),
     6    mptval2(5,1),mptval2(5,2),mptval2(5,1)+mptval2(5,2)
c
c     Quadrupole moments
c
      rsquar = mptval(1,1)+mptval(4,1)+mptval(6,1)
      mptval(1,1)=(3.0d0*mptval(1,1)-rsquar)/2.0d0
      mptval(4,1)=(3.0d0*mptval(4,1)-rsquar)/2.0d0
      mptval(6,1)=(3.0d0*mptval(6,1)-rsquar)/2.0d0
      mptval(2,1)= 3.0d0*mptval(2,1)/2.0d0
      mptval(3,1)= 3.0d0*mptval(3,1)/2.0d0
      mptval(5,1)= 3.0d0*mptval(5,1)/2.0d0
      rsquar = mptval(1,2)+mptval(4,2)+mptval(6,2)
      mptval(1,2)=(3.0d0*mptval(1,2)-rsquar)/2.0d0
      mptval(4,2)=(3.0d0*mptval(4,2)-rsquar)/2.0d0
      mptval(6,2)=(3.0d0*mptval(6,2)-rsquar)/2.0d0
      mptval(2,2)= 3.0d0*mptval(2,2)/2.0d0
      mptval(3,2)= 3.0d0*mptval(3,2)/2.0d0
      mptval(5,2)= 3.0d0*mptval(5,2)/2.0d0
      write(luout,3997)
      write(luout,3996) mptval(1,1),mptval(1,2),mptval(1,1)+mptval(1,2),
     1                  mptval(4,1),mptval(4,2),mptval(4,1)+mptval(4,2),
     2                  mptval(6,1),mptval(6,2),mptval(6,1)+mptval(6,2),
     3                  mptval(2,1),mptval(2,2),mptval(2,1)+mptval(2,2),
     4                  mptval(3,1),mptval(3,2),mptval(3,1)+mptval(3,2),
     5                  mptval(5,1),mptval(5,2),mptval(5,1)+mptval(5,2)
c
      do j = 1, 2
        do i = 1, 6
          mptval2(i,j) = mptval(i,j)*buck
        enddo
      enddo
c
      write(luout,3994)
      write(luout,3996) 
     1    mptval2(1,1),mptval2(1,2),mptval2(1,1)+mptval2(1,2),
     2    mptval2(4,1),mptval2(4,2),mptval2(4,1)+mptval2(4,2),
     3    mptval2(6,1),mptval2(6,2),mptval2(6,1)+mptval2(6,2),
     4    mptval2(2,1),mptval2(2,2),mptval2(2,1)+mptval2(2,2),
     5    mptval2(3,1),mptval2(3,2),mptval2(3,1)+mptval2(3,2),
     6    mptval2(5,1),mptval2(5,2),mptval2(5,1)+mptval2(5,2)
      write(luout,3993)
c
c     Done, goto cleanup and sync with other nodes
c
      goto 999 
c
c     Octupole output
c
c     Ordering in mptval (from defNxyz):
c     xxx, xxy, xxz, yyx, xyz, zzx, yyy, zzy, zzz
c
c     Third moments
c
  103 write(luout,9999)
      write(luout,9996) 
     1              mptval(1,1),mptval(1,2),mptval(1,1)+mptval(1,2),
     2              mptval(7,1),mptval(7,2),mptval(7,1)+mptval(7,2),
     3              mptval(10,1),mptval(10,2),mptval(10,1)+mptval(10,2),
     4              mptval(2,1),mptval(2,2),mptval(2,1)+mptval(2,2),
     5              mptval(3,1),mptval(3,2),mptval(3,1)+mptval(3,2),
     6              mptval(4,1),mptval(4,2),mptval(4,1)+mptval(4,2),
     7              mptval(8,1),mptval(8,2),mptval(8,1)+mptval(8,2),
     8              mptval(6,1),mptval(6,2),mptval(6,1)+mptval(6,2),
     9              mptval(9,1),mptval(9,2),mptval(9,1)+mptval(9,2),
     1              mptval(5,1),mptval(5,2),mptval(5,1)+mptval(5,2)
c
      do j = 1, 2
        do i = 1, 10
          mptval2(i,j) = mptval(i,j)*octo
        enddo
      enddo
c
      write(luout,9995)
      write(luout,9996)
     1        mptval2(1,1),mptval2(1,2),mptval2(1,1)+mptval2(1,2),
     2        mptval2(7,1),mptval2(7,2),mptval2(7,1)+mptval2(7,2),
     3        mptval2(10,1),mptval2(10,2),mptval2(10,1)+mptval2(10,2),
     4        mptval2(2,1),mptval2(2,2),mptval2(2,1)+mptval2(2,2),
     5        mptval2(3,1),mptval2(3,2),mptval2(3,1)+mptval2(3,2),
     6        mptval2(4,1),mptval2(4,2),mptval2(4,1)+mptval2(4,2),
     7        mptval2(8,1),mptval2(8,2),mptval2(8,1)+mptval2(8,2),
     8        mptval2(6,1),mptval2(6,2),mptval2(6,1)+mptval2(6,2),
     9        mptval2(9,1),mptval2(9,2),mptval2(9,1)+mptval2(9,2),
     1        mptval2(5,1),mptval2(5,2),mptval2(5,1)+mptval2(5,2)
c
c     Octupole moments
c
      do i=1,2
         dum1=mptval(1,i)
         dum2=mptval(4,i)
         mptval(1,i)=mptval(1,i)-3.0d0*(mptval(4,i)+mptval(6,i))/2.0d0
         mptval(4,i)=(4.0d0*mptval(4,i)-dum1-mptval(6,i))/2.0d0
         mptval(6,i)=(4.0d0*mptval(6,i)-dum1-dum2)/2.0d0
         dum1=mptval(7,i)
         dum2=mptval(2,i)
         mptval(7,i)=mptval(7,i)-3.0d0*(mptval(2,i)+mptval(9,i))/2.0d0
         mptval(2,i)=(4.0d0*mptval(2,i)-dum1-mptval(9,i))/2.0d0
         mptval(9,i)=(4.0d0*mptval(9,i)-dum1-dum2)/2.0d0
         dum1=mptval(10,i)
         dum2=mptval(3,i)
         mptval(10,i)=mptval(10,i)-3.0d0*(mptval(3,i)+mptval(8,i))/2.0d0
         mptval(3,i)=(4.0d0*mptval(3,i)-dum1-mptval(8,i))/2.0d0
         mptval(8,i)=(4.0d0*mptval(8,i)-dum1-dum2)/2.0d0
         mptval(5,i)=5.0d0*mptval(5,i)/2.0d0
      enddo
      write(luout,9997)
      write(luout,9996) 
     1              mptval(1,1),mptval(1,2),mptval(1,1)+mptval(1,2),
     2              mptval(7,1),mptval(7,2),mptval(7,1)+mptval(7,2),
     3              mptval(10,1),mptval(10,2),mptval(10,1)+mptval(10,2),
     4              mptval(2,1),mptval(2,2),mptval(2,1)+mptval(2,2),
     5              mptval(3,1),mptval(3,2),mptval(3,1)+mptval(3,2),
     6              mptval(4,1),mptval(4,2),mptval(4,1)+mptval(4,2),
     7              mptval(8,1),mptval(8,2),mptval(8,1)+mptval(8,2),
     8              mptval(6,1),mptval(6,2),mptval(6,1)+mptval(6,2),
     9              mptval(9,1),mptval(9,2),mptval(9,1)+mptval(9,2),
     1              mptval(5,1),mptval(5,2),mptval(5,1)+mptval(5,2)
c
      do j = 1, 2
        do i = 1, 10
          mptval2(i,j) = mptval(i,j)*octo
        enddo
      enddo
c
      write(luout,9994)
      write(luout,9996)
     1        mptval2(1,1),mptval2(1,2),mptval2(1,1)+mptval2(1,2),
     2        mptval2(7,1),mptval2(7,2),mptval2(7,1)+mptval2(7,2),
     3        mptval2(10,1),mptval2(10,2),mptval2(10,1)+mptval2(10,2),
     4        mptval2(2,1),mptval2(2,2),mptval2(2,1)+mptval2(2,2),
     5        mptval2(3,1),mptval2(3,2),mptval2(3,1)+mptval2(3,2),
     6        mptval2(4,1),mptval2(4,2),mptval2(4,1)+mptval2(4,2),
     7        mptval2(8,1),mptval2(8,2),mptval2(8,1)+mptval2(8,2),
     8        mptval2(6,1),mptval2(6,2),mptval2(6,1)+mptval2(6,2),
     9        mptval2(9,1),mptval2(9,2),mptval2(9,1)+mptval2(9,2),
     1        mptval2(5,1),mptval2(5,2),mptval2(5,1)+mptval2(5,2)
      write(luout,9993)
c
c     ----- release MA memory blocks -----
c
  999 if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_mtpole, ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_mtpole, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
c
c     Synchronize all nodes (as only node 0 was writing)
c     Reset rtdb access to parallel
c
 1000 call ga_sync()
      status = rtdb_parallel(.true.)
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('mtpole: ga_destroy failed g_dens',0,GA_ERR)
      enddo
c
c     Terminate integrals
c
      call schwarz_tidy()
      call int_terminate()
c
      return
c
c     Formatting dipole
c
 2001 format(/,10x,13(1h-),/,10x,'Dipole Moment',/,10x,13(1h-))
 2996 format(/,2x,' Dipole moment',f20.10,' Debye(s)',
     2 /,12x,' DMX',f20.10,' DMXEFC',f20.10,
     3 /,12x,' DMY',f20.10,' DMYEFC',f20.10,
     4 /,12x,' DMZ',f20.10,' DMZEFC',f20.10,
     5 /,2x,' -EFC- dipole ',f20.10,' DEBYE(S)',
     6 /,2x,' Total dipole ',f20.10,' DEBYE(S)')
 2995 format(/,2x,' Dipole moment',f20.10,' A.U.',
     2 /,12x,' DMX',f20.10,' DMXEFC',f20.10,
     3 /,12x,' DMY',f20.10,' DMYEFC',f20.10,
     4 /,12x,' DMZ',f20.10,' DMZEFC',f20.10,
     5 /, 2x,' -EFC- dipole ',f20.10,' A.U.',
     6 /, 2x,' Total dipole ',f20.10,' A.U.')
 2994 format(/,' 1 a.u. = 2.541766 Debyes ')
c
c     Formatting quadrupole
c
 3001 format(/,10x,17(1h-),/,10x,'Quadrupole Moment',/,10x,17(1h-))
 3999 format(/,2x,' Second moments in atomic units')
 3997 format(/,2x,' Quadrupole moments in atomic units')
 3995 format(/,2x,' Second moments in buckingham(s)')
 3994 format(/,2x,' Quadrupole moments in buckingham(s)')
 3996 format(/,2x,' Component','  Electronic+nuclear',4x,
     2      ' Point charges',12x,' Total',/,2x,74(1h-),
     2 /,2x,'    XX    ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    YY    ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    ZZ    ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    XY    ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    XZ    ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    YZ    ',f20.10,2x,f20.10,2x,f20.10)
 3993 format(/,' 1 a.u. = 1.344911 Buckinghams ',
     1                 '= 1.344911 10**(-26) esu*cm**2 ')
 3992 format(/,' < R**2 > = ',f10.6,' a.u. ',
     1 ' ( 1 a.u. = 0.280023 10**(-16) cm**2 ) ',/,
     2 ' ( also called diamagnetic susceptibility ) ')
c
c     Formatted printing octupole moments
c
 9001 format(/,10x,15(1h-),/,10x,'Octupole Moment',/,10x,15(1h-))
 9999 format(/,2x,' Third moments in atomic units')
 9997 format(/,2x,' Octupole moments in atomic units')
 9995 format(/,2x,' Third moments in 10**(-34) esu*cm**3')
 9994 format(/,2x,' Octupole moments in 10**(-34) esu*cm**3')
 9996 format(/,2x,' Component','  Electronic+nuclear',4x,
     2      ' Point charges',12x,' Total',/,2x,74(1h-),
     2 /,2x,'    XXX   ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    YYY   ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    ZZZ   ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    XXY   ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    XXZ   ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    YYX   ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    YYZ   ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    ZZX   ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    ZZY   ',f20.10,2x,f20.10,2x,f20.10,
     2 /,2x,'    XYZ   ',f20.10,2x,f20.10,2x,f20.10)
 9993 format(/,' 1 a.u. = 0.711688 10**(-34) esu*cm**3 ')
      end
c
      subroutine prp_cent(rtdb,geom,cm)
c
c     Routine calculates coordinates of the expansion center
c     depending on what kind of expansion the user chooses
c     1: center of charge (default)
c     2: center of mass
c     3: origin
c     4: arbitrary point
c
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
c
      integer rtdb           ! [input] rtdb handle
      integer geom           ! [input] geometry handle
      double precision cm(3) ! [output] expansion center
c
      integer i, i_expan
      double precision scale
c
c     Read in point of expansion
c
      if (.not. rtdb_get(rtdb, 'prop:center', mt_int, 1, i_expan))
     $  i_expan = 1    ! default of center of charge
c
      if (i_expan.eq.4) then   ! arbitrary point of expansion
        if (.not. rtdb_get(rtdb, 'prop:center_val',mt_dbl,3,cm))
     $  call errquit('prp_cent: rtdb_get center_val failed',555,
     &       RTDB_ERR)
c
c get scale factor based on geometry UNITS keyword
c
        if (.not. geom_get_user_scale(geom,scale))
     $    call errquit('prp_cent: trouble getting units scale',555,
     &       GEOM_ERR)
        do i = 1, 3             ! scale according to geometry units
          cm(i) = cm(i)*scale
        end do
        if(ga_nodeid().eq.0) write(luout,5994) cm(1),cm(2),cm(3)
c
      elseif (i_expan.eq.3) then   ! origin
        do i = 1, 3
          cm(i) = 0.d0
        end do
        if (ga_nodeid().eq.0) write(luout,5995) cm(1),cm(2),cm(3)
c
      elseif (i_expan.eq.2) then   ! center of mass
        if (.not.geom_center_of_mass(geom,cm)) call
     &     errquit('prp_cent: geom_center_of_mass failed',0,GEOM_ERR)
        if (ga_nodeid().eq.0) write(luout,5996) cm(1),cm(2),cm(3)
c
      elseif (i_expan.eq.1) then   ! center of charge
        if (.not.geom_center_of_charge(geom,cm)) call
     &     errquit('prp_cent: geom_center_of_charge failed',0,GEOM_ERR)
        if (ga_nodeid().eq.0) write(luout,5997) cm(1),cm(2),cm(3)
      else
        call errquit('prp_cent: improper center value',555, INPUT_ERR)
      endif
c
      return
 5994 format(/,' Center of expansion (in au) is the arbitrary point',
     1 /,8x,' X = ',f15.7,' Y = ',f15.7,' Z = ',f15.7)
 5995 format(/,' Center of expansion (in au) is the origin',
     1 /,8x,' X = ',f15.7,' Y = ',f15.7,' Z = ',f15.7)
 5996 format(/,' Center of mass (in au) is the expansion point',
     1 /,8x,' X = ',f15.7,' Y = ',f15.7,' Z = ',f15.7)
 5997 format(/,' Center of charge (in au) is the expansion point',
     1 /,8x,' X = ',f15.7,' Y = ',f15.7,' Z = ',f15.7)
      end
c
c
c
c
c
c     Jeff wrote a new version which returns mtpval - in atomic units - and
c     prints less stuff, for use in the TCE.
c
      subroutine hnd_mtpole2(rtdb,basis,geom,mpole,mtpval)
c
c $Id$
c
c     This routine calculates the multipole moments around
c     choosen center of expansion
c
c     mpole = 1 : Dipole
c     mpole = 2 : Quadrupole
c     mpole = 3 : Octupole
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "cosmo.fh"
c
      integer rtdb    ! [input] rtdb handle
      integer basis   ! [input] basis handle
      integer geom    ! [input] geometry handle
      integer mpole   ! [input] order of multipole
      integer max_pole, maxcomp
      parameter (max_pole = 3, maxcomp = (max_pole+1)*(max_pole+2)/2)
c
      double precision cm(3), mtpval(maxcomp,3)
      double precision cx, cy, cz
      double precision dipol, dipefc, diptot, rsquar, dum1, dum2
      logical status
      integer Nxyz(3)
      integer ncomp, nat, iat, icomp, i
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt
      integer nefc, l_efcc, k_efcc, l_efcz, k_efcz
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      character*3 scftyp
      character*16 at_tag
      if (ga_nodeid().eq.0.and.mpole.eq.1) write(luout,2001)
      if (ga_nodeid().eq.0.and.mpole.eq.2) write(luout,3001)
      if (ga_nodeid().eq.0.and.mpole.eq.3) write(luout,9001)
c
c     Get density matrix
c
      call hnd_prp_get_dens(rtdb,geom,basis,g_dens,ndens,scftyp,
     &                      nclosed,nopen,nvirt)
c
c     Determine number of components and zero mtpval
c
      ncomp = (mpole+1)*(mpole+2)/2
c
c     It would appear dcopy screws things up be zeroing the wrong area.
c     This bothers me but I do not care now that it is fixed.
c      call dcopy(2*maxcomp,0.0d0,0,mtpval,1)
c
      do icomp = 1, ncomp
         mtpval(icomp,1) = 0.0d0
         mtpval(icomp,2) = 0.0d0
         mtpval(icomp,3) = 0.0d0
      enddo
c
c     Get center of expansion
c
      call prp_cent(rtdb,geom,cm)
c
c     Electronic contribution
c
      call hnd_mtpcon(basis,geom,g_dens(ndens),mtpval,mpole,cm)
c
c     Only node 0 writing to rtdb
c
      status = rtdb_parallel(.false.)
c
      if (ga_nodeid().gt.0) goto 1000
c
c     Allocate some local memory
c
      if (.not.geom_ncent(geom,nat)) call
     &   errquit('hnd_efgmap: geom_ncent failed',911,GEOM_ERR)
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_mtpole2: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_mtpole2: ma failed',911,MA_ERR)
c
      do iat=1,nat
         if (.not.geom_cent_get(geom,iat,at_tag,
     &       dbl_mb(k_xyzpt+3*(iat-1)),dbl_mb(k_zanpt+iat-1))) call
     &       errquit('hnd_mtpole2: geom_cent_get failed',911,GEOM_ERR)
      enddo ! iat
c
c     For quadrupole case, need electronic contribution for
c     diamagnetic susceptibility
c
      if (mpole.eq.2) rsquar=mtpval(1,1)+mtpval(4,1)+mtpval(6,1)
c
c     Nuclear contribution
c
      do icomp = 1, ncomp
        mtpval(icomp,1) = -mtpval(icomp,1)
        call getNxyz(mpole,icomp,Nxyz)
        do iat = 1, nat
           cx = 1d0
           if(Nxyz(1).ne.0)
     &     cx = (dbl_mb(k_xyzpt  +3*(iat-1))-cm(1))**Nxyz(1)
           cy = 1d0
           if(Nxyz(2).ne.0)
     &     cy = (dbl_mb(k_xyzpt+1+3*(iat-1))-cm(2))**Nxyz(2)
           cz = 1d0
           if(Nxyz(3).ne.0)
     &     cz = (dbl_mb(k_xyzpt+2+3*(iat-1))-cm(3))**Nxyz(3)
           mtpval(icomp,1) = mtpval(icomp,1) +
     &                       dbl_mb(k_zanpt+iat-1) * cx * cy * cz
        enddo ! iat
      enddo ! icomp
c
c     Print output for
c     mpole = 1 : Dipole
c     mpole = 2 : Quadrupole
c     mpole = 3 : Octupole
c
      goto (101,102,103) mpole
c
c     Dipole output
c
  101 continue
c
      do icomp = 1, ncomp
         mtpval(icomp,1) = -1.0*mtpval(icomp,1) 
      enddo
c
      goto 999 
c
c     Quadrupole output
c
c     Ordering in mtpval (from defNxyz):
c     xx, xy, xz, yy, yz, zz
c
  102 continue
c
c     Quadrupole moments
c
      rsquar = mtpval(1,1)+mtpval(4,1)+mtpval(6,1)
      mtpval(1,1)=(3.0d0*mtpval(1,1)-rsquar)/2.0d0
      mtpval(4,1)=(3.0d0*mtpval(4,1)-rsquar)/2.0d0
      mtpval(6,1)=(3.0d0*mtpval(6,1)-rsquar)/2.0d0
      mtpval(2,1)= 3.0d0*mtpval(2,1)/2.0d0
      mtpval(3,1)= 3.0d0*mtpval(3,1)/2.0d0
      mtpval(5,1)= 3.0d0*mtpval(5,1)/2.0d0
      rsquar = mtpval(1,2)+mtpval(4,2)+mtpval(6,2)
      mtpval(1,2)=(3.0d0*mtpval(1,2)-rsquar)/2.0d0
      mtpval(4,2)=(3.0d0*mtpval(4,2)-rsquar)/2.0d0
      mtpval(6,2)=(3.0d0*mtpval(6,2)-rsquar)/2.0d0
      mtpval(2,2)= 3.0d0*mtpval(2,2)/2.0d0
      mtpval(3,2)= 3.0d0*mtpval(3,2)/2.0d0
      mtpval(5,2)= 3.0d0*mtpval(5,2)/2.0d0
c
      goto 999 
c
c     Octupole output
c
c     Ordering in mtpval (from defNxyz):
c     xxx, xxy, xxz, yyx, xyz, zzx, yyy, zzy, zzz
c
  103 continue
c
c     Octupole moments
c
      do i=1,2
         dum1=mtpval(1,i)
         dum2=mtpval(4,i)
         mtpval(1,i)=mtpval(1,i)-3.0d0*(mtpval(4,i)+mtpval(6,i))/2.0d0
         mtpval(4,i)=(4.0d0*mtpval(4,i)-dum1-mtpval(6,i))/2.0d0
         mtpval(6,i)=(4.0d0*mtpval(6,1)-dum1-dum2)/2.0d0
         dum1=mtpval(7,i)
         dum2=mtpval(2,i)
         mtpval(7,i)=mtpval(7,i)-3.0d0*(mtpval(2,i)+mtpval(9,i))/2.0d0
         mtpval(2,i)=(4.0d0*mtpval(2,i)-dum1-mtpval(9,i))/2.0d0
         mtpval(9,i)=(4.0d0*mtpval(9,i)-dum1-dum2)/2.0d0
         dum1=mtpval(10,i)
         dum2=mtpval(3,i)
         mtpval(10,i)=mtpval(10,i)-3.0d0*(mtpval(3,i)+mtpval(8,i))/2.0d0
         mtpval(3,i)=(4.0d0*mtpval(3,i)-dum1-mtpval(8,i))/2.0d0
         mtpval(8,i)=(4.0d0*mtpval(8,i)-dum1-dum2)/2.0d0
         mtpval(5,i)=5.0d0*mtpval(5,i)/2.0d0
      enddo
c
c     ----- release MA memory blocks -----
c
  999 if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_mtpole2, ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_mtpole2, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
c
c     Synchronize all nodes (as only node 0 was writing)
c     Reset rtdb access to parallel
c
 1000 call ga_sync()
      status = rtdb_parallel(.true.)
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('mtpole: ga_destroy failed g_dens',0,GA_ERR)
      enddo
c
      return
c
 2001 format(/,10x,13(1h-),/,10x,'Dipole Moment',/,10x,13(1h-))
 3001 format(/,10x,17(1h-),/,10x,'Quadrupole Moment',/,10x,17(1h-))
 9001 format(/,10x,15(1h-),/,10x,'Octupole Moment',/,10x,15(1h-))
      end

      subroutine hnd_mtpole_nuclear(rtdb,geom,mpole,mtpval)
c
c $Id$
c
c     This routine calculates the multipole moments around
c     choosen center of expansion
c
c     mpole = 1 : Dipole
c     mpole = 2 : Quadrupole
c     mpole = 3 : Octupole
c
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "cosmo.fh"
c
      integer rtdb    ! [input] rtdb handle
      integer geom    ! [input] geometry handle
      integer mpole   ! [input] order of multipole
      integer max_pole, maxcomp
      parameter (max_pole = 3, maxcomp = (max_pole+1)*(max_pole+2)/2)
c
      double precision cm(3), mtpval(maxcomp,3)
      double precision cx, cy, cz
      double precision dipol, dipefc, diptot, rsquar, dum1, dum2
      logical status
      integer Nxyz(3)
      integer ncomp, nat, iat, icomp, i
      integer l_xyzpt, k_xyzpt, l_zanpt, k_zanpt
      integer nefc, l_efcc, k_efcc, l_efcz, k_efcz
      integer g_dens(3),ndens,nclosed(2),nopen(2),nvirt(2)
      character*3 scftyp
      character*16 at_tag
      if (ga_nodeid().eq.0.and.mpole.eq.1) write(luout,2001)
      if (ga_nodeid().eq.0.and.mpole.eq.2) write(luout,3001)
      if (ga_nodeid().eq.0.and.mpole.eq.3) write(luout,9001)
c
c     Determine number of components and zero mtpval
c
      ncomp = (mpole+1)*(mpole+2)/2
c
c     It would appear dcopy screws things up be zeroing the wrong area.
c     This bothers me but I do not care now that it is fixed.
c      call dcopy(2*maxcomp,0.0d0,0,mtpval,1)
c
      do icomp = 1, ncomp
         mtpval(icomp,1) = 0.0d0
         mtpval(icomp,2) = 0.0d0
         mtpval(icomp,3) = 0.0d0
      enddo
c
c     Get center of expansion
c
      call prp_cent(rtdb,geom,cm)
c
c     Only node 0 writing to rtdb
c
      status = rtdb_parallel(.false.)
c
      if (ga_nodeid().gt.0) goto 1000
c
c     Allocate some local memory
c
      if (.not.geom_ncent(geom,nat)) call
     &   errquit('hnd_efgmap: geom_ncent failed',911,GEOM_ERR)
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_mtpole2: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_mtpole2: ma failed',911,MA_ERR)
c
      do iat=1,nat
         if (.not.geom_cent_get(geom,iat,at_tag,
     &       dbl_mb(k_xyzpt+3*(iat-1)),dbl_mb(k_zanpt+iat-1))) call
     &       errquit('hnd_mtpole2: geom_cent_get failed',911,GEOM_ERR)
      enddo ! iat
c
c     For quadrupole case, need electronic contribution for
c     diamagnetic susceptibility
c
      if (mpole.eq.2) rsquar=mtpval(1,1)+mtpval(4,1)+mtpval(6,1)
c
c     Nuclear contribution
c
      do icomp = 1, ncomp
        mtpval(icomp,1) = -mtpval(icomp,1)
        call getNxyz(mpole,icomp,Nxyz)
        do iat = 1, nat
           cx = 1d0
           if(Nxyz(1).ne.0)
     &     cx = (dbl_mb(k_xyzpt  +3*(iat-1))-cm(1))**Nxyz(1)
           cy = 1d0
           if(Nxyz(2).ne.0)
     &     cy = (dbl_mb(k_xyzpt+1+3*(iat-1))-cm(2))**Nxyz(2)
           cz = 1d0
           if(Nxyz(3).ne.0)
     &     cz = (dbl_mb(k_xyzpt+2+3*(iat-1))-cm(3))**Nxyz(3)
           mtpval(icomp,1) = mtpval(icomp,1) +
     &                       dbl_mb(k_zanpt+iat-1) * cx * cy * cz
        enddo ! iat
      enddo ! icomp
c
c     Print output for
c     mpole = 1 : Dipole
c     mpole = 2 : Quadrupole
c     mpole = 3 : Octupole
c
      goto (101,102,103) mpole
c
c     Dipole output
c
  101 continue
c
      do icomp = 1, ncomp
         mtpval(icomp,1) = -1.0*mtpval(icomp,1) 
      enddo
c
      goto 999 
c
c     Quadrupole output
c
c     Ordering in mtpval (from defNxyz):
c     xx, xy, xz, yy, yz, zz
c
  102 continue
c
c     Quadrupole moments
c
      rsquar = mtpval(1,1)+mtpval(4,1)+mtpval(6,1)
      mtpval(1,1)=(3.0d0*mtpval(1,1)-rsquar)/2.0d0
      mtpval(4,1)=(3.0d0*mtpval(4,1)-rsquar)/2.0d0
      mtpval(6,1)=(3.0d0*mtpval(6,1)-rsquar)/2.0d0
      mtpval(2,1)= 3.0d0*mtpval(2,1)/2.0d0
      mtpval(3,1)= 3.0d0*mtpval(3,1)/2.0d0
      mtpval(5,1)= 3.0d0*mtpval(5,1)/2.0d0
      rsquar = mtpval(1,2)+mtpval(4,2)+mtpval(6,2)
      mtpval(1,2)=(3.0d0*mtpval(1,2)-rsquar)/2.0d0
      mtpval(4,2)=(3.0d0*mtpval(4,2)-rsquar)/2.0d0
      mtpval(6,2)=(3.0d0*mtpval(6,2)-rsquar)/2.0d0
      mtpval(2,2)= 3.0d0*mtpval(2,2)/2.0d0
      mtpval(3,2)= 3.0d0*mtpval(3,2)/2.0d0
      mtpval(5,2)= 3.0d0*mtpval(5,2)/2.0d0
c
      goto 999 
c
c     Octupole output
c
c     Ordering in mtpval (from defNxyz):
c     xxx, xxy, xxz, yyx, xyz, zzx, yyy, zzy, zzz
c
  103 continue
c
c     Octupole moments
c
      do i=1,2
         dum1=mtpval(1,i)
         dum2=mtpval(4,i)
         mtpval(1,i)=mtpval(1,i)-3.0d0*(mtpval(4,i)+mtpval(6,i))/2.0d0
         mtpval(4,i)=(4.0d0*mtpval(4,i)-dum1-mtpval(6,i))/2.0d0
         mtpval(6,i)=(4.0d0*mtpval(6,1)-dum1-dum2)/2.0d0
         dum1=mtpval(7,i)
         dum2=mtpval(2,i)
         mtpval(7,i)=mtpval(7,i)-3.0d0*(mtpval(2,i)+mtpval(9,i))/2.0d0
         mtpval(2,i)=(4.0d0*mtpval(2,i)-dum1-mtpval(9,i))/2.0d0
         mtpval(9,i)=(4.0d0*mtpval(9,i)-dum1-dum2)/2.0d0
         dum1=mtpval(10,i)
         dum2=mtpval(3,i)
         mtpval(10,i)=mtpval(10,i)-3.0d0*(mtpval(3,i)+mtpval(8,i))/2.0d0
         mtpval(3,i)=(4.0d0*mtpval(3,i)-dum1-mtpval(8,i))/2.0d0
         mtpval(8,i)=(4.0d0*mtpval(8,i)-dum1-dum2)/2.0d0
         mtpval(5,i)=5.0d0*mtpval(5,i)/2.0d0
      enddo
c
c     ----- release MA memory blocks -----
c
  999 if (.not.ma_pop_stack(l_zanpt)) call errquit
     &   ('hnd_mtpole2, ma_pop_stack of l_zanpt failed',911,MA_ERR)
      if (.not.ma_pop_stack(l_xyzpt)) call errquit
     &   ('hnd_mtpole2, ma_pop_stack of l_xyzpt failed',911,MA_ERR)
c
c     Synchronize all nodes (as only node 0 was writing)
c     Reset rtdb access to parallel
c
 1000 call ga_sync()
      status = rtdb_parallel(.true.)
c
      do i = 1, ndens
         if (.not.ga_destroy(g_dens(i))) call
     &       errquit('mtpole: ga_destroy failed g_dens',0,GA_ERR)
      enddo
c
      return
c
 2001 format(/,10x,13(1h-),/,10x,'Dipole Moment',/,10x,13(1h-))
 3001 format(/,10x,17(1h-),/,10x,'Quadrupole Moment',/,10x,17(1h-))
 9001 format(/,10x,15(1h-),/,10x,'Octupole Moment',/,10x,15(1h-))
      end

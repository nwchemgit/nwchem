c =============== Fredy Aquino's routines ======= START
c Another version of grid_zoracorr
      subroutine grid_zoracorr_FA(nqpts,qxyz,qwght,
     &                            natoms,g_dens,amat)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "zora.fh"
#include "geom.fh"
c
      integer nqpts,iatom
      integer g_dens(2),igrid,natoms,npol
      integer closegridpts(nqpts)
      double precision qxyz(*),qwght(*)   ! [ in  ]
      double precision amat(nqpts,ipol)   ! [ out ]
      double precision amat_coul(nqpts,ipol)
      double precision amat_nucl(nqpts)
      double precision amat_Qnucl(nqpts) ! Added by FA
      double precision tol,rx,ry,rz,dist
      double precision nucCharge, nucCoords(3)
      character*16 tags(natoms)
      logical lSuccess
      external fn0,fn1,fn2,fn3,fn4
      external gridQpqPotential1
c
c     == preliminaries ==
      do igrid = 1,nqpts
        amat(igrid,1) = 0.d0            
        amat_coul(igrid,1) = 0.d0            
        if (ipol.gt.1) then
           amat(igrid,2) = 0.d0            
           amat_coul(igrid,2) = 0.d0            
        end if 
        amat_nucl(igrid) = 0.d0
        amat_Qnucl(igrid) = 0.d0 ! Added by FA
        closegridpts(igrid) = 0
      end do
c
c     == calculate the hartree potential on a supplied list of points ==
      tol = 1d-8
      call potential_list(ao_bas_han, g_dens(1), nqpts,
     &     qxyz, amat_coul(1,1), tol)
    
      if (ipol.gt.1) then
        call potential_list(ao_bas_han, g_dens(2), nqpts,
     &     qxyz, amat_coul(1,2), tol)
      end if
c
c     == calculate the total nuclear potential on the grid ==
      call gridNuclearPotential(geom,natoms,nqpts,qxyz,qwght,
     &                          closegridpts,amat_nucl)
      if ((zora_calc_type.eq.3).or.(zora_calc_type.eq.4)) then
c     == calculate Quadrupole potential on the grid ==
       call gridQpqPotential1(nqpts,qxyz,amat_Qnucl,
     &                        amat_Qnucl)     
      end if
c
c     == assemble zora correction ==
      if      (zora_calc_type.eq.0) then ! pure kinetic test
       call get_amat(fn0,nqpts,qwght,
     &               amat_coul,amat_nucl,amat_Qnucl,
     &               amat,closegridpts)
      else if (zora_calc_type.eq.1) then ! zora correction
       call get_amat(fn1,nqpts,qwght,
     &               amat_coul,amat_nucl,amat_Qnucl,
     &               amat,closegridpts)
      else if (zora_calc_type.eq.2) then ! zora energy scaling
       call get_amat(fn2,nqpts,qwght,
     &               amat_coul,amat_nucl,amat_Qnucl,
     &               amat,closegridpts)
      else if (zora_calc_type.eq.3) then ! zora EFG 
       call get_amat(fn3,nqpts,qwght,
     &               amat_coul,amat_nucl,amat_Qnucl,
     &               amat,closegridpts)
      else if (zora_calc_type.eq.4) then ! num  EFG 
       call get_amat(fn4,nqpts,qwght,
     &               amat_coul,amat_nucl,amat_Qnucl,
     &               amat,closegridpts)
      endif
      return
      end

      subroutine get_amat(fcn,nqpts,qwght,
     &                    amat_coul,amat_nucl,amat_Qnucl,
     &                    amat,closegridpts)
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "zora.fh"
#include "geom.fh"
      integer nqpts,igrid
      integer closegridpts(*)
      double precision qwght(*)         ![in]
      double precision amat(nqpts,ipol)
      double precision amat_coul(nqpts,ipol)
      double precision amat_nucl(nqpts),amat_Qnucl(nqpts)
      double precision valfn,totPot
      external fcn

       do igrid = 1,nqpts
        if (ipol.gt.1) then
         totPot = -amat_coul(igrid,1)-amat_coul(igrid,2) 
     &            + amat_nucl(igrid)
        else
         totPot = -amat_coul(igrid,1)+amat_nucl(igrid)
        end if
        valfn=0.0
        if (igrid.ne.closegridpts(igrid)) then   
          valfn=fcn(totPot,amat_Qnucl(igrid))
        endif
        amat(igrid,1)=valfn*qwght(igrid)
        if (ipol.gt.1) amat(igrid,2) = valfn*qwght(igrid)
       end do
      return
      end
       
      double precision function fn0(totPot,Qnucl)
       double precision toPot,Qnucl
       fn0=0.5d0
      return
      end

      double precision function fn1(totPot,Qnucl)
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "zora.fh"
#include "geom.fh"
       double precision totPot,Qnucl
       double precision clight_au2
       clight_au2 = clight_au*clight_au
       fn1=totPot/(4.d0*clight_au2-2.d0*totPot)
      return
      end

      double precision function fn2(totPot,Qnucl)
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "zora.fh"
#include "geom.fh"
       double precision totPot,Qnucl
       double precision clight_au2
       double precision denomFac
       clight_au2 = clight_au*clight_au
       denomFac = (2.d0*clight_au2-totPot)
       fn2=clight_au2/denomFac/denomFac
      return
      end

      double precision function fn3(totPot,Qnucl)
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "zora.fh"
#include "geom.fh"
       double precision totPot,Qnucl
       double precision clight_au2
       double precision denomFac
       clight_au2 = clight_au*clight_au
       denomFac = (2.d0*clight_au2-totPot)
       fn3=clight_au2/denomFac/denomFac*Qnucl  
      return
      end

      double precision function fn4(totPot,Qnucl)
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "zora.fh"
#include "geom.fh"
       double precision totPot,Qnucl
       fn4=Qnucl  
      return
      end
czora...
czora...Calculate the Quadrupole potential on the grid
czora...
      subroutine gridQpqPotential1(nqpts,qxyz,
     &                             amat_Qnucl,
     &                             closegridpts)
c    Purpose: Evaluates Q_pq(\vec{r},\vec{r}')
c    Input: 
c           zora_Qpq = 1  -> Evaluates Qxx
c                    = 2  -> Evaluates Qyy
c                    = 3  -> Evaluates Qzz
c                    = 4  -> Evaluates Qxy
c                    = 5  -> Evaluates Qxz
c                    = 6  -> Evaluates Qyz
c      ===> zora_Qpq, defined in zora.fh
c           nqpts    , number of grid points
c           qxyz     , grid points
c           xyz_EFGcoords, Quadrupole potential is evaluated
c                          in this point (\vec{r})
c      ===> xyz_EFGcoords, defined in zora.fh 
c    Output:
c          amat_Qnucl, Quadrupole potential ev. in the grid
c                      for integration purpose (\vec{r}')
c
      implicit none
c
#include "global.fh" 
#include "stdio.fh"
#include "zora.fh"

c
      integer igrid,nqpts
      double precision qxyz(3,nqpts)
      double precision rx,ry,rz,dist
      double precision amat_Qnucl(nqpts),dist5
      integer closegridpts(*)
      external fxx,fyy,fzz,fxy,fxz,fyz

c     == loop over the grid points ==
       if      (zora_Qpq.eq.1) then      ! Qxx
        call ev_Qpot(fxx,nqpts,qxyz,closegridpts,amat_Qnucl)
       else if (zora_Qpq.eq.2) then      ! Qyy
        call ev_Qpot(fyy,nqpts,qxyz,closegridpts,amat_Qnucl)
       else if (zora_Qpq.eq.3) then      ! Qzz
        call ev_Qpot(fzz,nqpts,qxyz,closegridpts,amat_Qnucl)
       else if (zora_Qpq.eq.4) then      ! Qxy
        call ev_Qpot(fxy,nqpts,qxyz,closegridpts,amat_Qnucl)
       else if (zora_Qpq.eq.5) then      ! Qxz
        call ev_Qpot(fxz,nqpts,qxyz,closegridpts,amat_Qnucl)
       else if (zora_Qpq.eq.6) then      ! Qyz
        call ev_Qpot(fyz,nqpts,qxyz,closegridpts,amat_Qnucl)
       endif
      return
      end 

      subroutine ev_Qpot(fcn,nqpts,qxyz,
     &                   closegridpts,amat_Qnucl)
#include "stdio.fh"
#include "zora.fh"
c
      integer igrid,nqpts
      double precision qxyz(3,nqpts)
      double precision rx,ry,rz,dist
      double precision amat_Qnucl(nqpts),dist5
      integer closegridpts(*)
      external fcn

       do igrid = 1,nqpts
       amat_Qnucl(igrid) = 0.d0
c     == distance from the grid points to given xyz_EFGcoords() == 
       rx = xyz_EFGcoords(1) - qxyz(1,igrid)
       ry = xyz_EFGcoords(2) - qxyz(2,igrid)
       rz = xyz_EFGcoords(3) - qxyz(3,igrid)
       dist = dsqrt(rx*rx + ry*ry + rz*rz)
        if (dist.gt.zoracutoff_EFG) then
          dist5=dist*dist*dist*dist*dist
          amat_Qnucl(igrid)=fcn(rx,ry,rz,dist5)
        else
          closegridpts(igrid) = igrid
        end if
       end do ! end-grid
      return
      end

      double precision function fxx(rx,ry,rz,dist5)
      double precision rx,ry,rz,dist5
      fxx=(2.d0*rx*rx-(ry*ry+rz*rz))/dist5
      return
      end

      double precision function fyy(rx,ry,rz,dist5)
      double precision rx,ry,rz,dist5
      fyy=(2.d0*ry*ry-(rx*rx+rz*rz))/dist5
      return
      end

      double precision function fzz(rx,ry,rz,dist5)
      double precision rx,ry,rz,dist5
      fzz=(2.d0*rz*rz-(rx*rx+ry*ry))/dist5
      return
      end

      double precision function fxy(rx,ry,rz,dist5)
      double precision rx,ry,rz,dist5
      fxy=(3.d0*rx*ry)/dist5
      return
      end

      double precision function fxz(rx,ry,rz,dist5)
      double precision rx,ry,rz,dist5
      fxz=(3.d0*rx*rz)/dist5
      return
      end

      double precision function fyz(rx,ry,rz,dist5)
      double precision rx,ry,rz,dist5
      fyz=(3.d0*ry*rz)/dist5
      return
      end
c =============== Fredy Aquino's routines ======= END
c
c $Id$

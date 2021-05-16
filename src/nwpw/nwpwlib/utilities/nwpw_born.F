
*     *******************************************
*     *                                         *
*     *            nwpw_born_radius             *
*     *                                         *
*     *******************************************
      real*8 function nwpw_born_radius(ii,nion,rion,vradii)
      implicit none
      integer ii,nion
      real*8 rion(3,nion),vradii(nion)

*     **** local variables ****
      integer jj
      real*8  dist,L,U,invL,invU,bornr1,delta
  
      bornr1 = 1.0d0/vradii(ii)
      do jj=1,nion
         dist = dsqrt( (rion(1,ii)-rion(1,jj))**2
     >               + (rion(2,ii)-rion(2,jj))**2
     >               + (rion(3,ii)-rion(3,jj))**2)
         if (dist.gt.0.1d0) then
            if ((dist+vradii(jj)).le.vradii(ii)) then
               L = 1.0d0
            else if (((dist-vradii(jj)).le.vradii(ii)).and.
     >               (vradii(ii).le.(dist+vradii(jj)))) then
               L = vradii(jj)
            else
               L = dist - vradii(jj)
            end if
            invL = 1.0d0/L

            if ((dist+vradii(jj)).le.vradii(ii)) then
               U = 1.0d0
            else
               U = dist + vradii(jj)
            end if
            invU = 1.0d0/U

            if ((U.gt.0.0d0).and.(L.gt.0.0d0)) then
               delta = -0.5d0*( (invL - invU)
     >                        + 0.25d0*dist*(invU**2 - invL**2)
     >                        + 0.50d0/dist*dlog(L/U)
     >                        + 0.25d0*vradii(jj)**2
     >                          /dist*(invL**2 - invU**2) )
               !if (delta.lt.0.0d0) bornr1 = bornr1 + delta
               bornr1 = bornr1 + delta
            end if
         end if
      end do
      nwpw_born_radius = 1.0d0/bornr1
      return
      end 


*     *******************************************
*     *                                         *
*     *             nwpw_born_init              *
*     *                                         *
*     *******************************************
      subroutine nwpw_born_init(rtdb0)
      implicit none
      integer rtdb0

#include "bafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** local variables ****
      integer  MASTER
      parameter(MASTER=0)
      logical  omaster,oprint,value
      integer  taskid
      integer  nion,ii,jj,iat
      real*8   rr

      integer mxelm
      parameter(mxelm = 102)
      real*8 vdwr(mxelm)

c      data vdwr /
c     1   0.80,0.49,0.00,0.00,0.00,1.65,1.55,1.50,1.50,0.00,
c     2   2.30,1.70,2.05,2.10,1.85,1.80,1.80,0.00,2.80,2.75,
c     3   0.00,0.00,1.20,0.00,0.00,0.00,2.70,0.00,0.00,0.00,
c     4   0.00,0.00,0.00,1.90,1.90,0.00,0.00,0.00,0.00,1.55,
c     5   0.00,1.64,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     6   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     7   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     8   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     9   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,
c     1   0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.65,
c     2   0.00,0.00,0.00/


c     radii for heavy elements: 1.17*1.9

      data (vdwr(iat),iat=1,mxelm)
     1 / 1.300,1.638,1.404,1.053,2.0475,2.00,
     2   1.830,1.720,1.720,1.8018,1.755,1.638,
     3   1.404,2.457,2.106,2.160,2.05,2.223,
     4   2.223,2.223,2.223,2.223,2.223,2.223,
     5   2.223,2.223,2.223,2.223,2.223,2.223,
     6   2.223,2.223,2.223,2.223,2.160,2.223,
     7   2.223,2.223,2.223,2.223,2.223,2.223,
     8   2.223,2.223,2.223,2.223,2.223,2.223,
     9   2.223,2.223,2.223,2.223,2.320,2.223,
     1   2.223,2.223,2.223,2.223,2.223,2.223,
     2   2.223,2.223,2.223,2.223,2.223,2.223,
     3   2.223,2.223,2.223,2.223,2.223,2.223,
     4   2.223,2.223,2.223,2.223,2.223,2.223,
     5   2.223,2.223,2.223,2.223,2.223,2.223,
     6   2.223,2.223,2.223,2.223,2.223,2.223,
     7   2.223,2.223,2.223,2.223,2.223,2.223,
     7   2.223,2.223,2.223,2.223,2.223,2.223/



*     **** nwpw_born common blocks ****
      logical born_on
      integer bradii(2),vradii(2),qborn(2),rtdb
      real*8  dielec
      common /nwpw_born_blk/ bradii,vradii,qborn,dielec,rtdb,born_on


*     ***** external functions *****
      logical  control_print
      external control_print
      integer  Parallel_threadid
      external Parallel_threadid
      integer  ion_nion,ion_rion_ptr,ion_katm
      external ion_nion,ion_rion_ptr,ion_katm
      real*8   ion_q,nwpw_born_radius
      external ion_q,nwpw_born_radius
      character*4 ion_atom
      external    ion_atom

c
      call Parallel_taskid(taskid)
      omaster =  (taskid.eq.MASTER).and.(Parallel_threadid().eq.MASTER)
      oprint = omaster.and.control_print(print_debug)


      if (.not.btdb_get(rtdb,'nwpw:born',mt_log,1,born_on)) 
     >   born_on = .false.

      if (born_on) then
         if (.not.btdb_get(rtdb,'nwpw:born_dielec',mt_dbl,1,dielec)) 
     >      dielec = 78.4d0

         nion = ion_nion()
         value = BA_alloc_get(mt_dbl,nion,'bradii',bradii(2),bradii(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,nion,'vradii',vradii(2),vradii(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,nion,'qborn',qborn(2),qborn(1))
         if (.not. value) 
     >     call errquit('nwpw_born_init:out of heap memory',0,MA_ERR)

         if (.not.btdb_get(rtdb,'nwpw:born_vradii',mt_dbl,
     >                            nion,dbl_mb(vradii(1)))) then
            do ii=1,nion
               dbl_mb(vradii(1)+ii-1) = vdwr(nint(ion_q(ii)))/0.529177d0
            end do
         end if
         do ii=1,nion
            rr =  dbl_mb(vradii(1)+ii-1)
            if (rr.lt.1.0d-3) rr = 1.17d0*1.90d0/0.529177d0
            dbl_mb(vradii(1)+ii-1) = rr
         end do

         if (.not.btdb_get(rtdb,'nwpw:born_bradii',mt_dbl,
     >                            nion,dbl_mb(bradii(1)))) then
            do ii=1,nion
               dbl_mb(bradii(1)+ii-1) 
     >            = nwpw_born_radius(ii,nion,
     >                               dbl_mb(ion_rion_ptr()),
     >                               dbl_mb(vradii(1)))
            end do
         end if

         if (oprint) then
            write(luout,*)
            write(luout,*) "initializing generalized Born radii"
            write(luout,*) "-----------------------------------"
            write(luout,*) " - G.D. Hawkins, C.R. Cramer, ",
     >                     "D.G. Truhlar (1995) ",
     >                     "Pairwise solute descreening of solute ",
     >                     "charge from a dielectric medium, ",
     >                     "Chem. Phys. Lett., vol. 246, pages 122-129."
            write(luout,*)
            do ii=1,nion
               write(luout,101) ion_atom(ion_katm(ii)),ii,
     >                          dbl_mb(bradii(1)+ii-1)*0.529177d0,
     >                          dbl_mb(vradii(1)+ii-1)*0.529177d0
            end do
            write(luout,*)
         end if

      endif

  101 format(1x,"Born Radius ",A4," - a(",I1,") = ",F6.3,
     >       " Angstroms (1.17*vdw radius=",F6.3," Angstroms)")
         
      return
      end


*     *******************************************
*     *                                         *
*     *             nwpw_born_end               *
*     *                                         *
*     *******************************************
      subroutine nwpw_born_end()
      implicit none
      
#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"

*     **** nwpw_born common blocks ****
      logical born_on
      integer bradii(2),vradii(2),qborn(2),rtdb
      real*8  dielec
      common /nwpw_born_blk/ bradii,vradii,qborn,dielec,rtdb,born_on

*     **** local variables ****
      logical value
      integer nion

*     **** external functions ****
      integer  ion_nion
      external ion_nion

c     *** Clean up born stuff, if necessary ***
      if (born_on) then
         nion = ion_nion()
         value = btdb_put(rtdb,'nwpw:born_bradii',mt_dbl,
     >                           nion,dbl_mb(bradii(1))) 
         value = value.and.
     >           btdb_put(rtdb,'nwpw:born_vradii',mt_dbl,
     >                           nion,dbl_mb(vradii(1))) 
         if (.not.value) 
     >      call errquit('nwpw_born_end:writing radii',0,MA_ERR)

         born_on = .false.
         value = BA_free_heap(qborn(2))
         value = value.and.BA_free_heap(vradii(2))
         value = value.and.BA_free_heap(bradii(2))
         if (.not.value) 
     >      call errquit('nwpw_born_end:heap memory',0,MA_ERR)
      endif
      return
      end

*     *******************************************
*     *                                         *
*     *              nwpw_born_on               *
*     *                                         *
*     *******************************************
      logical function nwpw_born_on()
      implicit none

*     **** nwpw_born common blocks ****
      logical born_on
      integer bradii(2),vradii(2),qborn(2),rtdb
      real*8  dielec
      common /nwpw_born_blk/ bradii,vradii,qborn,dielec,rtdb,born_on

      nwpw_born_on = born_on
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_born_print              *
*     *                                         *
*     *******************************************
      subroutine nwpw_born_print(nga,nion_qm,qgaus)
      implicit none
      integer nga,nion_qm
      real*8 qgaus(nga,nion_qm)

#include "bafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"

*     **** nwpw_born common blocks ****
      logical born_on
      integer bradii(2),vradii(2),qborn(2),rtdb
      real*8  dielec
      common /nwpw_born_blk/ bradii,vradii,qborn,dielec,rtdb,born_on

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical oprint,omaster
      integer ii,nion,j,indx
      real*8  Gsolv,sum

*     **** external functions ****
      character*4 ion_atom
      external    ion_atom
      integer     ion_nion,ion_katm,ion_rion_ptr
      external    ion_nion,ion_katm,ion_rion_ptr
      integer     ion_katm_qm
      external    ion_katm_qm
      real*8      nwpw_born_energy,psp_zv
      external    nwpw_born_energy,psp_zv
      integer     Parallel_threadid
      external    Parallel_threadid
      logical     control_print
      external    control_print

      call Parallel_taskid(taskid)
      omaster =  (taskid.eq.MASTER).and.(Parallel_threadid().eq.MASTER)
      oprint = omaster.and.control_print(print_medium)

      if (born_on) then
         nion  = ion_nion()
         do ii=1,nion_qm
            sum = 0.0d0
            do j=1,nga
              sum = sum + qgaus(j,ii)
            end do
            dbl_mb(qborn(1)+ii-1) = psp_zv(ion_katm_qm(ii))-sum
         end do
         do ii=nion_qm+1,nion
            dbl_mb(qborn(1)+ii-1) = psp_zv(ion_katm(ii))
         end do

         Gsolv = nwpw_born_energy(nion,
     >                            dbl_mb(ion_rion_ptr()),
     >                            dbl_mb(bradii(1)),
     >                            dbl_mb(qborn(1)),dielec)
         if (oprint) then
            write(luout,*)
            write(luout,*) "Generalized Born Solvation"
            write(luout,*) "--------------------------"
            write(luout,*) " - Radii defined by G.D. Hawkins,",
     >                     " C.R. Cramer, D.G. Truhlar (1995) ",
     >                     "Pairwise" 
            write(luout,*) "   solute descreening of solute ",
     >                     "charge from a dielectric medium, ",
     >                     "Chem. Phys. Lett.,"
            write(luout,*) "   vol. 246, pages 122-129."
            write(luout,*)
            write(luout,'(2x,A,F8.2)') "Dielectric constant -eps- = ",
     >                              dielec
            write(luout,*)
            do ii=1,nion
               write(luout,101) ion_atom(ion_katm(ii)),ii,
     >                          dbl_mb(bradii(1)+ii-1)*0.529177d0,
     >                          dbl_mb(vradii(1)+ii-1)*0.529177d0,
     >                          dbl_mb(qborn(1)+ii-1)
            end do
            write(luout,*)
            write(luout,102) Gsolv,Gsolv*27.2116d0*23.06d0
            write(luout,*)
         end if

      end if

  101 format(2x,"Generalized Born radius ",A4," - a(",I1,") = ",F8.3,
     >       " Angstroms (1.17*vdw radius=",F8.3,
     >       ") - atom charge = ",F8.3)
  102 format(2x,"Gsolvation(Born) = ",E14.6," (",F8.3," kcal/mol)")
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_born_energy             *
*     *                                         *
*     *******************************************
      real*8 function nwpw_born_energy(nion,rion,bradii,q,dielec)
      implicit none
      integer nion
      real*8  rion(3,nion),bradii(nion),q(nion)
      real*8  dielec

*     **** local variables ****
      integer MASTER,taskid,np
      parameter (MASTER=0)
      integer ii,jj,itask
      real*8  Gsolv,screen,C,f,dist2

      call Parallel_np(np)
      call Parallel_taskid(taskid)


      Gsolv = 0.0d0
      screen = (1.0d0 - 1.0d0/dielec)
      
      itask = 0
!$OMP MASTER
      do ii=1,nion
         do jj=1,nion
            if (itask.eq.taskid) then
               dist2 = ((rion(1,ii)-rion(1,jj))**2
     >                + (rion(2,ii)-rion(2,jj))**2
     >                + (rion(3,ii)-rion(3,jj))**2)
               C = dexp(-0.25d0*dist2/(bradii(ii)*bradii(jj)))
               f = dsqrt(dist2 + bradii(ii)*bradii(jj)*C)
               Gsolv = Gsolv - 0.5d0*screen*q(ii)*q(jj)/f
            end if
            itask = mod(itask+1,np)
         end do
      end do
!$OMP END MASTER
      call Parallel_SumAll(Gsolv)

      nwpw_born_energy = Gsolv
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_born_dVdq               *
*     *                                         *
*     *******************************************
      subroutine nwpw_born_dVdq(nion,rion,bradii,q,dielec,u)
      implicit none
      integer nion
      real*8  rion(3,nion),bradii(nion),q(nion)
      real*8  dielec
      real*8  u(nion)

*     **** local variables ****
      integer MASTER,taskid,np
      parameter (MASTER=0)
      integer ii,jj,itask
      real*8  Gsolv,screen,C,f,dist2

      call Parallel_np(np)
      call Parallel_taskid(taskid)


      !Gsolv = 0.0d0
      screen = (1.0d0 - 1.0d0/dielec)

      call Parallel_shared_vector_zero(.true.,nion,u)

      itask = 0
!$OMP MASTER
      do ii=1,nion
         do jj=1,nion
            if (itask.eq.taskid) then
               dist2 = ((rion(1,ii)-rion(1,jj))**2
     >                + (rion(2,ii)-rion(2,jj))**2
     >                + (rion(3,ii)-rion(3,jj))**2)
               C = dexp(-0.25d0*dist2/(bradii(ii)*bradii(jj)))
               f = dsqrt(dist2 + bradii(ii)*bradii(jj)*C)
               u(ii) = u(ii) - 0.5d0*screen*q(jj)/f
               u(jj) = u(jj) - 0.5d0*screen*q(jj)/f
               !Gsolv = Gsolv - 0.5d0*screen*q(ii)*q(jj)/f
            end if
            itask = mod(itask+1,np)
         end do
      end do
!$OMP END MASTER
      call Parallel_Vector_SumAll(nion,u)

      return
      end

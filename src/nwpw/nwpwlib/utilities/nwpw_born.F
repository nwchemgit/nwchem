
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
            L=0.0d0
            if ((dist+vradii(jj)).le.vradii(ii)) then
               L = 1.0d0
            else if (((dist-vradii(jj)).le.vradii(ii)).and.
     >               (vradii(ii).le.(dist+vradii(jj)))) then
               L = vradii(jj)
            else if (vradii(ii).le.(dist+vradii(jj))) then
               L = dist - vradii(jj)
            end if

            U=0.0d0
            if ((dist+vradii(jj)).le.vradii(ii)) then
               U = 1.0d0
            else if (vradii(ii).lt.(dist+vradii(jj))) then
               U = dist + vradii(jj)
            end if

            if ((U.gt.0.0d0).and.(L.gt.0.0d0)) then
               invL = 1.0d0/L
               invU = 1.0d0/U
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
      if (bornr1.lt.1.0d-6) bornr1 = 1.0d-6
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
      logical born_on,born_relax
      integer uborn(2),qborn(2)
      integer bradii(2),vradii(2),rtdb
      real*8  dielec,rcut,gg0
      common /nwpw_born_blk/ uborn,qborn,bradii,vradii,dielec,rcut,gg0,
     >                       rtdb,born_on,born_relax


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
      oprint = omaster.and.control_print(print_medium)


      if (.not.btdb_get(rtdb,'nwpw:born_relax',mt_log,1,born_relax)) 
     >   born_relax = .false.

      if (.not.btdb_get(rtdb,'nwpw:born',mt_log,1,born_on)) 
     >   born_on = .false.


      if (born_on) then
         if (.not.btdb_get(rtdb,'nwpw:born_dielec',mt_dbl,1,dielec)) 
     >      dielec = 78.4d0

         if (.not.btdb_get(rtdb,'nwpw:born_rcut',mt_dbl,1,rcut)) 
     >      rcut = 0.1d0

         if (.not.btdb_get(rtdb,'nwpw:born_gg0',mt_dbl,1,gg0)) 
     >      gg0 = 1.00d0

         nion  = ion_nion()
         value = BA_alloc_get(mt_dbl,nion,'bradii',bradii(2),bradii(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,nion,'vradii',vradii(2),vradii(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,nion,'qborn',qborn(2),qborn(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,nion,'uborn',uborn(2),uborn(1))
         if (.not.value) 
     >     call errquit('nwpw_born_init:out of heap memory',0,MA_ERR)

         if (.not.btdb_get(rtdb,'nwpw:born_vradii',mt_dbl,
     >                            nion,dbl_mb(vradii(1)))) then
            do ii=1,nion
               dbl_mb(vradii(1)+ii-1) = vdwr(nint(ion_q(ii)))/0.529177d0
            end do
         end if
         do ii=1,nion
            rr = dbl_mb(vradii(1)+ii-1)
            if (rr.lt.1.0d-3) rr = 1.17d0*1.90d0/0.529177d0
            dbl_mb(vradii(1)+ii-1) = rr
         end do

         if (.not.btdb_get(rtdb,'nwpw:born_bradii',mt_dbl,
     >                            nion,dbl_mb(bradii(1)))) then
            do ii=1,nion
               dbl_mb(bradii(1)+ii-1) = nwpw_born_radius(ii,nion,
     >                                      dbl_mb(ion_rion_ptr()),
     >                                      dbl_mb(vradii(1)))
            end do
         end if

         if (oprint) then
            write(luout,*)
            write(luout,*) "extended Born solvation model:"
            write(luout,'(5x,3A)') "G.D. Hawkins, C.R. Cramer, ",
     >                     "D.G. Truhlar (1995) ",
     >                     "Pairwise solute descreening of solute"
            write(luout,'(5x,2A)') "  charge from a dielectric medium,",
     >                    " Chem. Phys. Lett., vol. 246, pages 122-129."
            write(luout,*)
            write(luout,'(5x,A,F11.6)') "dielectric constant = ",dielec
            write(luout,'(5x,A,F11.6)') "dielectric rcut     = ",rcut
            write(luout,'(5x,A,F11.6)') "dielectric gg0      = ",gg0
            if (born_relax) then
               write(luout,'(5x,A)') "self-consistent solvation"
            else
               write(luout,'(5x,A)') "perturbative solvation"
            end if
            write(luout,'(5x,A)') "generalized Born radii:"
            do ii=1,nion
               write(luout,101) ion_atom(ion_katm(ii)),ii,
     >                          dbl_mb(bradii(1)+ii-1)*0.529177d0,
     >                          dbl_mb(vradii(1)+ii-1)*0.529177d0
            end do
            write(luout,*)
         end if

      endif

  101 format(7x,"- Born Radius ",A4," - a(",I4,") = ",F10.3,
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
      logical born_on,born_relax
      integer uborn(2),qborn(2)
      integer bradii(2),vradii(2),rtdb
      real*8  dielec,rcut,gg0
      common /nwpw_born_blk/ uborn,qborn,bradii,vradii,dielec,rcut,gg0,
     >                       rtdb,born_on,born_relax

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
         value = BA_free_heap(uborn(2))
         value = value.and.BA_free_heap(qborn(2))
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
      logical born_on,born_relax
      integer uborn(2),qborn(2)
      integer bradii(2),vradii(2),rtdb
      real*8  dielec,rcut,gg0
      common /nwpw_born_blk/ uborn,qborn,bradii,vradii,dielec,rcut,gg0,
     >                       rtdb,born_on,born_relax

      nwpw_born_on = born_on
      return
      end

*     *******************************************
*     *                                         *
*     *              nwpw_born_relax            *
*     *                                         *
*     *******************************************
      logical function nwpw_born_relax()
      implicit none

*     **** nwpw_born common blocks ****
      logical born_on,born_relax
      integer uborn(2),qborn(2)
      integer bradii(2),vradii(2),rtdb
      real*8  dielec,rcut,gg0
      common /nwpw_born_blk/ uborn,qborn,bradii,vradii,dielec,rcut,gg0,
     >                       rtdb,born_on,born_relax

      nwpw_born_relax = born_relax
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_born_Qprint             *
*     *                                         *
*     *******************************************
      subroutine nwpw_born_Qprint(nga,nion_qm,qgaus)
      implicit none
      integer nga,nion_qm
      real*8 qgaus(nga,nion_qm)

#include "bafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"


*     **** nwpw_born common blocks ****
      logical born_on,born_relax
      integer uborn(2),qborn(2)
      integer bradii(2),vradii(2),rtdb
      real*8  dielec,rcut,gg0
      common /nwpw_born_blk/ uborn,qborn,bradii,vradii,dielec,rcut,gg0,
     >                       rtdb,born_on,born_relax

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
      real*8      nwpw_born_energy0,psp_zv
      external    nwpw_born_energy0,psp_zv
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

         Gsolv = nwpw_born_energy0(nion,
     >                            dbl_mb(ion_rion_ptr()),
     >                            dbl_mb(bradii(1)),
     >                            dbl_mb(qborn(1)),dielec,rcut,gg0)
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
     >                                 dielec
            write(luout,'(2x,A,F8.2)') "Dielectric rcut           = ",
     >                                 rcut
            write(luout,'(2x,A,F8.2)') "Dielectric gg0            = ",
     >                                 gg0
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
*     *            nwpw_born_screen             *
*     *                                         *
*     *******************************************
      real*8 function nwpw_born_screen()
      implicit none

*     **** nwpw_born common blocks ****
      logical born_on,born_relax
      integer uborn(2),qborn(2)
      integer bradii(2),vradii(2),rtdb
      real*8  dielec,rcut,gg0
      common /nwpw_born_blk/ uborn,qborn,bradii,vradii,dielec,rcut,gg0,
     >                       rtdb,born_on,born_relax

      nwpw_born_screen = (1.0d0 - 1.0d0/dielec)
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_born_energy             *
*     *                                         *
*     *******************************************
      real*8 function nwpw_born_energy()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_born common blocks ****
      logical born_on,born_relax
      integer uborn(2),qborn(2)
      integer bradii(2),vradii(2),rtdb
      real*8  dielec,rcut,gg0
      common /nwpw_born_blk/ uborn,qborn,bradii,vradii,dielec,rcut,gg0,
     >                       rtdb,born_on,born_relax

*     **** external functions ****
      integer  ion_rion_ptr,ion_nion_qm
      external ion_rion_ptr,ion_nion_qm
      real*8   nwpw_born_energy0
      external nwpw_born_energy0

      nwpw_born_energy = nwpw_born_energy0(ion_nion_qm(),
     >                            dbl_mb(ion_rion_ptr()),
     >                            dbl_mb(bradii(1)),
     >                            dbl_mb(qborn(1)),dielec,rcut,gg0)
      return
      end

      real*8 function nwpw_born_energy0(nion,rion,bradii,q,
     >                                  dielec,rcut,gg0)
      implicit none
      integer nion
      real*8  rion(3,nion),bradii(nion),q(nion)
      real*8  dielec,rcut,gg0

*     **** local variables ****
      integer MASTER,taskid,np
      parameter (MASTER=0)
      integer ii,jj,itask
      real*8  Gsolv,screen,C,f,dist2,gg

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
               if (ii.eq.jj) then
                  gg = gg0
               else
                  gg = erf(dsqrt(dist2)/rcut)
               end if
               Gsolv = Gsolv - 0.5d0*screen*q(ii)*q(jj)*gg/f

            end if
            itask = mod(itask+1,np)
         end do
      end do
!$OMP END MASTER
      call Parallel_SumAll(Gsolv)

c      write(*,*) "EE dielec=",dielec
c      do ii=1,nion
c         write(*,*) "EE ii,q,Gsolv=",ii,q(ii),Gsolv
c      end do

      nwpw_born_energy0 = Gsolv
      return
      end



*     *******************************************
*     *                                         *
*     *            nwpw_born_fion               *
*     *                                         *
*     *******************************************
      subroutine nwpw_born_fion(fion)
      implicit none
      real*8 fion(3,*) 

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_born common blocks ****
      logical born_on,born_relax
      integer uborn(2),qborn(2)
      integer bradii(2),vradii(2),rtdb
      real*8  dielec,rcut,gg0
      common /nwpw_born_blk/ uborn,qborn,bradii,vradii,dielec,rcut,gg0,
     >                       rtdb,born_on,born_relax

*     **** external functions ****
      integer  ion_rion_ptr,ion_nion_qm
      external ion_rion_ptr,ion_nion_qm

      integer nion,ftmp(2)

*     *** push stack ftmp ***
      nion = ion_nion_qm()
      if (.not.BA_push_get(mt_dbl,3*nion,'ftmp',ftmp(2),ftmp(1)))
     >   call errquit('nwpw_born_fion:push stack',0,MA_ERR)

      call nwpw_born_fion0(nion,
     >                     dbl_mb(ion_rion_ptr()),
     >                     dbl_mb(bradii(1)),
     >                     dbl_mb(qborn(1)),dielec,rcut,dbl_mb(ftmp(1)))

      call DAXPY_OMP(3*nion,1.0d0,dbl_mb(ftmp(1)),1,fion,1)

*     *** pop stack ftmp ***
      if (.not.BA_pop_stack(ftmp(2)))
     >   call errquit('nwpw_born_fion:pop stack',0,MA_ERR)

      return
      end

*     **************************************
*     *                                    *
*     *          nwpw_born_gsolv           *
*     *                                    *
*     **************************************
*
*     Calculates gsolv = -0.5*screen*qi*qj*gg/f
*   
*     where f = sqrt(xx + bi*bj*C)
*           gg = erf(sqrt(xx)/rcut)
*           C = exp(-xx/(4*bi*bj))  
*

      real*8 function nwpw_born_gsolv(screen,qi,qj,bi,bj,xx,rcut)
      implicit none
      real*8 screen,qi,qj,bi,bj,xx,rcut
      real*8 C,f,gg

      gg = erf(dsqrt(xx)/rcut) 
      C = dexp(-0.25d0*xx/(bi*bj))
      f = dsqrt(xx + bi*bj*C)
      nwpw_born_gsolv  = -0.5d0*screen*qi*qj*gg/f
      return
      end

*     **************************************
*     *                                    *
*     *          nwpw_born_dgsolv          *
*     *                                    *
*     **************************************
*
*     Calculates dgsolv/dxx
*
*     dC/dxx = -1/(4*bi*bj)*exp(-xx/(4*bi*bj)) = -1/(4*bi*bj)*C
*     df/dxx = 0.5*1/f * (dxx/dxx + bi*bj*dC/dxx) =  0.5*1/f*(1-0.25*C)
*     dgg/dxx = 1/sqrt(xx*pi) * exp(-xx/rcut**2)/rcut
*     dsqrt(xx)/dxx = 0.5/xx
*
      real*8 function nwpw_born_dgsolv(screen,qi,qj,bi,bj,xx,rcut)
      implicit none
      real*8 screen,qi,qj,bi,bj,xx,rcut
      real*8 C,f,gsolv,gg,dgg


      gg  = erf(dsqrt(xx)/rcut)
      dgg = (1.0d0/dsqrt(xx*4.0d0*datan(1.0d0)))*exp(-xx/rcut**2)/rcut
      C = dexp(-0.25d0*xx/(bi*bj))
      f = dsqrt(xx + bi*bj*C)
      gsolv  = -0.5d0*screen*qi*qj*gg/f

      nwpw_born_dgsolv = -0.5d0*gsolv*(1.0-0.25d0*C)/f**2
     >                  - 0.5d0*screen*qi*qj*dgg/f
      return
      end

      subroutine nwpw_born_fion0(nion,rion,bradii,q,dielec,rcut,fion)
      implicit none
      integer nion
      real*8  rion(3,nion),bradii(nion),q(nion)
      real*8  dielec,rcut
      real*8  fion(3,nion)

*     **** local variables ****
      integer MASTER,taskid,np
      parameter (MASTER=0)
      integer ii,jj,itask
      real*8  Gsolv,dGsolv,screen,C,f,dist2,dx,dy,dz

      real*8   nwpw_born_dgsolv
      external nwpw_born_dgsolv

      call Parallel_np(np)
      call Parallel_taskid(taskid)
      
      screen = (1.0d0 - 1.0d0/dielec)
      call Parallel_shared_vector_zero(.true.,3*nion,fion)

      itask = 0
!$OMP MASTER
      do ii=1,nion
         do jj=1,nion
            if (itask.eq.taskid) then
               if (ii.ne.jj) then
                  dx = rion(1,ii)-rion(1,jj)
                  dy = rion(2,ii)-rion(2,jj)
                  dz = rion(3,ii)-rion(3,jj)
                  dist2 = dx*dx + dy*dy + dz*dz

                  dGsolv = nwpw_born_dgsolv(screen,q(ii),q(jj),
     >                                      bradii(ii),bradii(jj),
     >                                      dist2,rcut)
                  fion(1,ii) = fion(1,ii) - 2.0d0*dGsolv*dx
                  fion(2,ii) = fion(2,ii) - 2.0d0*dGsolv*dy
                  fion(3,ii) = fion(3,ii) - 2.0d0*dGsolv*dz

                  fion(1,jj) = fion(1,jj) + 2.0d0*dGsolv*dx
                  fion(2,jj) = fion(2,jj) + 2.0d0*dGsolv*dy
                  fion(3,jj) = fion(3,jj) + 2.0d0*dGsolv*dz
               end if
            end if
            itask = mod(itask+1,np)
         end do
      end do
!$OMP END MASTER
      call Parallel_Vector_SumAll(3*nion,fion)

      return
      end 






*     *******************************************
*     *                                         *
*     *            nwpw_born_dVdq               *
*     *                                         *
*     *******************************************
      subroutine  nwpw_born_dVdq(nion,q,u)
      implicit none
      integer nion
      real*8 q(nion),u(nion)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_born common blocks ****
      logical born_on,born_relax
      integer uborn(2),qborn(2)
      integer bradii(2),vradii(2),rtdb
      real*8  dielec,rcut,gg0
      common /nwpw_born_blk/ uborn,qborn,bradii,vradii,dielec,rcut,gg0,
     >                       rtdb,born_on,born_relax

      integer  ion_rion_ptr
      external ion_rion_ptr

      call nwpw_born_dVdq0(nion,
     >                     dbl_mb(ion_rion_ptr()),
     >                     dbl_mb(bradii(1)),
     >                     q,dielec,rcut,gg0,u)
      return
      end

      subroutine nwpw_born_dVdq0(nion,rion,bradii,q,dielec,rcut,gg0,u)
      implicit none
      integer nion
      real*8  rion(3,nion),bradii(nion),q(nion)
      real*8  dielec,rcut,gg0
      real*8  u(nion)

*     **** local variables ****
      integer MASTER,taskid,np
      parameter (MASTER=0)
      integer ii,jj,itask
      real*8  Gsolv,screen,C,f,dist2,gg

      call Parallel_np(np)
      call Parallel_taskid(taskid)


      Gsolv = 0.0d0
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
               if (ii.eq.jj) then
                  gg = gg0
               else
                  gg = erf(dsqrt(dist2)/rcut)
               end if
               u(ii) = u(ii) + 0.5d0*screen*q(jj)*gg/f
               u(jj) = u(jj) + 0.5d0*screen*q(ii)*gg/f
               Gsolv = Gsolv - 0.5d0*screen*q(ii)*q(jj)*gg/f
            end if
            itask = mod(itask+1,np)
         end do
      end do
!$OMP END MASTER
      call Parallel_Vector_SumAll(nion,u)

c      write(*,*) "             ","         dielec,Gsolv=",dielec,Gsolv
c      do ii=1,nion
c         write(*,*) "             ","         ii,q,u=",ii,q(ii),u(ii)
c      end do

      return
      end

*     *******************************************
*     *                                         *
*     *              nwpw_born_u_ptr            *
*     *                                         *
*     *******************************************
      integer function nwpw_born_u_ptr()
      implicit none

*     **** nwpw_born common blocks ****
      logical born_on,born_relax
      integer uborn(2),qborn(2)
      integer bradii(2),vradii(2),rtdb
      real*8  dielec,rcut,gg0
      common /nwpw_born_blk/ uborn,qborn,bradii,vradii,dielec,rcut,gg0,
     >                       rtdb,born_on,born_relax

      nwpw_born_u_ptr = uborn(1)
      return
      end


*     *******************************************
*     *                                         *
*     *              nwpw_born_q_ptr            *
*     *                                         *
*     *******************************************
      integer function nwpw_born_q_ptr()
      implicit none

*     **** nwpw_born common blocks ****
      logical born_on,born_relax
      integer uborn(2),qborn(2)
      integer bradii(2),vradii(2),rtdb
      real*8  dielec,rcut,gg0
      common /nwpw_born_blk/ uborn,qborn,bradii,vradii,dielec,rcut,gg0,
     >                       rtdb,born_on,born_relax

      nwpw_born_q_ptr = qborn(1)
      return
      end



*     *******************************************
*     *                                         *
*     *            nwpw_born_print              *
*     *                                         *
*     *******************************************
      subroutine nwpw_born_print(iunit,egas,esol)
      implicit none
      integer iunit
      real*8 egas,esol

         write(iunit,911)
         if (dabs(egas).gt.0.d0) then
            write(iunit,912) egas
            write(iunit,913) esol
            write(iunit,914) (egas-esol),
     >                       (egas-esol)*27.2116d0*23.06d0
         else
            write(iunit,915)
            write(iunit,913) esol
         endif

c
c    Standard cosmo output
c
 911  format(
     $/5x,' extended Born solvation results',/,
     $5x, ' -------------------------------')
 912  format(5x,' gas phase energy                 : ',e19.10)
 913  format(5x,' sol phase energy                 : ',e19.10)
 914  format(5x,' (electrostatic) solvation energy : ',e19.10,
     $     ' (',f8.3,' kcal/mol)'                       )
 915  format('      skipped: no gas phase energy')
c
c     Custom cosmo-smd output
c
 819  format(
     $/5x,' extended Born-SMD solvation results',/,
     $5x,' -----------------------------------',/
     $5x,' Reference for the SMD model:',/,
     $5x,' Marenich, A. V.; Cramer, C. J.; Truhlar, D. G.',
     $' J. Phys. Chem. B 2009, 113, 6378',/)
c 820  format(1x,
c     $'            internal energy in gas <Psi(g)|H|Psi(g)> = ',
c     $ f22.12)
c 821  format(1x,
c     $'        internal energy in solvent <Psi(s)|H|Psi(s)> = ',
c     $ f22.12)
c 822  format(1x,
c     $'                              delta internal energy  = ',
c     $ f22.12,' (',f8.2,' kcal/mol)')
c 823  format(1x,
c     $'  total free energy in solvent <Psi(s)|H+V/2|Psi(s)> = ',
c     $ f22.12)
c 824  format(1x,
c     $'polarization energy contribution <Psi(s)|V/2|Psi(s)> = ',
c     $ f22.12,' (',f8.2,' kcal/mol)')

 825  format(5x,
     $' total free energy in solvent including G(SMD-CDS) : ',
     $ e19.10)
 826  format(5x,
     $' G(SMD-CDS) energy contribution                    : ',
     $ e19.10,
     $     ' (',f9.3,' kcal/mol)')
 827  format(5x,
     $' 1 M fixed-concentration free energy of solvation  : ',
     $ e19.10,' (',f9.3,' kcal/mol)'                       )


      return
      end

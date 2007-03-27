c
c $Id: pspw_LJ.F,v 1.12 2007-03-27 02:02:50 bylaska Exp $
c

*     *************************
*     *                       *
*     *     pspw_LJ_init	  *
*     *                       *
*     *************************
      subroutine pspw_LJ_init(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "pspw_LJ.fh"
#include "errquit.fh"

*     **** local variables
      logical value
      integer taskid
      integer MASTER
      parameter(MASTER=0)

      integer i,j
      integer nkatm
      real*8  sigma_epsilon(2)
      character*30 rtdbname

*     **** external functions ****
      character*2 ion_atom
      integer     ion_nkatm
      integer     ion_katm
      external    ion_atom
      external    ion_nkatm
      external    ion_katm
      

      nkatm = ion_nkatm()

*     **** allocate LJ parameters ****
      value = rtdb_parallel(.true.)
      value = MA_alloc_get(mt_dbl,nkatm,
     >                 'sigma_all',sigma_all(2),sigma_all(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,nkatm,
     >                 'epsilon_ion',epsilon_all(2),epsilon_all(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)


      
*     **** read in LJ parameters ***
      do i=1,nkatm
         rtdbname = 'pspw_LJ_param_ion:'//ion_atom(i)
         j = index(rtdbname,' ')-1
         value = rtdb_get(rtdb,rtdbname(1:j),
     >                  mt_dbl,2,sigma_epsilon)
         if (.not. value) call dcopy(2,0.0d0,0,sigma_epsilon,1)

        dbl_mb(sigma_all(1)+i-1)  =sigma_epsilon(1)/(0.529177d0)
        dbl_mb(epsilon_all(1)+i-1)=sigma_epsilon(2)/(23.06d0*27.2116d0)
      end do
      
      
      call Parallel_taskid(taskid)
      if (taskid.eq.MASTER) then
          write(*,*)
          write(*,*) 'Ion LJ Parameters:'
          do i=1,nkatm
            write(*,'(A3,4x,A12,E10.3,A9,A12,E10.3,A9)') 
     >         ion_atom(i),
     >         'sigma:',  
     >         dbl_mb(sigma_all(1)+i-1)*(0.529177d0), 
     >         'Angstrom',
     >         'epsilon:',
     >         dbl_mb(epsilon_all(1)+ i-1)*(23.06d0*27.2116d0),
     >         'kcal/mol'
          end do
         write(*,*)
      end if

      return
      end

*     *************************
*     *                       *
*     *     pspw_LJ_end	      *
*     *                       *
*     *************************
      subroutine pspw_LJ_end()
      implicit none

#include "mafdecls.fh"
#include "pspw_LJ.fh"
#include "errquit.fh"

      logical value

      value =           MA_free_heap(sigma_all(2))
      value = value.and.MA_free_heap(epsilon_all(2))
      if (.not.value) call errquit('pspw_LJ_end: error MA_free_heap',
     >                             0,MA_ERR)
      return
      end


*     ***************************
*     *				*
*     *    pspw_LJ_param_Ion	*
*     *				*
*     ***************************
      subroutine pspw_LJ_param_Ion(i,s,e)
      implicit none
      integer i
      real*8  s,e

#include "mafdecls.fh"
#include "pspw_LJ.fh"

*     **** local variables ****
      integer ka

*     **** external functions ****
      integer  ion_katm
      external ion_katm

      ka = ion_katm(i)
      s  = dbl_mb(sigma_all(1)  + ka-1)
      e  = dbl_mb(epsilon_all(1)+ ka-1)
      return
      end


c     *************************************
c     *                                   *
c     *           pspw_LJ_E               *
c     *                                   *
c     *************************************
      real*8 function pspw_LJ_E(nion,nion_qm,katm,
     >                          nfrag,indx_frag_start,size_frag,
     >                          nshl3d,rcell,
     >                          rion)
      implicit none
      integer nion,nion_qm
      integer katm(*)
      integer nfrag
      integer indx_frag_start(*),size_frag(*) 
      integer nshl3d
      real*8  rcell(nshl3d,3)
      real*8  rion(3,*)

#include "mafdecls.fh"
#include "pspw_LJ.fh"

*     **** local variables ****
      integer i,j
      integer w1,a,k1,kk1,n1
      integer w2,b,k2,kk2
      real*8  e1,s1,e2,s2,E

*     **** external functions ****
      real*8   pspw_VLJ_E_periodic,pspw_VLJ_E_periodic_self
      external pspw_VLJ_E_periodic,pspw_VLJ_E_periodic_self
      

      E = 0.0d0

c     **** QM/MM LJ energy ****
      do i=1,nion_qm
         e1 = dbl_mb(epsilon_all(1)+katm(i)-1)
         s1 = dbl_mb(sigma_all(1)  +katm(i)-1)
         do j = nion_qm+1,nion
            e2 = dbl_mb(epsilon_all(1)+katm(j)-1)
            if ((e1.gt.1.0d-11).and.(e2.gt.1.0d-11)) then
            s2 = dbl_mb(sigma_all(1)+katm(j)-1)
            E = E + pspw_VLJ_E_periodic(e1,s1,rion(1,i),
     >                                  e2,s2,rion(1,j),
     >                                  nshl3d,rcell)
            end if
         end do
      end do

c     **** MM/MM LJ energy ****
      do w1 = 1,nfrag-1
      do w2 = w1+1,nfrag
         k1 = indx_frag_start(w1)
         k2 = indx_frag_start(w2)
         kk1 = k1
         do a=1,size_frag(w1)
            kk2 = k2
            do b=1,size_frag(w2)
               e1 = dbl_mb(epsilon_all(1)+katm(kk1)-1)
               e2 = dbl_mb(epsilon_all(1)+katm(kk2)-1)
               if ((e1.gt.1.0d-11).and.(e2.gt.1.0d-11)) then
               s1 = dbl_mb(sigma_all(1)+katm(kk1)-1)
               s2 = dbl_mb(sigma_all(1)+katm(kk2)-1)
               E = E + pspw_VLJ_E_periodic(e1,s1,rion(1,kk1),
     >                                     e2,s2,rion(1,kk2),
     >                                     nshl3d,rcell)
               end if
               kk2 = kk2 + 1
            end do
            kk1 = kk1 + 1
         end do
      end do
      end do

c     **** MM/MM LJ self image energy ****
      if (nshl3d.gt.1) then
      do w1 = 1,nfrag
         k1 = indx_frag_start(w1)
         n1 = size_frag(w1)
         kk1 = k1
         do a=1,n1
            kk2 = k1
            do b=1,n1
               e1 = dbl_mb(epsilon_all(1)+katm(kk1)-1)
               e2 = dbl_mb(epsilon_all(1)+katm(kk2)-1)
               if ((e1.gt.1.0d-11).and.(e2.gt.1.0d-11)) then
               s1 = dbl_mb(sigma_all(1)+katm(kk1)-1)
               s2 = dbl_mb(sigma_all(1)+katm(kk2)-1)
               E = E + pspw_VLJ_E_periodic_self(e1,s1,rion(1,kk1),
     >                                          e2,s2,rion(1,kk2),
     >                                          nshl3d,rcell)
               end if
               kk2 = kk2 + 1
            end do
            kk1 = kk1 + 1
         end do
      end do
      end if

      pspw_LJ_E = E
      return
      end

c     *********************************************
c     *                                           *
c     *              pspw_VLJ_E_periodic          *
c     *                                           *
c     *********************************************

      real*8 function pspw_VLJ_E_periodic(e1,s1,r1,e2,s2,r2,
     >                                    nshl3d,rcell)
      implicit none
      real*8 e1,s1,r1(3)
      real*8 e2,s2,r2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  sigma12,epsilon12
      real*8  u,u6,u12,E

      E         = 0.0d0
      sigma12   = dsqrt(s1*s2)
      epsilon12 = dsqrt(e1*e2)
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=1,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         u = (sigma12/r)
         u6  = u**6
         u12 = u6**2
         E = E + 4.0d0*epsilon12*(u12-u6)
      end do

      pspw_VLJ_E_periodic = E
      return
      end



c     *********************************************
c     *                                           *
c     *              pspw_VLJ_E_periodic_self     *
c     *                                           *
c     *********************************************

      real*8 function pspw_VLJ_E_periodic_self(e1,s1,r1,e2,s2,r2,
     >                                    nshl3d,rcell)
      implicit none
      real*8 e1,s1,r1(3)
      real*8 e2,s2,r2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  sigma12,epsilon12
      real*8  u,u6,u12,E

      E         = 0.0d0
      sigma12   = dsqrt(s1*s2)
      epsilon12 = dsqrt(e1*e2)
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=2,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         u = (sigma12/r)
         u6  = u**6
         u12 = u6**2
         E = E + 4.0d0*epsilon12*(u12-u6)
      end do

      pspw_VLJ_E_periodic_self = E
      return
      end



c     *************************************
c     *                                   *
c     *           pspw_LJ_fion            *
c     *                                   *
c     *************************************
      subroutine pspw_LJ_fion(nion,nion_qm,katm,
     >                        nfrag,indx_frag_start,size_frag,
     >                        nshl3d,rcell,
     >                        rion,fion)
      implicit none
      integer nion,nion_qm
      integer katm(*)
      integer nfrag
      integer indx_frag_start(*),size_frag(*) 
      integer nshl3d
      real*8  rcell(nshl3d,3)
      real*8  rion(3,*)
      real*8  fion(3,*)

#include "mafdecls.fh"
#include "pspw_LJ.fh"

*     **** local variables ****
      integer i,j
      integer w1,a,k1,kk1,n1
      integer w2,b,k2,kk2
      real*8  e1,s1,e2,s2


c     **** QM/MM LJ energy ****
      do i=1,nion_qm
         e1 = dbl_mb(epsilon_all(1)+katm(i)-1)
         s1 = dbl_mb(sigma_all(1)  +katm(i)-1)
         do j = nion_qm+1,nion
            e2 = dbl_mb(epsilon_all(1)+katm(j)-1)
            if ((e1.gt.1.0d-11).and.(e2.gt.1.0d-11)) then
            s2 = dbl_mb(sigma_all(1)+katm(j)-1)
            call pspw_VLJ_fion_periodic(e1,s1,rion(1,i),fion(1,i),
     >                                  e2,s2,rion(1,j),fion(1,j),
     >                                  nshl3d,rcell)
            end if
         end do
      end do

c     **** MM/MM LJ energy ****
      do w1 = 1,nfrag-1
      do w2 = w1+1,nfrag
         k1 = indx_frag_start(w1)
         k2 = indx_frag_start(w2)
         kk1 = k1
         do a=1,size_frag(w1)
            kk2 = k2
            do b=1,size_frag(w2)
              e1 = dbl_mb(epsilon_all(1)+katm(kk1)-1)
              e2 = dbl_mb(epsilon_all(1)+katm(kk2)-1)
              if ((e1.gt.1.0d-11).and.(e2.gt.1.0d-11)) then
              s1 = dbl_mb(sigma_all(1)+katm(kk1)-1)
              s2 = dbl_mb(sigma_all(1)+katm(kk2)-1)
              call pspw_VLJ_fion_periodic(e1,s1,rion(1,kk1),fion(1,kk1),
     >                                    e2,s2,rion(1,kk2),fion(1,kk2),
     >                                    nshl3d,rcell)
              end if
              kk2 = kk2 + 1
            end do
            kk1 = kk1 + 1
         end do
      end do
      end do

c     **** MM/MM LJ self image energy ****
      if (nshl3d.gt.1) then
      do w1 = 1,nfrag
         k1 = indx_frag_start(w1)
         n1 = size_frag(w1)
         kk1 = k1
         do a=1,n1
            kk2 = k1
            do b=1,n1
               e1 = dbl_mb(epsilon_all(1)+katm(kk1)-1)
               e2 = dbl_mb(epsilon_all(1)+katm(kk2)-1)
               if ((e1.gt.1.0d-11).and.(e2.gt.1.0d-11)) then
               s1 = dbl_mb(sigma_all(1)+katm(kk1)-1)
               s2 = dbl_mb(sigma_all(1)+katm(kk2)-1)
               call pspw_VLJ_fion_periodic_self(e1,s1,rion(1,kk1),
     >                                                fion(1,kk1),
     >                                          e2,s2,rion(1,kk2),
     >                                                fion(1,kk2),
     >                                          nshl3d,rcell)
               end if
               kk2 = kk2 + 1
            end do
            kk1 = kk1 + 1
         end do
      end do
      end if


      return
      end


c     *********************************************
c     *                                           *
c     *              pspw_VLJ_fion_periodic       *
c     *                                           *
c     *********************************************

      subroutine pspw_VLJ_fion_periodic(e1,s1,r1,f1,
     >                                  e2,s2,r2,f2,
     >                                  nshl3d,rcell)
      implicit none
      real*8 e1,s1,r1(3),f1(3)
      real*8 e2,s2,r2(3),f2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  sigma12,epsilon12
      real*8  u,u6,u12,dVLJ

      sigma12   = dsqrt(s1*s2)
      epsilon12 = dsqrt(e1*e2)
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=1,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         u = (sigma12/r)
         u6  = u**6
         u12 = u6**2
         dVLJ = -(4.0d0*epsilon12/r)*(12.0d0*u12-6.0d0*u6)
         f1(1) = f1(1) - (x/r)*dVLJ
         f1(2) = f1(2) - (y/r)*dVLJ
         f1(3) = f1(3) - (z/r)*dVLJ
         f2(1) = f2(1) + (x/r)*dVLJ
         f2(2) = f2(2) + (y/r)*dVLJ
         f2(3) = f2(3) + (z/r)*dVLJ
      end do

      return
      end



c     *********************************************
c     *                                           *
c     *         pspw_VLJ_fion_periodic_self       *
c     *                                           *
c     *********************************************

      subroutine pspw_VLJ_fion_periodic_self(e1,s1,r1,f1,
     >                                  e2,s2,r2,f2,
     >                                  nshl3d,rcell)
      implicit none
      real*8 e1,s1,r1(3),f1(3)
      real*8 e2,s2,r2(3),f2(3)
      integer nshl3d
      real*8  rcell(nshl3d,3)

*     **** local variables ****
      integer l
      real*8  dx,dy,dz
      real*8  x,y,z,r
      real*8  sigma12,epsilon12
      real*8  u,u6,u12,dVLJ

      sigma12   = dsqrt(s1*s2)
      epsilon12 = dsqrt(e1*e2)
      dx = r1(1) - r2(1)
      dy = r1(2) - r2(2)
      dz = r1(3) - r2(3)
      do l=2,nshl3d
         x = dx + rcell(l,1)
         y = dy + rcell(l,2)
         z = dz + rcell(l,3)
         r = dsqrt(x**2 + y**2 + z**2)
         u = (sigma12/r)
         u6  = u**6
         u12 = u6**2
         dVLJ = -(4.0d0*epsilon12/r)*(12.0d0*u12-6.0d0*u6)
         f1(1) = f1(1) - (x/r)*dVLJ
         f1(2) = f1(2) - (y/r)*dVLJ
         f1(3) = f1(3) - (z/r)*dVLJ
         f2(1) = f2(1) + (x/r)*dVLJ
         f2(2) = f2(2) + (y/r)*dVLJ
         f2(3) = f2(3) + (z/r)*dVLJ
      end do

      return
      end



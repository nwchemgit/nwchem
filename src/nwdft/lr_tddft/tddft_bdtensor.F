      subroutine tddft_bdtensor(tmom, eexcit, el_transdip, magtd_giao)

c     calculate the Buckingham-Dunn rotatory strength tensor
c     for an excitation

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"

c     input:

      double precision tmom(20), eexcit, el_transdip(3), magtd_giao(3)

c     local variables:

      double precision trace, qmat(3,3), rmat(3,3), tmat(3,3)
      integer i,j, k, l, m, n

c     Levi Civita tensor:
      integer LCTensor(3,3,3)
      data LCTensor/27 * 0/ ! nonzero elements are initialized below
c     Define translation table for quadrupole incices in
c     packed storage in tmom(i) below
c     XX=5,XY=YX=6, XZ=ZX=7,YY=8,YZ=ZY=9,ZZ=10
      integer qindex(3,3)
      data qindex(1,1),qindex(1,2),qindex(1,3)/5,6,7/
      data qindex(2,1),qindex(2,2),qindex(2,3)/6,8,9/
      data qindex(3,1),qindex(3,2),qindex(3,3)/7,9,10/

c     caution: speed of light and rotatory strength conversion factor
c     are hardwired here:
      double precision cl, rau2cgs
      parameter (cl = 137.0359895d0, rau2cgs = 235.726327d0)

      double precision zero, one, two, three, half, third
      parameter (zero = 0d0, one=1d0, two=2d0, three=3d0,
     &   half=one/two,
     &   third=one/three)

      double precision const

      logical nodezero

c     ==================================================================

      nodezero=(ga_nodeid().eq.0)

c     initialize nonzero LC Tensor elements here. Doing so in a data
c     statement caused gfortran to choke
      
      LCtensor(1,2,3) = 1
      LCtensor(2,3,1) = 1
      LCtensor(3,1,2) = 1      
      LCtensor(2,1,3) = -1
      LCtensor(1,3,2) = -1
      LCtensor(3,2,1) = -1        

c     step 1: traceless transition quadrupole moment
            
      do i = 1,3
        do j = 1,3
          qmat(i,j) = three * tmom(qindex(i,j))
          if (j.ne.i)
     &       qmat(j,i) = three * tmom(qindex(j,i))
        end do
      end do

c     we assume that this transition dipole matrix is symmetric
      
      trace = zero
      do i = 1,3
        trace = trace +  tmom(qindex(i,i))
      end do
      
      do i = 1,3
        qmat(i,i) = qmat(i,i) - trace
      end do                
      
      do i = 1,3
        do j = 1,3
         qmat(i,j) = half * qmat(i,j)
        end do
      end do
      
      trace = zero
      do i = 1,3
        trace = trace + qmat(i,i)
      end do
      
      if (nodezero) then
        write (luout,*) 'traceless quadrupole moment tensor'
        call output(qmat,1,3,1,3,3,3,1)
        write (luout,'(1X,a,f15.5)') 'Trace = ', trace
      end if
      
c     the SOS equation for the D-Q polarizability
c     has a factor of the excitation energy,
c     just like other real polarizabilities. So we need to add that
c     here

      do i = 1,3
        do j = 1,3
         qmat(i,j) = qmat(i,j) * eexcit
        end do
      end do

c     step 2: the symmetrized rotatory strength tensor
      
c     rmat(:,:) = zero  ! comment out implicit
      do i = 1,3
        do j = 1,3
          rmat(i,j) = zero
        end do
      end do
c      
      do i = 1,3
        do j = 1,i
          rmat(i,j) =
     &       half * el_transdip(i)*magtd_giao(j) +
     &       half * el_transdip(j)*magtd_giao(i) 
          if (j.ne.i) rmat(j,i) = rmat(i,j)
        end do
      end do
      
      trace = zero
      do i = 1,3
        trace = trace + rmat(i,i)
      end do
      
      if (nodezero) then
        write (luout,*) 'Symmetrized Rotatory strength tensor'
        call output(rmat,1,3,1,3,3,3,1)
        write (luout,'(1X,a,f15.5)') 'Trace = ', trace
      end if
      
c     step 3: combine rotatory strength and quadrupole terms

c      rmat = 0 ! debug
      
      do i = 1,3
        do j = 1,3
          do k = 1,3
            do l = 1,3
              rmat(i,j) = rmat(i,j) 
     &           + half * third * LCtensor(i,k,l)
     &           * el_transdip(k) * qmat(l,j) 
     &           + half * third * LCtensor(j,k,l)
     &           * el_transdip(k) * qmat(l,i) 
              
            end do
          end do
        end do
      end do
      
c      rmat(:,:) = rmat(:,:) * two * rau2cgs
      
      trace = zero
      do i = 1,3
        trace = trace + rmat(i,i)
      end do
      
      if (nodezero) then
        write (luout,'(/5x,a/5x,a)')
     &     'Symmetrized Buckingham-Dunn Rotatory strength tensor',
     &     '(see JCP 125 (2006), 034102 and references therein)'
        call output(rmat,1,3,1,3,3,3,1)
        write (luout,'(5X,a,f15.5)') 'Trace = ', trace
      end if

c     also output the orientation-dependent rotatory strength

      trace = zero
      do i = 1,3
        trace = trace + rmat(i,i)
      end do

c     rmat(:,:) = -half * rmat(:,:) ! comment out implicit
      do i = 1,3
        do j = 1,3
          rmat(i,j) = -half * rmat(i,j)
        end do
      end do
c
      do i = 1,3
        rmat(i,i) = rmat(i,i) + half * trace
      end do
c
      trace = zero
      do i = 1,3
        trace = trace + rmat(i,i)
      end do
c      
      if (nodezero) then
        write (luout,'(/5x,a)')
     &     'Rotatory strength tensor'
        call output(rmat,1,3,1,3,3,3,1)
        write (luout,'(5X,a,f15.5)') 'Trace = ', trace
      end if

c     step 4: redo steps 2 and 3 for the antisymmetric part
c     of the BD tensor. 


c$$$      rmat(:,:) = zero
c$$$      
c$$$      do i = 1,3
c$$$        do j = 1,3
c$$$          rmat(i,j) =
c$$$     &         half * el_transdip(i)*magtd_giao(j) 
c$$$     &       - half * el_transdip(j)*magtd_giao(i) 
c$$$        end do
c$$$      end do
c$$$      
c$$$      trace = zero
c$$$      do i = 1,3
c$$$        trace = trace + rmat(i,i)
c$$$      end do
c$$$      
c$$$      if (nodezero) then
c$$$        write (luout,*) 'Antisymmetrized Rotatory strength tensor'
c$$$        call output(rmat,1,3,1,3,3,3,1)
c$$$        write (luout,'(1X,a,f15.5)') 'Trace = ', trace
c$$$      end if
c$$$
c$$$
c$$$      tmat(:,:) = zero
c$$$
c$$$      do i = 1,3
c$$$        do j = 1,3
c$$$          do k = 1,3
c$$$            do l = 1,3
c$$$              tmat(i,j) = tmat(i,j) 
c$$$     &           - half * third * LCtensor(i,k,l)
c$$$     &           * el_transdip(k) * qmat(l,j) 
c$$$     &           + half * third * LCtensor(j,k,l)
c$$$     &           * el_transdip(k) * qmat(l,i) 
c$$$              
c$$$            end do
c$$$          end do
c$$$        end do
c$$$      end do
c$$$
c$$$      trace = zero
c$$$      do i = 1,3
c$$$        trace = trace + tmat(i,i)
c$$$      end do
c$$$      
c$$$      if (nodezero) then
c$$$        write (luout,*) 'Antisymmetrized e-q tensor'
c$$$        call output(tmat,1,3,1,3,3,3,1)
c$$$        write (luout,'(1X,a,f15.5)') 'Trace = ', trace
c$$$      end if
c$$$
c$$$      rmat = rmat + tmat
c$$$
c$$$      trace = zero
c$$$      do i = 1,3
c$$$        trace = trace + rmat(i,i)
c$$$      end do
c$$$      
c$$$      if (nodezero) then
c$$$        write (luout,*) 'Antisymmetrized BD tensor in au'
c$$$        call output(rmat,1,3,1,3,3,3,1)
c$$$        write (luout,'(1X,a,f15.5)') 'Trace = ', trace
c$$$      end if
c$$$      
c$$$      rmat(:,:) = rmat(:,:) * two * rau2cgs
c$$$      
c$$$      trace = zero
c$$$      do i = 1,3
c$$$        trace = trace + rmat(i,i)
c$$$      end do
c$$$      
c$$$      if (nodezero) then
c$$$        write (luout,'(/5x,a)')
c$$$     &     'Antisymmetrized Buckingham-Dunn Rotatory strength tensor'
c$$$        call output(rmat,1,3,1,3,3,3,1)
c$$$        write (luout,'(5X,a,f15.5)') 'Trace = ', trace
c$$$      end if


c     step 5: calculate the Graham-Raab formula

      tmat(:,:) = zero
      
      do i = 1,3
        do j = 1,3
          tmat(i,j) = el_transdip(i)*magtd_giao(j) 
        end do
      end do
      
      rmat(:,:) = zero

      do i = 1,3                ! alpha
        do j = 1,3              ! beta

          const = 0.25d0

          do k = 1,3            ! gamma
            do l = 1,3          ! delta 
              do m = 1,3        ! epsilon
                rmat(i,j) = rmat(i,j) +
     &             const*LCTensor(i,k,l)*LCTensor(j,k,m)*tmat(l,m) +
     &             const*LCTensor(j,k,l)*LCTensor(i,k,m)*tmat(l,m)
              end do
            end do
          end do


          const = - 0.25d0 / 6d0

          do k = 1,3            ! gamma
            do l = 1,3          ! delta 
                rmat(i,j) = rmat(i,j) +
     &             const*el_transdip(k)*qmat(l,j) *
     &             (LCTensor(i,k,l) - LCTensor(i,l,k) )

                rmat(i,j) = rmat(i,j) +
     &             const*el_transdip(k)*qmat(l,i) *
     &             (LCTensor(j,k,l) - LCTensor(j,l,k) )
              end do
          end do 

        end do                  ! j
        end do                  ! i

        trace = zero
        do i = 1,3
          trace = trace + rmat(i,i)
        end do
        
        if (nodezero) then
          write (luout,'(/5x,a)')
     &       'Symmetric Graham-Raab tensor'
          call output(rmat,1,3,1,3,3,3,1)
          write (luout,'(5X,a,f15.5)') 'Trace = ', trace
        end if

      rmat(:,:) = zero

      do i = 1,3                ! alpha
        do j = 1,3              ! beta

          const = 0.25d0

          do k = 1,3            ! gamma
            do l = 1,3          ! delta 
              do m = 1,3        ! epsilon
                rmat(i,j) = rmat(i,j) +
     &             const*LCTensor(i,k,l)*LCTensor(j,k,m)*tmat(l,m) -
     &             const*LCTensor(j,k,l)*LCTensor(i,k,m)*tmat(l,m)
              end do
            end do
          end do


          const = - 0.25d0 / 6d0

          do k = 1,3            ! gamma
            do l = 1,3          ! delta 
                rmat(i,j) = rmat(i,j) +
     &             const*el_transdip(k)*qmat(l,j) *
     &             (LCTensor(i,k,l) - LCTensor(i,l,k) )

                rmat(i,j) = rmat(i,j) -
     &             const*el_transdip(k)*qmat(l,i) *
     &             (LCTensor(j,k,l) - LCTensor(j,l,k) )
              end do
            end do 
            
          end do                ! j
        end do                  ! i

        trace = zero
        do i = 1,3
          trace = trace + rmat(i,i)
        end do
        
        if (nodezero) then
          write (luout,'(/5x,a)')
     &       'Antisymmetric Graham-Raab tensor'
          call output(rmat,1,3,1,3,3,3,1)
          write (luout,'(5X,a,f15.5)') 'Trace = ', trace
        end if

      return

c     =================================================================

      end
c $Id$

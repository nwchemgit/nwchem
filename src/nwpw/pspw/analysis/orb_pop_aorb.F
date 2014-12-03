c
c $Id$
c
 
*     ***********************************
*     *					*
*     *        Orb_pop_aorb		*
*     *					*
*     ***********************************

      subroutine Orb_pop_aorb(flag,ispin,ne,npack1,nemax,c2,
     >                        lmax,a,sum)
      implicit none
      integer    iunit,flag
      integer    ispin,ne(2)
      integer    npack1,nemax
      complex*16 c2(npack1,*)

      integer lmax(*)
      real*8  a(36,nemax,*)
      real*8  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer spin
      integer i,k,l,n,ll,l1,l2,jj
      integer llmax,nfft3d
      real*8  w,pi,twopi
      integer n1(2),n2(2),q,pj,taskid_j,nnn

*     *** MA Arrays ****
      integer phi(2)


*     **** external functions ****
      integer  ion_nion_qm,ion_katm_qm,aorbs_get_basis_number
      external ion_nion_qm,ion_katm_qm,aorbs_get_basis_number

      call D3dB_nfft3d(1,nfft3d)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)


      value = BA_push_get(mt_dcpl, nfft3d,'phi',phi(2),phi(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


*     -------------------
*     POPULATION ANALYSIS
*     -------------------
      pi = 4.0d0*datan(1.0d0)
      twopi = 2.0d0*pi

      nnn = 36*nemax*ion_nion_qm()
      call dcopy(nnn,0.0d0,0,A,1)
      call Parallel2d_taskid_j(taskid_j)

      call phafac()
      do i=1,ion_nion_qm()


        DO spin=1,ispin
          do n=n1(spin),n2(spin)
            if (flag.eq.0) then
               call Dneall_ntoqp(n,q,pj)
            else
               q = n
               pj=taskid_j
            end if
            if (pj.eq.taskid_j) then

* S-COMPONENT
            jj=aorbs_get_basis_number(i,1)
            call aorbs_aorb(jj,dcpl_mb(phi(1)))
            call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
            A(1,N,I)=w


* P-COMPONENT
            IF(lmax(ion_katm_qm(i)).GE.1) THEN

*             **** px component ****
              jj=aorbs_get_basis_number(i,2)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(2,N,I)=w

*             **** pz component ****
              jj=aorbs_get_basis_number(i,3)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(3,N,I)=w

*             **** py component ****
              jj=aorbs_get_basis_number(i,4)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(4,N,I)=w
            ENDIF

* D-COMPONENT
            IF(lmax(ion_katm_qm(I)).GE.2) THEN

*             **** dx2-y2 component ****
              jj=aorbs_get_basis_number(i,5)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(5,N,I)=w

*             **** dzx component ****
              jj=aorbs_get_basis_number(i,6)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(6,N,I)=w

*             **** d3z2-1 component ****
              jj=aorbs_get_basis_number(i,7)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(7,N,I)=w

*             **** dyz component ****
              jj=aorbs_get_basis_number(i,8)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(8,N,I)=w

*             **** dxy component ****
              jj=aorbs_get_basis_number(i,9)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(9,N,I)=w
            ENDIF 

* F-COMPONENT
            IF(lmax(ion_katm_qm(I)).GE.3) THEN

*             **** fx(x2-3y2) component ****
              jj=aorbs_get_basis_number(i,10)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(10,N,I)=w

*             **** fz(x2-y2) component ****
              jj=aorbs_get_basis_number(i,11)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(11,N,I)=w

*             **** fx(5z2-1) component ****
              jj=aorbs_get_basis_number(i,12)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(12,N,I)=w

*             **** fz(5z2-3) component ****
              jj=aorbs_get_basis_number(i,13)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(13,N,I)=w

*             **** fy(5z2-1) component ****
              jj=aorbs_get_basis_number(i,14)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(14,N,I)=w

*             **** fxyz component ****
              jj=aorbs_get_basis_number(i,15)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(15,N,I)=w

*             **** fy(3x2-y2) component ****
              jj=aorbs_get_basis_number(i,16)
              call aorbs_aorb(jj,dcpl_mb(phi(1)))
              call Pack_cc_idot(1,dcpl_mb(phi(1)),c2(1,q),w)
              A(16,N,I)=w
            ENDIF

            end if
          end do
        end do
      end do
      if (flag.eq.0) then
         call Parallel_Vector_SumAll(nnn,A)
      else
         call D3dB_Vector_SumAll(nnn,A)
      end if
 
      DO spin=1,ispin
        DO n=n1(spin),n2(spin)
          sum(n)=0.0d0

          DO i=1,ion_nion_qm()
            llmax=(lmax(ion_katm_qm(I))+1)**2
            DO ll=1,llmax
              sum(n)=sum(n)+A(ll,n,i)**2
            END DO
          END DO

          DO I=1,ion_nion_qm()
            llmax=(lmax(ion_katm_qm(i))+1)**2
            DO ll=1,llmax
              A(ll,n,i)=A(ll,n,i)/dsqrt(sum(n))
            END DO
          END DO

        END DO
      END DO

      value = BA_pop_stack(phi(2))

      return 
      end
 

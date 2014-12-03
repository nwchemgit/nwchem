c
c $Id$
c
 
*     ***********************************
*     *					*
*     *        Orb_pop_borb		*
*     *					*
*     ***********************************
      subroutine Orb_pop_borb(flag,nbq,ispin,ne,
     >                        npack1,nemax,c2,
     >                        lmax,a,sum)
      implicit none
      integer    flag
      integer    nbq,ispin,ne(2)
      integer    npack1,nemax
      complex*16 c2(npack1,*)

      integer lmax(*)
      complex*16 a(36,nemax,*)
      real*8  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer spin
      integer ii,k,l,n,ll,l1,l2,jj
      integer llmax,nfft3d
      real*8  pi,twopi
      integer n1(2),n2(2),q,nnn
      complex*16 w

*     *** MA Arrays ****
      integer phi(2)

*     **** external functions ****
      integer  ion_nion_qm,ion_katm_qm,borbs_get_basis_number
      external ion_nion_qm,ion_katm_qm,borbs_get_basis_number
      real*8   ddot
      external ddot

      call C3dB_nfft3d(1,nfft3d)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)

      if (.not.BA_push_get(mt_dcpl, nfft3d,'phi',phi(2),phi(1)))
     >   call errquit('Orb_pop_borb:out of stack memory',0,MA_ERR)

*     -------------------
*     POPULATION ANALYSIS
*     -------------------
      pi = 4.0d0*datan(1.0d0)
      twopi = 2.0d0*pi

      nnn = 2*36*nemax*ion_nion_qm()
      call dcopy(nnn,0.0d0,0,A,1)

      call cphafac()
      do ii=1,ion_nion_qm()
           do spin=1,ispin
           do n=n1(spin),n2(spin)


* S-COMPONENT
            jj=borbs_get_basis_number(ii,1)
            call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
            call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
            A(1,n,ii)=w


* P-COMPONENT
            if(lmax(ion_katm_qm(ii)).GE.1) then

*             **** px component ****
              jj=borbs_get_basis_number(ii,2)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(2,n,ii)=w

*             **** pz component ****
              jj=borbs_get_basis_number(ii,3)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(3,n,ii)=w

*             **** py component ****
              jj=borbs_get_basis_number(ii,4)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(4,n,ii)=w
            end if

* D-COMPONENT
            if(lmax(ion_katm_qm(ii)).GE.2) then

*             **** dx2-y2 component ****
              jj=borbs_get_basis_number(ii,5)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(5,n,ii)=w

*             **** dzx component ****
              jj=borbs_get_basis_number(ii,6)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(6,n,ii)=w

*             **** d3z2-1 component ****
              jj=borbs_get_basis_number(ii,7)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(7,n,ii)=w

*             **** dyz component ****
              jj=borbs_get_basis_number(ii,8)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(8,n,ii)=w

*             **** dxy component ****
              jj=borbs_get_basis_number(ii,9)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(9,n,ii)=w
            end if

* F-COMPONENT
            if(lmax(ion_katm_qm(ii)).GE.3) then

*             **** fx(x2-3y2) component ****
              jj=borbs_get_basis_number(ii,10)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(10,n,ii)=w

*             **** fz(x2-y2) component ****
              jj=borbs_get_basis_number(ii,11)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(11,n,ii)=w

*             **** fx(5z2-1) component ****
              jj=borbs_get_basis_number(ii,12)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(12,n,ii)=w

*             **** fz(5z2-3) component ****
              jj=borbs_get_basis_number(ii,13)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(13,n,ii)=w

*             **** fy(5z2-1) component ****
              jj=borbs_get_basis_number(ii,14)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(14,n,ii)=w

*             **** fxyz component ****
              jj=borbs_get_basis_number(ii,15)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(15,n,ii)=w

*             **** fy(3x2-y2) component ****
              jj=borbs_get_basis_number(ii,16)
              call borbs_borb(nbq,jj,dcpl_mb(phi(1)))
              call Cram_cc_izdot(nbq,dcpl_mb(phi(1)),c2(1,n),w)
              A(16,n,ii)=w
            end if

           end do
           end do
      end do
      call C3dB_Vector_SumAll(nnn,A)
 
      do spin=1,ispin
        do n=n1(spin),n2(spin)
          sum(n)=0.0d0

          do ii=1,ion_nion_qm()
            llmax=(lmax(ion_katm_qm(ii))+1)**2
            sum(n) = sum(n) + ddot(2*llmax,A(1,n,ii),1,A(1,n,ii),1)
c            do ll=1,llmax
c              sum(n)=sum(n)+A(ll,n,ii)**2
c            end do
          end do

          do ii=1,ion_nion_qm()
            llmax=(lmax(ion_katm_qm(ii))+1)**2
            do ll=1,llmax
              A(ll,n,ii)=A(ll,n,ii)/dsqrt(sum(n))
            end do
          end do

        end do
      end do

      if (.not.BA_pop_stack(phi(2)))
     > call errquit('Orb_pop_borb:popping stack',0, MA_ERR)


      return 
      end
 

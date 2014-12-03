c
c $Id$
c
 
*     ***********************************
*     *					*
*     *        Orb_pop_kawai		*
*     *					*
*     ***********************************

      subroutine Orb_pop_kawai(flag,ispin,ne,npack1,nemax,c2,
     >                        lmax,b0,a,sum)
      implicit none
      integer    iunit,flag
      integer    ispin,ne(2)
      integer    npack1,nemax
      complex*16 c2(npack1,*)

      integer lmax(*)
      real*8  b0(0:5,*)
      real*8  a(36,nemax,*)
      real*8  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer spin
      integer i,k,l,n,ll,l1,l2
      integer llmax,nfft3d
      real*8  omega,w,pi,twopi
      real*8  f0,f1,f2,f3,f21,f22,f23,f24,f25
      integer n1(2),n2(2),q,pj,nnn,taskid_j

*     *** MA Arrays ****
      integer exr(2),cw(2)
      integer Gx(2),Gy(2),Gz(2)


*     **** external functions ****
      integer  G_indx,ion_nion_qm,ion_katm_qm
      external G_indx,ion_nion_qm,ion_katm_qm
      real*8   lattice_omega
      external lattice_omega

      call Parallel2d_taskid_j(taskid_j)
      call D3dB_nfft3d(1,nfft3d)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)


      value = BA_push_get(mt_dcpl, npack1,'cw',cw(2),cw(1))
      value = value.and.
     >        BA_push_get(mt_dcpl, npack1,'exr',exr(2),exr(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

*     **** define Gx,Gy and Gz in packed space ****
      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(1,dbl_mb(Gx(1)))
      call Pack_t_pack(1,dbl_mb(Gy(1)))
      call Pack_t_pack(1,dbl_mb(Gz(1)))


*     -------------------
*     POPULATION ANALYSIS
*     -------------------
      omega = lattice_omega()
      pi = 4.0d0*datan(1.0d0)
      twopi = 2.0d0*PI
      F0=dsqrt(4.0d0*PI/omega)
      F1=F0*dsqrt(3.0d0)
      F2=F0*dsqrt(15.0d0)
      F3=F0*dsqrt(105.0d0)

      F21=dsqrt(PI/5.0d0/omega)/3.0d0
      F22=dsqrt(2.d0*PI/15.d0/omega)
      F23=2.0d0*F22
      F24=F23
      F25=F23

      nnn = 36*nemax*ion_nion_qm()
      call dcopy(nnn,0.0d0,0,A,1)

      call phafac()
      do i=1,ion_nion_qm()

*       **** get the structure factor ****
        call strfac_pack(1,i,dcpl_mb(exr(1)))
c        call Pack_c_pack(1,dcpl_mb(exr(1)))


        DO spin=1,ispin
          do n=n1(spin),n2(spin)
          if (flag.eq.0) then
             call Dneall_ntoqp(n,q,pj)
          else
             q = n
             pj = taskid_j
          end if
          if (pj.eq.taskid_j) then



* S-COMPONENT
            call Pack_cc_idot(1,dcpl_mb(exr(1)),c2(1,q),w)
            A(1,N,I)=F0/b0(0,ion_katm_qm(I))*w


* P-COMPONENT
            IF(lmax(ion_katm_qm(i)).GE.1) THEN

*             **** px component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(exr(1)+K-1)
     >                            *DCMPLX(0.0d0,-dbl_mb(gx(1)+K-1))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(2,N,I)=F1/b0(1,ion_katm_qm(I))*w

*             **** pz component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >                            *DCMPLX(0.0d0,-dbl_mb(gz(1)+K-1))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(3,N,I)=F1/b0(1,ion_katm_qm(I))*w

*             **** py component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >                            *DCMPLX(0.0d0,-dbl_mb(gy(1)+K-1))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(4,N,I)=F1/b0(1,ion_katm_qm(I))*w
            ENDIF


* D-COMPONENT
            IF(lmax(ion_katm_qm(I)).GE.2) THEN

*             **** dx2-y2 component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *( dbl_mb(gx(1)+k-1)**2
     >             - dbl_mb(gy(1)+k-1)**2)/(2.0d0)
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(5,N,I)=F2/b0(2,ion_katm_qm(I))*w

*             **** dzx component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *( dbl_mb(gz(1)+k-1)*dbl_mb(gx(1)+k-1))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(6,N,I)=F2/b0(2,ion_katm_qm(I))*w

*             **** d3z2-1 component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *(3.0d0*dbl_mb(gz(1)+k-1)**2 - 1.0d0)
     >            /(2.0d0*dsqrt(3.0d0))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(7,N,I)=F2/b0(2,ion_katm_qm(I))*w

*             **** dyz component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *dbl_mb(gy(1)+k-1)*dbl_mb(gz(1)+k-1)
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(8,N,I)=F2/b0(2,ion_katm_qm(I))*w

*             **** dxy component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *( dbl_mb(gx(1)+k-1)*dbl_mb(gy(1)+k-1))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(9,N,I)=F2/b0(2,ion_katm_qm(I))*w
            ENDIF 

* F-COMPONENT
            IF(lmax(ion_katm_qm(I)).GE.3) THEN

*             **** fx(x2-3y2) component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *dcmplx(0.0d0, dbl_mb(gx(1)+k-1)
     >            *( 4.0d0*dbl_mb(gx(1)+k-1)**2 
     >             - 3.0d0*(1.0d0-dbl_mb(gz(1)+k-1)**2)))
     >            /(dsqrt(24.0d0))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(10,N,I)=F3/b0(3,ion_katm_qm(I))*w

*             **** fz(x2-y2) component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *dcmplx(0.0d0,dbl_mb(gz(1)+k-1)
     >            *(dbl_mb(gx(1)+k-1)**2 - dbl_mb(gy(1)+k-1)**2))
     >            /(2.0d0)
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(11,N,I)=F3/b0(3,ion_katm_qm(I))*w

*             **** fx(5z2-1) component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *dcmplx(0.0d0,dbl_mb(gx(1)+k-1)
     >            *(5.0d0*dbl_mb(gz(1)+k-1)**2 - 1.0d0))
     >            /(dsqrt(40.0d0))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(12,N,I)=F3/b0(3,ion_katm_qm(I))*w

*             **** fz(5z2-3) component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *dcmplx(0.0d0,dbl_mb(gz(1)+k-1)
     >            *(5.0d0*dbl_mb(gz(1)+k-1)**2 - 3.0d0))
     >            /(dsqrt(60.0d0))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(13,N,I)=F3/b0(3,ion_katm_qm(I))*w

*             **** fy(5z2-1) component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *dcmplx(0.0d0,dbl_mb(gy(1)+k-1)
     >            *(3.0d0*dbl_mb(gz(1)+k-1)**2 - 1.0d0))
     >            /(dsqrt(40.0d0))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(14,N,I)=F3/b0(3,ion_katm_qm(I))*w

*             **** fxyz component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *dcmplx(0.0d0,dbl_mb(gx(1)+k-1)
     >                    *dbl_mb(gy(1)+k-1)*dbl_mb(gz(1)+k-1))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(15,N,I)=F3/b0(3,ion_katm_qm(I))*w

*             **** fy(3x2-y2) component ****
              do k=1,npack1
                dcpl_mb(CW(1)+k-1)=dcpl_mb(EXR(1)+k-1)
     >            *dcmplx(0.0d0,dbl_mb(gy(1)+k-1)
     >            *(3.0d0*(1.0d0-dbl_mb(gz(1)+k-1)**2) 
     >             -4.0d0*dbl_mb(gy(1)+k-1)**2))
     >            /(dsqrt(24.0d0))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(16,N,I)=F3/b0(3,ion_katm_qm(I))*w

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
            llmax=(lmax(ion_katm_qm(i))+1)**2
            DO ll=1,llmax
              sum(n)=sum(n)+A(ll,n,i)**2
            END DO
          END DO

          DO i=1,ion_nion_qm()
            llmax=(lmax(ion_katm_qm(I))+1)**2
            DO ll=1,llmax
              A(ll,n,i)=A(ll,n,i)/dsqrt(sum(n))
            END DO
          END DO

        END DO
      END DO

      value = BA_pop_stack(Gz(2))
      value = BA_pop_stack(Gy(2))
      value = BA_pop_stack(Gx(2))
      value = BA_pop_stack(exr(2))
      value = BA_pop_stack(cw(2))

      return 
      end
 

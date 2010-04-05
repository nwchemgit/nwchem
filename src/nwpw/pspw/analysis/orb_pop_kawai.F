c
c $Id: orb_pop_kawai.F,v 1.2 2008-11-03 06:06:13 bylaska Exp $
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
      real*8  B0(0:5,*)
      real*8  a(36,nemax,*)
      real*8  sum(*)

#include "mafdecls.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer spin
      integer i,k,l,n,ll,l1,l2
      integer llmax,nfft3d
      real*8  omega,w,pi,twopi
      real*8  f0,f1,f2,f21,f22,f23,f24,f25
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


      value = MA_push_get(mt_dcpl, npack1,'cw',cw(2),cw(1))
      value = value.and.
     >        MA_push_get(mt_dcpl, nfft3d,'exr',exr(2),exr(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
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
        call strfac(i,dcpl_mb(exr(1)))
        call Pack_c_pack(1,dcpl_mb(exr(1)))


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
            A(1,N,I)=F0/B0(0,ion_katm_qm(I))*w


* P-COMPONENT
            IF(lmax(ion_katm_qm(i)).GE.1) THEN

*             **** px component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >                            *DCMPLX(0.0d0,-dbl_mb(gx(1)+K-1))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(2,N,I)=F1/B0(1,ion_katm_qm(I))*w

*             **** py component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >                            *DCMPLX(0.0d0,-dbl_mb(gy(1)+K-1))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(3,N,I)=F1/B0(1,ion_katm_qm(I))*w

*             **** pz component ****
              do k=1,npack1
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >                            *DCMPLX(0.0d0,-dbl_mb(gz(1)+K-1))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(4,N,I)=F1/B0(1,ion_katm_qm(I))*w
            ENDIF

* D-COMPONENT
            IF(lmax(ion_katm_qm(I)).GE.2) THEN

*             **** dzz component ****
              do k=1,npack1
c               CW(K)=EXR(K)*(3.0d0*g(K,3)**2-1.0)/(2.0d0*dsqrt(3.0d0)
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *(3.0d0*dbl_mb(gz(1)+k-1)**2 - 1.0d0)
     >            /(2.0d0*dsqrt(3.0d0))
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(5,N,I)=F2/B0(2,ion_katm_qm(I))*w

*             **** dx2-y2 component ****
              do k=1,npack1
c               CW(K)=EXR(K)*(g(K,1)**2-g(K,2)**2)/(2.0d0)
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *( dbl_mb(gx(1)+k-1)**2
     >             - dbl_mb(gy(1)+k-1)**2)/(2.0d0)
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(6,N,I)=F2/B0(2,ion_katm_qm(I))*w

*             **** dxy component ****
              do k=1,npack1
c               CW(K)=EXR(K)*g(K,1)*g(K,2)
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *( dbl_mb(gx(1)+k-1)*dbl_mb(gy(1)+k-1)**2)
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(7,N,I)=F2/B0(2,ion_katm_qm(I))*w

*             **** dyz component ****
              do k=1,npack1
c               CW(K)=EXR(K)*g(K,2)*g(K,3)
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *( dbl_mb(gy(1)+k-1)*dbl_mb(gz(1)+k-1)**2)
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(8,N,I)=F2/B0(2,ion_katm_qm(I))*w

*             **** dzx component ****
              do k=1,npack1
c               CW(K)=EXR(K)*g(K,3)*g(K,1)
                dcpl_mb(CW(1)+K-1)=dcpl_mb(EXR(1)+K-1)
     >            *( dbl_mb(gz(1)+k-1)*dbl_mb(gx(1)+k-1)**2)
              end do
              call Pack_cc_idot(1,dcpl_mb(cw(1)),c2(1,q),w)
              A(9,N,I)=F2/B0(2,ion_katm_qm(I))*w
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
 
      DO SPIN=1,ISPIN
        DO N=N1(SPIN),N2(SPIN)
          SUM(N)=0.0d0

          DO I=1,ion_nion_qm()
            LLMAX=(lmax(ion_katm_qm(I))+1)**2
            DO LL=1,LLMAX
              SUM(N)=SUM(N)+A(LL,N,I)**2
            END DO
          END DO

          DO I=1,ion_nion_qm()
            LLMAX=(lmax(ion_katm_qm(I))+1)**2
            DO LL=1,LLMAX
              A(LL,N,I)=A(LL,N,I)/DSQRT(SUM(N))
            END DO
          END DO

        END DO
      END DO

      value = MA_pop_stack(Gz(2))
      value = MA_pop_stack(Gy(2))
      value = MA_pop_stack(Gx(2))
      value = MA_pop_stack(exr(2))
      value = MA_pop_stack(cw(2))

      return 
      end
 

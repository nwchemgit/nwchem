*     ***********************************
*     *             			*
*     *          c_OrthoCheck		*
*     *             			*
*     ***********************************

      subroutine c_OrthoCheck(ispin,ne,nbrill,psi1)
      implicit none 
      integer ispin,ne(2),nbrill
      double complex psi1(*)


*    *** local variables ***

      integer MASTER,taskid
      parameter(MASTER=0)
      integer  nb,ms,npack1,i,j,ii,jj
      integer n1(2),n2(2),neall,indx
      real*8  error
      complex*16 w
      character*255 full_filename

      call Parallel_taskid(taskid)
      call Cram_max_npack(npack1)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)
      neall = ne(1)+ne(2)

*     **** produce CHECK FILE ****
      call util_file_name('ORTHOCHECK',
     >                             .true.,
     >                             .false.,
     >                             full_filename)
      if (taskid.eq.MASTER) then
         open(unit=17,file=full_filename,form='formatted')
      end if

     
*     **** check orthonormality ****
      if (taskid.eq.MASTER) then
         write(17,1350)
         write(17,*) ne,n1,n2
      end if

      do nb=1,nbrill
      error = 0.0d0
      indx = (nb-1)*npack1*neall
      do ms=1,ispin
         do i=n1(ms),n2(ms)
            ii = i-n1(ms)+1
            do j=i,n2(ms)
               jj = j-n1(ms)+1
               
               call Cram_cc_zdot(nb,
     >                            psi1(1+(i-1)*npack1+indx),
     >                            psi1(1+(j-1)*npack1+indx),
     >                            w)

               if (i.eq.j) then
                   error = dabs(1.0d0-dble(w)) + dabs(dimag(w))
               else
                   error = dabs(dble(w)) + dabs(dimag(w))
               end if
             
               if (taskid.eq.MASTER) then
                  write(17,1360) nb,ms,ii,jj,w
               end if
            end do
         end do
      end do
      if (taskid.eq.MASTER) write(17,1370) error
      if (taskid.eq.MASTER) write(17,1380)
      end do
*     **** produce CHECK FILE ****
      if (taskid.eq.MASTER) then
         close(17)
      end if


      return
 1350 FORMAT(/'******** band structure orthonormality **********')
 1360 FORMAT(2I3,2I3,'(',2E18.7,')')
 1370 FORMAT('ERROR = ',E18.7)
 1380 FORMAT(/'*************************************************')
      end


*     ***********************************
*     *             			*
*     *       Stfl_g_OrthoCheckMake	*
*     *             			*
*     ***********************************

      subroutine Stfl_g_OrthoCheckMake(ispin,ne,npack1,neall,nbrill,A)
      implicit none
      integer ispin,ne(2)
      integer npack1,neall,nbrill
      complex*16 A(npack1,neall,nbrill)

#include "mafdecls.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      logical value,gramschmidt
      integer nb,ms
      real*8 sum1,sum2,sum3
      integer sm2(2),sm3(2)


*     **** Ortho Check ****
      value = MA_push_get(mt_dbl,2*nbrill,
     >                     'sm2',sm2(2),sm2(1))
      value = value.and.MA_push_get(mt_dbl,2*nbrill,     
     >                     'sm3',sm3(2),sm3(1))
      if (.not. value) call errquit('out of stack memory',0,0)


      call Stfl_gg_trace(ispin,ne,npack1,neall,nbrill,
     >                        A,
     >                        A,
     >                        dbl_mb(sm2(1)))

      gramschmidt = .false.
      do nb=1,nbrill
      do ms=1,ispin
         sum1 = dble(ne(ms))
         sum2 = dbl_mb(sm2(1)+(ms-1)+2*(nb-1))
         if (dabs(sum2-sum1).gt.1.0d-10) then
           gramschmidt = .true.
           call Stfl_g_MakeOrtho(nb,npack1,ne(ms),
     >                           A(1,1+(ms-1)*ne(1),nb))
         end if
      end do
      end do

      if (gramschmidt) then
      call Stfl_gg_trace(ispin,ne,npack1,neall,nbrill,
     >                        A,
     >                        A,
     >                        dbl_mb(sm3(1)))

      call Parallel_taskid(taskid)
      do nb=1,nbrill
      do ms=1,ispin
         sum1 = dble(ne(ms))
         sum2 = dbl_mb(sm2(1)+(ms-1)+(nb-1)*2)
         sum3 = dbl_mb(sm3(1)+(ms-1)+(nb-1)*2)
         if (dabs(sum2-sum1).gt.1.0d-10) then
           if (taskid.eq.MASTER) then
            write(*,1400) ms,nb
            write(*,1401) dabs(sum2-sum1),sum2
            write(*,1402) dabs(sum3-sum1),sum3
           end if
         end if
      end do
      end do
      end if

      value = MA_pop_stack(sm3(2))
      value = value.and.MA_pop_stack(sm2(2))
      if (.not. value)
     > call errquit('cpsi_initialize:popping stack memory',0,0)

 1400 FORMAT('Warning - Orthogonalization performed, spin,zone:',I2,I4)
 1401 FORMAT(8x,'- error(before)=', E14.6,' (',E14.6,')')
 1402 FORMAT(8x,'- error(after)= ', E14.6,' (',E14.6,')')
      return
      end 



*     ***********************************
*     *             			*
*     *          Stfl_MakeOrtho		*
*     *             			*
*     ***********************************

      subroutine Stfl_g_MakeOrtho(nb,nfft3d,ne,psi)
      implicit none 
      integer    nb,nfft3d,ne
      complex*16 psi(nfft3d,ne)

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
      integer j,k
      real*8  w
      complex*16 wc

*     **** Gram-Schmidt ****
      do k=1,ne
         call Cram_cc_dot(nb,psi(1,k),psi(1,k),w)
         if (w.gt.1.0d-9) then
            w = 1.0d0/dsqrt(w)
         else !should never happen!!!!      
            call Parallel_taskid(taskid)
            if (taskid.eq.MASTER) then
             write(*,*) "Stfl_g_MakeOrtho: this should not happen"
             write(*,*) 
     >       "Stfl_g_MakeOrtho: <psi_k|psi_k> = 0.0, k=",k,w
            end if
            w = 1.0d0
         end if
         call Cram_c_SMul(nb,w,psi(1,k),psi(1,k))

         do j=k+1,ne
            call Cram_cc_zdot(nb,psi(1,k),psi(1,j),wc)
            wc = -wc
            call Cram_cc_zaxpy(nb,wc,psi(1,k),psi(1,j))
            !call Cram_cc_zdot(nb,psi(1,k),psi(1,j),wc)
         end do
      end do

      return
      end




*
* $Id: orthocheck.f,v 1.1 2003-04-07 21:12:13 windus Exp $
*

*     ***********************************
*     *             			*
*     *          OrthoCheck		*
*     *             			*
*     ***********************************

      subroutine OrthoCheck(ispin,ne,psi1)
      implicit none 
      integer ispin,ne(2)
      double complex psi1(*)


*    *** local variables ***

      integer MASTER,taskid
      parameter(MASTER=0)
      integer ms,npack1,i,j,ii,jj
      integer n1(2),n2(2)
      real*8  w,error
      character*255 full_filename

      call Parallel_taskid(taskid)
      call Pack_npack(1,npack1)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)

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

      error = 0.0d0
      do ms=1,ispin
         do i=n1(ms),n2(ms)
            ii = i-n1(ms)+1
            do j=i,n2(ms)
               jj = j-n1(ms)+1
               call Pack_cc_dot(1,psi1(1+(i-1)*npack1),
     >                            psi1(1+(j-1)*npack1),
     >                            w)

               if (i.eq.j) then
                   error = dabs(1.0d0-w)
               else
                   error = dabs(w)
               end if
             
               if (taskid.eq.MASTER) then
                  write(17,1360) ms,ii,jj,w
               end if
            end do
         end do
      end do
      if (taskid.eq.MASTER) write(17,1370) error
      if (taskid.eq.MASTER) write(17,1380)

*     **** produce CHECK FILE ****
      if (taskid.eq.MASTER) then
         close(17)
      end if


      return
 1350 FORMAT(/'******** orthonormality **********')
 1360 FORMAT(I3,2I3,E18.7)
 1370 FORMAT('ERROR = ',E18.7)
 1380 FORMAT(/'**********************************')
      end



*     ***********************************
*     *             			*
*     *          Grsm_MakeOrtho		*
*     *             			*
*     ***********************************

      subroutine Grsm_g_MakeOrtho(npack,ne,psi)
      implicit none 
      integer npack,ne
      double complex psi(npack,ne)

*     **** local variables ****
      integer j,k
      real*8  w

      do k=1,ne
         call Pack_cc_dot(1,psi(1,k),psi(1,k),w)
         w = 1.0d0/dsqrt(w)
         call Pack_c_SMul(1,w,psi(1,k),psi(1,k))

         do j=k+1,ne
            call Pack_cc_dot(1,psi(1,k),psi(1,j),w)
            w = -w
            call Pack_cc_daxpy(1,w,psi(1,k),psi(1,j))
         end do
      end do

      return
      end



*     ***********************************
*     *             			*
*     *          OrthoCheck_geo		*
*     *             			*
*     ***********************************

      subroutine OrthoCheck_geo(ispin,ne,psi1)
      implicit none 
      integer ispin,ne(2)
      double complex psi1(*)


*    *** local variables ***

      integer MASTER,taskid
      parameter(MASTER=0)
      integer ms,npack1,i,j,ii,jj
      integer n1(2),n2(2)
      real*8  w,error
      character*255 full_filename

      call Parallel_taskid(taskid)
      call Pack_npack(1,npack1)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)



     
*     **** check orthonormality ****
      if (taskid.eq.MASTER) then
         write(*,1350)
         write(*,*) ne,n1,n2
      end if

      error = 0.0d0
      do ms=1,ispin
         do i=n1(ms),n2(ms)
            ii = i-n1(ms)+1
            do j=i,n2(ms)
               jj = j-n1(ms)+1
               call Pack_cc_dot(1,psi1(1+(i-1)*npack1),
     >                            psi1(1+(j-1)*npack1),
     >                            w)

               if (i.eq.j) then
                   error = dabs(1.0d0-w)
               else
                   error = dabs(w)
               end if
             
               if (taskid.eq.MASTER) then
                  write(*,1360) ms,ii,jj,w
               end if
            end do
         end do
      end do
      if (taskid.eq.MASTER) write(*,1370) error
      if (taskid.eq.MASTER) write(*,1380)


      return
 1350 FORMAT(/'******** geo orthonormality **********')
 1360 FORMAT(I3,2I3,E18.7)
 1370 FORMAT('ERROR = ',E18.7)
 1380 FORMAT(/'**********************************')
      end


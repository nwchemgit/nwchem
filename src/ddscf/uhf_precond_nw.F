c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++ FA precond routines ++++++++++++++++++ START
      subroutine uhf_hessv_precond3(g_x, 
     &                              omega)
c
c     Author : Fredy W. Aquino, Northwestern University
c     Date   : 03-24-12
c Note.- slight modification of uhf_precond()
c --> Experimental (not published yet)

      implicit none
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cscf.fh"
#include "cuhf.fh"
      integer g_x
      double precision omega
c      
      integer nvir, nocc(2), iset
      integer i, j, ioff, ibase
      integer l_diag, k_diag, l_x, k_x
      integer g_f(2), ivec, nvec, gtype, vlen
      double precision diag, denominator, shift, dnrm
      integer ilo(2), ihi(2)
c
      logical oprint, olprint

      diag(i) = dbl_mb(k_diag + i - 1)
c     
      g_f(1) = cuhf_g_falpha
      g_f(2) = cuhf_g_fbeta
      nocc(1) = nalpha
      nocc(2) = nbeta

      oprint = util_print('precond',print_high)
      olprint = oprint .and. (ga_nodeid().eq.0)

      if (.not.ma_push_get(MT_DBL,nmo,'uhf: tmp',l_diag,k_diag))
     $     call errquit('uhf_precond: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nmo,'uhf: tmp',l_x,k_x))
     $     call errquit('uhf_hdiag_scale: cannot allocate',0, MA_ERR)

      call ga_sync()
      call ga_inquire(g_x, gtype, vlen, nvec)

c      shift = lshift - solveshift + 4.0d0*omega
c      shift = lshift + 4.0d0*omega
      shift = lshift + 2.0d0*omega
      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,a,": in g_x = ",i4,f24.8)')
     +      __FILE__,ivec,dnrm
          endif
        enddo
      endif

c           if (ga_nodeid().eq.0)
c     &      write(*,*) '------ g_x-BEF-scaling-uhf---START'
c            call ga_print(g_x)
c           if (ga_nodeid().eq.0)
c     &      write(*,*) '------ g_x-BEF-scaling-uhf---END'

      ibase = 1
      do iset = 1, 2
         nvir = nmo - nocc(iset)
         if (nvir .gt. 0) then
            call ga_get_diagonal(g_f(iset), dbl_mb(k_diag))
            do ivec = 1, nvec
               do i = ga_nodeid()+1, nocc(iset), ga_nnodes()
                  ioff = (i-1)*nvir + ibase
                  call ga_get(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $                 dbl_mb(k_x),1)
                  do j=1,nvir
                     denominator = 2.d0 * (diag(j+nocc(iset)) - diag(i))
     $                    + shift
                     if (denominator .lt. 0.1d0) denominator = 0.1d0
                     dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
                  end do
                  call ga_put(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $                 dbl_mb(k_x),1)
               end do
            enddo
         endif
         ibase = ibase + nocc(1)*(nmo-nocc(1))
      end do

c           if (ga_nodeid().eq.0)
c     &      write(*,*) '------ g_x-AFT-scaling-uhf---START'
c            call ga_print(g_x)
c           if (ga_nodeid().eq.0)
c     &      write(*,*) '------ g_x-AFT-scaling-uhf---END'

      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,a,": out g_x = ",i4,f24.8)')
     +      __FILE__,ivec,dnrm
          endif
        enddo
      endif

      if (.not. ma_pop_stack(l_x)) call errquit('uhf:pop x',0, MA_ERR)
      if (.not. ma_pop_stack(l_diag)) call errquit('uhf:pop ',0, MA_ERR)
      call ga_sync()

      end

      subroutine uhf_hessv_precond3_damp(
     &                              g_x, 
     &                              g_x_im,
     &                              omega,
     &                              gamwidth)
c
c     Author : Fredy W. Aquino, Northwestern University
c     Date   : 03-24-12
c Note.- modification of uhf_precond() following
c        rohf_hessv_precond3_damp()
c --> Experimental (not published yet)

      implicit none
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cscf.fh"
#include "cuhf.fh"
      integer g_x,g_x_im
      double precision omega,    ! [input] frequency
     &                 gamwidth, ! [input]
     &                 gamwidth1
c      
      integer nvir, nocc(2), iset
      integer i, j, ioff, ibase
      integer l_diag, k_diag, 
     &        l_x   , k_x,
     &        l_x_im, k_x_im
      integer g_f(2), ivec, nvec, gtype, vlen
      double precision diag, denominator, shift, dnrm
      integer ilo(2), ihi(2)
      logical oprint, olprint

      double precision tmpre,tmpim,
     &                 factor_re,factor_im,
     &                 val


      external scale_gx_jth

      diag(i) = dbl_mb(k_diag + i - 1)
     
      g_f(1) = cuhf_g_falpha
      g_f(2) = cuhf_g_fbeta
      nocc(1) = nalpha
      nocc(2) = nbeta

      oprint = util_print('precond',print_high)
      olprint = oprint .and. (ga_nodeid().eq.0)

      if (.not.ma_push_get(MT_DBL,nmo,'uhf: diag',l_diag,k_diag))
     $     call errquit('uhf_precond_damp: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nmo,'uhf: x',l_x,k_x))
     $     call errquit('uhf_precond_damp: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nmo,'uhf: x_im',l_x_im,k_x_im))
     $     call errquit('uhf_precond_damp: cannot allocate',0, MA_ERR)

      call ga_sync()
      call ga_inquire(g_x, gtype, vlen, nvec)

      gamwidth1= 2d0 * gamwidth
      shift    = lshift + 2.0d0*omega

c      if (ga_nodeid().eq.0) then
c       write(*,1) omega,gamwidth,gamwidth1,shift
c 1     format('(omega,gamwidth,gamwidth1,shift)=(',
c     &        f15.8,',',f15.8,',',f15.8,',',f15.8,')')
c      endif

      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,a,": in g_x = ",i4,f24.8)')
     +      __FILE__,ivec,dnrm
          endif
        enddo
      endif

      ibase = 1
      do iset = 1, 2
         nvir = nmo - nocc(iset)
         if (nvir .gt. 0) then
            call ga_get_diagonal(g_f(iset), dbl_mb(k_diag))
            do ivec = 1, nvec
               do i = ga_nodeid()+1, nocc(iset), ga_nnodes()
                  ioff = (i-1)*nvir + ibase
                  call ga_get(g_x   , ioff, ioff+nvir-1, ivec, ivec,
     $                        dbl_mb(k_x),1)
                  call ga_get(g_x_im, ioff, ioff+nvir-1, ivec, ivec,
     $                        dbl_mb(k_x_im),1)

c                  if (ga_nodeid().eq.0) then
c                   write(*,2) ioff,ioff+nvir-1,ivec
c  2                format('(ioff,ioff+nvir-1,ivec)=(',
c     &                    i4,',',i4,',',i4,')')
c                  endif

                  do j=1,nvir
                     denominator = 2.d0 * (diag(j+nocc(iset)) - diag(i))
     $                             + shift

c                     if (ga_nodeid().eq.0) then
c                      write(*,3) iset,ivec,i,j,
c     &                           dbl_mb(k_x  +j-1),dbl_mb(k_x_im+j-1),
c     &                           denominator,diag(j+nocc(iset)),diag(i),
c     &                           shift
c  3                   format('FA-BEF-scl:(iset,ivec,i,j,',
c     &                       'xj,xj_im,denominator,enj,eni,shift)=(',
c     &                       i4,',',i4,',',i4,',',i4,',',
c     &                       f15.8,',',f15.8,',',f15.8,',',
c     &                       f15.8,',',f15.8,',',f15.8,')')         
c                     endif
 
                    if (denominator .lt. 0.1d0) denominator = 0.1d0
c ========================================START
                     val= denominator**2 + gamwidth1**2        
                     factor_re =  denominator
                     factor_im = -gamwidth1  
                     tmpre = dbl_mb(k_x   +j-1)   
                     tmpim = dbl_mb(k_x_im+j-1)   
                     dbl_mb(k_x   +j-1)=tmpre*factor_re+tmpim*factor_im
                     dbl_mb(k_x_im+j-1)=tmpim*factor_re+tmpre*factor_im
                     dbl_mb(k_x   +j-1)=dbl_mb(k_x   +j-1)/val
                     dbl_mb(k_x_im+j-1)=dbl_mb(k_x_im+j-1)/val
c==========================================END

c                     if (ga_nodeid().eq.0) then
c                      write(*,33) iset,ivec,i,j,
c     &                            dbl_mb(k_x  +j-1),dbl_mb(k_x_im+j-1)
c 33                   format('FA-AFT-scl:(iset,ivec,i,j,xj,xj_im)=(',
c     &                       i4,',',i4,',',i4,',',i4,',',
c     &                       f15.8,',',f15.8,')')         
c                     endif
                  end do ! end-loop-j
                  call ga_put(g_x   , ioff, ioff+nvir-1, ivec, ivec,
     $                        dbl_mb(k_x)   ,1)
                  call ga_put(g_x_im, ioff, ioff+nvir-1, ivec, ivec,
     $                        dbl_mb(k_x_im),1)
               end do
            enddo
         endif
         ibase = ibase + nocc(1)*(nmo-nocc(1))
      end do

      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,a,": out g_x = ",i4,f24.8)')
     +      __FILE__,ivec,dnrm
          endif
        enddo
      endif
      if (.not. ma_pop_stack(l_x_im)) 
     &   call errquit('uhf:pop x-im',0, MA_ERR)
      if (.not. ma_pop_stack(l_x))    
     &   call errquit('uhf:pop x-re',0, MA_ERR)
      if (.not. ma_pop_stack(l_diag)) 
     &   call errquit('uhf:pop ' ,0, MA_ERR)

      call ga_sync()

      end

      subroutine uhf_hessv_precond3_damp_anti(
     &                              g_x, 
     &                              g_x_im,
     &                              omega,
     &                              gamwidth)
c
c     Author : Fredy W. Aquino, Northwestern University
c     Date   : 03-24-12
c
c Note.- modification of uhf_precond() following
c        rohf_hessv_precond3_damp()
c --> Experimental (not published yet)

      implicit none
#include "stdio.fh"
#include "util.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cscf.fh"
#include "cuhf.fh"
      integer g_x,g_x_im
      double precision omega,    ! [input] frequency
     &                 gamwidth, ! [input]
     &                 gamwidth1
c      
      integer nvir, nocc(2), iset
      integer i, j, ioff, ibase
      integer l_diag, k_diag, 
     &        l_x   , k_x,
     &        l_x_im, k_x_im
      integer g_f(2), ivec, nvec, gtype, vlen
      double precision diag, denominator, shift, dnrm
      integer ilo(2), ihi(2)
      logical oprint, olprint

      double precision tmpre,tmpim,
     &                 factor_re,factor_im,
     &                 val,val1,rval
      external scale_gx_jth

      diag(i) = dbl_mb(k_diag + i - 1)
     
      g_f(1) = cuhf_g_falpha
      g_f(2) = cuhf_g_fbeta
      nocc(1) = nalpha
      nocc(2) = nbeta

      oprint = util_print('precond',print_high)
      olprint = oprint .and. (ga_nodeid().eq.0)

      if (.not.ma_push_get(MT_DBL,nmo,'uhf: diag',l_diag,k_diag))
     $     call errquit('uhf_precond_damp: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nmo,'uhf: x',l_x,k_x))
     $     call errquit('uhf_precond_damp: cannot allocate',0, MA_ERR)
      if (.not.ma_push_get(MT_DBL,nmo,'uhf: x_im',l_x_im,k_x_im))
     $     call errquit('uhf_precond_damp: cannot allocate',0, MA_ERR)

      call ga_sync()
      call ga_inquire(g_x, gtype, vlen, nvec)

      gamwidth1= 2d0 * gamwidth
      shift    = lshift + 2.0d0*omega
      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,a,": in g_x = ",i4,f24.8)')
     +      __FILE__,ivec,dnrm
          endif
        enddo
      endif

      ibase = 1
      do iset = 1, 2
         nvir = nmo - nocc(iset)
         if (nvir .gt. 0) then
            call ga_get_diagonal(g_f(iset), dbl_mb(k_diag))
            do ivec = 1, nvec
               do i = ga_nodeid()+1, nocc(iset), ga_nnodes()
                  ioff = (i-1)*nvir + ibase
                  call ga_get(g_x   , ioff, ioff+nvir-1, ivec, ivec,
     $                        dbl_mb(k_x),1)
                  call ga_get(g_x_im, ioff, ioff+nvir-1, ivec, ivec,
     $                        dbl_mb(k_x_im),1)

                  do j=1,nvir
                     denominator = 2.d0 * (diag(j+nocc(iset)) - diag(i))
     $                             + shift
                     if (denominator .lt. 0.1d0) denominator = 0.1d0
c =======Do anti-transformation ========================START
                     val =(denominator**2 + gamwidth1**2)
                     val1=(denominator**2 - gamwidth1**2)
                     rval=val/val1
                     factor_re = denominator
                     factor_im = gamwidth1   
                     tmpre = dbl_mb(k_x   +j-1)   
                     tmpim = dbl_mb(k_x_im+j-1)   
                     dbl_mb(k_x   +j-1)=tmpre*factor_re+tmpim*factor_im
                     dbl_mb(k_x_im+j-1)=tmpim*factor_re+tmpre*factor_im
                     dbl_mb(k_x   +j-1)=dbl_mb(k_x   +j-1)*rval
                     dbl_mb(k_x_im+j-1)=dbl_mb(k_x_im+j-1)*rval
c =======Do anti-transformation ========================END
                  end do ! end-loop-j
                  call ga_put(g_x   , ioff, ioff+nvir-1, ivec, ivec,
     $                        dbl_mb(k_x)   ,1)
                  call ga_put(g_x_im, ioff, ioff+nvir-1, ivec, ivec,
     $                        dbl_mb(k_x_im),1)
               end do ! end-loop-i
            enddo ! end-loop-ivec
         endif
         ibase = ibase + nocc(1)*(nmo-nocc(1))
      enddo ! end-loop-iset

      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,a,": out g_x = ",i4,f24.8)')
     +      __FILE__,ivec,dnrm
          endif
        enddo
      endif
      if (.not. ma_pop_stack(l_x_im)) 
     &   call errquit('uhf:pop x-im',0, MA_ERR)
      if (.not. ma_pop_stack(l_x))    
     &   call errquit('uhf:pop x-re',0, MA_ERR)
      if (.not. ma_pop_stack(l_diag)) 
     &   call errquit('uhf:pop ' ,0, MA_ERR)

      call ga_sync()

      end
c++++++++++++++++++++ FA precond routines ++++++++++++++++++ END
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

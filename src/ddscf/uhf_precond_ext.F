c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++ FA precond routines ++++++++++++++++++ START
      subroutine uhf_hessv_precond3(g_x, 
     &                              omega)
c
c     Author : Fredy W. Aquino, Northwestern University
c     Date   : 03-24-12
c Note.- slight modification of uhf_precond()

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
      logical oprint, olprint, warning

      double precision small, sign
      parameter (small = 1d-5)

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

      shift = lshift + 2.0d0*omega ! jochen 02/2016

c      shift = lshift + omega ! FA-03-13-14

      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,a,": in g_x = ",i4,f24.8)')
     +           'uhf_hessv_precond3',ivec,dnrm
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
                  call ga_get(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $                 dbl_mb(k_x),1)
                  do j=1,nvir
                     denominator = 2.d0 * (diag(j+nocc(iset)) - diag(i))
     $                    + shift

c                   denominator = (diag(j+nocc(iset))-diag(i))+shift ! FA-03-13-14

c                    if (denominator .lt. 0.1d0) denominator = 0.1d0

                     if (abs(denominator).lt.small) then
                       sign = 1d0
                       if (denominator.lt.0d0) sign = -1d0
                       denominator = small * sign
                       if (ga_nodeid().eq.0 .and. .not.warning)
     &                   write (LuOut,9000) 
                       warning = .true.
                     endif
                     
                     dbl_mb(k_x+j-1) = dbl_mb(k_x+j-1) / denominator
c                   dbl_mb(k_x+j-1)=2.0d0* dbl_mb(k_x+j-1)/denominator ! FA-03-13-14
                  end do
                  call ga_put(g_x, ioff, ioff+nvir-1, ivec, ivec,
     $                 dbl_mb(k_x),1)
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
     +      'uhf_hessv_precond3',ivec,dnrm
          endif
        enddo
      endif

      if (.not. ma_pop_stack(l_x)) call errquit('uhf:pop x',0, MA_ERR)
      if (.not. ma_pop_stack(l_diag)) call errquit('uhf:pop ',0, MA_ERR)
      call ga_sync()

      return

 9000 format (/' *** WARNING: SMALL OCC-VIRT ORBITAL ENERGY DIFFERENCE',
     &  ' DETECTED'/' (Including frequency factor for response)'/
     &  ' This is going to lead to diverging respose properties'/
     &  ' I have adjusted the denominator to a value of 1E-5')
      
      end

c     =================================================================

      subroutine uhf_hessv_precond3_damp(
     &                              g_x, 
     &                              g_x_im,
     &                              omega,
     &                              gamwidth)
c
c     Author : Fredy W. Aquino, Northwestern University
c     Date   : 03-24-12
c Note.- modification of uhf_precond() following
c     rohf_hessv_precond3_damp()
      
c     jochen 02/16: applied a similar change here as in the
c     closed-shell counterpart in order to get clean convergence

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
     &                 gamwidth  ! [input]
c      
      integer nvir, nocc(2), iset
      integer i, j, ioff, ibase
      integer l_diag, k_diag, 
     &        l_x   , k_x,
     &        l_x_im, k_x_im
      integer g_f(2), ivec, nvec, gtype, vlen
      double precision diag, denominator, shift, dnrm
      integer ilo(2), ihi(2)
      logical oprint, olprint, warning

      double precision tmpre,tmpim,
     &                 factor_re,factor_im,
     &                 den, gam

      double precision small, sign
      parameter (small = 1d-5)

      external scale_gx_jth

      diag(i) = dbl_mb(k_diag + i - 1)

c     =================================================================

      warning = .false.
     
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

c     jochen 02/16: changed such that omega and gam
c     include factors of 2
      
      shift    = lshift + 2.0d0 * omega
      gam = 2.0d0 * gamwidth

      if (oprint) then
        do ivec = 1, nvec
          ilo(1) = 1
          ilo(2) = ivec
          ihi(1) = vlen
          ihi(2) = ivec
          call nga_normf_patch(g_x,ilo,ihi,dnrm)
          if (olprint) then
            write(LuOut,'(1x,a,": in g_x = ",i4,f24.8)')
     +           'uhf_hessv_precond3_damp',ivec,dnrm
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
                    
c                   denominator = (diag(j+nocc(iset))-diag(i))+shift
                    denominator =
     &                2.0d0 * (diag(j+nocc(iset))-diag(i))+shift
                    
c                    if (denominator .lt. 0.1d0) denominator = 0.1d0

                    if (abs(denominator).lt.small) then
                      sign = 1d0
                      if (denominator.lt.0d0) sign = -1d0
                      denominator = small * sign
                      if (ga_nodeid().eq.0 .and. .not.warning)
     &                  write (LuOut,9000) 
                      warning = .true.
                    endif
                    
                     den= denominator**2 + gam**2        
c                    factor_re =  2.0d0*denominator/den
c                    factor_im =  2.0d0*gamwidth/den
                     
                     factor_re =  denominator/den
                     factor_im =  gamwidth/den
                     
                     tmpre = dbl_mb(k_x   +j-1)   
                     tmpim = dbl_mb(k_x_im+j-1)   
                     dbl_mb(k_x   +j-1)=tmpre*factor_re+tmpim*factor_im
                     dbl_mb(k_x_im+j-1)=tmpim*factor_re-tmpre*factor_im 
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
     +           'uhf_hessv_precond3_damp',ivec,dnrm
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

      return

 9000 format (/' *** WARNING: SMALL OCC-VIRT ORBITAL ENERGY DIFFERENCE',
     &  ' DETECTED'/' (Including frequency factor for response)'/
     &  ' This is going to lead to diverging respose properties'/
     &  ' I have adjusted the denominator to a value of 1E-5')


      end
c++++++++++++++++++++ FA precond routines ++++++++++++++++++ END
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c $Id$

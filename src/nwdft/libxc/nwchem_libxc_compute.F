      subroutine nwchem_libxc_compute(nq,ex,ec,qwght,rho,delrho,ttau,
     $                                laprho,Amat,Amat2,Amat3,Cmat,
     $                                Cmat2,Cmat3,
     $                                Mmat,Lmat,func,
     $                                ksgrad,kske,kslap,ldew,
     $                                do_2nd, do_3rd)
      use, intrinsic :: iso_c_binding

#ifdef USE_LIBXC
      use xc_f03_lib_m
#endif
      implicit none

#include "libxc.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "cdft.fh"
#include "dft2drv.fh"
#include "dft3drv.fh"

      logical :: ksgrad, kske, kslap, ldew, do_2nd, do_3rd

      integer :: q,nq
      integer :: ifunc
      integer :: size1, size2
      integer :: lexc,kexc
      integer :: lrho,krho,lvrho,kvrho
      integer :: lsigma,ksigma,lvsigma,kvsigma
      integer :: ltau,ktau,lvtau,kvtau
      integer :: llapl,klapl,lvlapl,kvlapl

      integer :: kv2rho2,lv2rho2
      integer :: kv2sig2,lv2sig2
      integer :: kv2rhosig,lv2rhosig

      integer :: kv3rho3,lv3rho3
      integer :: kv3sig3,lv3sig3
      integer :: kv3rho2sig,lv3rho2sig
      integer :: kv3rhosig2,lv3rhosig2

      integer(c_int) :: polarized

      double precision :: ex,ec,func(nq),qwght(nq)
      double precision :: rho(nq,(ipol*(ipol+1))/2)
      double precision :: delrho(nq,3,ipol)
      double precision :: laprho(nq,ipol)
      double precision :: ttau(nq,ipol)

      double precision :: Amat(nq,ipol)
      double precision :: Cmat(nq,3)
      double precision :: Mmat(nq,ipol)
      double precision :: Lmat(nq,ipol)

      double precision :: Amat2(nq,NCOL_AMAT2)
      double precision :: Cmat2(nq,NCOL_CMAT2)

      double precision :: Amat3(nq,NCOL_AMAT3)
      double precision :: Cmat3(nq,NCOL_CMAT3)

#ifdef USE_LIBXC
      double precision :: fac 

      double precision, external :: ddot
      logical gga,mgga,dolap
      logical,external :: nwchem_libxc_family

#if LIBXC_NP_SIZE == 4
      integer*4 nqs
#elif LIBXC_NP_SIZE == 8
      integer(c_size_t) :: nqs
#else
#error "unexpected value for LIBXC_NP_SIZE"
#endif
      type(xc_f03_func_t) :: xcfunc      

      gga = .false.
      mgga = .false.
      dolap = .false.

      nqs = nq
      size1 = nq*ipol
      size2 = nq*(ipol*(ipol+1))/2

      if (.not.(ma_push_get(mt_dbl,size1,'rhoval',lrho,krho)))
     $  call errquit(" could not allocate",0,ma_err) 
      if (.not.(ma_push_get(mt_dbl,size1,'vrho',lvrho,kvrho)))
     $  call errquit(" could not allocate",0,ma_err) 
      if (.not.(ma_push_get(mt_dbl,nq,'exc',lexc,kexc)))
     $  call errquit(" could not allocate",0,ma_err)

      if (ksgrad) then
        if (.not.(ma_push_get(mt_dbl,size2,'sigma',lsigma,ksigma)))
     $    call errquit(" could not allocate",0,ma_err) 
        if (.not.(ma_push_get(mt_dbl,size2,'vsigma',lvsigma,kvsigma)))
     $    call errquit(" could not allocate",0,ma_err) 
      endif

      if (kske.or.kslap) then
          if (.not.(ma_push_get(mt_dbl,size2,'tau',ltau,ktau)))
     $      call errquit(" could not allocate",0,ma_err) 
          if (.not.(ma_push_get(mt_dbl,size2,'vtau',lvtau,kvtau)))
     $      call errquit(" could not allocate",0,ma_err) 
          if (.not.(ma_push_get(mt_dbl,size2,'lapl',llapl,klapl)))
     $      call errquit(" could not allocate",0,ma_err) 
          if (.not.(ma_push_get(mt_dbl,size2,'vlapl',lvlapl,kvlapl)))
     $      call errquit(" could not allocate",0,ma_err) 
      endif

      if (do_2nd .or. do_3rd) then
        if (.not.(ma_push_get(mt_dbl,3*nq,'v2rho2',lv2rho2,kv2rho2)))
     $    call errquit(" could not allocate v2rho2",0,ma_err)
        if (ksgrad) then
          if (.not.(ma_push_get(mt_dbl,6*nq,'v2rhosig',lv2rhosig,
     $                                                 kv2rhosig)))
     $      call errquit(" could not allocate v2rhosig",0,ma_err)
          if (.not.(ma_push_get(mt_dbl,6*nq,'v2sig2',lv2sig2,kv2sig2)))
     $      call errquit(" could not allocate v2sig2",0,ma_err)
        endif
      endif

      if (do_3rd) then
        if (.not.(ma_push_get(mt_dbl,4*nq,'v3rho3',lv3rho3,kv3rho3)))
     $    call errquit(" could not allocate v3rho3",0,ma_err)
        if (ksgrad) then
          if (.not.(ma_push_get(mt_dbl,9*nq,'v3rho2sig',lv3rho2sig,
     $                                                   kv3rho2sig)))
     $      call errquit(" could not allocate v3rho2sig",0,ma_err)
          if (.not.(ma_push_get(mt_dbl,12*nq,'v3rhosig2',lv3rhosig2,
     $                                                   kv3rhosig2)))
     $      call errquit(" could not allocate v3rhosig2",0,ma_err)
          if (.not.(ma_push_get(mt_dbl,10*nq,'v3sig3',lv3sig3,kv3sig3)))
     $      call errquit(" could not allocate v3sig3",0,ma_err)     
        endif
      endif

      if (kske) call dcopy(nq,ttau,1,dbl_mb(ktau),ipol)
      if (kslap) call dcopy(nq,laprho,1,dbl_mb(klapl),ipol)

      if (ipol.eq.1) then
        polarized = xc_unpolarized

        call dcopy(nq,rho,1,dbl_mb(krho),1)
        if (ksgrad) then
           call nwchem_libxc_util1(nq,dbl_mb(ksigma),delrho)
        endif
      else
        polarized = xc_polarized
        call dcopy(nq,rho(1,2),1,dbl_mb(krho),2)
        call dcopy(nq,rho(1,3),1,dbl_mb(krho+1),2)
        if (ksgrad) then
           call nwchem_libxc_util2(nq,dbl_mb(ksigma),delrho)
        endif
        if (kske) call dcopy(nq,ttau(1,2),1,dbl_mb(ktau+1),ipol)
        if (kslap) call dcopy(nq,laprho(1,2),1,dbl_mb(klapl+1),ipol)
      endif

      do ifunc=1,libxc_nfuncs

        fac = libxc_facts(ifunc)

        call xc_f03_func_init(xcfunc,libxc_funcs(ifunc),polarized)
        call xc_f03_func_set_dens_threshold(xcfunc, tol_rho)
#if XC_MAJOR_VERSION > 4
        call xc_f03_func_set_sigma_threshold(xcfunc, tol_rho)
        call xc_f03_func_set_zeta_threshold(xcfunc, 1d-10)
#endif
        select case(libxc_family(ifunc))
#if XC_MAJOR_VERSION > 4
        case (XC_FAMILY_LDA, XC_FAMILY_HYB_LDA)
#else
        case (XC_FAMILY_LDA)
#endif
          if ((.not.do_2nd) .and. (.not.do_3rd)) then
            call xc_f03_lda_exc_vxc(xcfunc,nqs,
     $            dbl_mb(krho),dbl_mb(kexc),dbl_mb(kvrho))
          elseif (.not.do_3rd) then
#if XC_MAJOR_VERSION > 4
            call xc_f03_lda_exc_vxc_fxc(xcfunc,nqs,dbl_mb(krho),
     $            dbl_mb(kexc),dbl_mb(kvrho),dbl_mb(kv2rho2))
#else
            call xc_f03_lda_exc_vxc(xcfunc,nqs,dbl_mb(krho),
     $                              dbl_mb(kexc),dbl_mb(kvrho))
            call xc_f03_lda_fxc(xcfunc,nqs,dbl_mb(krho),dbl_mb(kv2rho2))
#endif
          else
#if XC_MAJOR_VERSION > 4
            call xc_f03_lda_exc_vxc_fxc_kxc(xcfunc,nqs,dbl_mb(krho),
     $            dbl_mb(kexc),dbl_mb(kvrho),dbl_mb(kv2rho2),
     $            dbl_mb(kv3rho3))
#else
            call xc_f03_lda_exc_vxc(xcfunc,nqs,dbl_mb(krho),
     $                              dbl_mb(kexc),dbl_mb(kvrho))
            call xc_f03_lda_fxc(xcfunc,nqs,dbl_mb(krho),dbl_mb(kv2rho2))
            call xc_f03_lda_kxc(xcfunc,nqs,dbl_mb(krho),dbl_mb(kv3rho3))
#endif
          endif


        case (XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)

          gga = .true.

          if (iand(libxc_flags(ifunc),xc_flags_have_exc).eq.
     $                                xc_flags_have_exc) then

            if ((.not.do_2nd).and.(.not.do_3rd)) then
              call xc_f03_gga_exc_vxc(xcfunc,nqs,
     $                                dbl_mb(krho),dbl_mb(ksigma),
     $                                dbl_mb(kexc),
     $                                dbl_mb(kvrho),dbl_mb(kvsigma))
            elseif (.not.do_3rd) then
#if XC_MAJOR_VERSION > 4
              call xc_f03_gga_exc_vxc_fxc(xcfunc,nqs,
     $             dbl_mb(krho),dbl_mb(ksigma),dbl_mb(kexc),
     $             dbl_mb(kvrho),dbl_mb(kvsigma),dbl_mb(kv2rho2),
     $             dbl_mb(kv2rhosig),dbl_mb(kv2sig2))
#else
              call xc_f03_gga_exc_vxc(xcfunc,nqs,
     $                                dbl_mb(krho),dbl_mb(ksigma),
     $                                dbl_mb(kexc),
     $                                dbl_mb(kvrho),dbl_mb(kvsigma))
              call xc_f03_gga_fxc(xcfunc,nqs,
     $                                dbl_mb(krho),dbl_mb(ksigma),
     $                                dbl_mb(kv2rho2),dbl_mb(kv2rhosig),
     $                                dbl_mb(kv2sig2))
#endif
            else
#if XC_MAJOR_VERSION > 4
              call xc_f03_gga_exc_vxc_fxc_kxc(xcfunc,nqs,
     $             dbl_mb(krho),dbl_mb(ksigma),dbl_mb(kexc),
     $             dbl_mb(kvrho),dbl_mb(kvsigma),dbl_mb(kv2rho2),
     $             dbl_mb(kv2rhosig),dbl_mb(kv2sig2),
     $             dbl_mb(kv3rho3),dbl_mb(kv3rho2sig),
     $             dbl_mb(kv3rhosig2),dbl_mb(kv3sig3))
#else
              call xc_f03_gga_exc_vxc(xcfunc,nqs,
     $                                dbl_mb(krho),dbl_mb(ksigma),
     $                                dbl_mb(kexc),
     $                                dbl_mb(kvrho),dbl_mb(kvsigma))
              call xc_f03_gga_fxc(xcfunc,nqs,
     $                            dbl_mb(krho),dbl_mb(ksigma),
     $                            dbl_mb(kv2rho2),dbl_mb(kv2rhosig),
     $                            dbl_mb(kv2sig2))
              call xc_f03_gga_kxc(xcfunc,nqs,
     $                            dbl_mb(krho),dbl_mb(ksigma),
     $                            dbl_mb(kv3rho3),dbl_mb(kv3rho2sig),
     $                            dbl_mb(kv3rhosig2),dbl_mb(kv3sig3))
#endif
            endif
          else
            if ((.not.do_2nd).and.(.not.do_3rd)) then
              call xc_f03_gga_vxc(xcfunc,nqs,
     $                            dbl_mb(krho),dbl_mb(ksigma),
     $                            dbl_mb(kvrho),dbl_mb(kvsigma))
            elseif (.not.do_3rd) then
#if XC_MAJOR_VERSION > 4
              call xc_f03_gga_vxc_fxc(xcfunc,nqs,
     $             dbl_mb(krho),dbl_mb(ksigma),
     $             dbl_mb(kvrho),dbl_mb(kvsigma),dbl_mb(kv2rho2),
     $             dbl_mb(kv2rhosig),dbl_mb(kv2sig2))
#else
              call xc_f03_gga_vxc(xcfunc,nqs,
     $                            dbl_mb(krho),dbl_mb(ksigma),
     $                            dbl_mb(kvrho),dbl_mb(kvsigma))
              call xc_f03_gga_fxc(xcfunc,nqs,
     $                            dbl_mb(krho),dbl_mb(ksigma),
     $                            dbl_mb(kv2rho2),dbl_mb(kv2rhosig),
     $                            dbl_mb(kv2sig2))
#endif
            else
#if XC_MAJOR_VERSION > 4
              call xc_f03_gga_vxc_fxc_kxc(xcfunc,nqs,
     $             dbl_mb(krho),dbl_mb(ksigma),
     $             dbl_mb(kvrho),dbl_mb(kvsigma),dbl_mb(kv2rho2),
     $             dbl_mb(kv2rhosig),dbl_mb(kv2sig2),
     $             dbl_mb(kv3rho3),dbl_mb(kv3rho2sig),
     $             dbl_mb(kv3rhosig2),dbl_mb(kv3sig3))
#else
              call xc_f03_gga_vxc(xcfunc,nqs,
     $                            dbl_mb(krho),dbl_mb(ksigma),
     $                            dbl_mb(kvrho),dbl_mb(kvsigma))
              call xc_f03_gga_fxc(xcfunc,nqs,
     $                            dbl_mb(krho),dbl_mb(ksigma),
     $                            dbl_mb(kv2rho2),dbl_mb(kv2rhosig),
     $                            dbl_mb(kv2sig2))
              call xc_f03_gga_kxc(xcfunc,nqs,
     $                            dbl_mb(krho),dbl_mb(ksigma),
     $                            dbl_mb(kv3rho3),dbl_mb(kv3rho2sig),
     $                            dbl_mb(kv3rhosig2),dbl_mb(kv3sig3))
#endif
            endif
            call dfill(nq,0d0,dbl_mb(kexc),1)
          endif

        case (XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA)
          gga = .true.
          mgga = .true.
          dolap = iand(libxc_flags(ifunc),xc_flags_needs_laplacian).eq.
     $            xc_flags_needs_laplacian
          
          if (iand(libxc_flags(ifunc),xc_flags_have_exc).eq.
     $                                xc_flags_have_exc) then

            call xc_f03_mgga_exc_vxc(xcfunc,nqs,
     $       dbl_mb(krho),dbl_mb(ksigma),dbl_mb(klapl),dbl_mb(ktau),
     $       dbl_mb(kexc),
     $       dbl_mb(kvrho),dbl_mb(kvsigma),dbl_mb(kvlapl),dbl_mb(kvtau))
          else
            call xc_f03_mgga_vxc(xcfunc,nqs,
     $       dbl_mb(krho),dbl_mb(ksigma),dbl_mb(klapl),dbl_mb(ktau),
     $       dbl_mb(kvrho),dbl_mb(kvsigma),dbl_mb(kvlapl),dbl_mb(kvtau))
            call dfill(nq,0d0,dbl_mb(kexc),1)
          endif

        end select

        call xc_f03_func_end(xcfunc)
        call nwchem_libxc_util3(ldew,nq,fac,
     e       dbl_mb(kexc),func,qwght,rho)

        if (libxc_kind(ifunc).eq.xc_exchange .or. 
     $      libxc_kind(ifunc).eq.xc_exchange_correlation) then
          ex = ex + fac*ddot(nq,dbl_mb(kexc),1,rho,1)
        elseif(libxc_kind(ifunc).eq.xc_correlation) then
          ec = ec + fac*ddot(nq,dbl_mb(kexc),1,rho,1)
        endif

        call daxpy(nq,fac,dbl_mb(kvrho),ipol,Amat,1)
        if (mgga) call daxpy(nq,0.5d0*fac,dbl_mb(kvtau),ipol,Mmat,1)
        if (dolap) call daxpy(nq,fac,dbl_mb(kvlapl),ipol,Lmat,1)
        
        if ((polarized.eq.xc_unpolarized) .and. (gga)) then
          call daxpy(nq,fac,dbl_mb(kvsigma),1,Cmat(1,D1_GAA),1)
          call daxpy(nq,2d0*fac,dbl_mb(kvsigma),1,Cmat(1,D1_GAB),1)
        elseif (polarized.eq.xc_polarized) then
          call daxpy(nq,fac,dbl_mb(kvrho+1),2,Amat(1,2),1)
          if (mgga) call daxpy(nq,0.5d0*fac,dbl_mb(kvtau+1),2,
     $                         Mmat(1,2),1)
          if (dolap) call daxpy(nq,fac,dbl_mb(kvlapl+1),2,Lmat(1,2),1)
          if (gga) then
            call daxpy(nq,fac,dbl_mb(kvsigma),3,Cmat(1,D1_GAA),1)
            call daxpy(nq,fac,dbl_mb(kvsigma+1),3,Cmat(1,D1_GAB),1)
            call daxpy(nq,fac,dbl_mb(kvsigma+2),3,Cmat(1,D1_GBB),1)
          endif
        endif

        if (do_2nd .or. do_3rd) then
          if (polarized.eq.xc_unpolarized) then
            call daxpy(nq,fac,dbl_mb(kv2rho2),1,Amat2(1,D2_RA_RA),1)
            call daxpy(nq,fac,dbl_mb(kv2rho2),1,Amat2(1,D2_RA_RB),1)
            if (gga) then
              call daxpy(nq,fac,dbl_mb(kv2rhosig),1,
     $                   Cmat2(1,D2_RA_GAA),1)
              call daxpy(nq,2d0*fac,dbl_mb(kv2rhosig),1,
     $                   Cmat2(1,D2_RA_GAB),1)
              call daxpy(nq,fac,dbl_mb(kv2rhosig),1,
     $                   Cmat2(1,D2_RA_GBB),1)
              call daxpy(nq,fac,dbl_mb(kv2sig2),1,
     $                   Cmat2(1,D2_GAA_GAA),1)
              call daxpy(nq,2d0*fac,dbl_mb(kv2sig2),1,
     $                   Cmat2(1,D2_GAA_GAB),1)
              call daxpy(nq,fac,dbl_mb(kv2sig2),1,
     $                   Cmat2(1,D2_GAA_GBB),1)
              call daxpy(nq,4d0*fac,dbl_mb(kv2sig2),1,
     $                   Cmat2(1,D2_GAB_GAB),1)
            endif
          else
            call daxpy(nq,fac,dbl_mb(kv2rho2),3,Amat2(1,D2_RA_RA),1)
            call daxpy(nq,fac,dbl_mb(kv2rho2+1),3,Amat2(1,D2_RA_RB),1)
            call daxpy(nq,fac,dbl_mb(kv2rho2+2),3,Amat2(1,D2_RB_RB),1)
            if (gga) then
              call daxpy(nq,fac,dbl_mb(kv2rhosig),6,
     $                   Cmat2(1,D2_RA_GAA),1)
              call daxpy(nq,fac,dbl_mb(kv2rhosig+1),6,
     $                   Cmat2(1,D2_RA_GAB),1)
              call daxpy(nq,fac,dbl_mb(kv2rhosig+2),6,
     $                   Cmat2(1,D2_RA_GBB),1)
              call daxpy(nq,fac,dbl_mb(kv2rhosig+3),6,
     $                   Cmat2(1,D2_RB_GAA),1)
              call daxpy(nq,fac,dbl_mb(kv2rhosig+4),6,
     $                   Cmat2(1,D2_RB_GAB),1)
              call daxpy(nq,fac,dbl_mb(kv2rhosig+5),6,
     $                   Cmat2(1,D2_RB_GBB),1)

              call daxpy(nq,fac,dbl_mb(kv2sig2),6,
     $                   Cmat2(1,D2_GAA_GAA),1)
              call daxpy(nq,fac,dbl_mb(kv2sig2+1),6,
     $                   Cmat2(1,D2_GAA_GAB),1)
              call daxpy(nq,fac,dbl_mb(kv2sig2+2),6,
     $                   Cmat2(1,D2_GAA_GBB),1)
              call daxpy(nq,fac,dbl_mb(kv2sig2+3),6,
     $                   Cmat2(1,D2_GAB_GAB),1)
              call daxpy(nq,fac,dbl_mb(kv2sig2+4),6,
     $                   Cmat2(1,D2_GAB_GBB),1)
              call daxpy(nq,fac,dbl_mb(kv2sig2+5),6,
     $                   Cmat2(1,D2_GBB_GBB),1)
            endif
          endif
        endif

        if (do_3rd) then
          if (polarized.eq.xc_unpolarized) then
           call daxpy(nq,fac,dbl_mb(kv3rho3),1,Amat3(1,D3_RA_RA_RA),1)
           call daxpy(nq,fac,dbl_mb(kv3rho3),1,Amat3(1,D3_RA_RA_RB),1)
           call daxpy(nq,fac,dbl_mb(kv3rho3),1,Amat3(1,D3_RA_RB_RB),1)
           if (gga) then
            call daxpy(nq,fac,dbl_mb(kv3rho2sig),1,
     $                 Cmat3(1,D3_RA_RA_GAA),1)
            call daxpy(nq,2d0*fac,dbl_mb(kv3rho2sig),1,
     $                 Cmat3(1,D3_RA_RA_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3rho2sig),1,
     $                 Cmat3(1,D3_RA_RA_GBB),1)

            call daxpy(nq,fac,dbl_mb(kv3rho2sig),1,
     $                 Cmat3(1,D3_RA_RB_GAA),1)
            call daxpy(nq,2d0*fac,dbl_mb(kv3rho2sig),1,
     $                 Cmat3(1,D3_RA_RB_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3rho2sig),1,
     $                 Cmat3(1,D3_RA_RB_GBB),1)

            call daxpy(nq,fac,dbl_mb(kv3rhosig2),1,
     $                 Cmat3(1,D3_RA_GAA_GAA),1)
            call daxpy(nq,2d0*fac,dbl_mb(kv3rhosig2),1,
     $                 Cmat3(1,D3_RA_GAA_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2),1,
     $                 Cmat3(1,D3_RA_GAA_GBB),1)
            call daxpy(nq,4d0*fac,dbl_mb(kv3rhosig2),1,
     $                 Cmat3(1,D3_RA_GAB_GAB),1)
            call daxpy(nq,2d0*fac,dbl_mb(kv3rhosig2),1,
     $                 Cmat3(1,D3_RA_GAB_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2),1,
     $                 Cmat3(1,D3_RA_GBB_GBB),1)

            call daxpy(nq,fac,dbl_mb(kv3sig3),1,
     $                 Cmat3(1,D3_GAA_GAA_GAA),1)
            call daxpy(nq,2d0*fac,dbl_mb(kv3sig3),1,
     $                 Cmat3(1,D3_GAA_GAA_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3sig3),1,
     $                 Cmat3(1,D3_GAA_GAA_GBB),1)
            call daxpy(nq,4d0*fac,dbl_mb(kv3sig3),1,
     $                 Cmat3(1,D3_GAA_GAB_GAB),1)
            call daxpy(nq,2d0*fac,dbl_mb(kv3sig3),1,
     $                 Cmat3(1,D3_GAA_GAB_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3sig3),1,
     $                 Cmat3(1,D3_GAA_GBB_GBB),1)
            call daxpy(nq,8d0*fac,dbl_mb(kv3sig3),1,
     $                 Cmat3(1,D3_GAB_GAB_GAB),1)
           endif
          elseif (polarized.eq.xc_polarized) then
           call daxpy(nq,fac,dbl_mb(kv3rho3),4,Amat3(1,D3_RA_RA_RA),1)
           call daxpy(nq,fac,dbl_mb(kv3rho3+1),4,Amat3(1,D3_RA_RA_RB),1)
           call daxpy(nq,fac,dbl_mb(kv3rho3+2),4,Amat3(1,D3_RA_RB_RB),1)
           call daxpy(nq,fac,dbl_mb(kv3rho3+3),4,Amat3(1,D3_RB_RB_RB),1)
           if (gga) then
            call daxpy(nq,fac,dbl_mb(kv3rho2sig),9,
     $                 Cmat3(1,D3_RA_RA_GAA),1)
            call daxpy(nq,fac,dbl_mb(kv3rho2sig+1),9,
     $                 Cmat3(1,D3_RA_RA_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3rho2sig+2),9,
     $                 Cmat3(1,D3_RA_RA_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3rho2sig+3),9,
     $                 Cmat3(1,D3_RA_RB_GAA),1)
            call daxpy(nq,fac,dbl_mb(kv3rho2sig+4),9,
     $                 Cmat3(1,D3_RA_RB_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3rho2sig+5),9,
     $                 Cmat3(1,D3_RA_RB_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3rho2sig+6),9,
     $                 Cmat3(1,D3_RB_RB_GAA),1)
            call daxpy(nq,fac,dbl_mb(kv3rho2sig+7),9,
     $                 Cmat3(1,D3_RB_RB_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3rho2sig+8),9,
     $                 Cmat3(1,D3_RB_RB_GBB),1)

            call daxpy(nq,fac,dbl_mb(kv3rhosig2),12,
     $                 Cmat3(1,D3_RA_GAA_GAA),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2+1),12,
     $                 Cmat3(1,D3_RA_GAA_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2+2),12,
     $                 Cmat3(1,D3_RA_GAA_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2+3),12,
     $                 Cmat3(1,D3_RA_GAB_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2+4),12,
     $                 Cmat3(1,D3_RA_GAB_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2+5),12,
     $                 Cmat3(1,D3_RA_GBB_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2+6),12,
     $                 Cmat3(1,D3_RB_GAA_GAA),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2+7),12,
     $                 Cmat3(1,D3_RB_GAA_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2+8),12,
     $                 Cmat3(1,D3_RB_GAA_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2+9),12,
     $                 Cmat3(1,D3_RB_GAB_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2+10),12,
     $                 Cmat3(1,D3_RB_GAB_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3rhosig2+11),12,
     $                 Cmat3(1,D3_RB_GBB_GBB),1)

            call daxpy(nq,fac,dbl_mb(kv3sig3),10,
     $                 Cmat3(1,D3_GAA_GAA_GAA),1)
            call daxpy(nq,fac,dbl_mb(kv3sig3+1),10,
     $                 Cmat3(1,D3_GAA_GAA_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3sig3+2),10,
     $                 Cmat3(1,D3_GAA_GAA_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3sig3+3),10,
     $                 Cmat3(1,D3_GAA_GAB_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3sig3+4),10,
     $                 Cmat3(1,D3_GAA_GAB_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3sig3+5),10,
     $                 Cmat3(1,D3_GAA_GBB_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3sig3+6),10,
     $                 Cmat3(1,D3_GAB_GAB_GAB),1)
            call daxpy(nq,fac,dbl_mb(kv3sig3+7),10,
     $                 Cmat3(1,D3_GAB_GAB_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3sig3+8),10,
     $                 Cmat3(1,D3_GAB_GBB_GBB),1)
            call daxpy(nq,fac,dbl_mb(kv3sig3+9),10,
     $                 Cmat3(1,D3_GBB_GBB_GBB),1)
           endif
          endif
        endif

      enddo

      if (.not.ma_chop_stack(lrho))
     $    call errquit(" could not chop stack",0,ma_err)
#endif
      end subroutine
      subroutine nwchem_libxc_util1(nq,sigma,delrho)
      implicit none
      integer nq
      double precision sigma(*)
      double precision delrho(nq,3,*)
c
      integer q
c
      do q=1,nq
         sigma(q) = delrho(q,1,1)**2 +
     $        delrho(q,2,1)**2 +
     $        delrho(q,3,1)**2
      enddo
      end subroutine
      subroutine nwchem_libxc_util2(nq,sigma,delrho)
      implicit none
      integer nq
      double precision sigma(3,*)
      double precision delrho(nq,3,*)
c
      integer q
c
      do q=1,nq
         sigma(1,q) = delrho(q,1,1)**2 +
     $        delrho(q,2,1)**2 +
     $        delrho(q,3,1)**2
         sigma(2,q) = delrho(q,1,1)*delrho(q,1,2) +
     $        delrho(q,2,1)*delrho(q,2,2) +
     $        delrho(q,3,1)*delrho(q,3,2)
         sigma(3,q) = delrho(q,1,2)**2 +
     $        delrho(q,2,2)**2 +
     $        delrho(q,3,2)**2
      enddo
      end subroutine
      subroutine nwchem_libxc_util3(ldew,nq,fac,
     e     exc,func,qwght,rho)
      implicit none
      logical ldew
      integer nq
      double precision fac
      double precision exc(nq),qwght(nq),
     r     rho(nq,*),func(nq)
c
      integer q
c
      if(ldew) then
         do q=1,nq
            func(q) = func(q) + fac*exc(q)*rho(q,1)
         enddo
      endif
      do q=1,nq
         exc(q) = exc(q)*qwght(q)
      enddo
      end subroutine

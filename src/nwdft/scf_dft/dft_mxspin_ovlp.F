      Subroutine dft_mxspin_ovlp(basis, noc, g_alpha, g_beta, g_tmp)

C$Id: dft_mxspin_ovlp.F,v 1.7 2002-01-04 18:10:05 d3k958 Exp $
      Implicit none
      integer basis
      integer g_s              ! overlap
      integer g_alpha          ! alpha eigenvecs [input]
      integer g_beta           ! beta eigenvecs [input]
      integer g_tmp            ! scratch space

#include "bas.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "stdio.fh"
c
c     local
c
      integer me,nproc
c
      integer i,j,jbig,n,ichunks,nbe,nend
      integer k_tmpr1, l_tmpr1, k_tmpr2, l_tmpr2, k_tmpi1, l_tmpi1
      integer g_ss, g_vt, g_u, g_t, k_vals, l_vals, noc(2), g_alphaT
      integer nalp, g_ualpha, k_unp, l_unp
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      double precision prodbig, prodtmp
      double precision alp_thresh, eval
c
      me=ga_nodeid()
      nproc=ga_nnodes()
c     
      if ( .not. bas_numbf(basis,nbf) ) 
     &     call errquit(' exiting in max_ovl',1)
c
c      call gacrea(geom,basis,nbf,nbf,'ovl ',g_s,'atom')
      g_s = ga_create_atom_blocked(geom, basis, 'AO ovl')
c
      if(.not.MA_Push_Get(MT_Dbl,nbf,'real_tmp1',l_tmpr1, k_tmpr1))
     &     call errquit('dft_mxspin_ovlp: cannot allocate real_tmp1',0)
      if(.not.MA_Push_Get(MT_Dbl,nbf,'real_tmp2',l_tmpr2, k_tmpr2))
     &     call errquit('dft_mxspin_ovlp: cannot allocate real_tmp2',0)
      if(.not.MA_Push_Get(MT_Int,nbf,'int_tmp1',l_tmpi1, k_tmpi1))
     &     call errquit('dft_mxspin_ovlp: cannot allocate int_tmp1',0)
c
      call ga_zero(g_s)
      call int_1e_ga(basis,basis,g_s,'overlap',.false.)
c
c     Compute matrix mult (C_alpha)T * S * C_beta = S`
c     
      call ga_dgemm('T','N',nbf,nbf,nbf,1.d0,g_alpha,g_s,0.d0,g_tmp)
      call ga_dgemm('N','N',nbf,nbf,nbf,1.d0,g_tmp,g_beta,0.d0,g_s)
c
c
c-----------------------------------------------------------------
c Corresponding Orbital Tansformation
c
c COT diagonalizes Sab*SabT (the S' calculated above). 
c After COT, the overlap is maximized and it is clear which
c alpha orbitals do not have beta partners. For alpha orbitals
c with partners, the overlap is 1.0 (or very nearly). For
c un-partnered orbitals, the overlap is zero.
c See J. Chem. Phys. (1967) 47, 1936.
c
c
c calculate Sab*SabT
c 
      if(.not.ga_create(mt_dbl,noc(1),noc(1),'SS',0,0,g_ss))
     $        call errquit('ga_create failed', g_ss)
c
      call ga_dgemm('N','T',noc(1),noc(1),noc(2),1.d0,g_s,g_s,0.d0,g_ss)
c
c create arrays needed for SVD
c
      if(.not.ga_create(mt_dbl,noc(1),noc(1),'u',0,0,g_u))
     $        call errquit('ga_create failed', g_u)
c
      if(.not.ga_create(mt_dbl,noc(1),noc(1),'vt',0,0,g_vt))
     $        call errquit('ga_create failed', g_vt)
c
      if(.not.MA_Push_Get(MT_dbl,noc(1),'values',l_vals, k_vals))
     &     call errquit('dft_mxspin_ovlp: cannot allocate values',0)
c
c perform SVD on Sab*SabT to determine unpaired alpha MO's
c
      call ga_svd_seq(g_ss,g_u,g_vt,dbl_mb(k_vals))   
c
      if (.not. ga_destroy(g_vt)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_vt', 0)
      if (.not. ga_destroy(g_ss)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_ss', 0)
c
c      write(6,*)'SVD eigenvalues'
c      write(6,*) (dbl_mb(k_vals+i),i=0,noc(1)-1)
c
c calculate transformed alpha vectors, alphaT
c
      if(.not.ga_create(mt_dbl,noc(1),noc(1),'t',0,0,g_t))
     $        call errquit('ga_create failed', g_t)
c
      call ga_zero(g_t)
c
      do i = 1,noc(1)
       call ga_put(g_t,i,i,i,i,dbl_mb(k_vals+i-1),k_vals)
      enddo
c
c      call ga_print(g_t)
c
      if(.not.ga_create(mt_dbl,nbf,noc(1),'alphaT',0,0,g_alphaT))
     $        call errquit('ga_create failed', g_alphaT)
c
      call ga_dgemm('N','N',nbf,noc(1),noc(1),1.d0
     &     ,g_alpha,g_u,0.d0,g_alphaT)
c
      if (.not. ga_destroy(g_u)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_u', 0)
c
c      call ga_print(g_alphaT)
c
      if (.not. ga_destroy(g_t)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_t', 0)
c
c     create array containing only alpha MO's which don't
c     have beta partners
c
      nalp = 0
      alp_thresh = 1.0d-10
c
      if(.not.MA_Push_Get(MT_int,noc(1),'unpaired',l_unp, k_unp))
     &     call errquit('dft_mxspin_ovlp: cannot allocate values',0)
c
      do i = 1,noc(1)
        eval = dbl_mb(k_vals+i-1) 
        if (dabs(eval).LT.alp_thresh) then
          nalp = nalp + 1
          int_mb(k_unp+i-1) = 1
        else
          int_mb(k_unp+i-1) = 0    
        endif
      enddo
c
c      write(6,*) 'paired/unpaired alpha orbitals'
c      write(6,*) (int_mb(k_unp+i-1),i=1,noc(1))
c
      if(.not.ga_create(mt_dbl,nbf,nalp,'unp alphaT',0,0,g_ualpha))
     $        call errquit('ga_create failed', g_ualpha)
c
      nalp=0
      do i = 1,noc(1)
       if(int_mb(k_unp+i-1).EQ.1) then
         nalp = nalp + 1
         call ga_copy_patch('N',g_alphaT,1,nbf,i,i
     &         ,g_ualpha,1,nbf,nalp,nalp)
        endif
      enddo
c
      if (.not. ga_destroy(g_alphaT)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_alphaT', 0)
c
c      call ga_print(g_ualpha)                 
c
      if(.not.MA_Pop_Stack(l_unp))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0)
c  
      if(.not.MA_Pop_Stack(l_vals))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0)
c
c  resulting orbitals are printed after the alpha/beta overlap 
c
c-----------------------------------------------------------------
c
      if(me.eq.0) then
        write(LuOut, 9996)
c        call ga_print(g_s)
      endif
      if(me.eq.0) then
        jbig = 1         ! take care of compiler warnings
        do i = 1, nbf
c
c         get row of g_s 
c
          call ga_get(g_s,i,i,1,nbf,DBL_MB(k_tmpr1),1)
          prodbig=0.0d0
          do j = 1, nbf
            prodtmp = abs(dbl_mb(k_tmpr1+j-1))
            if(prodtmp.gt.prodbig) then
              prodbig = prodtmp
              jbig = j
            endif
          enddo
          dbl_mb(k_tmpr2+i-1) = prodbig
          int_mb(k_tmpi1+i-1) = jbig
        enddo
        ichunks = nbf/10
        if(ichunks*10.lt.nbf)ichunks = ichunks +1
        do i = 1, ichunks
          nbe = 10*(i-1) + 1
          nend = nbe + 9
          if(nend.gt.nbf)nend = nbf
          write(LuOut,9997)(n,n=nbe,nend)
          write(LuOut,9998)(int_mb(k_tmpi1+n-1),n=nbe,nend)
          write(LuOut,9999)(dbl_mb(k_tmpr2+n-1),n=nbe,nend)
        enddo
      endif
 9996 format(/,1x,'  alpha - beta orbital overlaps ',/,
     &         1x,'  ----------------------------- ',/)
 9997 format(/,1x,'  alpha ',10(1x,i5,1x))
 9998 format(  1x,'   beta ',10(1x,i5,1x))
 9999 format(  1x,'overlap ',10(f7.3),/)
c
c
c -------------------------------------------------------------
c  print alpha orbitals without beta partners which resulted
c  from COT above
c
       call movecs_print_anal(basis, 1, nalp, 0.10d0, g_ualpha,
     & 'Alpha Orbitals without Beta Partners (after COT) ', 
     &   .false., 0.0 ,.false., 0 , .false., 0 )
c
      if (.not. ga_destroy(g_ualpha)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_ualpha', 0)
c
c -------------------------------------------------------------
c      call gadest(g_s)
      if (.not. ga_destroy(g_s)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_s', 0)
c
      if(.not.MA_Pop_Stack(l_tmpi1))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0)
      if(.not.MA_Pop_Stack(l_tmpr2))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0)
      if(.not.MA_Pop_Stack(l_tmpr1))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0)
c
      return
      end

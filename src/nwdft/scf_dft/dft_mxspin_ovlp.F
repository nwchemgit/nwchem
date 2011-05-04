      Subroutine dft_mxspin_ovlp(nao,nmo,basis,g_alpha,g_beta,g_tmp)

C$Id$
      Implicit none
#include "errquit.fh"
      integer basis
      integer g_s              ! overlap
      integer g_alpha          ! alpha eigenvecs [input]
      integer g_beta           ! beta eigenvecs [input]
      integer g_tmp            ! scratch space
      integer nao              ! # of basis functions
      integer nmo              ! # of molecular orbitals

#include "bas.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "util.fh"   
c
c     local
c
      integer me,nproc
c
      integer i,j,jbig,n,ichunks,nbe,nend
      integer k_tmpr1, l_tmpr1, k_tmpr2, l_tmpr2, k_tmpi1, l_tmpi1
      integer g_ss, g_vt, g_u, g_t, k_vals, l_vals, g_alphaT
      integer nalp, g_ualpha, k_unp, l_unp, l_part, k_part
      integer nct, l_non, k_non, ind
c     
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked
c
      double precision prodbig, prodtmp
      double precision alp_thresh, eval, ovlmax, ovl, absovl
c
      me=ga_nodeid()
      nproc=ga_nnodes()
c     
      g_s = ga_create_atom_blocked(geom, basis, 'AO ovl')
c
      if(.not.MA_Push_Get(MT_Dbl,nao,'real_tmp1',l_tmpr1, k_tmpr1))
     &     call errquit('dft_mxspin_ovlp: cannot allocate real_tmp1',0,
     &       MA_ERR)
      if(.not.MA_Push_Get(MT_Dbl,nao,'real_tmp2',l_tmpr2, k_tmpr2))
     &     call errquit('dft_mxspin_ovlp: cannot allocate real_tmp2',0,
     &       MA_ERR)
      if(.not.MA_Push_Get(MT_Int,nao,'int_tmp1',l_tmpi1, k_tmpi1))
     &     call errquit('dft_mxspin_ovlp: cannot allocate int_tmp1',0,
     &       MA_ERR)
c
      call ga_sync()
      call ga_zero(g_s)
      call int_1e_ga(basis,basis,g_s,'overlap',.false.)
c
c     Compute matrix mult (C_alpha)T * S * C_beta = S`
c     
      call ga_dgemm('T','N',nmo,nao,nao,1.d0,g_alpha,g_s,0.d0,g_tmp)
      call ga_dgemm('N','N',nmo,nmo,nao,1.d0,g_tmp,g_beta,0.d0,g_s)
c
c
      if(me.eq.0) then
        write(LuOut, 9996)
c        call ga_print(g_s)
      endif
      if(me.eq.0) then
        jbig = 1         ! take care of compiler warnings
        do i = 1, nmo
c
c         get row of g_s 
c
          call ga_get(g_s,i,i,1,nmo,DBL_MB(k_tmpr1),1)
          prodbig=0.0d0
          do j = 1, nmo
            prodtmp = abs(dbl_mb(k_tmpr1+j-1))
            if(prodtmp.gt.prodbig) then
              prodbig = prodtmp
              jbig = j
            endif
          enddo
          dbl_mb(k_tmpr2+i-1) = prodbig
          int_mb(k_tmpi1+i-1) = jbig
        enddo
        ichunks = nmo/10
        if(ichunks*10.lt.nbf)ichunks = ichunks +1
        do i = 1, ichunks
          nbe = 10*(i-1) + 1
          nend = nbe + 9
          if(nend.gt.nmo)nend = nmo
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
c Print some useful information
c about the un-paired alpha orbitals.
c
      if(.not.util_print('alpha partner info', print_high)) goto 2001
      if(noc(2).eq.0) goto 2001
c
c Find alpha orbital which best overlaps with each beta
c orbital, then the remaining alpha orbitals are the un-partnered
c ones:
c
      if(.not.MA_Push_Get(MT_int,noc(2),'partners',l_part, k_part))
     &    call errquit('dft_mxspin_ovlp: cannot allocate values',0,
     &       MA_ERR)
c
c find alpha partners for the beta orbitals:
c
      do j = 1,noc(2)
        ovlmax = 0.0d0
      do i = 1,noc(1)
        call ga_get(g_s,i,i,j,j,ovl,1)
        absovl = dabs(ovl)
        if (absovl.gt.ovlmax) then
          ovlmax = absovl
          int_mb(k_part+j-1) = i
        endif
      enddo !i
      enddo !j
c
      if(me.eq.0) then
      write(luout,*)
      write(luout,*)'ALPHA/BETA PARTNERS'
      endif
c
      call ga_sync()
      if(me.eq.0) then
        ichunks = noc(2)/10
        if(ichunks*10.lt.noc(2))ichunks = ichunks +1
        do i = 1, ichunks
          nbe = 10*(i-1) + 1
          nend = nbe + 9
          if(nend.gt.noc(2))nend = noc(2)
          write(LuOut,9997)(int_mb(k_part+n-1),n=nbe,nend) 
          write(LuOut,9998)(n,n=nbe,nend)
        enddo
        write(luout,*)
      endif
c
c determine which alpha orbitals are left over:
c
      call ga_sync()
      if(.not.MA_Push_Get(MT_int,noc(1),'non',l_non, k_non))
     &    call errquit('dft_mxspin_ovlp: cannot allocate values',0,
     &       MA_ERR)
c
      nct=1
      do 160 i = 1,noc(1)
      do 150 j = 1,noc(2)
        if(int_mb(k_part+j-1).EQ.i) goto 160
c          if(nct.le.noc(1)-noc(2)) then
            int_mb(k_non+nct-1) = i
c          endif 
  150 continue    
        nct = nct + 1
  160 continue
c
      nct=nct-1
c
      if(me.eq.0) then
      if (nct.Ge.1) then
         write(luout,9995) nct, (int_mb(k_non+n-1),n=1,nct)
      endif 
      endif
c
 9995 format(/,1X,'THERE ARE ',i3,' UN-PARTNERED ALPHA ORBITALS :',20I4)  
c
c print the un-partnered alpha orbitals
c
      if (nct.GE.1) then  !only do the rest if there are un-partnered orbitals
c
      call ga_sync()
c
      call movecs_print_anal(basis,int_mb(k_non),int_mb(k_non)
     & ,0.15d0,g_alpha,'Alpha Orbitals without Beta Partners',
     &   .false., 0.0 ,.false., 0 , .false., 0 )
c
      if (nct.GE.2) then
      do i = 2,nct
      ind = int_mb(k_non+i-1)
      call movecs_print_anal(basis,ind,ind
     & ,0.15d0,g_alpha,' ',
     &   .false., 0.0 ,.false., 0 , .false., 0 )
      enddo
      endif
c
c Overlap Diagonalization
c
c SVD diagonalizes Sab*SabT (Sab is the S' calculated above). 
c Following diagonalization, it is clearer which
c alpha orbitals do not have beta partners. For alpha orbitals
c with partners, the overlap is very nearly 1.0. For
c un-partnered orbitals, the overlap is zero.
c See J. Chem. Phys. (1967) 47, 1936.
c
c
c if the SVD vectors and overlap aren't desired in the output,
c don't bother calculating them.
c
c calculate Sab*SabT
c 
      call ga_sync()
      if(.not.ga_create(mt_dbl,noc(1),noc(1),'SS',0,0,g_ss))
     $        call errquit('ga_create failed', g_ss, GA_ERR)
c
      call ga_dgemm('N','T',noc(1),noc(1),noc(2),1.d0,g_s,g_s,0.d0,g_ss)
c
c      call ga_print(g_ss)
c create arrays needed for SVD
c
      if(.not.ga_create(mt_dbl,noc(1),noc(1),'u',0,0,g_u))
     $        call errquit('ga_create failed', g_u, GA_ERR)
c
      if(.not.ga_create(mt_dbl,noc(1),noc(1),'vt',0,0,g_vt))
     $        call errquit('ga_create failed', g_vt, GA_ERR)
c
      if(.not.MA_Push_Get(MT_dbl,noc(1),'values',l_vals, k_vals))
     &     call errquit('dft_mxspin_ovlp: cannot allocate values',0,
     &       MA_ERR)
      call ga_sync()
c
c perform SVD on Sab*SabT to determine unpaired alpha MO's
c
      call ga_svd_seq(g_ss,g_u,g_vt,dbl_mb(k_vals))   
c
c      call ga_print(g_u)
c      call ga_print(g_vt)
c
      if (.not. ga_destroy(g_vt)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_vt', 0, GA_ERR)
      if (.not. ga_destroy(g_ss)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_ss', 0, GA_ERR)
c
      call ga_sync()
c      write(luout,*)'SVD eigenvalues'
c      write(luout,*) (dbl_mb(k_vals+i),i=0,noc(1)-1)
c      call ga_sync()
c
c calculate transformed alpha vectors, alphaT
c
      if(.not.ga_create(mt_dbl,noc(1),noc(1),'t',0,0,g_t))
     $        call errquit('ga_create failed', g_t, GA_ERR)
c
      call ga_zero(g_t)
c
      call ga_sync()
      do i = 1,noc(1)
       call ga_put(g_t,i,i,i,i,dbl_mb(k_vals+i-1),k_vals)
      enddo
      call ga_sync()
c
c      call ga_print(g_t)
c
      if(.not.ga_create(mt_dbl,nao,noc(1),'alphaT',0,0,g_alphaT))
     $        call errquit('ga_create failed', g_alphaT, GA_ERR)
c
      call ga_dgemm('N','N',nao,noc(1),noc(1),1.d0
     &     ,g_alpha,g_u,0.d0,g_alphaT)
c
c      call ga_print(g_alpha)
      if (.not. ga_destroy(g_u)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_u', 0, GA_ERR)
c
c      call ga_print(g_alphaT)
c
      if (.not. ga_destroy(g_t)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_t', 0, GA_ERR)
c
c     create array containing only alpha MO's which don't
c     have beta partners
c
      call ga_sync()
      nalp = 0
      alp_thresh = 1.0d-10
c
      if(.not.MA_Push_Get(MT_int,noc(1),'unpaired',l_unp, k_unp))
     &     call errquit('dft_mxspin_ovlp: cannot allocate values',0,
     &       MA_ERR)
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
c
c      write(luout,*) 'paired/unpaired alpha orbitals'
c      write(luout,*) (int_mb(k_unp+i-1),i=1,noc(1))
c
      call ga_sync()
      if(.not.ga_create(mt_dbl,nao,max(nalp,1),'unp alphaT',
     $                  0,0,g_ualpha))
     $        call errquit('ga_create failed', g_ualpha, GA_ERR)
c
      nalp=0
      do i = 1,noc(1)
       if(int_mb(k_unp+i-1).EQ.1) then
         nalp = nalp + 1
         call ga_copy_patch('N',g_alphaT,1,nbf,i,i
     &         ,g_ualpha,1,nao,nalp,nalp)
        endif
      enddo
      call ga_sync()

c      call ga_print(g_ualpha)                 
c
c associate SVD transformed alpha orbitals with the original
c alpha orbitals:
c
      if(me.eq.0) then
      write(luout,9989)
     &        '=================================================='
      write(luout,9989)
     &     '  Performing Singular Value Decomposition (SVD)   '
      write(luout,9989)
     &     'to diagonalize and maximize the alpha/beta overlap'
      write(luout,9989)
     &     '       See J. Chem. Phys. (1967) 47, 1936.        '
      write(luout,9989)
      write(luout,9989)
     &     'NOTE: The vector numbering of the SVD transformed '
      write(luout,9989)
     &     'orbitals is different from the original orbitals. '
      write(luout,9989)
     &     '=================================================='
      write(luout,9989)
      endif
 9989 format(13x,a51)
c
c      if(me.eq.0) then
c      if (nalp.GT.1) then
c        write(luout,9990) nalp
c      endif 
c      endif
c 9990 format(/,18x,'THERE ARE',i3,1x,'UN-PARTNERED ALPHA ORBITALS')
c
       call movecs_print_anal(basis, 1, nalp, 0.15d0, g_ualpha,
     & 'Alpha Orb. w/o Beta Partners (after maxim. alpha/beta overlap)',
     &   .false., 0.0 ,.false., 0 , .false., 0 )
c
c print the SVD eigenvalues
c 
      call ga_sync()
      if(me.eq.0) then
      write(LuOut,9994)
        ichunks = noc(1)/10
        if(ichunks*10.lt.noc(1))ichunks = ichunks +1
        do i = 1, ichunks
          nbe = 10*(i-1) + 1
          nend = nbe + 9
          if(nend.gt.noc(1)) nend = noc(1)
          write(LuOut,9997)(n,n=nbe,nend)
          write(LuOut,9998)(n,n=nbe,nend)
          write(LuOut,9999)(dbl_mb(k_vals+n-1),n=nbe,nend)
        enddo
        write(luout,*)
      endif
c
 9994 format(/,1x,'  SVD maximized alpha - beta orbital overlaps ',/,
     &         1x,'  ------------------------------------------- ',/)
c
c
      if (.not. ga_destroy(g_ualpha)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_ualpha', 0, GA_ERR)
c
      if (.not. ga_destroy(g_alphaT)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_alphaT', 0, GA_ERR)
c
      if(.not.MA_Pop_Stack(l_unp))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0, MA_ERR)
c
      if(.not.MA_Pop_Stack(l_vals))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0, MA_ERR)
c
      endif !if na>nb
c
      if(.not.MA_Pop_Stack(l_non))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0, MA_ERR)
c
      if(.not.MA_Pop_Stack(l_part))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0, MA_ERR)
c
c---------------------------------------------------------------
2001  continue
      if (.not. ga_destroy(g_s)) call errquit
     &   ('dft_mxspin_ovlp: could not destroy g_s', 0, GA_ERR)
c
      if(.not.MA_chop_Stack(l_tmpr1))
     & call errquit('dft_mxspin_ovlp: cannot pop stack',0, MA_ERR)
c
      return
      end

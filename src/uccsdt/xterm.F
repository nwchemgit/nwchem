      subroutine uccsdt_makex(g_omega,d_x,spini,spinj)
c
c     x(k,l,i,j) = Sum(u,v) t(u,v,i,j)C(u,k,spini)C(v,l,spinj)
c     where u & v are SO's (no spin), k,b,i,j are MO's with spin labels
c
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
      integer d_x          ! [output] DRA handle for X(k,l,i,j)
      integer g_omega      ! [input] GA handle for Omega(u,v,i,j)
      integer spini, spinj ! [input] Spins of i and j
      integer spink, spinl ! Spins of k and l
      integer g_tmp        ! Temporary GA t(u,l,i,j) & t(u,v,i,j)
      integer g_x          ! Temporary GA x(k,l,i,j)
      integer k_buf,l_buf  ! Temporary MA for tranformation
      integer k_tmp,l_tmp  ! Temporary MA for tranformation
      integer k_c,l_c      ! Temporary MA for MO coeffs
      integer lenij,lenul,lenuv,lenu,lenv,lenl
      integer lenik,lenjl,leni,lenk,leniu
      integer maxlenuv
      integer symij,symj,symi,symul,syml
      integer symuv,symu,symv,symjl,symik,symk
      integer joff(nw_max_nbf_small,0:7)
      integer loff(nw_max_nbf_small,0:7)
      integer symvoff(0:7,0:7)
      integer symloff(0:7,0:7)
      integer ij,jl,ik,iu
      integer ijlo,ijhi,jllo,jlhi
      integer i,j,l
      integer llo,lhi,vlo,klo,ilo,ihi,ulo,uhi
      integer dummy
c
      spink = spini
      spinl = spinj
c
c     Addressing & sizes
c
      lenij = 0
      do symij = 0,7
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         joff(j,symij) = lenij
         lenij = lenij + no_sym(symi,spini)
        end do
       end do
      end do
c
c     Addressing & sizes
c
      do symul = 0,7
       lenul = 0
       do syml = 0,7
        symu = ieor(symul,syml)
        do l = o_sym(1,syml,spinl),o_sym(2,syml,spinl)
         loff(l,symul) = lenul
         lenul = lenul + bf_per_ir(symu)
        end do
       end do
      end do
c
c     Addressing & sizes
c
      do symul = 0,7
       lenul = 0
       do syml = 0,7
        symu = ieor(symul,syml)
        symloff(syml,symul) = lenul
        lenul = bf_per_ir(symu)*no_sym(syml,spinl)
       end do
      end do
c
c     Addressing & sizes
c
      do symuv = 0,7
       lenuv = 0
       do symv = 0,7
        symu = ieor(symuv,symv)
        symvoff(symv,symuv) = lenuv
        lenuv = lenuv + bf_per_ir(symu)*bf_per_ir(symv)
       end do
       if (lenuv .gt. maxlenuv) maxlenuv = lenuv
      end do
c
c     Transformed MO coefficients
c
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c)
      if (.not. ma_push_get(mt_dbl, nbf*no(spinl), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*no(spinl))
      call ga_get(g_hole(spinl),1,nbf,1,no(spinl),
     $ dbl_mb(k_c),nbf)
c
c     Allocate temporary GA
c
      if (.not.ga_create(mt_dbl,maxlenuv,lenij,'sopair', 
     $ maxlenuv,-1,g_tmp)) call errquit
     $ ('uccsdt_makex: room for tmp?',maxlenuv*lenij)
c
c     Data parallel transform Omega(u,v,i,j)C(v,l) = tmp(u,l,i,j)
c     (tmp was allocated big enough to do this).
c
      call ga_distribution(g_omega,ga_nodeid(),ijlo,ijhi,dummy,dummy)
      ij = 0
      do symij = 0,7
       symuv = symij
       do symj = 0,7
        symi = ieor(symij,symj)
        do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
         do i = o_sym(1,symi,spini),o_sym(2,symi,spini)
          ij = ij + 1
          if (ij.ge.ijlo .and. ij.le.ijhi) then
           do syml = 0,7
            symu = ieor(symij,syml)
            symv = syml
            lenl = no_sym(syml,spinl)
            lenu = bf_per_ir(symu)
            lenv = bf_per_ir(symv)
            lenul = lenu * lenl
            lenuv = lenu * lenv
            llo = o_sym(1,syml,spinl)
            vlo = bf_per_ir_cum(symv) + 1
            if (lenuv.gt.0 .and. lenl.gt.0) then
             if (.not. ma_push_get(mt_dbl, lenuv, 'buf',
     $        l_buf, k_buf)) 
     $       call errquit('ma? lenul',lenul)
             if (.not. ma_push_get(mt_dbl, lenul, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? lenuv',lenuv)
c
c     (uv) -> (ul) for l (virtual) & v (SO) of same irrep 
c
             call ga_get(g_omega,1,lenuv,ij,ij,
     $        dbl_mb(k_buf),lenu)
             call dgemm('n','n',lenu,lenl,lenv,
     $        1.0d0,dbl_mb(k_buf),lenu,
     $        dbl_mb(k_c),nbf,
     $        0.0d0,dbl_mb(k_tmp),lenu)
             call ga_put(g_tmp,
     $        symloff(syml,symul)+1,symloff(syml,symul)+lenul,
     $        ij,ij,dbl_mb(k_tmp),lenu)
             if (.not. ma_pop_stack(l_tmp))
     $        call errquit('ma_pop?',l_tmp)
             if (.not. ma_pop_stack(l_buf))
     $        call errquit('ma_pop?',l_buf)
            end if
           end do
          end if
         end do
        end do
       end do
      end do
c
c     Perform quarter transformation from tmp(u,l,i,j) to t(k,l,i,j)
c     t(i,k,symk,j,l,syml,symjl)
c
      if (.not. ma_push_get(mt_dbl, nbf*no(spink), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*no(spink))
      call ga_get(g_hole(spink),1,nbf,1,no(spink),
     $ dbl_mb(k_c),nbf)
      do symjl = 0,7
       symik = symjl
       do syml = 0,7
        symj = ieor(symjl,syml)
c
c     Allocate X block
c
        llo = o_sym(1,syml,spinl)
        lhi = o_sym(2,syml,spinl)
        lenik = 0
        do symk = 0,7
         symi = ieor(symik,symk)
         lenik = lenik + no_sym(symi,spini) * no_sym(symk,spink)
        end do
        lenjl = (lhi-llo+1)*no_sym(symj,spinj)
        if (.not.ga_create(mt_dbl,lenik,lenjl,'l', 
     $   lenik,-1,g_x)) call errquit
     $   ('uccsdt_makex: room for x?',lenik*lenjl)
c
c     Data-parallel transformation of k
c
        call ga_distribution(g_tmp,ga_nodeid(),jllo,jlhi,dummy,dummy)
        jl = 0
        do l = llo, lhi
         do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
          jl = jl + 1
          if (jl.ge.jllo .and. jl.le.jlhi) then
           ik = 0
           iu = 0
c     Compare the next line with the comments in Amplitude.F
c     I think the comments are not accurate
           do symk = 0, 7
            symi = ieor(symik,symk)
            symij = ieor(symi,symj)
            symu = symk
            lenk = no_sym(symk,spink)
            leni = no_sym(symi,spini)
            lenu = bf_per_ir(symu)
            lenik = leni*lenk
            leniu = leni*lenu
            klo = o_sym(1,symk,spink)
            ilo = o_sym(1,symi,spini)
            ihi = o_sym(2,symi,spini)
            ulo = bf_per_ir_cum(symu) + 1
            uhi = ulo + lenu - 1
            if (lenik.gt.0 .and. lenu.gt.0) then
             if (.not. ma_push_get(mt_dbl, leniu, 'buf',
     $        l_buf, k_buf)) 
     $        call errquit('ma? lenia',leniu)
             if (.not. ma_push_get(mt_dbl, lenik, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? leniu',lenik)
c
c     (iu) -> (ik) for k (occupied) & u (SO) of same irrep 
c
             ij = joff(j,symij)
             do i = ilo,ihi
              ij = ij + 1
              call ga_get(g_tmp,
     $         loff(l,symul)+1,loff(l,symul)+lenu,
     $         ij,ij,dbl_mb(k_buf+lenu*(i-ilo)),lenu)
             end do
             call dgemm('t','t',lenk,leni,lenu,
     $        1.0d0,dbl_mb(k_c),nbf,
     $        dbl_mb(k_tmp),leni,
     $        0.0d0,dbl_mb(k_buf),lenu)
             call ga_put(g_x,ik+1,ik+lenik,jl,jl,
     $        dbl_mb(k_buf),leni)
             ik = ik + lenik
             if (.not. ma_pop_stack(l_tmp))
     $        call errquit('ma_pop?',l_tmp)
             if (.not. ma_pop_stack(l_buf))
     $        call errquit('ma_pop?',l_buf)
            end if
           end do
          end if
         end do
        end do
c
c     Write one symmetry block of X to disk
c
        if (.not. uccsdt_ampfile_write_t2(d_x, 
     $   spini, spink, spinj, spinl, symik, 
     $   llo, lhi, g_x))
     $   call errquit('write_t2 failed',0)
        if (.not. ga_destroy(g_x)) call errquit('GA 999',0)
       end do
      end do
      if (.not. ga_destroy(g_tmp)) call errquit('GA 999',0)
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c)
c
      return
      end

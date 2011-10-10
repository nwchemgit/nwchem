      subroutine uccsdt_makex(g_omega,g_x_big,spini,spinj)
csh   subroutine uccsdt_makex(g_omega,d_x,spini,spinj)
c
c     x(k,l,i,j) = Sum(u,v) t(u,v,i,j)C(u,k,spini)C(v,l,spinj)
c     where u & v are SO's (no spin), k,b,i,j are MO's with spin labels
c
      implicit none
#include "mafdecls.fh"
#include "cuccsdtP.fh"
#include "amplitudes.fh"
#include "global.fh"
      integer g_x_big      ! [output] DRA handle for X(k,l,i,j)
      integer g_omega      ! [input] GA handle for Omega(u,v,i,j)
      integer spini, spinj ! [input] Spins of i and j
      integer spink, spinl ! Spins of k and l
      integer g_tmp        ! Temporary GA x(u,l,i,j) & x(u,v,i,j)
      integer g_x          ! Temporary GA x(k,l,i,j)
      integer k_buf,l_buf  ! Temporary MA for tranformation
      integer k_tmp,l_tmp  ! Temporary MA for tranformation
      integer k_c,l_c      ! Temporary MA for MO coeffs
      integer lenij,lenul,lenuv,lenu,lenv,lenl
      integer lenik,lenjl,leni,lenk,leniu
      integer maxlenuv
      integer symij,symj,symi,symul,syml
      integer symuv,symu,symv,symjl,symik,symk
      integer joff(nbf,0:7)
      integer loff(nbf,0:7)
      integer symvoff(0:7,0:7)
      integer symloff(0:7,0:7)
      integer ij,jl,ik
      integer ijlo,ijhi,jllo,jlhi
      integer i,j,k,l
      integer llo,lhi,vlo,klo,ilo,ihi,ulo,uhi
      integer dummy
csh
      logical odebug
      odebug = .false.
csh
c
csh
      if (odebug) then
      write(*,*) 'no = ',no(1),no(2)
      write(*,*) 'nv = ',nv(1),nv(2)
      write(*,*) 'no_sym = ',(no_sym(i,1),i=0,7),(no_sym(i,2),i=0,7)
      write(*,*) 'nv_sym = ',(nv_sym(i,1),i=0,7),(nv_sym(i,2),i=0,7)
      write(*,*) 'o_sym(a) = ',(o_sym(1,i,1),i=0,7),(o_sym(2,i,1),i=0,7)
      write(*,*) 'o_sym(b) = ',(o_sym(1,i,2),i=0,7),(o_sym(2,i,2),i=0,7)
      write(*,*) 'v_sym(a) = ',(v_sym(1,i,1),i=0,7),(v_sym(2,i,1),i=0,7)
      write(*,*) 'v_sym(b) = ',(v_sym(1,i,2),i=0,7),(v_sym(2,i,2),i=0,7)
      write(*,*) 'bf_per_ir = ',(bf_per_ir(i),i=0,7)
      endif
csh
      spink = spini
      spinl = spinj
csh ***** We have to create an empty g_b_big
      dummy = 0
      do symjl = 0,7
      do syml = 0,7
      symj = ieor(symjl,syml)
      do l = 1,no_sym(syml,spinl)
      do j = 1,no_sym(symj,spinj)
      do symk = 0,7
      symi = ieor(symjl,symk)
      do k = 1,no_sym(symk,spink)
      do i = 1,no_sym(symi,spini)
      dummy = dummy + 1
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      if (odebug) write(*,*) 'Size of g_x_big = ',dummy
      if (.not.ga_create(mt_dbl,dummy,1,'x_big',
     $ -1,1,g_x_big)) call errquit
     $ ('uccsdt_makex: room for x_big?',dummy)
csh
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
      if (odebug) write(*,*) 'lenij = ',lenij
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
      if (odebug) write(*,*) 'maxlenuv = ',maxlenuv
c
c     Transformed MO coefficients
c
      if (.not. ma_push_get(mt_dbl, nbf*no(spinl), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*no(spinl),0)
      call ga_get(g_part(spinl),1,nbf,1,no(spinl),
     $ dbl_mb(k_c),nbf)
c
c     Allocate temporary GA
c
      if (.not.ga_create(mt_dbl,maxlenuv,lenij,'X tmp', 
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
       symul = symij
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
csh
            if (odebug) write(*,*) 'symij,symj,j,i,syml = ',
     $                  symij,symj,j,i,syml
csh
            if (lenuv.gt.0 .and. lenl.gt.0) then
csh
             if (odebug) write(*,*) 'lenul,lenuv = ',lenul,lenuv
csh
             if (.not. ma_push_get(mt_dbl, lenuv, 'buf',
     $        l_buf, k_buf)) 
     $       call errquit('ma? lenuv',lenuv,0)
             if (.not. ma_push_get(mt_dbl, lenul, 'tmp',
     $        l_tmp, k_tmp)) 
     $        call errquit('ma? lenul',lenul,0)
c
c     (uv) -> (ul) for l (occupied) & v (SO) of same irrep 
c
             call ga_get(g_omega,1,lenuv,ij,ij,
     $        dbl_mb(k_buf),lenu)
csh
             if (odebug) call ma_print(dbl_mb(k_buf),lenu,lenv,'UV')
csh
             call dgemm('n','n',lenu,lenl,lenv,
     $        1.0d0,dbl_mb(k_buf),lenu,
     $        dbl_mb(k_c),nbf,
     $        0.0d0,dbl_mb(k_tmp),lenu)
csh
             if (odebug) call ma_print(dbl_mb(k_tmp),lenu,lenl,'UL')
             if (odebug) write(*,*) 'ga_put in the range ... ',
     $       symloff(syml,symul)+1,'~',symloff(syml,symul)+lenul,
     $       ij,'~',ij
csh
             call ga_put(g_tmp,
     $        symloff(syml,symul)+1,symloff(syml,symul)+lenul,
     $        ij,ij,dbl_mb(k_tmp),lenu)
             if (.not. ma_pop_stack(l_tmp))
     $        call errquit('ma_pop?',l_tmp,0)
             if (.not. ma_pop_stack(l_buf))
     $        call errquit('ma_pop?',l_buf,0)
            end if
           end do
          end if
         end do
        end do
       end do
      end do
csh
      if (odebug) write(*,*) '**************************'
      if (odebug) write(*,*) 'transformation v -> l done'
      if (odebug) write(*,*) '**************************'
csh
c
c     Perform quarter transformation from tmp(u,l,i,j) to t(k,l,i,j)
c     t(i,k,symk,j,l,syml,symjl)
c
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c,0)
      if (.not. ma_push_get(mt_dbl, nbf*no(spink), 'c',
     $ l_c, k_c)) 
     $ call errquit('ma? nbf*no',nbf*no(spink),0)
      call ga_get(g_part(spink),1,nbf,1,no(spink),
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
        if ((lenik.gt.0).and.(lenjl.gt.0)) then
csh
         if (odebug) write(*,*) 'lenik,lenjl = ',lenik,lenjl
csh
         if (.not.ga_create(mt_dbl,lenik,lenjl,'l', 
     $    lenik,-1,g_x)) call errquit
     $    ('uccsdt_makex: room for x?',lenik*lenjl)
c
c     Data-parallel transformation of k
c
         call ga_distribution(g_x,ga_nodeid(),jllo,jlhi,dummy,dummy)
         jl = 0
         do l = llo, lhi
          do j = o_sym(1,symj,spinj),o_sym(2,symj,spinj)
           jl = jl + 1
           if (jl.ge.jllo .and. jl.le.jlhi) then
            ik = 0
c     Compare the next line with the comments in Amplitude.F
c     I think the comments are not accurate
            do symk = 0, 7
             symi = ieor(symik,symk)
             symij = ieor(symi,symj)
             symu = symk
             symul = ieor(symu,syml)
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
             if (lenu.gt.0) then
              if (.not. ma_push_get(mt_dbl, leniu, 'buf',
     $         l_buf, k_buf)) 
     $         call errquit('ma? leniu',leniu,0)
              if (.not. ma_push_get(mt_dbl, lenik, 'tmp',
     $         l_tmp, k_tmp)) 
     $         call errquit('ma? lenik',lenik,0)
c
c     (iu) -> (ik) for k (occupied) & u (SO) of same irrep 
c
              ij = joff(j,symij)
              do i = ilo,ihi
               ij = ij + 1
csh
               if (odebug) write(*,*) 'i = ',i
               if (odebug) write(*,*) 'get from range ... ',
     $          loff(l,symul)+1,'~',loff(l,symul)+lenu,
     $          ij,'~',ij,
     $          ' to k_buf+',lenu*(i-ilo)
csh
               call ga_get(g_tmp,
     $          loff(l,symul)+1,loff(l,symul)+lenu,
     $          ij,ij,dbl_mb(k_buf+lenu*(i-ilo)),lenu)
              end do
csh
              if (odebug) call ma_print(dbl_mb(k_buf),lenu,leni,'X(UI)')
csh
              call dgemm('t','n',leni,lenk,lenu,
     $         1.0d0,dbl_mb(k_buf),lenu,
     $         dbl_mb(k_c),nbf,
     $         0.0d0,dbl_mb(k_tmp),leni)
csh
              if (odebug) call ma_print(dbl_mb(k_tmp),leni,lenk,'X(IK)')
              if (odebug) write(*,*) 'writing to range ... ',
     $        ik+1,'~',ik+lenik,jl,'~',jl
csh
              call ga_put(g_x,ik+1,ik+lenik,jl,jl,
     $         dbl_mb(k_tmp),leni)
              ik = ik + lenik
              if (.not. ma_pop_stack(l_tmp))
     $         call errquit('ma_pop?',l_tmp,0)
              if (.not. ma_pop_stack(l_buf))
     $         call errquit('ma_pop?',l_buf,0)
             end if
            end do
           end if
          end do
         end do
c
c     Write one symmetry block of X to disk
c
         if (spini.eq.spinj) call ga_scale(g_x,2.0d0)
         if (odebug) call ga_print(g_x)
         call sh_put_x(g_x,syml,symjl,spini,spinj,g_x_big)
csh      if (.not. uccsdt_ampfile_write_t2(d_x, 
csh  $    spini, spink, spinj, spinl, symik, 
csh  $    llo, lhi, g_x))
csh  $    call errquit('write_t2 failed',0)
         if (.not. ga_destroy(g_x)) call errquit('GA 999',0,0)
        end if
       end do
      end do
      if (.not. ga_destroy(g_tmp)) call errquit('GA 999',0,0)
      if (.not. ma_pop_stack(l_c))
     $ call errquit('ma_pop?',l_c,0)
c
      return
      end
c $Id$

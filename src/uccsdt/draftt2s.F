      subroutine uccsdt_triples_amplitudes(d_amp, spina, spinc)
      implicit none
#include "errquit.fh"
#include "cuccsdtP.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "amplitudes.fh"
      integer ind, list,
     &     max_lenia, symkc, lenkc, l_t2, k_t2, symia, 
     &     a, alo, ahi, syma, asub, asublo, asubhi, asubdim, spina,
     &     b, blo, bhi, symb, bsub, bsublo, bsubhi, bsubdim, spinb,
     &     c, clo, chi, symc, csub, csublo, csubhi, csubdim, spinc, 
     &     e, elo, ehi, syme, esub, esublo, esubhi, esubdim, spine, 
     &     i, ilo, ihi, symi, isub, isublo, isubhi, isubdim, spini, 
     &     j, jlo, jhi, symj, jsub, jsublo, jsubhi, jsubdim, spinj, 
     &     k, klo, khi, symk, ksub, ksublo, ksubhi, ksubdim, spink, 
     &     m, mlo, mhi, symm, msub, msublo, msubhi, msubdim, spinm, 
     $     nproc, me, ei, cik, ma, kbc, mb, ia, ek, aki, l, lenia,
     $     symib, symie, symke, symmc, lenke, kc, max_lenke, 
     $     ptr, d_amp, g_t2, kac, offset, actual_lenkc, actual_lenia,
     $     mc, ica, actual_lenib, ib, bki, icb,
     $     max_lenje, symje, lenje, symmb, aji, ej, iab, bji

      double precision buf(1000)
c     
c     Offset maps a 4-d array into a 1-d array.  It is used to look up
c     ptr(ind,i,j,k,l) where ptr is dimensioned ptr(2,dim1,dim2,dim3,dim4) 
c     with dim1-4 = listinfo(6-9,list)
c     
      offset(ind,i,j,k,l,list) = 
     $     ind-1 + 
     $     2*(i-1 + 
     $     listinfo(6,list)*(j-1 + 
     $     listinfo(7,list)*(k-1 + 
     $     listinfo(8,list)*(l-1))))
c     
c     t(ia,kc)
c     
      spink = spinc
      spini = spina
      spinb = spina
      spine = spina
c     
c     get local memory of size max "lenkc"
c     
      max_lenia = 0
      do symia = 0,7
        lenia = ov_len(symia,spini,spina)
        if (lenia.gt.max_lenia)max_lenia = lenia
      enddo
      if(.not.ma_push_get(mt_dbl,max_lenia,'t2',l_t2,k_t2))
     $     call errquit('t2s: t2?',max_lenia, MA_ERR)
c     
      me = ga_nodeid()
      nproc = ga_nnodes()
c     
c     loop over symkc
c     
      alo = asuper(1)
      ahi = asuper(2)
      blo = bsuper(1)
      bhi = bsuper(2)
      clo = csuper(1)
      chi = csuper(2)
      ilo = nc(spini) + 1
      ihi = nc(spini) + no(spini)
      do symkc = 0,7
        symia = symkc
        symib = symkc
        symie = symkc
        lenia = ov_len(symia,spini,spina)
        actual_lenkc = ov_off(chi+1,symkc,spink,spinc) - 
     $       ov_off(clo,symkc,spink,spinc)
        if(actual_lenkc.gt.0.and.lenia.gt.0) then
          if (.not.uccsdt_ampfile_read_t2(d_amp, 
     $         spini, spina, spink, spinc, symkc, clo, chi, 
     $         g_t2, .true., 'column'))
     $         call errquit('amp_read_t2: reading t2 failed', d_amp,
     &       DISK_ERR)
          kc = 1
          do c = clo, chi
            csub = cblock_inv(c)
            csublo = cblock(1,csub)
            csubdim = cblock(2,csub) - cblock(1,csub) + 1
            symc = cblock(3,csub)
            symk = ieor(symkc,symc)
            do k = o_sym(1,symk,spink), o_sym(2,symk,spink)
              IF (MOD(KC,NPROC).eq.ME) THEN
                ksub = oblock_inv(k,spink)
                ksublo = oblock(1,ksub,spink)
                ksubdim = oblock(2,ksub,spink) - 
     $               oblock(1,ksub,spink) + 1
                call ga_get(g_t2,1,lenia,kc,kc,dbl_mb(k_t2),1)
c               if (.not.ma_verify_allocator_stuff())
c               $               call errquit(' after ga ',0)
c               
c               14.  t(e,c,i,k) spin(e)=spin(a) mixed = t(a,c,i,k)
c               
                do i = ilo, ihi
                  isub = oblock_inv(i,spini)
                  isublo = oblock(1,isub,spini)
                  isubdim = oblock(2,isub,spini) - 
     $                 oblock(1,isub,spini) + 1
                  symi = oblock(3,isub,spini)
                  syme = ieor(symie,symi)
                  esublo = v_sym(1,syme,spine)
                  esubhi = v_sym(2,syme,spine)
                  esubdim = esubhi - esublo + 1
                  if (esubdim.gt.0)then
                    cik = (c-csublo + csubdim*(i-isublo + 
     $                   isubdim*(k-ksublo)))
                    list = 14
                    ptr = int_mb(listinfo(2,list) +
     $                   offset(1,1,csub,isub,ksub,list))
                    ptr = ptr + esubdim*cik
                    ei = k_t2 + ov_off(esublo,symie,spini,spine) + 
     $                   i - o_sym(1,symi,spini)
                    call dfill(1000, 0.0d0, buf, 1)
                    do e = 1, esubdim
                      buf(e) = dbl_mb(ei+(e-1)*no_sym(symi,spini))
                    enddo
                    call ga_put(listinfo(5,list),ptr,
     $                   ptr+esubdim-1,1,1,buf,1)
                  endif
                enddo
c               
c               20.  t(m,k,a,c) spin(m)=spin(a) mixed = t(i,a,k,c)  m=i
c               
                do a = alo, ahi
                  asub = ablock_inv(a)
                  asublo = ablock(1,asub)
                  asubdim = ablock(2,asub) - ablock(1,asub) + 1
                  syma = ablock(3,asub)
                  symm = ieor(symia,syma)
                  spinm = spini
                  msublo = o_sym(1,symm,spinm)
                  msubhi = o_sym(2,symm,spinm)
                  msubdim = msubhi - msublo + 1
                  if (msubdim.gt.0)then
                    kac = (k-ksublo + ksubdim*(a-asublo + 
     $                   asubdim*(c-csublo)))
                    list = 20
                    ptr = int_mb(listinfo(2,list) +
     $                   offset(1,1,ksub,asub,csub,list))
                    ptr = ptr + msubdim*kac
                    ma = k_t2 + ov_off(a,symia,spini,spina)
                    call ga_put(listinfo(5,list),ptr,
     $                   ptr+msubdim-1,1,1,dbl_mb(ma),1)
                  endif
                enddo
c               
c               21.  t(m,k,b,c) spin(m)=spin(b) mixed = t(i,a,k,c)  m=i
c               
                do b = blo, bhi
                  bsub = bblock_inv(b)
                  bsublo = bblock(1,bsub)
                  bsubdim = bblock(2,bsub) - bblock(1,bsub) + 1
                  symb = bblock(3,bsub)
                  symm = ieor(symib,symb)
                  spinm = spini
                  msublo = o_sym(1,symm,spinm)
                  msubhi = o_sym(2,symm,spinm)
                  msubdim = msubhi - msublo + 1
                  if (msubdim.gt.0)then
                    kbc = (k-ksublo + ksubdim*(b-bsublo + 
     $                   bsubdim*(c-csublo)))
                    list = 21
                    ptr = int_mb(listinfo(2,list) +
     $                   offset(1,1,ksub,bsub,csub,list))
                    ptr = ptr + msubdim*kbc
                    mb = k_t2 + ov_off(b,symib,spini,spinb)
                    call ga_put(listinfo(5,list),ptr,
     $                   ptr+msubdim-1,1,1,dbl_mb(mb),1)
                  endif
                enddo
              ENDIF             ! end parallel work
              kc = kc + 1
            enddo
          enddo
          if (.not. ga_destroy(g_t2))
     $         call errquit('t2s: ga_destroy?',1, GA_ERR)
        endif
      enddo
      if (.not. ma_pop_stack(l_t2))
     $     call errquit('t2s: ma_pop_stack?',1, MA_ERR)
c     
c     get local memory of size max "lenkc"
c     
      spine = spink
      max_lenke = 0
      do symke = 0,7
        lenke = ov_len(symke,spink,spine)
        if (lenke.gt.max_lenke)max_lenke = lenke
      enddo
      if(.not.ma_push_get(mt_dbl,max_lenke,'t2',l_t2,k_t2))
     $     call errquit('t2s: t2?',max_lenke, MA_ERR)
c     
      klo = nc(spink) + 1
      khi = nc(spink) + no(spink)
      do symia = 0,7
        symke = symia
        symmc = symia
        lenke = ov_len(symke,spink,spine)
        actual_lenia = ov_off(ahi+1,symia,spini,spina) - 
     $       ov_off(alo,symia,spini,spina)
        if(actual_lenia.gt.0.and.lenke.gt.0) then
          if (.not.uccsdt_ampfile_read_t2(d_amp, 
     $         spink, spine, spini, spina, symia, alo, ahi, 
     $         g_t2, .true., 'column'))
     $         call errquit('amp_read_t2: reading t2 failed', d_amp,
     &       DISK_ERR)
          ia = 1
          do a = alo, ahi
            asub = ablock_inv(a)
            asublo = ablock(1,asub)
            asubdim = ablock(2,asub) - ablock(1,asub) + 1
            syma = ablock(3,asub)
            symi = ieor(symia,syma)
            do i = o_sym(1,symi,spini), o_sym(2,symi,spini)
              IF (MOD(IA,NPROC).eq.ME) THEN
                isub = oblock_inv(i,spini)
                isublo = oblock(1,isub,spini)
                isubdim = oblock(2,isub,spini) - 
     $               oblock(1,isub,spini) + 1
                call ga_get(g_t2,1,lenke,ia,ia,dbl_mb(k_t2),1)
c               if (.not.ma_verify_allocator_stuff())
c               $               call errquit(' after ga ',0)
c               
c               17.  t(e,a,k,i) spin(e)=spin(k) mixed = t(k,e,i,a)
c               
                do k = klo, khi
                  ksub = oblock_inv(k,spink)
                  ksublo = oblock(1,ksub,spink)
                  ksubdim = oblock(2,ksub,spink) - 
     $                 oblock(1,ksub,spink) + 1
                  symk = oblock(3,ksub,spink)
                  syme = ieor(symke,symk)
                  esublo = v_sym(1,syme,spine)
                  esubhi = v_sym(2,syme,spine)
                  esubdim = esubhi - esublo + 1
                  if (esubdim.gt.0)then
                    aki = (a-asublo + asubdim*(k-ksublo + 
     $                   ksubdim*(i-isublo)))
                    list = 17
                    ptr = int_mb(listinfo(2,list) +
     $                   offset(1,1,asub,ksub,isub,list))
                    ptr = ptr + esubdim*aki
                    ek = k_t2 + ov_off(esublo,symke,spink,spine) + 
     $                   k - o_sym(1,symk,spink)
                    call dfill(1000, 0.0d0, buf, 1)
                    do e = 1, esubdim
                      buf(e) = dbl_mb(ek+(e-1)*no_sym(symk,spink))
                    enddo
                    call ga_put(listinfo(5,list),ptr,
     $                   ptr+esubdim-1,1,1,buf,1)
                  endif
                enddo
c               
c               22.  t(m,i,c,a) spin(m)=spin(c) mixed = t(m,c,i,a)
c               
                do c = clo, chi
                  csub = cblock_inv(c)
                  csublo = cblock(1,csub)
                  csubdim = cblock(2,csub) - cblock(1,csub) + 1
                  symc = cblock(3,csub)
                  symm = ieor(symmc,symc)
                  spinm = spinc
                  msublo = o_sym(1,symm,spinm)
                  msubhi = o_sym(2,symm,spinm)
                  msubdim = msubhi - msublo + 1
                  if (msubdim.gt.0)then
                    ica = (i-isublo + isubdim*(c-csublo + 
     $                   csubdim*(a-asublo)))
                    list = 22
                    ptr = int_mb(listinfo(2,list) +
     $                   offset(1,1,isub,csub,asub,list))
                    ptr = ptr + msubdim*ica
                    mc = k_t2 + ov_off(c,symmc,spinm,spinc)
                    call ga_put(listinfo(5,list),ptr,
     $                   ptr+msubdim-1,1,1,dbl_mb(mc),1)
                  endif
                enddo
              ENDIF             ! end parallel work
              ia = ia + 1
            enddo
          enddo
          if (.not. ga_destroy(g_t2))
     $         call errquit('t2s: ga_destroy?',1, GA_ERR)
        endif
      enddo
c
      do symib = 0,7
        symke = symib
        symmc = symib
        lenke = ov_len(symke,spink,spine)
        actual_lenib = ov_off(bhi+1,symib,spini,spinb) - 
     $       ov_off(blo,symib,spini,spinb)
        if(actual_lenib.gt.0.and.lenke.gt.0) then
          if (.not.uccsdt_ampfile_read_t2(d_amp, 
     $         spink, spine, spini, spinb, symib, blo, bhi, 
     $         g_t2, .true., 'column'))
     $         call errquit('amp_read_t2: reading t2 failed', d_amp,
     &       DISK_ERR)
          ib = 1
          do b = blo, bhi
            bsub = bblock_inv(b)
            bsublo = bblock(1,bsub)
            bsubdim = bblock(2,bsub) - bblock(1,bsub) + 1
            symb = bblock(3,bsub)
            symi = ieor(symib,symb)
            do i = o_sym(1,symi,spini), o_sym(2,symi,spini)
              IF (MOD(IB,NPROC).eq.ME) THEN
                isub = oblock_inv(i,spini)
                isublo = oblock(1,isub,spini)
                isubdim = oblock(2,isub,spini) - 
     $               oblock(1,isub,spini) + 1
                call ga_get(g_t2,1,lenke,ib,ib,dbl_mb(k_t2),1)
c               if (.not.ma_verify_allocator_stuff())
c               $               call errquit(' after ga ',0)
c               
c               18.  t(e,b,k,i) spin(e)=spin(k) mixed = t(k,e,i,b)
c               
                do k = klo, khi
                  ksub = oblock_inv(k,spink)
                  ksublo = oblock(1,ksub,spink)
                  ksubdim = oblock(2,ksub,spink) - 
     $                 oblock(1,ksub,spink) + 1
                  symk = oblock(3,ksub,spink)
                  syme = ieor(symke,symk)
                  esublo = v_sym(1,syme,spine)
                  esubhi = v_sym(2,syme,spine)
                  esubdim = esubhi - esublo + 1
                  if (esubdim.gt.0)then
                    bki = (b-bsublo + bsubdim*(k-ksublo + 
     $                   ksubdim*(i-isublo)))
                    list = 18
                    ptr = int_mb(listinfo(2,list) +
     $                   offset(1,1,bsub,ksub,isub,list))
                    ptr = ptr + esubdim*bki
                    ek = k_t2 + ov_off(esublo,symke,spink,spine) + 
     $                   k - o_sym(1,symk,spink)
                    call dfill(1000, 0.0d0, buf, 1)
                    do e = 1, esubdim
                      buf(e) = dbl_mb(ek+(e-1)*no_sym(symk,spink))
                    enddo
                    call ga_put(listinfo(5,list),ptr,
     $                   ptr+esubdim-1,1,1,buf,1)
                  endif
                enddo
c               
c               23.  t(m,i,c,b) spin(m)=spin(c) mixed = t(m,c,i,b)
c               
                do c = clo, chi
                  csub = cblock_inv(c)
                  csublo = cblock(1,csub)
                  csubdim = cblock(2,csub) - cblock(1,csub) + 1
                  symc = cblock(3,csub)
                  symm = ieor(symmc,symc)
                  spinm = spinc
                  msublo = o_sym(1,symm,spinm)
                  msubhi = o_sym(2,symm,spinm)
                  msubdim = msubhi - msublo + 1
                  if (msubdim.gt.0)then
                    icb = (i-isublo + isubdim*(c-csublo + 
     $                   csubdim*(b-bsublo)))
                    list = 23
                    ptr = int_mb(listinfo(2,list) +
     $                   offset(1,1,isub,csub,bsub,list))
                    ptr = ptr + msubdim*icb
                    mc = k_t2 + ov_off(c,symmc,spinm,spinc)
                    call ga_put(listinfo(5,list),ptr,
     $                   ptr+msubdim-1,1,1,dbl_mb(mc),1)
                  endif
                enddo
              ENDIF             ! end parallel work
              ib = ib + 1
            enddo
          enddo
          if (.not. ga_destroy(g_t2))
     $         call errquit('t2s: ga_destroy?',1, GA_ERR)
        endif
      enddo

      if (.not. ma_pop_stack(l_t2))
     $     call errquit('t2s: ma_pop_stack?',1, MA_ERR)

c     
c     get local memory of size max "lenkc"
c     
      spinj = spini
      spine = spinj
      max_lenje = 0
      do symje = 0,7
        lenje = ov_len(symje,spinj,spine)
        if (lenje.gt.max_lenje)max_lenje = lenje
      enddo
      if(.not.ma_push_get(mt_dbl,max_lenje,'t2',l_t2,k_t2))
     $     call errquit('t2s: t2?',max_lenje, MA_ERR)
c     
      jlo = nc(spinj) + 1
      jhi = nc(spinj) + no(spinj)
      do symia = 0,7
        symje = symia
        symmb = symia
        lenje = ov_len(symje,spinj,spine)
        actual_lenia = ov_off(ahi+1,symia,spini,spina) - 
     $       ov_off(alo,symia,spini,spina)
        if(actual_lenia.gt.0.and.lenje.gt.0) then
          if (.not.uccsdt_ampfile_read_t2(d_amp, 
     $         spinj, spine, spini, spina, symia, alo, ahi, 
     $         g_t2, .true., 'column'))
     $         call errquit('amp_read_t2: reading t2 failed', d_amp,
     &       DISK_ERR)
          ia = 1
          do a = alo, ahi
            asub = ablock_inv(a)
            asublo = ablock(1,asub)
            asubdim = ablock(2,asub) - ablock(1,asub) + 1
            syma = ablock(3,asub)
            symi = ieor(symia,syma)
            do i = o_sym(1,symi,spini), o_sym(2,symi,spini)
              IF (MOD(IA,NPROC).eq.ME) THEN
                isub = oblock_inv(i,spini)
                isublo = oblock(1,isub,spini)
                isubdim = oblock(2,isub,spini) - 
     $               oblock(1,isub,spini) + 1
                call ga_get(g_t2,1,lenje,ia,ia,dbl_mb(k_t2),1)
c               if (.not.ma_verify_allocator_stuff())
c               $               call errquit(' after ga ',0)
c               
c               15.  t(e,a,j,i) spin(e)=spin(j) pure = t(j,e,i,a)
c               
                do j = jlo, jhi
                  jsub = oblock_inv(j,spinj)
                  jsublo = oblock(1,jsub,spinj)
                  jsubdim = oblock(2,jsub,spinj) - 
     $                 oblock(1,jsub,spinj) + 1
                  symj = oblock(3,jsub,spinj)
                  syme = ieor(symje,symj)
                  esublo = v_sym(1,syme,spine)
                  esubhi = v_sym(2,syme,spine)
                  esubdim = esubhi - esublo + 1
                  if (esubdim.gt.0)then
                    aji = (a-asublo + asubdim*(j-jsublo + 
     $                   jsubdim*(i-isublo)))
                    list = 15
                    ptr = int_mb(listinfo(2,list) +
     $                   offset(1,1,asub,jsub,isub,list))
                    ptr = ptr + esubdim*aji
                    ej = k_t2 + ov_off(esublo,symje,spinj,spine) + 
     $                   j - o_sym(1,symj,spinj)
                    call dfill(1000, 0.0d0, buf, 1)
                    do e = 1, esubdim
                      buf(e) = dbl_mb(ej+(e-1)*no_sym(symj,spinj))
                    enddo
                    call ga_put(listinfo(5,list),ptr,
     $                   ptr+esubdim-1,1,1,buf,1)
                  endif
                enddo
c               
c               19.  t(m,i,a,b) pure spin = -t(m,b,i,a)
c               
                do b = blo, bhi
                  bsub = bblock_inv(b)
                  bsublo = bblock(1,bsub)
                  bsubdim = bblock(2,bsub) - bblock(1,bsub) + 1
                  symb = bblock(3,bsub)
                  symm = ieor(symmb,symb)
                  spinm = spinb
                  msublo = o_sym(1,symm,spinm)
                  msubhi = o_sym(2,symm,spinm)
                  msubdim = msubhi - msublo + 1
                  if (msubdim.gt.0)then
                    iab = (i-isublo + isubdim*(a-asublo + 
     $                   asubdim*(b-bsublo)))
                    list = 19
                    ptr = int_mb(listinfo(2,list) +
     $                   offset(1,1,isub,asub,bsub,list))
                    ptr = ptr + msubdim*iab
                    mb = k_t2 + ov_off(b,symmb,spinm,spinb)
                    call dscal(msubdim,-1.0d0,dbl_mb(mb),1)
                    call ga_put(listinfo(5,list),ptr,
     $                   ptr+msubdim-1,1,1,dbl_mb(mb),1)
                    call dscal(msubdim,-1.0d0,dbl_mb(mb),1)
                  endif
                enddo
              ENDIF             ! end parallel work
              ia = ia + 1
            enddo
          enddo
          if (.not. ga_destroy(g_t2))
     $         call errquit('t2s: ga_destroy?',1, GA_ERR)
        endif
      enddo
      do symib = 0,7
        symje = symib
        symmb = symib
        lenje = ov_len(symje,spinj,spine)
        actual_lenib = ov_off(bhi+1,symib,spini,spinb) - 
     $       ov_off(blo,symib,spini,spinb)
        if(actual_lenib.gt.0.and.lenje.gt.0) then
          if (.not.uccsdt_ampfile_read_t2(d_amp, 
     $         spinj, spine, spini, spinb, symib, blo, bhi, 
     $         g_t2, .true., 'column'))
     $         call errquit('amp_read_t2: reading t2 failed', d_amp,
     &       DISK_ERR)
          ib = 1
          do b = blo, bhi
            bsub = bblock_inv(b)
            bsublo = bblock(1,bsub)
            bsubdim = bblock(2,bsub) - bblock(1,bsub) + 1
            symb = bblock(3,bsub)
            symi = ieor(symib,symb)
            do i = o_sym(1,symi,spini), o_sym(2,symi,spini)
              IF (MOD(IB,NPROC).eq.ME) THEN
                isub = oblock_inv(i,spini)
                isublo = oblock(1,isub,spini)
                isubdim = oblock(2,isub,spini) - 
     $               oblock(1,isub,spini) + 1
                call ga_get(g_t2,1,lenje,ib,ib,dbl_mb(k_t2),1)
c               if (.not.ma_verify_allocator_stuff())
c               $               call errquit(' after ga ',0)
c               
c               16.  t(e,b,j,i) spin(e)=spin(j) pure = t(j,e,i,b)
c               
                do j = jlo, jhi
                  jsub = oblock_inv(j,spinj)
                  jsublo = oblock(1,jsub,spinj)
                  jsubdim = oblock(2,jsub,spinj) - 
     $                 oblock(1,jsub,spinj) + 1
                  symj = oblock(3,jsub,spinj)
                  syme = ieor(symje,symj)
                  esublo = v_sym(1,syme,spine)
                  esubhi = v_sym(2,syme,spine)
                  esubdim = esubhi - esublo + 1
                  if (esubdim.gt.0)then
                    bji = (b-bsublo + bsubdim*(j-jsublo + 
     $                   jsubdim*(i-isublo)))
                    list = 16
                    ptr = int_mb(listinfo(2,list) +
     $                   offset(1,1,bsub,jsub,isub,list))
                    ptr = ptr + esubdim*bji
                    ej = k_t2 + ov_off(esublo,symje,spinj,spine) + 
     $                   j - o_sym(1,symj,spinj)
                    call dfill(1000, 0.0d0, buf, 1)
                    do e = 1, esubdim
                      buf(e) = dbl_mb(ej+(e-1)*no_sym(symj,spinj))
                    enddo
                    call ga_put(listinfo(5,list),ptr,
     $                   ptr+esubdim-1,1,1,buf,1)
                  endif
                enddo
              ENDIF             ! end parallel work
              ib = ib + 1
            enddo
          enddo
          if (.not. ga_destroy(g_t2))
     $         call errquit('t2s: ga_destroy?',1, GA_ERR)
        endif
      enddo
      if (.not. ma_pop_stack(l_t2))
     $     call errquit('t2s: ma_pop_stack?',1, MA_ERR)
c
      end

c $Id$

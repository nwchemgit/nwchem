      subroutine soints_block(basis, nbf,
     $     ishuqlo, ishuqhi, ioff, irange, 
     $     jshuqlo, jshuqhi, joff, jrange,
     $     kshuqlo, kshuqhi, koff, krange,
     $     lshuqlo, lshuqhi, loff, lrange,
     $     tol, g)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "schwarz.fh"
#include "nwc_const.fh"
#include "cbassym.fh"
      integer basis
      integer nbf
      integer ishuqlo, ishuqhi, irange(2,0:7), ioff(nbf)
      integer jshuqlo, jshuqhi, jrange(2,0:7), joff(nbf,0:7)
      integer kshuqlo, kshuqhi, krange(2,0:7), koff(nbf,0:7)
      integer lshuqlo, lshuqhi, lrange(2,0:7), loff(nbf,0:7)
      logical odirac
      double precision tol
      double precision g(*)
c
c     Assumes that int_init(), schwarz_init(), sym_abelian_init(),
c     have all been called and that soints_block_info() was used
c     to generate the addressing arrays *off() and *range().
c
c     For SOs arising from unique AO shells in the given ranges,
c     compute the 2-electron integrals in the SO basis using
c     the symmetry-blocked, or second, labelling for the SOs.
c     In this scheme SOs are blocked by irrep and within the
c     irrep are in the natural order of symmetry unique basis 
c     functions. See sym_abelian() for  more info on the labelling.
c
c     g(i,j,k,l) = (ij|kl)
c
c     The integrals non-zero by symmetry are stored densely as if
c     g() were declared as g(i,j,symj,k,l,syml,symkl) corresponding
c     to the below loops.  The array can also be thought as
c     being in the order g(ij,kl,symkl) where ij and kl run
c     over the pairs with symmetry kl.  The offset to the first 
c     element of a given symkl is symkloff(symkl).
c     
c     ind = 0
c     do symkl = 0, nir-1
c     .  symij = symkl
c     .  do syml = 0, nir-1
c     .     symk = ieor(syml,symkl)
c     .     do l = SOs of syml from l-shell-block
c     .        do k = SOs of symk from k-shell-block
c     .           do symj = 0, nir-1
c     .              symi = ieor(symij,symj)
c     .              do j = SOs of symj from j-shell-block
c     .                 do i = SOs of symi from i-shell-block
c     .                    ind = ind + 1
c     .                    g(ind) = (ij|kl) or <ij|kl>
c
c     The block of integrals may be addressed using arrays generated 
c     by soints_block_info()
c     (this will only work for indices non-zero by symmetry)
c     
c     ijkl = ioff(i) + joff(j,symij) + koff(k,symkl) 
c     .      + loff(l,symkl)
c     integral = g(ijkl)
c
c     Also useful,
c
c     irange(1,symi) = first orbital i of symi in i shell block 
c     irange(2,symi) = last  orbital i of symi in i shell block 
c     ditto for j, k, l
c     lenij(symij) = no. of ij pairs in the block with symmetry ij
c     lenkl(symkl) = ditto  kl
c     symkloff(symkl) = first integral in the symkl block.
c
      integer xishuq, xjshuq, xkshuq, xlshuq
      integer ishuq, jshuq, kshuq, lshuq, lopuq, kopuq, jopuq
      integer ish, jsh, ksh, lsh, jdim, kdim, ldim
      integer i, j, k, l, i2, j2, k2, l2, ijkl, ind, iuq
*      integer ngot
      integer symi, symj, symk, syml, symij, symkl
      integer juq, kuq, luq, jnw, knw, lnw, jshnw, kshnw, lshnw, nsh
      integer ilo, ihi, jlo, jhi, klo, khi, llo, lhi
      integer nopuq
      integer l_q4, l_scr, l_eri, l_lab, l_q, l_ao
      integer k_q4, k_scr, k_eri, k_lab, k_q, k_ao
      integer k_iq, k_jq, k_kq, k_lq
      integer k_ilab, k_jlab, k_klab, k_llab
      integer leneri, lenscr, lenao
      integer maxq, nq, idimuq
      double precision smax, skl, gijkl, c(8,8)
      integer nops, geom
      integer sym_number_ops
      external sym_number_ops
      integer mapi(nw_max_nbf_small), mapj(nw_max_nbf_small), 
     $     mapk(nw_max_nbf_small), mapl(nw_max_nbf_small)
c
      if (.not. bas_geom(basis,geom)) call errquit
     $     ('soints_block: basis?', 0, BASIS_ERR)
      if (.not. bas_numcont(basis,nsh)) call errquit
     $     ('soints_block: basis?', 0, BASIS_ERR)
c
c     AO integrals are stored with symmetry unique i, and all symmetry 
c     related j, k, l in natural order.  map{ijkl} maps from NWChem AO 
c     labels into this order offset and modified to address the array.
c
      call ifill(nbf,99999999,mapi,1)
      call ifill(nbf,99999999,mapj,1)
      call ifill(nbf,99999999,mapk,1)
      call ifill(nbf,99999999,mapl,1)
      jdim = sh_uq_bf(2,jshuqhi) - sh_uq_bf(1,jshuqlo) + 1
      kdim = sh_uq_bf(2,kshuqhi) - sh_uq_bf(1,kshuqlo) + 1
      ldim = sh_uq_bf(2,lshuqhi) - sh_uq_bf(1,lshuqlo) + 1
      idimuq = 0
      do xishuq = ishuqlo, ishuqhi
         ish = sh_uq(xishuq)
         if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $        call errquit('soints_block: cn2bfr?', ish, BASIS_ERR)
         do i = ilo,ihi
            mapi(i) = idimuq + i - ilo
         end do
         idimuq = idimuq + ihi - ilo + 1
      end do
      lenao = idimuq*jdim*kdim*ldim
c
      jlo = sh_uq_bf(1,jshuqlo)
      jhi = sh_uq_bf(2,jshuqhi)
      do j = jlo, jhi
         mapj(bf_nat(2,j)) = (j-jlo)*idimuq
      end do
      klo = sh_uq_bf(1,kshuqlo)
      khi = sh_uq_bf(2,kshuqhi)
      do k = klo, khi
         mapk(bf_nat(2,k)) = (k-klo)*(idimuq*jdim)
      end do
      llo = sh_uq_bf(1,lshuqlo)
      lhi = sh_uq_bf(2,lshuqhi)
      do l = llo, lhi
         mapl(bf_nat(2,l)) = (l-llo)*(idimuq*jdim*kdim)
      end do
c     
c     Allocate memory for the integrals.  Cannot do them all at once since
c     Texas internal memory allocation fails if we give it more than
c     10000 quartets at once.
c     
      nops = sym_number_ops(geom) + 1
      maxq = 5000
*      maxq = (ishuqhi-ishuqlo+1)*(jshuqhi-jshuqlo+1)*
*     $     (kshuqhi-kshuqlo+1)*(lshuqhi-lshuqlo+1)*(nops**3)
c
      call intb_mem_2e4c(leneri, lenscr)
      leneri = max(leneri,100*(6**4)) ! Room for 100 (dd|dd)
      lenscr = lenscr + lenscr/2
C
      if (.not. ma_push_get(mt_dbl, maxq, 'q4', l_q4, k_q4))
     $     call errquit('soints_block: ma q4?',maxq, MA_ERR)
      if (.not. ma_push_get(mt_int, 4*maxq, 'quartets', l_q, k_q))
     $     call errquit('soints_block: ma q?', 4*maxq, MA_ERR)
      if (.not. ma_push_get(mt_int, 4*leneri, 'labels', l_lab, k_lab))
     $     call errquit('soints_block: ma lab?', 4*leneri, MA_ERR)
      if (.not. ma_push_get(mt_dbl, leneri, 'eri', l_eri, k_eri))
     $     call errquit('soints_block: ma eri?', leneri, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lenscr, 'scr', l_scr, k_scr))
     $     call errquit('soints_block: ma scr?', lenscr, MA_ERR)
      if (.not. ma_push_get(mt_dbl, lenao, 'ao', l_ao, k_ao))
     $     call errquit('soints_block: ma ao?',lenao, MA_ERR)
      k_iq = k_q          
      k_jq = k_q + maxq   
      k_kq = k_q + maxq*2 
      k_lq = k_q + maxq*3 
      k_ilab = k_lab       
      k_jlab = k_lab + leneri  
      k_klab = k_lab + leneri*2
      k_llab = k_lab + leneri*3
      call dfill(lenao, 0d0, dbl_mb(k_ao), 1)
c
c     Loop thru symmetry unique shells, then symmetry related shells jkl
c
      nq = 0
      smax = schwarz_max()
      do xlshuq = lshuqlo, lshuqhi
         lshuq = sh_uq(xlshuq)
         do xkshuq = kshuqlo, kshuqhi
            kshuq = sh_uq(xkshuq)
            do xjshuq = jshuqlo, jshuqhi
               jshuq = sh_uq(xjshuq)
               do xishuq = ishuqlo, ishuqhi
                  ishuq = sh_uq(xishuq)
                  ish = ishuq   ! Don't use related ish
                  do lopuq = 1, sh_n_uq_op(lshuq)
                     lsh = sh_op_map(sh_uq_op(lopuq,lshuq),lshuq)
                     do kopuq = 1, sh_n_uq_op(kshuq)
                        ksh = sh_op_map(sh_uq_op(kopuq,kshuq),kshuq)
                        skl = schwarz_shell(ksh,lsh)
                        if (skl*smax.gt.tol) then
                           do jopuq = 1, sh_n_uq_op(jshuq)
                              jsh=sh_op_map(sh_uq_op(jopuq,jshuq),jshuq)
                              if (skl*schwarz_shell(ish,jsh).gt.tol)then

                                 if ( ish.le.0.or.ish.gt.nsh.or.
     $                                jsh.le.0.or.jsh.gt.nsh.or.
     $                                ksh.le.0.or.ksh.gt.nsh.or.
     $                                lsh.le.0.or.lsh.gt.nsh) then
                                    write(6,*) ish,jsh,ksh,lsh
                                    call errquit('bad sh',0, BASIS_ERR)
                                 end if
                                 int_mb(k_iq + nq) = ish
                                 int_mb(k_jq + nq) = jsh
                                 int_mb(k_kq + nq) = ksh
                                 int_mb(k_lq + nq) = lsh
                                 dbl_mb(k_q4 + nq) = 1.0d0
                                 nq = nq + 1
c
                                 if (nq .eq. maxq) then
                                    call soints_aoints(basis, nq,
     $                                   int_mb(k_iq), int_mb(k_jq), 
     $                                   int_mb(k_kq), int_mb(k_lq), 
     $                                   dbl_mb(k_q4), tol, 
     $                                   int_mb(k_ilab), int_mb(k_jlab),
     $                                   int_mb(k_klab), int_mb(k_llab),
     $                                   dbl_mb(k_eri), leneri, 
     $                                   dbl_mb(k_scr), lenscr,
     $                                   mapi, mapj, mapk, mapl, 
     $                                   dbl_mb(k_ao), lenao)
                                    nq = 0
                                 end if
                              end if
                           end do
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
      if (nq .gt. 0) call soints_aoints(basis, nq,
     $     int_mb(k_iq), int_mb(k_jq), 
     $     int_mb(k_kq), int_mb(k_lq), 
     $     dbl_mb(k_q4), tol, 
     $     int_mb(k_ilab), int_mb(k_jlab),
     $     int_mb(k_klab), int_mb(k_llab),
     $     dbl_mb(k_eri), leneri, 
     $     dbl_mb(k_scr), lenscr,
     $     mapi, mapj, mapk, mapl, 
     $     dbl_mb(k_ao), lenao)
c
c$$$ 789  format(i6,4x,4i5,f12.6)
c$$$      write(6,*) ' The AO integrals'
c$$$      ind = 0
c$$$      ngot = 0
c$$$      do l = sh_uq_bf(1,lshuqlo), sh_uq_bf(2,lshuqhi)
c$$$         do k = sh_uq_bf(1,kshuqlo), sh_uq_bf(2,kshuqhi)
c$$$            do j = sh_uq_bf(1,jshuqlo), sh_uq_bf(2,jshuqhi)
c$$$               do iuq = sh_uq_bf(3,ishuqlo), sh_uq_bf(4,ishuqhi)
c$$$                  gijkl = dbl_mb(k_ao+ind) ! The integral
c$$$                  ind = ind + 1
c$$$                  if (abs(gijkl).gt.1e-6) then
c$$$                      write(6,789) ind, iuq,j,k,l,gijkl
c$$$                      ngot = ngot + 1
c$$$                   end if
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      write(6,*) ngot
c
c     Now have all of the AO integrals.
c
c     For each of lsh, ksh, and jsh loop thru the generators 
c     (symmetry unique/lexically highest functions) in the block of 
c     unique shells and transform them to the SO basis.
c
      do lshuq = lshuqlo, lshuqhi ! Loop thru unique shells
         lshnw = sh_uq(lshuq)     ! NWChem order shell
         nopuq  = sh_n_uq_op(lshnw) ! No. of related functions
         do luq = sh_uq_bf(3,lshuq), sh_uq_bf(4,lshuq)
            lnw = bf_uq(luq)
            call soints_uq_to_so_matrix(luq, c)
c$$$            write(6,*) ' luq ', luq, ' nopuq ', nopuq, ' lnw ', lnw,
c$$$     $           ' mapl ', mapl(lnw)
c$$$            call output(c,1,nopuq,1,nopuq,8,8,1)
            call soints_inplace_mxm(idimuq*jdim*kdim, nopuq, 1,
     $           dbl_mb(k_ao+mapl(lnw)), ldim, c)
         end do
      end do
c$$$      write(6,*) ' After L'
c$$$      ind = 0
c$$$      ngot = 0
c$$$      do l = sh_uq_bf(1,lshuqlo), sh_uq_bf(2,lshuqhi)
c$$$         do k = sh_uq_bf(1,kshuqlo), sh_uq_bf(2,kshuqhi)
c$$$            do j = sh_uq_bf(1,jshuqlo), sh_uq_bf(2,jshuqhi)
c$$$               do iuq = sh_uq_bf(3,ishuqlo), sh_uq_bf(4,ishuqhi)
c$$$                  gijkl = dbl_mb(k_ao+ind) ! The integral
c$$$                  ind = ind + 1
c$$$                  if (abs(gijkl).gt.1e-6) then
c$$$                      write(6,789) ind, iuq,j,k,l,gijkl
c$$$                      ngot = ngot + 1
c$$$                   end if
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      write(6,*) ngot
      do kshuq = kshuqlo, kshuqhi ! Loop thru unique shells
         kshnw = sh_uq(kshuq)     ! NWChem order shell
         nopuq  = sh_n_uq_op(kshnw) ! No. of related functions
         do kuq = sh_uq_bf(3,kshuq), sh_uq_bf(4,kshuq)
            knw = bf_uq(kuq)
            call soints_uq_to_so_matrix(kuq, c)
            call soints_inplace_mxm(idimuq*jdim, nopuq, ldim,
     $           dbl_mb(k_ao+mapk(knw)), kdim, c)
         end do
      end do
c$$$      write(6,*) ' After K'
c$$$      ind = 0
c$$$      ngot = 0
c$$$      do l = sh_uq_bf(1,lshuqlo), sh_uq_bf(2,lshuqhi)
c$$$         do k = sh_uq_bf(1,kshuqlo), sh_uq_bf(2,kshuqhi)
c$$$            do j = sh_uq_bf(1,jshuqlo), sh_uq_bf(2,jshuqhi)
c$$$               do iuq = sh_uq_bf(3,ishuqlo), sh_uq_bf(4,ishuqhi)
c$$$                  gijkl = dbl_mb(k_ao+ind) ! The integral
c$$$                  ind = ind + 1
c$$$                  if (abs(gijkl).gt.1e-6) then
c$$$                      write(6,789) ind, iuq,j,k,l,gijkl
c$$$                      ngot = ngot + 1
c$$$                   end if
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      write(6,*) ngot
      do jshuq = jshuqlo, jshuqhi ! Loop thru unique shells
         jshnw = sh_uq(jshuq)     ! NWChem order shell
         nopuq  = sh_n_uq_op(jshnw) ! No. of related functions
         do juq = sh_uq_bf(3,jshuq), sh_uq_bf(4,jshuq)
            jnw = bf_uq(juq)
            call soints_uq_to_so_matrix(juq, c)
            call soints_inplace_mxm(idimuq, nopuq, kdim*ldim,
     $           dbl_mb(k_ao+mapj(jnw)), jdim, c)
         end do
      end do
c$$$      write(6,*) ' After K'
c$$$      ind = 0
c$$$      ngot = 0
c$$$      do l = sh_uq_bf(1,lshuqlo), sh_uq_bf(2,lshuqhi)
c$$$         do k = sh_uq_bf(1,kshuqlo), sh_uq_bf(2,kshuqhi)
c$$$            do j = sh_uq_bf(1,jshuqlo), sh_uq_bf(2,jshuqhi)
c$$$               do iuq = sh_uq_bf(3,ishuqlo), sh_uq_bf(4,ishuqhi)
c$$$                  gijkl = dbl_mb(k_ao+ind) ! The integral
c$$$                  ind = ind + 1
c$$$                  if (abs(gijkl).gt.1e-6) then
c$$$                      write(6,789) ind, iuq,j,k,l,gijkl
c$$$                      ngot = ngot + 1
c$$$                   end if
c$$$               end do
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      write(6,*) ngot
c
c     The transformation is complete.  Only the unique ish were
c     computed and they were not transformed, but the transformation
c     of the other three indices act as a projector for their
c     symmetry.
c
c     A range of unique shells will give rise to the same ranges
c     for both basis functions and SOs in natural order.
c
      ind = 0
      do l = sh_uq_bf(1,lshuqlo), sh_uq_bf(2,lshuqhi)
         syml = so_lab1(1,l)
         l2   = so_lab1(2,l)
         do k = sh_uq_bf(1,kshuqlo), sh_uq_bf(2,kshuqhi)
            symk = so_lab1(1,k)
            k2   = so_lab1(2,k)
            symkl= ieor(symk,syml)
            symij= symkl
            do j = sh_uq_bf(1,jshuqlo), sh_uq_bf(2,jshuqhi)
               symj = so_lab1(1,j)
               j2   = so_lab1(2,j)
               symi = ieor(symij,symj)
               do iuq = sh_uq_bf(3,ishuqlo), sh_uq_bf(4,ishuqhi)
                  i2  = so_uq_to_lab(2,symi,iuq)
c
                  gijkl = dbl_mb(k_ao+ind) ! The integral
                  ind = ind + 1
                  if (i2 .lt. 0) then
                     if (abs(gijkl) .gt. 1d-6) then
                        write(6,*) i2,j2,k2,l2,gijkl
                        call errquit('should be zero?',ijkl,
     &       UNKNOWN_ERR)
                     end if
                  else
                     ijkl = ioff(i2) + joff(j2,symij) + 
     $                       koff(k2,symkl) + loff(l2,symkl)
                     g(ijkl) = gijkl
                  end if
               end do
            end do
         end do
      end do
c
      if (.not. ma_chop_stack(l_q4)) call errquit
     $     ('soints_aoblock: ma corrupted?',0, MA_ERR)
c
      end
      subroutine soints_block_info(basis, nbf,
     $     ishuqlo, ishuqhi, ioff, irange, 
     $     jshuqlo, jshuqhi, joff, jrange,
     $     kshuqlo, kshuqhi, koff, krange,
     $     lshuqlo, lshuqhi, loff, lrange,
     $     lenij, lenkl, symkloff, leng)
      implicit none
#include "errquit.fh"
      integer basis
      integer nbf
      integer ishuqlo, ishuqhi, irange(2,0:7), ioff(nbf)
      integer jshuqlo, jshuqhi, jrange(2,0:7), joff(nbf,0:7)
      integer kshuqlo, kshuqhi, krange(2,0:7), koff(nbf,0:7)
      integer lshuqlo, lshuqhi, lrange(2,0:7), loff(nbf,0:7)
      integer lenij(0:7), lenkl(0:7), symkloff(0:7), leng
c
c     See the comments in soints_block() for detailed info.
c
c     lenij(symij) is the number of ij pairs for a given symmetry
c     lenkl(symkl) ................ kl ..........................
c     symkloff(symkl) = index of the first element in the symkl block.
c     leng returns the no. of integrals non-zero by symmetry.
c     
      integer symi, symj, symk, syml, symij, symkl
      integer i, j, k, l, ijkl, count
      integer leni, lenk
c
c     Generate the ranges of SOs in each block of shells
c
      call soints_shblock_info(basis, ishuqlo, ishuqhi, irange)
      call soints_shblock_info(basis, jshuqlo, jshuqhi, jrange)
      call soints_shblock_info(basis, kshuqlo, kshuqhi, krange)
      call soints_shblock_info(basis, lshuqlo, lshuqhi, lrange)
c
c     Generate the addressing arrays for the integrals.  Only entries
c     within the given blocks will be valid.
c
      call ifill(  nbf, 99999999, ioff, 1)
      call ifill(8*nbf, 99999999, joff, 1)
      call ifill(8*nbf, 99999999, koff, 1)
      call ifill(8*nbf, 99999999, loff, 1)
c
c     ioff(i) is just the index within its symmetry block
c
      do symi = 0, 7
         do i = irange(1,symi), irange(2,symi)
            ioff(i) = 1 + i - irange(1,symi)
         end do
      end do
c
c     For given (symij, symj) the matrix (i,j) is dense 
c     ... joff(j,symij) is the offset to the top of the
c     j'th column within the given symmetry, plus the size of
c     all the preceeding i symmetries.
c
      do symij = 0, 7
         lenij(symij) = 0
         do symj = 0, 7
            symi = ieor(symij,symj)
            do j = jrange(1,symj),jrange(2,symj)
               joff(j,symij) = lenij(symij)
               leni = max(0, irange(2,symi)-irange(1,symi)+1)
               lenij(symij) = lenij(symij) + leni
            end do
         end do
      end do
c
c     koff is offset within its symmetry block * lenij
c
      do symkl = 0,7
         symij = symkl
         do syml = 0, 7
            symk = ieor(symkl,syml)
            do k = krange(1,symk), krange(2,symk)
               koff(k,symkl) = (k - krange(1,symk))*lenij(symij)
            end do
         end do
      end do
c
c     loff as for joff but * lenij and accumulate all previous kl symmetries
c
      count = 0
      do symkl = 0, 7
         symij = symkl
         lenkl(symkl) = 0
         symkloff(symkl) = count + 1
         do syml = 0, 7
            symk = ieor(symkl,syml)
            do l = lrange(1,syml),lrange(2,syml)
               loff(l,symkl) = count
               lenk = max(0, krange(2,symk)-krange(1,symk)+1)
               lenkl(symkl) = lenkl(symkl) + lenk
               count = count + lenk*lenij(symij)
            end do
         end do
      end do
c
c     Verify the indexing ... comment this out when we are sure
c     it is working since it is O(N^4).
c
      count = 0
      do symkl = 0, 7
         symij = symkl
         do syml = 0, 7
            symk = ieor(symkl,syml)
            do l = lrange(1,syml),lrange(2,syml)
               do k = krange(1,symk),krange(2,symk)
                  do symj = 0, 7
                     symi = ieor(symij,symj)
                     do j = jrange(1,symj),jrange(2,symj)
                        do i = irange(1,symi),irange(2,symi)
                           count = count + 1
                           ijkl = ioff(i) + joff(j,symij) +
     $                          koff(k,symkl) + loff(l,symkl)
                           if (ijkl .ne. count) then
                              write(6,*) symi, symj, symk, syml,
     $                             symij, symkl
                              write(6,*) i,j,k,l,ioff(i),joff(j,symij),
     $                             koff(k,symkl),loff(l,symkl),
     $                             ijkl, count
                              call errquit(' soints_block_info: bad',0,
     &       UNKNOWN_ERR)
                           endif
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
c
      leng = count
c
      end
      subroutine soints_shblock_info(basis, ishuqlo, ishuqhi, irange)
      implicit none
#include "bas.fh"
#include "nwc_const.fh"
#include "cbassym.fh"
      integer basis, ishuqlo, ishuqhi ! [input]
      integer irange(2,0:7)     ! [output]
c
c     Given a range of symmetry unique AO shells determine
c     the range of SO functions within each irrep (in the
c     second labelling scheme).
c
c     Loop thru all of the unique basis functions.  Find the
c     first/last function in natural in the first/last shell. 
c     If there are no functions in an irrep, the range is set
c     to [0,-1].
c
      integer ir, lab1, lab1lo, lab1hi
c
      do ir = 0, 7
         irange(1,ir) = 0
         irange(2,ir) = -1
      end do
c
c     Determine the range of basis functions in natural order
c     ... the range of SOs in natual order (lab1) will be the same.
c
      lab1lo = sh_uq_bf(1,ishuqlo)
      lab1hi = sh_uq_bf(2,ishuqhi)
c
c     Loop forward thru the SOs to determine the first in each irrep
c     then backward to get the last, mapping into the symmetry
c     blocked order (lab2) at the same time.
c
      do ir = 0,7
         do lab1 = lab1lo,lab1hi
            if (so_lab1(1,lab1).eq.ir .and. irange(1,ir).eq.0) then
               irange(1,ir) = so_lab1(2,lab1)
               goto 10
            endif
         end do
 10      continue
      end do
      do ir = 0,7
         do lab1 = lab1hi,lab1lo,-1
            if (so_lab1(1,lab1).eq.ir .and. irange(2,ir).eq.-1) then
               irange(2,ir) = so_lab1(2,lab1)
               goto 20
            end if
         end do
 20      continue
      end do
c
c$$$      write(6,*) ' SO UQSHRANGE: ', ishuqlo, ishuqhi
c$$$      write(6,1) (irange(1,ir),irange(2,ir),ir=0,7)
c$$$ 1    format(8(2i4,4x))
c
      end
      subroutine soints_uq_to_so_matrix(luq, c)
      implicit none
#include "nwc_const.fh"
#include "cbassym.fh"
      integer luq               ! [input]
      double precision c(8,8)   ! [output]
c
c     For symmetry unique basis function number luq return
c     in C the matrix that transforms the symmetry related
c     AO functions into SO basis
c
c     |iruq> = sum(k=1..nopuq) C(iruq,luq) |k>
c
c     where iruq runs over the irreducible representations arising
c     from this generator (the actual irrep is ir=bf_so_ir(iruq,luq))
c
c     and k runs over the functions related by symmetry to the
c     unique basis function (their NWChem numbering is given 
c     as bf_op_map(bf_uq_op(opuq),bf_uq(luq)), opuq=1,nopuq. 
c
      integer nopuq, opuq, op, iruq, ir, lnw
      double precision scale, phase
c
      call dfill(8*8, 1d300, c, 1) ! Make unitialized numbers bad news
c
      lnw = bf_uq(luq)          ! NWChem label for unique function
      nopuq  = bf_n_uq_op(lnw)  ! No. of related functions
      scale = 1d0/dble(nopuq)
      do opuq = 1,nopuq
         op = bf_uq_op(opuq,lnw)
         phase = bf_phase(op,lnw)*scale
         do iruq = 1, nopuq
            ir = bf_so_ir(iruq,luq)
            c(iruq,opuq) = char_tab(op,ir)*phase
         end do
      end do
c
      end
      subroutine soints_inplace_mxm(n, nopuq, m, a, dim2, c)
      implicit none
      integer n, nopuq, m, dim2
      double precision a(n,dim2,m), c(8,8)
c
c     c() is a small matrix that transforms a set of symmetry
c     related functions from the AO to SO bases and a() is
c     a set of integrals that are to be transformed.
c
c     a(i,iruq,j) <= sum(opuq=1,nopuq) a(i,opuq,j)*c(iruq,opuq)
c
c     This routine does the transformation inplace via a 
c     small buffer.
c     
      integer maxbuf
      parameter (maxbuf = 64)
      double precision buf(maxbuf,8)
      integer ilo, ihi, idim, i, j, ii, opuq, iruq
c
      do j = 1, m
         do ilo = 1, n, maxbuf
            ihi = min(n,ilo+maxbuf-1)
            idim = ihi-ilo+1
            do iruq = 1, nopuq
               do i = 1, idim
                  buf(i,iruq) = 0d0
               end do
            end do
            do opuq = 1, nopuq
               do iruq = 1, nopuq
                  do i = ilo,ihi
                     ii = i-ilo+1 
                     buf(ii,iruq) = buf(ii,iruq) +
     $                    a(i,opuq,j)*c(iruq,opuq)
                  end do
               end do
            end do
            do iruq = 1, nopuq
               do i = ilo,ihi
                  ii = i-ilo+1 
                  a(i,iruq,j) = buf(ii,iruq)
               end do
            end do
         end do
      end do
c
      end
      subroutine soints_aoints(basis, nq, iq, jq, kq, lq,
     $        q4, tol, ilab, jlab, klab, llab, eri, leneri,
     $        scr, lenscr, mapi, mapj, mapk, mapl, ao, lenao)
      implicit none
#include "errquit.fh"
      integer basis
      integer nq
      integer iq(*), jq(*), kq(*), lq(*)
      integer ilab(*), jlab(*), klab(*), llab(*)
      integer mapi(*), mapj(*), mapk(*), mapl(*)
      integer leneri, lenscr, lenao
      double precision q4(*), tol, eri(*), scr(*), ao(*)
c     
      integer neri, ieri, ijkl
      logical more
      double precision block_eff
      logical intb_init4c, intb_2e4c
      external intb_init4c, intb_2e4c
c      
      if (nq .le. 0) return
c
      if (.not. intb_init4c(
     $     basis, iq, jq, basis, kq, lq,
     $     nq, q4, .true., 
     $     lenscr, scr, leneri, block_eff))
     $     call errquit('intb_init?',nq, INT_ERR)
c
 10   more = intb_2e4c(
     $     basis, iq, jq, basis, kq, lq,
     $     nq, q4, .true., tol, .false.,
     $     ilab, jlab, klab, llab, 
     $     eri, leneri, neri, lenscr, scr)
c     
c     Map the labels from NWChem order into symmetry natural order
c     
      do ieri = 1, neri
         ijkl = 1 + 
     $        mapi(ilab(ieri)) +
     $        mapj(jlab(ieri)) + 
     $        mapk(klab(ieri)) +
     $        mapl(llab(ieri))
         if (ijkl.lt.1 .or. ijkl.gt.lenao) call errquit('ijkl',ijkl,
     &       UNKNOWN_ERR)
         ao(ijkl) = eri(ieri)
      end do
c
      if (more) goto 10
c
      end
      subroutine soints_sovecs(nbf,nmo,c,order)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cbassym.fh"
      integer nbf, nmo, order
      double precision c(nbf,nmo)
c
c     On input C() contains a set of molecular orbital coefficients 
c     in the atomic orbital basis ... transform the AO index into the 
c     symmetry-adapted orbital basis.  
c
c     Order determines the output labelling
c
c     if (order == 1) then
c     .  The SOs are in natural order (not blocked).
c     else if (order == 2) then
c     .  The SOs are blocked by symmetry (second labelling scheme).
c     else
c     .  Error
c     endif
c
c     Assumes cbassym has been initialized.
c
      integer i, j, juq, nopuq, iruq, opuq
      double precision t(nw_max_nbf_small), s(8,8), sum
c
      if (order.lt.1 .or. order.gt.2) call errquit('what order?',0,
     &       INPUT_ERR)
c
      write(6,*) ' INPUT'
      call output(c,1,nbf,1,nmo,nbf,nmo,1)
c     
c     First reorder the AO functions into symmetry natural order,
c     and transform into SOs in natural order, then copy back
c     in blocked order. 
c
      do i = 1, nmo
         do j = 1, nbf
            t(bf_nat(1,j)) = c(j,i)
         end do
         j = 0
         do juq = 1, bf_n_uq
            nopuq = bf_n_uq_op(bf_uq(juq))
            call soints_uq_to_so_matrix(juq, s)
            do iruq = 1, nopuq
               sum = 0d0
               do opuq = 1, nopuq
                  sum = sum + s(iruq,opuq)*t(j+opuq)
               end do
               if (order .eq. 1) then
                  c(j+iruq,i) = sum
               else 
                  c(so_lab1(2,j+iruq),i) = sum
               endif
            end do
            j = j + nopuq
         end do
      end do
c
      write(6,*) ' OUTPUT in order ', order
      call output(c,1,nbf,1,nmo,nbf,nmo,1)
c
      end
      subroutine soints_uqshblock(uqlo, uqhi, maxnbf, blocks, nblocks)
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "cbassym.fh"
      integer uqlo, uqhi        ! [input]
      integer maxnbf            ! [input]
      integer blocks(2,*)       ! [output]
      integer nblocks           ! [output]
c
c     Group the given range of symmetry unqiue shells into blocks so
c     that there are no more than maxnbf functions arising from each
c     block.
c
c     nblocks returns the no. of blocks
c     blocks(1,*) = first symmetry unique shell in the block
c     blocks(2,*) = last ...
c
      integer ninblock, ishuq, ishuqlo, shdim, iblock
c
      ninblock = 0
      nblocks = 0
      ishuqlo = uqlo
      do ishuq = uqlo,uqhi
         shdim = sh_uq_bf(2,ishuq) - sh_uq_bf(1,ishuq) + 1
         if (shdim .gt. maxnbf) call errquit
     $        ('soints_uqshblock: shell block > maxnbf?', shdim,
     &       BASIS_ERR)
         if ((ninblock+shdim .gt. maxnbf) .or. (ishuq.eq.sh_n_uq)) then
            nblocks = nblocks + 1
            blocks(1,nblocks) = ishuqlo
            blocks(1,nblocks) = ishuq - 1
            if (ishuq .eq. sh_n_uq) blocks(1,nblocks) = ishuq
            ishuqlo = ishuq
            ninblock = 0
         else
            ninblock = ninblock + shdim
         end if
      end do
c
      write(6,*) ' UNIQUE shell blocking maxnbf = ', maxnbf
      do iblock = 1, nblocks
         write(6,1) iblock, blocks(1,iblock), blocks(2,iblock)
 1       format(1x,i5,5x,2i5)
      end do
c
      end

c $Id$

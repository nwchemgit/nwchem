      subroutine mp2_lai_uhf(rtdb,nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
     $     sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
     $     sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
     $     sym_lo_foa,sym_hi_foa,num_foa,
     $     sym_lo_fob,sym_hi_fob,num_fob,
     $     nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $     c_a,c_b,indx_v,nbf,ioff,nnbf,
     $     nva_lo_local,nva_hi_local,noa_lo,noa_hi,
     $     nob_lo,nob_hi,nva,nvb,
     $     noa,nob,t,ia_cb,ia_jk,ia_uv,tmp,tunita,tunitb,
     $     l_ai3_a,l_ai3_b,l_ai1_a,l_ai1_b,
     $     l_iy3_a, l_iy3_b,
     $     g_lai_a,g_lai_b,
     $     g_wai_a,g_wai_b,g_liy_a, g_liy_b,kunit)
*     
*     $Id$
*     
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "eaf.fh"
#include "rtdb.fh"
      integer rtdb, i
      integer noa_lo,noa_hi,nva_lo,nva_hi
      integer sym_lo_oa(0:*),sym_hi_oa(0:*)
      integer sym_lo_va(0:*),sym_hi_va(0:*)
      integer sym_lo_foa(0:7),sym_hi_foa(0:7)
      integer num_foa(0:7),num_oa(0:*),num_va(0:*)
      integer nob_lo,nob_hi,nvb_lo,nvb_hi
      integer sym_lo_ob(0:*),sym_hi_ob(0:*)
      integer sym_lo_vb(0:*),sym_hi_vb(0:*)
      integer sym_lo_fob(0:7),sym_hi_fob(0:7)
      integer num_fob(0:7),num_ob(0:*),num_vb(0:*)
      integer indx_v(*)
      integer nbf,nnbf,ioff,kunit
      integer noa,nva,nob,nvb
      integer nir
      integer a_a
      integer oseg_lo,oseg_hi
      double precision c_a(nbf,nbf),c_b(nbf,nbf)
      double precision zero,one,minusone,half
      integer syma,symi,symb,symj,symia,symc,symk,symy
      integer irs_a(*),irs_b(*)
      integer nva_lo_local,nva_hi_local
      double precision l_ai1_a(1:noa_hi,nva_lo:nva_hi)
      double precision l_ai1_b(1:nob_hi,nvb_lo:nvb_hi)
      double precision l_ai3_a(noa_lo:noa_hi,nva_lo:nva_hi)
      double precision l_ai3_b(nob_lo:nob_hi,nvb_lo:nvb_hi)
c     
      double precision l_iy3_a(noa_lo:noa_hi,1:*)
      double precision l_iy3_b(nob_lo:nob_hi,1:*)
c     
      double precision t(*),ia_cb(*),ia_jk(*),tmp(*),ia_uv(*)
      integer g_lai_a,g_lai_b,g_liy_a, g_liy_b
      integer tunita,tunitb
      double precision tunitptra,tunitptrb
      integer tcounta,tcountb
      integer msg_l_ai3_a_sum,msg_l_ai3_b_sum
      parameter(msg_l_ai3_a_sum=9992)
      parameter(msg_l_ai3_b_sum=9991)
      integer msg_l_ai1_a_sum,msg_l_ai1_b_sum
      parameter(msg_l_ai1_a_sum=9990)
      parameter(msg_l_ai1_b_sum=9989)
      integer i_batch
      integer g_wai_a,g_wai_b
c
      logical scs                         ! flag for SCS-MP2
      double precision fss,fos,mfss,mfos  ! scaling factors for SCS-MP2
c     
#include "bitops.fh"
c     
      zero=0.0d0
      one=1.0d0
      minusone=-1.0d0
      half=0.5d0
      tunitptra=1
      tunitptrb=1

c     flag and scaling factors for SCS-MP2

      if (.not. rtdb_get(rtdb, 'mp2:scs', mt_log, 1, scs))
     $  scs = .false.       ! flag for SCS gradient
      if(scs)then
        if (.not. rtdb_get(rtdb, 'mp2:fss', mt_dbl, 1, fss))
     $    fss = 1.0d0/3.0d0   ! same spin scaling factor for SCS-MP2
        if (.not. rtdb_get(rtdb, 'mp2:fos', mt_dbl, 1, fos))
     $    fos = 1.2d0         ! opposite spin scaling factor for SCS-MP2
      else
        if (.not. rtdb_get(rtdb, 'mp2:fss', mt_dbl, 1, fss))
     $    fss = one
        if (.not. rtdb_get(rtdb, 'mp2:fos', mt_dbl, 1, fos))
     $    fos = one
        if ((fss.ne.one).or.(fos.ne.one)) then
          if (ga_nodeid().eq.0) then
            write(6,*) "You set fss or fos but not scs; setting scs=T"
          endif
          scs = .true.
        endif
      endif
      mfss=-fss
      mfos=-fos

      if (noa_lo .gt. 1) then
         call dfill(noa*(noa_lo-1),zero,l_iy3_a,1)
         call dfill(nob*(nob_lo-1),zero,l_iy3_b,1)
      end if
      call dfill((noa*nva),zero,l_ai3_a,1)
      call dfill((nob*nvb),zero,l_ai3_b,1)
      call dfill((noa_hi*nva),zero,l_ai1_a,1)
      call dfill((nob_hi*nvb),zero,l_ai1_b,1)
c     
      do i=oseg_lo,oseg_hi 
         symi=irs_a(i)
         i_batch=i-oseg_lo+1
         do a_a=nva_lo_local,nva_hi_local
            syma=irs_a(a_a)
            symia=ieor(syma,symi)
            call dfill((nbf*nbf),zero,tmp,1)
            call moints_rdhfint(kunit,nva_lo_local,nva_hi_local,
     $           i_batch,a_a,nnbf,ioff,ia_uv)
            call scatter(nnbf,tmp,indx_v,ia_uv)
            call dfill((nbf*nbf),zero,ia_uv,1)
            call upper2square(nbf,tmp,ia_uv)
c     
c...  form pure spin t*(vv|vo) and t*(vo|oo) terms 
c     
            call mp2_read_tijab(nva_lo, nva_hi, irs_a, symia,
     $           num_oa, sym_hi_oa, sym_lo_oa, tunita, tunitptra, t)
c     
            tcounta=1
            do symb=0,nir-1
               symc=ieor(symia,symb)
               symj=symc
               if(num_va(symc).gt.0.and.num_va(symb).gt.0.
     $              and.num_oa(symj).gt.0)then
                  call dgemm('n','n',nbf,num_va(symc),nbf,one,
     $                 ia_uv,nbf,c_a(1,sym_lo_va(symc)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_va(symc),num_va(symb),
     $                 nbf,one,tmp,nbf,c_a(1,sym_lo_va(symb)),
     $                 nbf,zero,ia_cb,num_va(symc))
c     
c     Ljc <- sum(iab) Tjb,ia*(cb|ia)
c     
                  call dgemm('n','t',num_oa(symj),
     $                 num_va(symc),num_va(symb),
     $                 fss,t(tcounta),num_oa(symj),
     $                 ia_cb,num_va(symc),one,
     $                 l_ai3_a(sym_lo_oa(symj),sym_lo_va(symc)),noa)
               end if
               tcounta=tcounta+num_oa(symj)*num_va(symb)
            end do
c     
            if (noa_lo .gt. 1) then ! Frozen orbital piece
               tcounta=1
               do symb=0,nir-1
                  symy=ieor(symia,symb)
                  symj=symy
                  if(num_foa(symy).gt.0.and.num_va(symb).gt.0.
     $                 and.num_oa(symj).gt.0)then
                     call dgemm('n','n',nbf,num_foa(symy),nbf,one,
     $                    ia_uv,nbf,c_a(1,sym_lo_foa(symy)),nbf,
     $                    zero,tmp,nbf)
                     call dgemm('t','n',num_foa(symy),num_va(symb),
     $                    nbf,one,tmp,nbf,c_a(1,sym_lo_va(symb)),
     $                    nbf,zero,ia_cb,num_foa(symy))
c     
c     Ljy <- sum(iab) Tjb,ia*(yb|ia)
c     
                     call dgemm('n','t',num_oa(symj),
     $                    num_foa(symy),num_va(symb),
     $                    fss,t(tcounta),num_oa(symj),
     $                    ia_cb,num_foa(symy),one,
     $                    l_iy3_a(sym_lo_oa(symj),sym_lo_foa(symy)),noa)
                  end if
                  tcounta=tcounta+num_oa(symj)*num_va(symb)
               end do
            end if
c     
            tcounta=1
            do symb=0,nir-1
               symk=symb
               symj=ieor(symia,symk)
               if(num_oa(symj).gt.0.and.num_oa(symk).gt.0.
     $              and.num_va(symb).gt.0)then
                  call dgemm('n','n',nbf,num_oa(symj),nbf,one,
     $                 ia_uv,nbf,c_a(1,sym_lo_oa(symj)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_oa(symj),num_oa(symk),
     $                 nbf,one,tmp,nbf,c_a(1,sym_lo_oa(symk)),
     $                 nbf,zero,ia_jk,num_oa(symj))
c     
c     Lkb = sum(iaj) (ia|jk)*Tjb,ia
c     
                  call dgemm('t','n',num_oa(symk),
     $                 num_va(symb),num_oa(symj),
     $                 mfss,ia_jk,num_oa(symj),
     $                 t(tcounta),num_oa(symj),one,
     $                 l_ai1_a(sym_lo_oa(symk),sym_lo_va(symb)),noa_hi)
               end if
               tcounta=tcounta+num_oa(symj)*num_va(symb)
            end do
c     
            if (noa_lo .ne. 1) then ! Frozen orbital piece
               tcounta=1
               do symb=0,nir-1
                  symk=symb
                  symj=ieor(symia,symk)
                  if(num_oa(symj).gt.0.and.num_foa(symk).gt.0.
     $                 and.num_va(symb).gt.0)then
                     call dgemm('n','n',nbf,num_oa(symj),nbf,one, ! Done twice
     $                    ia_uv,nbf,c_a(1,sym_lo_oa(symj)),nbf,
     $                    zero,tmp,nbf)
                     call dgemm('t','n',num_oa(symj),num_foa(symk),
     $                    nbf,one,tmp,nbf,c_a(1,sym_lo_foa(symk)),
     $                    nbf,zero,ia_jk,num_oa(symj))
c     
c     Lyb = sum(iaj) (ia|jy)*Tjb,ia
c     
                     call dgemm('t','n',num_foa(symk),
     $                    num_va(symb),num_oa(symj),
     $                    mfss,ia_jk,num_oa(symj),
     $                    t(tcounta),num_oa(symj),one,
     $                    l_ai1_a(sym_lo_foa(symk),sym_lo_va(symb)),
     $                    noa_hi)
                  end if
                  tcounta=tcounta+num_oa(symj)*num_va(symb)
               end do
            end if
c     
c...  form mixed spin t*(vv|vo) and t*(vo|oo) terms
c     
            call mp2_read_tijab(nvb_lo, nvb_hi, irs_b, symia,
     $           num_ob, sym_hi_ob, sym_lo_ob, tunitb, tunitptrb, t)
c     
            tcountb=1
            do symb=0,nir-1
               symc=ieor(symia,symb)
               symj=symc
               if(num_vb(symc).gt.0.and.num_vb(symb).gt.0.
     $              and.num_ob(symj).gt.0)then
                  call dgemm('n','n',nbf,num_vb(symc),nbf,one,
     $                 ia_uv,nbf,c_b(1,sym_lo_vb(symc)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_vb(symc),num_vb(symb),
     $                 nbf,one,tmp,nbf,c_b(1,sym_lo_vb(symb)),
     $                 nbf,zero,ia_cb,num_vb(symc))
c     
c     Ljc <- sum(iab) Tjb,ia*(cb|ia)
c     
                  call dgemm('n','t',num_ob(symj),
     $                 num_vb(symc),num_vb(symb),
     $                 fos,t(tcountb),num_ob(symj),
     $                 ia_cb,num_vb(symc),one,
     $                 l_ai3_b(sym_lo_ob(symj),sym_lo_vb(symc)),nob)
               end if
               tcountb=tcountb+num_ob(symj)*num_vb(symb)
            end do
c     
            if (nob_lo .gt. 1) then ! Frozen orbital piece
               tcountb=1
               do symb=0,nir-1
                  symy=ieor(symia,symb)
                  symj=symy
                  if(num_fob(symy).gt.0.and.num_vb(symb).gt.0.
     $                 and.num_ob(symj).gt.0)then
                     call dgemm('n','n',nbf,num_fob(symy),nbf,one,
     $                    ia_uv,nbf,c_b(1,sym_lo_fob(symy)),nbf,
     $                    zero,tmp,nbf)
                     call dgemm('t','n',num_fob(symy),num_vb(symb),
     $                    nbf,one,tmp,nbf,c_b(1,sym_lo_vb(symb)),
     $                    nbf,zero,ia_cb,num_fob(symy))
c     
c     Ljy <- sum(iab) Tjb,ia*(yb|ia)
c     
                     call dgemm('n','t',num_ob(symj),
     $                    num_fob(symy),num_vb(symb),
     $                    fos,t(tcountb),num_ob(symj),
     $                    ia_cb,num_fob(symy),one,
     $                    l_iy3_b(sym_lo_ob(symj),sym_lo_fob(symy)),nob)
                  end if
                  tcountb=tcountb+num_ob(symj)*num_vb(symb)
               end do
            end if
c     
            tcountb=1
            do symb=0,nir-1
               symk=symb
               symj=ieor(symia,symk)
               if(num_ob(symj).gt.0.and.num_ob(symk).gt.0.
     $              and.num_vb(symb).gt.0)then
                  call dgemm('n','n',nbf,num_ob(symj),nbf,one,
     $                 ia_uv,nbf,c_b(1,sym_lo_ob(symj)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_ob(symj),num_ob(symk),
     $                 nbf,one,tmp,nbf,c_b(1,sym_lo_ob(symk)),
     $                 nbf,zero,ia_jk,num_ob(symj))
c     
c     Lkb = sum(iaj) (ia|kj)*Tjb,ia
c     
                  call dgemm('t','n',num_ob(symk),
     $                 num_vb(symb),num_ob(symj),
     $                 mfos,ia_jk,num_ob(symj),
     $                 t(tcountb),num_ob(symj),one,
     $                 l_ai1_b(sym_lo_ob(symk),sym_lo_vb(symb)),
     $                 nob_hi)
               end if
               tcountb=tcountb+num_ob(symj)*num_vb(symb)
            end do
c     
            if (nob_lo .ne. 1) then ! Frozen orbital piece
               tcountb=1
               do symb=0,nir-1
                  symk=symb
                  symj=ieor(symia,symk)
                  if(num_ob(symj).gt.0.and.num_fob(symk).gt.0.
     $                 and.num_vb(symb).gt.0)then
                     call dgemm('n','n',nbf,num_ob(symj),nbf,one, ! Done twice
     $                    ia_uv,nbf,c_b(1,sym_lo_ob(symj)),nbf,
     $                    zero,tmp,nbf)
                     call dgemm('t','n',num_ob(symj),num_fob(symk),
     $                    nbf,one,tmp,nbf,c_b(1,sym_lo_fob(symk)),
     $                    nbf,zero,ia_jk,num_ob(symj))
c     
c     Lyb = sum(iaj) (ia|jy)*Tjb,ia
c     
                     call dgemm('t','n',num_fob(symk),
     $                    num_vb(symb),num_ob(symj),
     $                    mfos,ia_jk,num_ob(symj),
     $                    t(tcountb),num_ob(symj),one,
     $                    l_ai1_b(sym_lo_fob(symk),sym_lo_vb(symb)),
     $                    nob_hi)
                  end if
                  tcountb=tcountb+num_ob(symj)*num_vb(symb)
               end do
            end if
         end do
      end do
c     
      call ga_dgop(msg_l_ai3_a_sum,l_ai3_a,noa*nva,'+')
      call ga_dgop(msg_l_ai3_b_sum,l_ai3_b,nob*nvb,'+')
      call ga_dgop(msg_l_ai1_a_sum,l_ai1_a,noa_hi*nva,'+')
      call ga_dgop(msg_l_ai1_b_sum,l_ai1_b,nob_hi*nvb,'+')
      if (noa_lo .gt. 1) then
         call ga_dgop(msg_l_ai3_a_sum,l_iy3_a,noa*(noa_lo-1),'+')
         call ga_dgop(msg_l_ai3_b_sum,l_iy3_b,nob*(nob_lo-1),'+')
      end if
c
c$$$      write(6,*) ' l_ai1_a '
c$$$      call output(l_ai1_a,1,noa_hi,1,nva,noa_hi,nva,1)
c$$$      write(6,*) ' l_ai1_b '
c$$$      call output(l_ai1_b,1,nob_hi,1,nvb,nob_hi,nvb,1)
c$$$c
c$$$      write(6,*) ' l_ai3_a '
c$$$      call output(l_ai3_a,1,noa,1,noa_lo-1,noa,noa_lo-1,1)
c$$$      write(6,*) ' l_ai3_b '
c$$$      call output(l_ai3_b,1,nob,1,nob_lo-1,nob,nob_lo-1,1)
c      
      call ga_sync()
      if(ga_nodeid().eq.0)then
         call ga_acc(g_lai_a,noa_lo,noa_hi,1,nva,l_ai3_a,noa,one)
         call ga_acc(g_lai_b,nob_lo,nob_hi,1,nvb,l_ai3_b,nob,one)
         call ga_acc(g_lai_a,1,noa_hi,1,nva,l_ai1_a,noa_hi,one)
         call ga_acc(g_lai_b,1,nob_hi,1,nvb,l_ai1_b,nob_hi,one)
         call ga_acc(g_wai_a,1,noa_hi,1,nva,l_ai1_a,noa_hi,minusone)
         call ga_acc(g_wai_b,1,nob_hi,1,nvb,l_ai1_b,nob_hi,minusone)
         if (noa_lo .gt. 1) then
            call ga_acc(g_liy_a,1,noa,1,noa_lo-1,l_iy3_a,noa,one)
            call ga_acc(g_liy_b,1,nob,1,nob_lo-1,l_iy3_b,nob,one)
         end if
      end if
c     
      end



      subroutine mp2_lai_uhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
     $     sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
     $     sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
     $     nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $     c_a,c_b,indx_v,nbf,ioff,nnbf,
     $     nva_lo_local,nva_hi_local,noa_lo,noa_hi,
     $     nob_lo,nob_hi,nva,nvb,
     $     noa,nob,t_a,t_b,ia_cb,ia_jk,ia_uv,tmp,tunita,tunitb,
     $     l_ai_a,l_ai_b,g_lai_a,g_lai_b)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer i,c,jk
      integer noa_lo,noa_hi,nva_lo,nva_hi
      integer sym_lo_oa(0:*),sym_hi_oa(0:*)
      integer sym_lo_va(0:*),sym_hi_va(0:*)
      integer num_oa(0:*),num_va(0:*)
      integer nob_lo,nob_hi,nvb_lo,nvb_hi
      integer sym_lo_ob(0:*),sym_hi_ob(0:*)
      integer sym_lo_vb(0:*),sym_hi_vb(0:*)
      integer num_ob(0:*),num_vb(0:*)
      integer indx_v(*)
      integer nbf,nnbf,ioff,kunit
      parameter(kunit=11)
      integer noa,nva,nob,nvb,counter
      integer nir
      integer a_a
      integer oseg_lo,oseg_hi
      double precision c_a(nbf,nbf),c_b(nbf,nbf)
      double precision zero,one,minusone,half
      integer syma,symi,symb,symj,symia,symc,symk,symik,symib
      integer irs_a(*),irs_b(*),b,j,k
      integer nva_lo_local,nva_hi_local
      double precision l_ai_a(*),l_ai_b(*)
      double precision t_a(*),t_b(*),ia_cb(*),ia_jk(*),tmp(*),ia_uv(*)
      integer g_lai_a,g_lai_b
      integer counta,countb,tunitptra,tunitptrb,tunita,tunitb,ierr
      integer tcounta,tcountb
      integer msg_l_ai_a_sum,msg_l_ai_b_sum
      parameter(msg_l_ai_a_sum=9992)
      parameter(msg_l_ai_b_sum=9991)
      integer i_batch
      integer laiacnt,laibcnt
c     
      logical do_3v,do_1v
c     
      do_3v=.true.
      do_1v=.true.
      zero=0.0d0
      one=1.0d0
      minusone=-1.0d0
      half=0.5d0
      tunitptra=1
      tunitptrb=1
      call dfill((noa*nva),zero,l_ai_a,1)
      call dfill((nob*nvb),zero,l_ai_b,1)
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
c...  form pure spin t*(vv|vo) and t*(vo|oo) terms 
            do symb=0,nir-1
               symc=ieor(symia,symb)
               if(num_va(symc).gt.0.and.num_va(symb).gt.0)then
                  call dgemm('n','n',nbf,num_va(symc),nbf,one,
     $                 ia_uv,nbf,c_a(1,sym_lo_va(symc)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_va(symc),num_va(symb),
     $                 nbf,one,tmp,nbf,c_a(1,sym_lo_va(symb)),
     $                 nbf,zero,ia_cb,num_va(symc))
               endif
            enddo
            do symj=0,nir-1
               symk=ieor(symia,symj)
               if(num_oa(symj).gt.0.and.num_oa(symk).gt.0)then
                  call dgemm('n','n',nbf,num_oa(symj),nbf,one,
     $                 ia_uv,nbf,c_a(1,sym_lo_oa(symj)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_oa(symk),num_oa(symj),
     $                 nbf,one,tmp,nbf,c_a(1,sym_lo_oa(symk)),
     $                 nbf,zero,ia_jk,num_oa(symk))
               endif
            enddo
            tcounta=1
            do b=nva_lo,nva_hi
               symb=irs_a(b)
               symk=ieor(symb,symia)
               if(num_oa(symk).ne.0)then
                  counta=sym_hi_oa(symk)-sym_lo_oa(symk)+1
                  call getwa(tunita,t_a(tcounta),
     $                 tunitptra,counta,ierr)
                  if(ierr.ne.0)call errquit(
     $                 'could not read t_a back',ierr)
                  tunitptra=tunitptra+counta
                  tcounta=tcounta+counta
               endif
            enddo
            if(do_3v)then
               tcounta=1
               laiacnt=1
               do symb=0,nir-1
                  symk=ieor(symia,symb)
                  symc=symk
                  if(num_va(symc).gt.0.and.num_oa(symk).gt.0.and.
     $                 num_va(symb).gt.0)then
                     call dgemm('n','t',num_oa(symk),
     $                    num_va(symc),num_va(symb),
     $                    one,t_a(tcounta),num_oa(symk),
     $                    ia_cb,num_va(symc),one,
     $                    l_ai_a(laiacnt),num_oa(symk))
                     tcounta=tcounta+num_oa(symk)*num_va(symc)
                     laiacnt=laiacnt+num_oa(symk)*num_va(symc)
                  endif
               enddo
            endif
            if(do_1v)then
               tcounta=1
               laiacnt=1
               do symk=0,nir-1
                  symj=ieor(symk,symia)
                  symb=symj
                  if(num_oa(symk).gt.0.and.num_va(symb).gt.0.and.
     $                 num_oa(symj).gt.0)then
                     call dgemm('t','n',num_oa(symj),
     $                    num_va(symb),num_oa(symk),
     $                    minusone,ia_jk,num_oa(symk),
     $                    t_a(tcounta),num_oa(symk),
     $                    one,l_ai_a(laiacnt),num_oa(symj))
                     tcounta=tcounta+num_oa(symj)*num_va(symb)
                     laiacnt=laiacnt+num_oa(symj)*num_va(symb)
                  endif
               enddo
            endif
c...  form mixed spin t*(vv|vo) and t*(vo|oo) terms
            do symb=0,nir-1
               symc=ieor(symia,symb)
               if(num_vb(symc).gt.0.and.num_vb(symb).gt.0)then
                  call dgemm('n','n',nbf,num_vb(symc),nbf,one,
     $                 ia_uv,nbf,c_b(1,sym_lo_vb(symc)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_vb(symc),num_vb(symb),
     $                 nbf,one,tmp,nbf,c_b(1,sym_lo_vb(symb)),
     $                 nbf,zero,ia_cb,num_vb(symc))
               endif
            enddo
            do symj=0,nir-1
               symk=ieor(symia,symj)
               if(num_ob(symk).gt.0.and.num_ob(symj).gt.0)then
                  call dgemm('n','n',nbf,num_ob(symk),nbf,one,
     $                 ia_uv,nbf,c_b(1,sym_lo_ob(symk)),nbf,
     $                 zero,tmp,nbf)
                  call dgemm('t','n',num_ob(symk),num_ob(symj),
     $                 nbf,one,tmp,nbf,c_b(1,sym_lo_ob(symj)),
     $                 nbf,zero,ia_jk,num_ob(symk))
               endif
            enddo
            tcountb=1
            do b=nvb_lo,nvb_hi
               symb=irs_b(b)
               symk=ieor(symb,symia)
               if(num_ob(symk).ne.0)then
                  countb=sym_hi_ob(symk)-sym_lo_ob(symk)+1
                  call getwa(tunitb,t_b(tcountb),
     $                 tunitptrb,countb,ierr)
                  if(ierr.ne.0)call errquit(
     $                 'could not read t_b back',ierr)
                  tunitptrb=tunitptrb+countb
                  tcountb=tcountb+countb
               endif
            enddo
            if(do_3v)then
               tcountb=1
               laibcnt=1
               do symb=0,nir-1
                  symk=ieor(symia,symb)
                  symc=symk
                  if(num_vb(symc).gt.0.and.num_ob(symk).gt.0.and.
     $                 num_ob(symb).gt.0)then
                     call dgemm('n','t',num_ob(symk),
     $                    num_vb(symc),num_vb(symb),
     $                    one,t_b(tcountb),num_ob(symk),
     $                    ia_cb,num_vb(symc),one,
     $                    l_ai_b(laibcnt),num_ob(symk))
                     tcountb=tcountb+num_ob(symk)*num_vb(symc)
                     laibcnt=laibcnt+num_ob(symk)*num_vb(symc)
                  endif
               enddo
            endif
            if(do_1v)then
               tcountb=1
               laibcnt=1
               do symk=0,nir-1
                  symj=ieor(symk,symia)
                  symb=symj
                  if(num_ob(symk).gt.0.and.num_vb(symb).gt.0.and.
     $                 num_ob(symj).gt.0)then
                     call dgemm('t','n',num_ob(symj),
     $                    num_vb(symb),num_ob(symk),
     $                    minusone,ia_jk,num_ob(symk),
     $                    t_b(tcountb),num_ob(symk),
     $                    one,l_ai_b(laibcnt),num_ob(symj))
                     tcountb=tcountb+num_ob(symj)*num_vb(symb)
                     laibcnt=laibcnt+num_ob(symj)*num_vb(symb)
                  endif
               enddo
            endif
         enddo
      enddo
      call ga_dgop(msg_l_ai_a_sum,l_ai_a,noa*nva,'+')
      call ga_dgop(msg_l_ai_b_sum,l_ai_b,nob*nvb,'+')
      call ga_sync()
      if(ga_nodeid().eq.0)then
         call ga_acc(g_lai_a,1,noa,1,nva,l_ai_a,noa,one)
         call ga_acc(g_lai_b,1,nob,1,nvb,l_ai_b,nob,one)
      endif
      end



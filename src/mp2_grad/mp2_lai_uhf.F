      subroutine mp2_lai_uhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
     $     sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,
     $     sym_lo_ob,sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,
     $     nir,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $     nva_lo_local,nva_hi_local,noa_lo_local,
     $     noa_hi_local,noa_lo,noa_hi,nob_lo,nob_hi,nva,nvb,
     $     noa,nob,t_a,t_b,ia_cb,ij_ka,vvvo_off_a,vvvo_off_b,
     $     g_vvvo_a,g_vvvo_b,g_iajk_a,g_iajk_b,ioff_vooo_a,
     $     ioff_vooo_b,tunita,tunitb,p_ij_a,p_ij_b,p_ab_a,p_ab_b,
     $     l_ai_a,l_ai_b,g_lai_a,g_lai_b,g_p_ij_a,g_p_ij_b,
     $     g_p_ab_a,g_p_ab_b)
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
      integer noa,nva,nob,nvb,counter
      integer nir
      integer sym_tot_c,sym_tot_k
      integer vvvo_off_a(nva*nva,0:nir-1),vvvo_off_b(nvb*nvb,0:nir-1)
      integer ioff_vooo_a(noa*noa,0:nir-1),ioff_vooo_b(nob*nob,0:nir-1)
      integer a_a
      integer oseg_lo,oseg_hi
      double precision zero,one,minusone,half
      integer syma,symi,symb,symj,symia,symc,symk,symik,symib
      integer irs_a(*),irs_b(*),b,j,k
      integer nva_lo_local,nva_hi_local,noa_lo_local,noa_hi_local,cb,jb
      double precision l_ai_a(noa,nva),l_ai_b(nob,nvb)
      double precision t_a(*),t_b(*),ia_cb(*),ij_ka(*)
      double precision p_ij_a(*),p_ij_b(*),p_ab_a(*),p_ab_b(*)
      integer g_vvvo_a,g_vvvo_b
      integer g_iajk_a,g_iajk_b
      integer g_lai_a,g_lai_b,i_batch
      integer counta,countb,tunitptra,tunitptrb,tunita,tunitb,ierr
      integer g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b
c     
      zero=0.0d0
      one=1.0d0
      minusone=-1.0d0
      half=0.5d0
      tunitptra=1
      tunitptrb=1
      call ga_zero(g_lai_a)
      call ga_zero(g_lai_b)
      do i=oseg_lo,oseg_hi 
         symi=irs_a(i)
         i_batch=i-oseg_lo+1
         do syma=0,nir-1
            counter=0
            do symb=0,nir-1
               do b=sym_lo_va(symb),sym_hi_va(symb)
                  symia=ieor(syma,symi)
                  symc=ieor(symia,symb)
                  vvvo_off_a(b,syma)=counter
                  counter=counter+num_va(symc)
               enddo
            enddo
         enddo
         do syma=0,nir-1
            counter=0
            do symb=0,nir-1
               do b=sym_lo_vb(symb),sym_hi_vb(symb)
                  symia=ieor(syma,symi)
                  symc=ieor(symia,symb)
                  vvvo_off_b(b,syma)=counter
                  counter=counter+num_vb(symc)
               enddo
            enddo
         enddo
         do syma=0,nir-1
            counter=0
            do symk=0,nir-1
               do k=sym_lo_oa(symk),sym_hi_oa(symk)
                  symia=ieor(syma,symi)
                  symj=ieor(symia,symk)
                  ioff_vooo_a(k,syma)=counter
                  counter=counter+num_oa(symj)
               enddo
            enddo
         enddo
         do syma=0,nir-1
            counter=0
            do symk=0,nir-1
               do k=sym_lo_ob(symk),sym_hi_ob(symk)
                  symia=ieor(syma,symi)
                  symj=ieor(symia,symk)
                  ioff_vooo_b(k,syma)=counter
                  counter=counter+num_ob(symj)
               enddo
            enddo
         enddo
         do a_a=nva_lo_local,nva_hi_local
            syma=irs_a(a_a)
            symia=ieor(syma,symi)
            do b=nva_lo,nva_hi
               symb=irs_a(b)
               symk=ieor(symb,symia)
               if(num_oa(symk).ne.0)then
                  counta=sym_hi_oa(symk)-sym_lo_oa(symk)+1
                  call getwa(tunita,t_a(sym_lo_oa(symk)),
     $                 tunitptra,counta,ierr)
                  if(ierr.ne.0)call errquit(
     $                 'could not read t_a back',ierr)
                  tunitptra=tunitptra+counta
               endif
            enddo
            do b=nvb_lo,nvb_hi
               symb=irs_b(b)
               symk=ieor(symb,symia)
               if(num_ob(symk).ne.0)then
                  countb=sym_hi_ob(symk)-sym_lo_ob(symk)+1
                  call getwa(tunitb,t_b(sym_lo_ob(symk)),
     $                 tunitptrb,countb,ierr)
                  if(ierr.ne.0)call errquit(
     $                 'could not read t_b back',ierr)
                  tunitptrb=tunitptrb+countb
               endif
            enddo
            do symb=0,nir-1
               symk=ieor(symia,symb)
               symc=ieor(symia,symb)
               sym_tot_k=ieor(symk,ieor(symia,symb))
               sym_tot_c=ieor(symc,ieor(symia,symb))
               if(num_va(symb).gt.0.and.num_oa(symk).gt.0.and.
     $              num_va(symc).gt.0)then
                  cb=1
                  do b=nva_lo,nva_hi
                     call ga_get(g_vvvo_a,vvvo_off_a(b,syma)+1,
     $                    vvvo_off_a(b,syma)+num_va(symc),
     $                    a_a-nva_lo+1,a_a-nva_lo+1,
     $                    ia_cb(cb),1)
                     cb=cb+num_va(symc)
                  enddo
                  call dgemm('n','t',num_oa(symk),
     $                 num_va(symc),num_va(symb),one,
     $                 t_a,num_oa(symk),ia_cb,num_va(symc),
     $                 half,l_ai_a,num_oa(symk))
                  do c=sym_lo_va(symc),sym_hi_va(symc)
                     call ga_acc(g_lai_a,sym_lo_oa(symk),
     $                    sym_hi_oa(symk),c-sym_lo_va(symc)+1,
     $                    c-sym_lo_va(symc)+1,
     $                    l_ai_a(1,c-nva_lo+1),1,one)
                  enddo
               endif
               if(num_vb(symb).gt.0.and.num_ob(symk).gt.0.and.
     $              num_vb(symc).gt.0)then
                  cb=1
                  do b=nvb_lo,nvb_hi
                     call ga_get(g_vvvo_b,vvvo_off_b(b,syma)+1,
     $                    vvvo_off_b(b,syma)+num_vb(symc),
     $                    a_a-nva_lo+1,a_a-nva_lo+1,
     $                    ia_cb(cb),1)
                     cb=cb+num_vb(symc)
                  enddo
                  call dgemm('n','t',num_ob(symk),
     $                 num_vb(symc),num_vb(symb),one,
     $                 t_b,num_ob(symk),ia_cb,num_vb(symc),
     $                 half,l_ai_b,num_ob(symk))
                  do c=sym_lo_vb(symc),sym_hi_vb(symc)
                     call ga_acc(g_lai_b,sym_lo_ob(symk),
     $                    sym_hi_ob(symk),c-sym_lo_vb(symc)+1,
     $                    c-sym_lo_vb(symc)+1,
     $                    l_ai_b(1,c-nvb_lo+1),1,1.0d0)
                  enddo
               endif
            enddo
            do symk=0,nir-1
               symj=ieor(symia,symk)
               symc=ieor(symia,symk)
               if(num_oa(symj).gt.0.and.num_oa(symk).gt.0.and.
     $              num_va(symc).gt.0)then
                  jk=1
                  do k=sym_lo_oa(symk),sym_hi_oa(symk)
                     call ga_get(g_iajk_a,ioff_vooo_a(k,syma)+1,
     $                    ioff_vooo_a(k,syma)+num_oa(symj),
     $                    a_a-nva_lo+1,a_a-nva_lo+1,
     $                    ij_ka(jk),1)
                     jk=jk+num_oa(symj)
                  enddo 
                  call dgemm('t','n',num_va(symc),
     $                 num_oa(symk),num_oa(symj),one,
     $                 ij_ka,num_oa(symj),t_a,num_oa(symj),
     $                 half,l_ai_a,num_va(symc))
                  do c=sym_lo_va(symc),sym_hi_va(symc)
                     call ga_acc(g_lai_a,sym_lo_oa(symk),
     $                    sym_hi_oa(symk),c-sym_lo_va(symc)+1,
     $                    c-sym_lo_va(symc)+1,
     $                    l_ai_a(1,c-nva_lo+1),1,one)
                  enddo
               endif
               if(num_ob(symj).gt.0.and.num_ob(symk).gt.0.and.
     $              num_vb(symc).gt.0)then
                  jk=1
                  do k=sym_lo_ob(symk),sym_hi_ob(symk)
                     call ga_get(g_iajk_b,ioff_vooo_b(k,syma)+1,
     $                    ioff_vooo_b(k,syma)+num_ob(symj),
     $                    a_a-nva_lo+1,a_a-nva_lo+1,
     $                    ij_ka(jk),1)
                     jk=jk+num_ob(symj)
                  enddo
                  call dgemm('t','n',num_vb(symc),
     $                 num_ob(symk),num_ob(symj),one,
     $                 ij_ka,num_ob(symj),t_b,num_ob(symj),
     $                 half,l_ai_b,num_vb(symc))
                  do c=sym_lo_vb(symc),sym_hi_vb(symc)
                     call ga_acc(g_lai_b,sym_lo_ob(symk),
     $                    sym_hi_ob(symk),c-sym_lo_vb(symc)+1,
     $                    c-sym_lo_vb(symc)+1,
     $                    l_ai_b(1,c-nvb_lo+1),1,1.0d0)
                  enddo
               endif
            enddo
         enddo
      enddo
c$$$      call ga_print(g_lai_a)
c$$$      call ga_print(g_lai_b)
      call ga_sync()
      end



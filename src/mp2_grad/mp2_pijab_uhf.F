      subroutine mp2_pijab_uhf(nva_lo,nva_hi,sym_lo_oa,sym_hi_oa,
     $     sym_lo_va,sym_hi_va,num_oa,nvb_lo,nvb_hi,sym_lo_ob,
     $     sym_hi_ob,num_ob,sym_lo_vb,sym_hi_vb,nir,num_va,num_vb,
     $     oseg_lo,oseg_hi,irs_a,irs_b,nva_lo_local,nva_hi_local,
     $     noa_lo,noa_hi,nob_lo,nob_hi,t_a,t_b,
     $     tunita,tunitb,p_ij_a,p_ij_b,p_ab_a,p_ab_b,eval_a,eval_b,
     $     g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer i
      integer noa_lo,noa_hi,nva_lo,nva_hi
      integer sym_lo_oa(0:*),sym_hi_oa(0:*)
      integer sym_lo_va(0:*),sym_hi_va(0:*)
      integer num_oa(0:*),num_va(0:*)
      integer nob_lo,nob_hi,nvb_lo,nvb_hi
      integer sym_lo_ob(0:*),sym_hi_ob(0:*)
      integer sym_lo_vb(0:*),sym_hi_vb(0:*)
      integer num_ob(0:*),num_vb(0:*)
      integer noa,nva,nob,nvb,counter
      integer a_a,tunita,tunitb
      integer oseg_lo,oseg_hi
      integer nir
      double precision zero,one,minusone,half,mhalf
      integer syma,symi,symb,symj,symia
      integer irs_a(*),irs_b(*),b,j
      integer nva_lo_local,nva_hi_local,jb
      double precision p_ij_a(*),p_ab_a(*),p_ij_b(*),p_ab_b(*)
      double precision t_a(*),t_b(*),te
      double precision denom,epsilon
      integer tunitptra,tunitptrb,counta,countb,ierr,length
      integer g_p_ij_a,g_p_ij_b,g_p_ab_a,g_p_ab_b
      integer msg_p_ij_a_sum,msg_p_ij_b_sum
      integer msg_p_ab_a_sum,msg_p_ab_b_sum
      parameter(msg_p_ij_a_sum=10000)
      parameter(msg_p_ij_b_sum=9999)
      parameter(msg_p_ab_a_sum=9998)
      parameter(msg_p_ab_b_sum=9997)
      double precision eval_a(*),eval_b(*)
      double precision eaa,eab
      integer tcounta,tcountb
      integer pijacnt,pijbcnt,pabacnt,pabbcnt
c     
      eaa=0.0d0
      eab=0.0d0
      zero=0.0d0
      one=1.0d0
      half=0.5d0
      minusone=-1.0d0
      mhalf=-0.5d0
      noa=noa_hi-noa_lo+1
      nva=nva_hi-nva_lo+1
      nob=nob_hi-nob_lo+1
      nvb=nvb_hi-nvb_lo+1
      tunitptra=1
      tunitptrb=1
      call dfill((noa*noa),zero,p_ij_a,1)
      call dfill((nva*nva),zero,p_ab_a,1)
      call dfill((nob*nob),zero,p_ij_b,1)
      call dfill((nvb*nvb),zero,p_ab_b,1)
      do i=oseg_lo,oseg_hi 
         symi=irs_a(i)
         do a_a=nva_lo_local,nva_hi_local 
            syma=irs_a(a_a)
            symia=ieor(syma,symi)
            tcounta=1
            do b=nva_lo,nva_hi
               symb=irs_a(b)
               symj=ieor(symb,symia)
               if(num_oa(symj).ne.0)then
                  counta=sym_hi_oa(symj)-sym_lo_oa(symj)+1
                  call getwa(tunita,t_a(tcounta),
     $                 tunitptra,counta,ierr)
                  if(ierr.ne.0)call errquit(
     $                 'problem reading t_a back in',ierr)
                  tunitptra=tunitptra+counta
                  do j=sym_lo_oa(symj),sym_hi_oa(symj)
                     denom=eval_a(a_a)+eval_a(b)
     $                    -eval_a(i)-eval_a(j)
                     eaa=eaa+t_a(tcounta+j-sym_lo_oa(symj))
     $                    *t_a(tcounta+j-sym_lo_oa(symj))*denom
                  enddo
                  tcounta=tcounta+counta
               endif
            enddo
            tcountb=1
            do b=nvb_lo,nvb_hi
               symb=irs_b(b)
               syma=irs_a(a_a)
               symj=ieor(symb,ieor(symi,syma))
               if(num_ob(symj).ne.0)then
                  countb=sym_hi_ob(symj)-sym_lo_ob(symj)+1
                  call getwa(tunitb,t_b(tcountb),
     $                 tunitptrb,countb,ierr)
                  if(ierr.ne.0)call errquit(
     $                 'problem reading t_b back in',ierr)
                  tunitptrb=tunitptrb+countb
                  do j=sym_lo_ob(symj),sym_hi_ob(symj)
                     denom=eval_a(a_a)+eval_b(b)
     $                    -eval_a(i)-eval_b(j)
                     eab=eab+t_b(tcountb+j-sym_lo_ob(symj))
     $                    *t_b(tcountb+j-sym_lo_ob(symj))*denom
                  enddo
                  tcountb=tcountb+countb
               endif
            enddo
            tcounta=1
            tcountb=1
            pijacnt=1
            pabacnt=1
            do symb=0,nir-1
               symj=ieor(symia,symb)
               if(num_oa(symj).gt.0.and.num_va(symb).gt.0)then
                  call dgemm('n','t',num_oa(symj),
     $                 num_oa(symj),num_va(symb),mhalf,
     $                 t_a(tcounta),num_oa(symj),
     $                 t_a(tcounta),num_oa(symj),
     $                 one,p_ij_a(pijacnt),num_oa(symj))
                  pijacnt=pijacnt+num_oa(symi)*num_oa(symj)
                  call dgemm('t','n',num_va(symb),
     $                 num_va(symb),num_oa(symj),half,
     $                 t_a(tcounta),num_oa(symj),
     $                 t_a(tcounta),num_oa(symj),
     $                 one,p_ab_a(pabacnt),num_va(symb))
                  tcounta=tcounta+num_oa(symj)*num_va(symb)
                  pabacnt=pabacnt+num_va(syma)*num_va(symb)
               endif
               pijbcnt=1
               pabbcnt=1
               if(num_ob(symj).gt.0.and.num_vb(symb).gt.0)then
                  call dgemm('n','t',num_ob(symj),
     $                 num_ob(symj),num_vb(symb),minusone,
     $                 t_b(tcountb),num_ob(symj),
     $                 t_b(tcountb),num_ob(symj),
     $                 one,p_ij_b(pijbcnt),num_ob(symj))
                  pijbcnt=pijbcnt+num_ob(symi)*num_ob(symj)
                  call dgemm('t','n',num_vb(symb),
     $                 num_vb(symb),num_ob(symj),one,
     $                 t_b(tcountb),num_ob(symj),
     $                 t_b(tcountb),num_ob(symj),
     $                 one,p_ab_b(pabbcnt),num_vb(symb))
                  tcountb=tcountb+num_ob(symj)*num_vb(symb)
                  pabbcnt=pabbcnt+num_vb(syma)*num_vb(symb)
               endif
            enddo
         enddo
      enddo
      call ga_dgop(1,eaa,1,'+')
      call ga_dgop(1,eab,1,'+')
      call ga_dgop(msg_p_ij_a_sum,p_ij_a,noa*noa,'+')
      call ga_dgop(msg_p_ab_a_sum,p_ab_a,nva*nva,'+')
      call ga_dgop(msg_p_ij_b_sum,p_ij_b,nob*nob,'+')
      call ga_dgop(msg_p_ab_b_sum,p_ab_b,nvb*nvb,'+')
      call ga_sync()
      if(ga_nodeid().eq.0)then
         call ga_acc(g_p_ij_a,1,noa,1,noa,p_ij_a,noa,one)
         call ga_acc(g_p_ij_b,1,nob,1,nob,p_ij_b,nob,one)
         call ga_acc(g_p_ab_a,1,nva,1,nva,p_ab_a,nva,one)
         call ga_acc(g_p_ab_b,1,nvb,1,nvb,p_ab_b,nvb,one)
      endif
      end
      
      
      
      
      
      



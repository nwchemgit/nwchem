      subroutine mp2_make_tuhf(nbf,noa_lo,noa_hi,nva_lo,nva_hi,
     $     sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,num_oa,num_va,
     $     nob_lo,nob_hi,nvb_lo,nvb_hi, sym_lo_ob,sym_hi_ob,
     $     sym_lo_vb,sym_hi_vb,num_ob,num_vb,   
     $     oseg_lo,oseg_hi,irs_a,c_a,irs_b,c_b,ia_uv,ia_jb,
     $     eval_a,eval_b,toff_a,toff_b,nva,noa,nvb,nob,nir,
     $     tmp,ioff,nnbf,indx_v,nva_lo_local,nva_hi_local,
     $     e2,task,tunita,tunitb)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
      character*8 task
      integer nbf,i
      integer noa_lo,noa_hi,nva_lo,nva_hi
      integer sym_lo_oa(0:*),sym_hi_oa(0:*)
      integer sym_lo_va(0:*),sym_hi_va(0:*)
      integer num_oa(0:*),num_va(0:*)
      integer nob_lo,nob_hi,nvb_lo,nvb_hi
      integer sym_lo_ob(0:*),sym_hi_ob(0:*)
      integer sym_lo_vb(0:*),sym_hi_vb(0:*)
      integer num_ob(0:*),num_vb(0:*)
      integer noa,nva,nob,nvb,counter,g_t_a,g_t_b
      integer a_a,tunita,tunitb,tunitptra,tunitptrb,counta,countb,ierr
      integer oseg_lo,oseg_hi
      integer nir
      integer toff_a(noa*nva,0:nir-1),toff_b(nob*nvb,0:nir-1)
      double precision eaa,eab
      double precision c_a(nbf,nbf),c_b(nbf,nbf)
      integer i_batch
      double precision e2,denom,t,minusone
      integer syma,symi,symb,symj,symia
      integer irs_a(*),irs_b(*),b,j
      integer indx_v(*)
      integer nnbf,ioff,kunit
      parameter(kunit=11)
      integer nva_lo_local,nva_hi_local,jb
      double precision tmp(*),ia_uv(*),ia_jb(*)
      double precision eval_a(*),eval_b(*)
c     
      eaa=0.0d0
      eab=0.0d0
      minusone=-1.0d0
      noa=noa_hi-noa_lo+1
      nva=nva_hi-nva_lo+1
      nob=nob_hi-nob_lo+1
      nvb=nvb_hi-nvb_lo+1
      if(.not.ga_create(mt_dbl,noa*nva,nva,'T_a',(noa*nva),1,g_t_a))
     $     call errquit('could not alloc ga for T_a',1)
      if(.not.ga_create(mt_dbl,nob*nvb,nva,'T_b',(nob*nvb),1,g_t_b))
     $     call errquit('could not alloc ga for T_b',1) 
      tunitptra=1
      tunitptrb=1
      do i=oseg_lo,oseg_hi 
         symi=irs_a(i)
         call ga_zero(g_t_a)
         call ga_zero(g_t_b)
         do syma=0,nir-1
            counter=0
            do symb=0,nir-1
               do b=sym_lo_va(symb),sym_hi_va(symb)
                  symia=ieor(syma,symi)
                  symj=ieor(symia,symb)
                  toff_a(b,syma)=counter
                  counter=counter+num_oa(symj)
               enddo
            enddo
         enddo
         do syma=0,nir-1
            counter=0
            do symb=0,nir-1
               do b=sym_lo_vb(symb),sym_hi_vb(symb)
                  symia=ieor(syma,symi)
                  symj=ieor(symia,symb)
                  toff_b(b,syma)=counter
                  counter=counter+num_ob(symj)
               enddo
            enddo
         enddo
         i_batch=i-oseg_lo+1
         do a_a=nva_lo_local,nva_hi_local 
            syma=irs_a(a_a)
            symi=irs_a(i)
            symia=ieor(syma,symi)
            call dfill((nbf*nbf),0.0d0,tmp,1)
            call moints_rdhfint(kunit,nva_lo_local,nva_hi_local,
     $           i_batch,a_a,nnbf,ioff,ia_uv)
            call scatter(nnbf,tmp,indx_v,ia_uv)
            call dfill((nbf*nbf),0.0d0,ia_uv,1)
            call upper2square(nbf,tmp,ia_uv)
            do symb=0,nir-1
               symj=ieor(symia,symb)
               if(num_oa(symj).gt.0.and.num_va(symb).gt.0)then
                  call dgemm('n','n',nbf,num_oa(symj),nbf,1.0d0,
     $                 ia_uv,nbf,c_a(1,sym_lo_oa(symj)),nbf,
     $                 0.0d0,tmp,nbf)
                  call dgemm('t','n',num_oa(symj),num_va(symb),
     $                 nbf,1.0d0,tmp,nbf,c_a(1,sym_lo_va(symb)),
     $                 nbf,0.0d0,ia_jb,num_oa(symj))
                  jb=1
                  do b=sym_lo_va(symb),sym_hi_va(symb)
                     call ga_acc(g_t_a,toff_a(b,syma)+1,
     $                    toff_a(b,syma)+num_oa(symj),
     $                    a_a-nva_lo+1,a_a-nva_lo+1,
     $                    ia_jb(jb),1,1.0d0)
                     call ga_acc(g_t_a,toff_a(a_a,symb)+1,
     $                    toff_a(a_a,symb)+num_oa(symj),
     $                    b-nva_lo+1,b-nva_lo+1,
     $                    ia_jb(jb),1,minusone)
                     jb=jb+num_oa(symj)
                  enddo
               endif
               if(num_ob(symj).gt.0.and.num_vb(symb).gt.0)then
                  call dgemm('n','n',nbf,num_ob(symj),nbf,1.0d0,
     $                 ia_uv,nbf,c_b(1,sym_lo_ob(symj)),nbf,
     $                 0.0d0,tmp,nbf)
                  call dgemm('t','n',num_ob(symj),num_vb(symb),
     $                 nbf,1.0d0,tmp,nbf,c_b(1,sym_lo_vb(symb)),
     $                 nbf,0.0d0,ia_jb,num_ob(symj))
                  jb=1
                  do b=sym_lo_vb(symb),sym_hi_vb(symb)
                     call ga_put(g_t_b,toff_b(b,syma)+1,
     $                    toff_b(b,syma)+num_ob(symj),
     $                    a_a-nva_lo+1,a_a-nva_lo+1,
     $                    ia_jb(jb),1)
                     jb=jb+num_ob(symj)
                  enddo
               endif
            enddo
         enddo
         call ga_sync()
         do a_a=nva_lo_local,nva_hi_local
            do b=nva_lo,nva_hi
               symb=irs_a(b)
               syma=irs_a(a_a)
               symj=ieor(symb,ieor(symi,syma))
               if(num_oa(symj).ne.0)then
                  call ga_get(g_t_a,toff_a(b,syma)+1,
     $                 toff_a(b,syma)+num_oa(symj),
     $                 a_a-nva_lo+1,a_a-nva_lo+1,
     $                 ia_jb(sym_lo_oa(symj)),1)
                  do j=sym_lo_oa(symj),sym_hi_oa(symj)
                     denom=eval_a(a_a)+eval_a(b)
     $                    -eval_a(i)-eval_a(j)
                     t=ia_jb(j)/denom
                     eaa=eaa+t*t*denom
                     ia_jb(j)=t
                  enddo
                  if(task.eq.'gradient')then
                     counta=sym_hi_oa(symj)-sym_lo_oa(symj)+1
                     call putwa(tunita,ia_jb(sym_lo_oa(symj)),
     $                    tunitptra,counta,ierr)
                     if(ierr.ne.0)call errquit(
     $                    'problem with writing t_a in make_t',ierr)
                     tunitptra=tunitptra+counta
                  endif
               endif
            enddo
            do b=nvb_lo,nvb_hi
               symb=irs_b(b)
               syma=irs_a(a_a)
               symj=ieor(symb,ieor(symi,syma))
               if(num_ob(symj).ne.0)then
                  call ga_get(g_t_b,toff_b(b,syma)+1,
     $                 toff_b(b,syma)+num_ob(symj),
     $                 a_a-nva_lo+1,a_a-nva_lo+1,
     $                 ia_jb(sym_lo_ob(symj)),1)
                  do j=sym_lo_ob(symj),sym_hi_ob(symj)
                     denom=eval_a(a_a)+eval_b(b)
     $                    -eval_a(i)-eval_b(j)
                     t=ia_jb(j)/denom
                     if(abs(t).lt.1d-12)t=0.0d0
                     eab=eab+t*t*denom
                     ia_jb(j)=t
                  enddo
                  if(task.eq.'gradient')then
                     countb=sym_hi_ob(symj)-sym_lo_ob(symj)+1
                     call putwa(tunitb,ia_jb(sym_lo_ob(symj)),
     $                    tunitptrb,countb,ierr)
                     if(ierr.ne.0)call errquit(
     $                    'problem with writing t_b in make_t',ierr)
                     tunitptrb=tunitptrb+countb
                  endif
               endif
            enddo
         enddo
      enddo
      call ga_dgop(1,eaa,1,'+')
      call ga_dgop(1,eab,1,'+')
      e2=e2-0.25d0*eaa-0.5d0*eab
      if(.not.ga_destroy(g_t_b))call errquit(
     $     'could not destroy t_b handle',1)
      if(.not.ga_destroy(g_t_a))call errquit(
     $     'could not destroy t_a handle',1)
      call ga_sync()
      end
      
      
      
      
      
      



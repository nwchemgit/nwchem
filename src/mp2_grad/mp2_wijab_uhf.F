      subroutine mp2_wijab_uhf(nir,nva_lo,nva_hi,
     $     sym_lo_oa,sym_hi_oa,sym_lo_va,sym_hi_va,sym_lo_vb,
     $     sym_hi_vb,nvb_lo,nvb_hi,sym_lo_ob,sym_hi_ob,num_oa,
     $     num_ob,num_va,num_vb,oseg_lo,oseg_hi,irs_a,irs_b,
     $     nva_lo_local,nva_hi_local,noa_lo,noa_hi,
     $     nob_lo,nob_hi,t_a,t_b,t_eps_a,t_eps_b,tunita,tunitb,
     $     w_ij_a,w_ab_a,w_ij_b,w_ab_b,p_ij_a,p_ab_a,p_ij_b,
     $     p_ab_b,eval_a,eval_b,trace)
      implicit none
      integer i,a_a
      integer noa_lo,noa_hi,nva_lo,nva_hi
      integer nob_lo,nob_hi,nvb_lo,nvb_hi
      integer sym_lo_oa(0:7),sym_hi_oa(0:7)
      integer sym_lo_va(0:7),sym_hi_va(0:7)
      integer sym_lo_ob(0:7),sym_hi_ob(0:7)
      integer sym_lo_vb(0:7),sym_hi_vb(0:7)
      integer noa,nva,nob,nvb
      integer num_oa(0:7),num_va(0:7)
      integer num_ob(0:7),num_vb(0:7)
      integer tunita,tunitb
      integer oseg_lo,oseg_hi
      integer nir
      integer g_iajb_a,g_iajb_b
      integer tcounta,tcountb
      double precision eval_a(*),eval_b(*),trace
      double precision eps_i_a,eps_a_a,eps_i_b,eps_a_b
      integer irs_a(*),irs_b(*)
      integer nva_lo_local,nva_hi_local
      double precision w_ij_a(*),w_ab_a(*),w_ij_b(*),w_ab_b(*)
      double precision p_ij_a(*),p_ab_a(*),p_ij_b(*),p_ab_b(*)
      double precision t_a(*),t_b(*),t_eps_a(*),t_eps_b(*)
      double precision zero,one,minusone,half,mhalf
      integer symi,syma,symia,symb,symj
      integer b,tunitptra,tunitptrb,counta,countb,ierr
      integer j,length_a,length_b
      double precision dsum
      double precision epsilon_a,epsilon_b,te
      integer msg_w_ij_a_sum,msg_w_ij_b_sum
      integer msg_w_ab_a_sum,msg_w_ab_b_sum
      parameter(msg_w_ij_a_sum=9996)
      parameter(msg_w_ij_b_sum=9995)
      parameter(msg_w_ab_a_sum=9994)
      parameter(msg_w_ab_b_sum=9993)
c     
      zero=0.0d0
      one=1.0d0
      minusone=-1.0d0
      half=0.5d0
      mhalf=-0.5d0
      noa=noa_hi-noa_lo+1
      nva=nva_hi-nva_lo+1
      nob=nob_hi-nob_lo+1
      nvb=nvb_hi-nvb_lo+1
      tunitptra=1
      tunitptrb=1
      call dfill((noa*noa),zero,w_ij_a,1)
      call dfill((nva*nva),zero,w_ab_a,1)
      call dfill((nob*nob),zero,w_ij_b,1)
      call dfill((nvb*nvb),zero,w_ab_b,1)
      do i=oseg_lo,oseg_hi 
         symi=irs_a(i)
         do a_a=nva_lo_local,nva_hi_local
            syma=irs_a(a_a)
            symia=ieor(syma,symi)
            tcounta=1
            tcountb=1
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
                  epsilon_a=eval_a(a_a)+eval_a(b)
     $                 -eval_a(i)-eval_a(j)
                  t_eps_a(tcounta+j-sym_lo_oa(symj))
     $                 =t_a(tcounta+j-sym_lo_oa(symj))*epsilon_a
               enddo
               tcounta=tcounta+counta
               endif
            enddo
            do b=nvb_lo,nvb_hi
               symb=irs_b(b)
               symj=ieor(symb,symia)
               if(num_ob(symj).ne.0)then
                  countb=sym_hi_ob(symj)-sym_lo_ob(symj)+1
                  call getwa(tunitb,t_b(tcountb),
     $                 tunitptrb,countb,ierr)
                  if(ierr.ne.0)call errquit(
     $                 'problem reading t_b back in',ierr)
                  tunitptrb=tunitptrb+countb
               do j=sym_lo_ob(symj),sym_hi_ob(symj)
                  epsilon_b=eval_a(a_a)+eval_b(b)
     $                 -eval_a(i)-eval_b(j)
                  t_eps_b(tcountb+j-sym_lo_ob(symj))
     $                 =t_b(tcountb+j-sym_lo_ob(symj))*epsilon_b
               enddo
               tcountb=tcountb+countb
               endif
            enddo
            tcounta=1
            tcountb=1
            do symb=0,nir-1
               symj=ieor(symia,symb)
               if(num_oa(symj).gt.0.and.num_va(symb).gt.0)then
                  call dgemm('n','t',num_oa(symj),
     $                 num_oa(symj),num_va(symb),one,
     $                 t_a(tcounta),num_oa(symj),t_eps_a(tcounta),
     $                 num_oa(symj),one,w_ij_a,num_oa(symj))
                  call dgemm('t','n',num_va(symb),
     $                 num_va(symb),num_oa(symj),one,
     $                 t_a(tcounta),num_oa(symj),t_eps_a(tcounta),
     $                 num_oa(symj),one,w_ab_a,num_va(symb))
                  tcounta=tcounta+num_oa(symj)*num_va(symb)
               endif
               if(num_ob(symj).gt.0.and.num_vb(symb).gt.0)then
                  call dgemm('n','t',num_ob(symj),
     $                 num_ob(symj),num_vb(symb),one,
     $                 t_b(tcountb),num_ob(symj),t_eps_b(tcountb),
     $                 num_ob(symj),one,w_ij_b,num_ob(symj))
                  call dgemm('t','n',num_vb(symb),
     $                 num_vb(symb),num_ob(symj),one,
     $                 t_b(tcountb),num_ob(symj),t_eps_b(tcountb),
     $                 num_ob(symj),one,w_ab_b,num_vb(symb))
                  tcountb=tcountb+num_ob(symj)*num_vb(symb)
               endif
            enddo
         enddo
      enddo
      call ga_dgop(msg_w_ij_a_sum,w_ij_a,noa*noa,'+')
      call ga_dgop(msg_w_ab_a_sum,w_ab_a,nva*nva,'+')
      call ga_dgop(msg_w_ij_b_sum,w_ij_b,nob*nob,'+')
      call ga_dgop(msg_w_ab_b_sum,w_ab_b,nvb*nvb,'+')
      trace=trace-0.25d0*dsum(noa,w_ij_a,noa+1)
     $     -0.5d0*dsum(nob,w_ij_b,nob+1)
      call ga_sync()
      end

      double precision function dsum(n,x,ix)
      double precision x(ix,n)
      dsum=0.0d0
      do i=1,n
         dsum=dsum+x(1,i)
      enddo
      end
      

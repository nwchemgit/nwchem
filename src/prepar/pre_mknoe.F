      logical function pre_mknoe(lfnout,lfntop,filtop,lfncmd,filcmd,
     + lfnnoe,filnoe,mnoe,lfnpmf,filpmf,slvnam)
c
c $Id: pre_mknoe.F,v 1.21 2001-08-07 18:10:13 d3j191 Exp $
c
      implicit none
c
#include "mafdecls.fh"
#include "util.fh"
c
      logical pre_topsiz,pre_toprd,pre_disres,pre_select,pre_equiv
      external pre_topsiz,pre_toprd,pre_disres,pre_select,pre_equiv
c
      integer lfnout,lfntop,lfncmd,lfnnoe,lfnpmf
      character*255 filtop,filcmd,filnoe,filpmf
      character*80 card,card2
c
      integer nsa,nwa,noe
      integer mat,msa,mwa,mnoe,msm,nsm,msb,nsb
      integer nnoe,nqu
      character*3 slvnam
c
      integer l_cwa,i_cwa,l_mas,i_mas,l_num,i_num,l_sat,i_sat
      integer l_csa,i_csa,l_sgm,i_sgm,l_sml,i_sml,l_qsa,i_qsa
      integer l_inoe,i_inoe,l_dnoe,i_dnoe,l_qu,i_qu
      integer l_qwa,i_qwa,l_wat,i_wat,l_sfr,i_sfr
      integer l_isb,i_isb,l_csb,i_csb,l_eq,i_eq
      integer l_temp,i_temp
c
      integer maxgrp,maxatm,maxpmf
      parameter(maxgrp=500)
      parameter(maxatm=500)
      parameter(maxpmf=500)
c
      integer igroup(maxgrp,maxatm),lgroup(maxgrp),mgroup(maxgrp)
      integer numpmf,ipmf(maxpmf,5),rpmf(maxpmf,4),iopt(maxpmf)
      integer npmfr,ipmfr,nsarep
c
      integer i
c
      npmfr=0
      ipmfr=0
c
      do 2 i=1,maxgrp
      lgroup(i)=0
      mgroup(i)=0
    2 continue
      numpmf=0
c
      if(util_print('where',print_debug)) then
      write(lfnout,1000)
 1000 format(' NOE FILE GENERATION')
      endif
c
      if(.not.pre_topsiz(lfntop,filtop,lfnout,mat,msa,mwa,msb,nqu))
     + call errquit('pre_topsiz failed',9999)
c
      if(util_print('where',print_high)) then
      write(lfnout,1001) filtop(1:index(filtop,' ')-1),mat,msa,mwa,msb
 1001 format(' Topology',t40,a,//,
     + ' Number of atom types',t40,i8,/,
     + ' Number of solute atoms',t40,i8,/,
     + ' Number of solvent atoms',t40,i8,/,
     + ' Number of solute bonds',t40,i8,/)
      endif
      nwa=0
      nsa=0
      nsb=0
      nnoe=0
c
c     allocate memory
c     ---------------
c
c     character*16 cwa(mwa)  : solvent atom names
c     character*16 csa(msa)  : solute atom names
c     integer isar(msa)      : solute atom types
c     integer isgm(msa)      : solute segment numbers
c     integer isfnd(msa)     : solute atom found flags
c     real*8 qsa(msa)        : solute atom charges
c     real*8 xs(3,msa)       : solute atom coordinates
c
      if(.not.ma_push_get(mt_dbl,nqu,'qu',l_qu,i_qu))
     + call errquit('Memory allocation failed for qu',9999)
      if(.not.ma_push_get(mt_int,2*mnoe,'inoe',l_inoe,i_inoe))
     + call errquit('Memory allocation failed for inoe',9999)
      if(.not.ma_push_get(mt_dbl,5*mnoe,'dnoe',l_dnoe,i_dnoe))
     + call errquit('Memory allocation failed for dnoe',9999)
      if(.not.ma_push_get(mt_int,mat,'anum',l_num,i_num))
     + call errquit('Memory allocation failed for anum',9999)
      if(.not.ma_push_get(mt_dbl,mat,'amass',l_mas,i_mas))
     + call errquit('Memory allocation failed for amass',9999)
      if(.not.ma_push_get(mt_byte,16*mwa,'cwa',l_cwa,i_cwa))
     + call errquit('Memory allocation failed for cwa',9999)
      if(.not.ma_push_get(mt_dbl,mwa,'qwa',l_qwa,i_qwa))
     + call errquit('Memory allocation failed for qwa',9999)
      if(.not.ma_push_get(mt_byte,16*msa,'csa',l_csa,i_csa))
     + call errquit('Memory allocation failed for csa',9999)
      if(.not.ma_push_get(mt_int,mwa,'wat',l_wat,i_wat))
     + call errquit('Memory allocation failed for wat',9999)
      if(.not.ma_push_get(mt_int,msa,'sat',l_sat,i_sat))
     + call errquit('Memory allocation failed for sat',9999)
      if(.not.ma_push_get(mt_int,msa,'sgm',l_sgm,i_sgm))
     + call errquit('Memory allocation failed for sgm',9999)
      if(.not.ma_push_get(mt_int,msa,'sfr',l_sfr,i_sfr))
     + call errquit('Memory allocation failed for sfr',9999)
      if(.not.ma_push_get(mt_int,msa,'sml',l_sml,i_sml))
     + call errquit('Memory allocation failed for sml',9999)
      if(.not.ma_push_get(mt_dbl,msa,'qsa',l_qsa,i_qsa))
     + call errquit('Memory allocation failed for qsa',9999)
      if(.not.ma_push_get(mt_int,2*msb,'idsb',l_isb,i_isb))
     + call errquit('Memory allocation failed for idsb',9999)
      if(.not.ma_push_get(mt_dbl,msb,'cdsb',l_csb,i_csb))
     + call errquit('Memory allocation failed for cdsb',9999)
c
      nwa=0
      nsa=0
      nsm=0
      noe=0
c
      if(util_print('where',print_debug)) then
      write(lfnout,1002) 
 1002 format(' Memory allocated')
      endif
c
c     read topology file
c     ------------------
c
      if(.not.pre_toprd(lfntop,filtop,lfnout,
     + int_mb(i_num),dbl_mb(i_mas),mat,
     + byte_mb(i_cwa),dbl_mb(i_qwa),mwa,nwa,
     + int_mb(i_wat),int_mb(i_sat),int_mb(i_sgm),int_mb(i_sml),
     + int_mb(i_sfr),
     + byte_mb(i_csa),dbl_mb(i_qsa),msa,nsa,nsm,int_mb(i_isb),
     + dbl_mb(i_csb),msb,nsb,dbl_mb(i_qu),nqu,slvnam))
     + call errquit('pre_toprd failed',9999)
      if(util_print('coordinates',print_default)) then
      write(lfnout,1003) filtop(1:index(filtop,' ')-1)
 1003 format(' Topology',t40,a,/)
      endif
c
      msm=nsm
c
      if(.not.ma_push_get(mt_int,4*msm,'eq',l_eq,i_eq))
     + call errquit('Memory allocation failed for eq',9999)
c
c     determine equivalent solute molecules
c     -------------------------------------
c
      if(.not.pre_equiv(byte_mb(i_csa),int_mb(i_sml),int_mb(i_sgm),
     + msa,nsa,int_mb(i_eq),msm,nsm))
     + call errquit('Solute molecule equivalent test failed',0)
c
c     open the command file
c     ---------------------
c
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old',err=999)
      rewind(lfncmd)
c
    1 continue
c
      read(lfncmd,1101,end=9,err=9997) card
 1101 format(a)
c
      if(util_print('where',print_high)) then
      write(lfnout,1099) card
 1099 format('command card ',a)
      endif
c
c     read distance restraints
c     ========================
c
      if(card(1:6).eq.'disres') then
      read(lfncmd,1101,end=9,err=9997) card2
      if(.not.pre_disres(card,card2,int_mb(i_sgm),byte_mb(i_csa),msa,
     + nsa,int_mb(i_inoe),dbl_mb(i_dnoe),mnoe,nnoe))
     + call errquit('pre_disres failed',9999)
      endif
c
c     read pmf replication
c     ====================
c
      if(card(1:7).eq.'rep_pmf') then
      read(card(8:22),'(i5,i10)') npmfr,nsarep
      ipmfr=numpmf
      endif
c
c     read group selection
c     ====================
c
      if(card(1:6).eq.'select') then
      if(.not.pre_select(card,int_mb(i_sgm),int_mb(i_sml),
     + byte_mb(i_csa),msa,
     + nsa,maxgrp,maxatm,igroup,lgroup,mgroup))
     + call errquit('pre_select failed',9999)
      endif
c
c     read pmf
c     ========
c
      if(card(1:3).eq.'pmf') then
      numpmf=numpmf+1
      if(numpmf.gt.maxpmf) call errquit('Increase maxpmf',maxpmf)
      ipmf(numpmf,1)=0
      if(card(5:8).eq.'dist') ipmf(numpmf,1)=1
      if(card(5:8).eq.'angl') ipmf(numpmf,1)=2
      if(card(5:8).eq.'tors') ipmf(numpmf,1)=3
      if(card(5:8).eq.'impr') ipmf(numpmf,1)=4
      if(card(5:8).eq.'line') ipmf(numpmf,1)=5
      if(card(5:8).eq.'plan') ipmf(numpmf,1)=6
      if(card(5:8).eq.'base') ipmf(numpmf,1)=7
      read(card(9:76),'(5i5,2f12.6,2e12.5)')
     + iopt(numpmf),(ipmf(numpmf,i),i=2,5),(rpmf(numpmf,i),i=1,4)
      endif
c
      if(card(1:3).ne.'end') goto 1
c
    9 continue
      close(unit=lfncmd,err=999)
c
      if(nnoe.gt.0) then
      call pre_wrtnoe(lfntop,filtop,mnoe,nnoe,
     + int_mb(i_inoe),dbl_mb(i_dnoe))
      write(lfnout,1004)
 1004 format(/,' Distance restraints appended to topology')
      endif
c
      if(numpmf.gt.0) then
      if(.not.ma_push_get(mt_int,msa,'temp',l_temp,i_temp))
     + call errquit('Memory allocation failed for temp',9999)
      call pre_wrtpmf(lfntop,filtop,numpmf,
     + maxpmf,ipmf,rpmf,iopt,maxgrp,maxatm,igroup,lgroup,mgroup,npmfr,
     + ipmfr,nsarep,msa,nsa,int_mb(i_sgm),int_mb(i_sml),byte_mb(i_csa),
     + int_mb(i_eq),msm,nsm,int_mb(i_temp))
      if(.not.ma_pop_stack(l_temp))
     + call errquit('Memory deallocation failed for temp',9999)
      write(lfnout,1005)
 1005 format(/,' Potential of mean force appended to topology')
      endif
c
  999 continue
c
c     deallocate memory
c     -----------------
c
      if(.not.ma_pop_stack(l_eq))
     + call errquit('Memory deallocation failed for eq',9999)
      if(.not.ma_pop_stack(l_csb))
     + call errquit('Memory deallocation failed for cdsb',9999)
      if(.not.ma_pop_stack(l_isb))
     + call errquit('Memory deallocation failed for idsb',9999)
      if(.not.ma_pop_stack(l_qsa))
     + call errquit('Memory deallocation failed for qsa',9999)
      if(.not.ma_pop_stack(l_sml))
     + call errquit('Memory deallocation failed for qsa',9999)
      if(.not.ma_pop_stack(l_sfr))
     + call errquit('Memory deallocation failed for sfr',9999)
      if(.not.ma_pop_stack(l_sgm))
     + call errquit('Memory deallocation failed for sgm',9999)
      if(.not.ma_pop_stack(l_sat))
     + call errquit('Memory deallocation failed for sat',9999)
      if(.not.ma_pop_stack(l_wat))
     + call errquit('Memory deallocation failed for wat',9999)
      if(.not.ma_pop_stack(l_csa))
     + call errquit('Memory deallocation failed for csa',9999)
      if(.not.ma_pop_stack(l_qwa))
     + call errquit('Memory deallocation failed for qwa',9999)
      if(.not.ma_pop_stack(l_cwa))
     + call errquit('Memory deallocation failed for cwa',9999)
      if(.not.ma_pop_stack(l_mas))
     + call errquit('Memory deallocation failed for amass',9999)
      if(.not.ma_pop_stack(l_num))
     + call errquit('Memory deallocation failed for anum',9999)
      if(.not.ma_pop_stack(l_dnoe))
     + call errquit('Memory deallocation failed for dnoe',9999)
      if(.not.ma_pop_stack(l_inoe))
     + call errquit('Memory deallocation failed for inoe',9999)
      if(.not.ma_pop_stack(l_qu))
     + call errquit('Memory deallocation failed for qu',9999)
c
      pre_mknoe=.true.
      return
c
 9997 continue
      write(lfnout,1017) filcmd(1:index(filcmd,' ')-1)
 1017 format(' Error reading commands',t40,a,/)
      pre_mknoe=.false.
      return
c
 9999 continue
      pre_mknoe=.false.
      return
      end

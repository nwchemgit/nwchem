      logical function argos_prep_mknoe(lfnout,lfntop,
     + filtop,lfncmd,filcmd,
     + lfnnoe,filnoe,mnoe,lfnpmf,filpmf,slvnam)
c
c $Id$
c
      implicit none
c
#include "mafdecls.fh"
#include "util.fh"
c
      logical argos_prep_topsiz,argos_prep_toprd,
     + argos_prep_disres,argos_prep_select,argos_prep_equiv
      external argos_prep_topsiz,argos_prep_toprd,
     + argos_prep_disres,argos_prep_select,argos_prep_equiv
c
      integer lfnout,lfntop,lfncmd,lfnnoe,lfnpmf
      character*255 filtop,filcmd,filnoe,filpmf
      character*255 card,card2
c
      integer nsa,nwa,noe
      integer mat,msa,mwa,mnoe,msm,nsm,msb,nsb
      integer nnoe,nqu,mseq,nseq
      character*3 slvnam
c
      integer l_cwa,i_cwa,l_mas,i_mas,l_num,i_num,l_sat,i_sat
      integer l_csa,i_csa,l_sgm,i_sgm,l_sml,i_sml,l_qsa,i_qsa
      integer l_inoe,i_inoe,l_dnoe,i_dnoe,l_qu,i_qu
      integer l_qwa,i_qwa,l_wat,i_wat,l_sfr,i_sfr
      integer l_isb,i_isb,l_csb,i_csb,l_eq,i_eq
      integer l_temp,i_temp,l_tempj,i_tempj
      integer l_lseq,i_lseq,i_ihop,l_ihop,l_istat,i_istat
c
      integer maxgrp,maxatm,maxpmf
      parameter(maxgrp=500)
      parameter(maxatm=2500)
      parameter(maxpmf=500)
c
      integer igroup(maxgrp,maxatm),lgroup(maxgrp),mgroup(maxgrp)
      integer numpmf,ipmf(maxpmf,5),iopt(maxpmf,3)
      real*8 rpmf(maxpmf,4)
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
      if(.not.argos_prep_topsiz(lfntop,filtop,lfnout,
     + mat,msa,mwa,msb,nqu,nseq))
     + call md_abort('argos_prep_topsiz failed',9999)
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
      mseq=nseq
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
     + call md_abort('Memory allocation failed for qu',9999)
      if(.not.ma_push_get(mt_int,2*mnoe,'inoe',l_inoe,i_inoe))
     + call md_abort('Memory allocation failed for inoe',9999)
      if(.not.ma_push_get(mt_dbl,5*mnoe,'dnoe',l_dnoe,i_dnoe))
     + call md_abort('Memory allocation failed for dnoe',9999)
      if(.not.ma_push_get(mt_int,mat,'anum',l_num,i_num))
     + call md_abort('Memory allocation failed for anum',9999)
      if(.not.ma_push_get(mt_dbl,mat,'amass',l_mas,i_mas))
     + call md_abort('Memory allocation failed for amass',9999)
      if(.not.ma_push_get(mt_byte,16*mwa,'cwa',l_cwa,i_cwa))
     + call md_abort('Memory allocation failed for cwa',9999)
      if(.not.ma_push_get(mt_dbl,mwa,'qwa',l_qwa,i_qwa))
     + call md_abort('Memory allocation failed for qwa',9999)
      if(.not.ma_push_get(mt_byte,16*msa,'csa',l_csa,i_csa))
     + call md_abort('Memory allocation failed for csa',9999)
      if(.not.ma_push_get(mt_int,mwa,'wat',l_wat,i_wat))
     + call md_abort('Memory allocation failed for wat',9999)
      if(.not.ma_push_get(mt_int,msa,'sat',l_sat,i_sat))
     + call md_abort('Memory allocation failed for sat',9999)
      if(.not.ma_push_get(mt_int,msa,'sgm',l_sgm,i_sgm))
     + call md_abort('Memory allocation failed for sgm',9999)
      if(.not.ma_push_get(mt_int,msa,'sfr',l_sfr,i_sfr))
     + call md_abort('Memory allocation failed for sfr',9999)
      if(.not.ma_push_get(mt_int,msa,'sml',l_sml,i_sml))
     + call md_abort('Memory allocation failed for sml',9999)
      if(.not.ma_push_get(mt_dbl,msa,'qsa',l_qsa,i_qsa))
     + call md_abort('Memory allocation failed for qsa',9999)
      if(.not.ma_push_get(mt_int,2*msb,'idsb',l_isb,i_isb))
     + call md_abort('Memory allocation failed for idsb',9999)
      if(.not.ma_push_get(mt_dbl,msb,'cdsb',l_csb,i_csb))
     + call md_abort('Memory allocation failed for cdsb',9999)
      if(.not.ma_push_get(mt_int,mseq,'lseq',l_lseq,i_lseq))
     + call md_abort('Memory allocation failed for lseq',9999)
      if(.not.ma_push_get(mt_int,msa,'ihop',l_ihop,i_ihop))
     + call md_abort('Memory allocation failed for ihop',9999)
      if(.not.ma_push_get(mt_int,msa,'istat',l_istat,i_istat))
     + call md_abort('Memory allocation failed for istat',9999)
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
      if(.not.argos_prep_toprd(lfntop,filtop,lfnout,
     + int_mb(i_num),dbl_mb(i_mas),mat,
     + byte_mb(i_cwa),dbl_mb(i_qwa),mwa,nwa,
     + int_mb(i_wat),int_mb(i_sat),int_mb(i_sgm),int_mb(i_sml),
     + int_mb(i_sfr),
     + byte_mb(i_csa),dbl_mb(i_qsa),msa,nsa,nsm,int_mb(i_isb),
     + dbl_mb(i_csb),msb,nsb,dbl_mb(i_qu),nqu,slvnam,
     + mseq,nseq,int_mb(i_lseq),int_mb(i_ihop),int_mb(i_istat)))
     + call md_abort('argos_prep_toprd failed',9999)
      if(util_print('coordinates',print_default)) then
      write(lfnout,1003) filtop(1:index(filtop,' ')-1)
 1003 format(' Topology',t40,a,/)
      endif
c
      msm=nsm
c
      if(.not.ma_push_get(mt_int,4*msm,'eq',l_eq,i_eq))
     + call md_abort('Memory allocation failed for eq',9999)
c
c     determine equivalent solute molecules
c     -------------------------------------
c
      if(.not.argos_prep_equiv(byte_mb(i_csa),int_mb(i_sml),
     + int_mb(i_sgm),
     + msa,nsa,int_mb(i_eq),msm,nsm))
     + call md_abort('Solute molecule equivalent test failed',0)
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
      if(.not.argos_prep_disres(card,card2,int_mb(i_sgm),
     + byte_mb(i_csa),msa,
     + nsa,int_mb(i_inoe),dbl_mb(i_dnoe),mnoe,nnoe))
     + call md_abort('argos_prep_disres failed',9999)
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
      if(.not.argos_prep_select(card,int_mb(i_sgm),int_mb(i_sml),
     + byte_mb(i_csa),msa,
     + nsa,maxgrp,maxatm,igroup,lgroup,mgroup))
     + call md_abort('argos_prep_select failed',9999)
      endif
c
c     read pmf
c     ========
c
      if(card(1:3).eq.'pmf') then
      numpmf=numpmf+1
      if(numpmf.gt.maxpmf) call md_abort('Increase maxpmf',maxpmf)
      ipmf(numpmf,1)=0
      if(card(5:7).eq.'con') ipmf(numpmf,1)=-1
      if(card(5:7).eq.'dis') ipmf(numpmf,1)=1
      if(card(5:7).eq.'ang') ipmf(numpmf,1)=2
      if(card(5:7).eq.'tor') ipmf(numpmf,1)=3
      if(card(5:7).eq.'imp') ipmf(numpmf,1)=4
      if(card(5:7).eq.'lin') ipmf(numpmf,1)=5
      if(card(5:7).eq.'pla') ipmf(numpmf,1)=6
      if(card(5:7).eq.'bas') ipmf(numpmf,1)=7
      if(card(5:7).eq.'zax') ipmf(numpmf,1)=8
      if(card(5:7).eq.'zdi') ipmf(numpmf,1)=9
      if(card(5:7).eq.'ZAX') ipmf(numpmf,1)=10
      read(card(8:90),'(7i5,2f12.6,2e12.5)')
     + (iopt(numpmf,i),i=1,3),(ipmf(numpmf,i),i=2,5),
     + (rpmf(numpmf,i),i=1,4)
      endif
c
      if(card(1:3).ne.'end') goto 1
c
    9 continue
      close(unit=lfncmd,err=999)
c
      if(nnoe.gt.0) then
      call argos_prep_wrtnoe(lfntop,filtop,mnoe,nnoe,
     + int_mb(i_inoe),dbl_mb(i_dnoe))
      write(lfnout,1004)
 1004 format(/,' Distance restraints appended to topology')
      endif
c
      if(numpmf.gt.0) then
      if(.not.ma_push_get(mt_int,msa,'temp',l_temp,i_temp))
     + call md_abort('Memory allocation failed for temp',9999)
      if(.not.ma_push_get(mt_int,msa,'tempj',l_tempj,i_tempj))
     + call md_abort('Memory allocation failed for tempj',9999)
      call argos_prep_wrtpmf(lfntop,filtop,numpmf,
     + maxpmf,ipmf,rpmf,iopt,maxgrp,maxatm,igroup,lgroup,mgroup,npmfr,
     + ipmfr,nsarep,msa,nsa,int_mb(i_sgm),int_mb(i_sml),byte_mb(i_csa),
     + int_mb(i_eq),msm,nsm,int_mb(i_temp),int_mb(i_tempj))
      if(.not.ma_pop_stack(l_tempj))
     + call md_abort('Memory deallocation failed for tempj',9999)
      if(.not.ma_pop_stack(l_temp))
     + call md_abort('Memory deallocation failed for temp',9999)
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
     + call md_abort('Memory deallocation failed for eq',9999)
      if(.not.ma_pop_stack(l_istat))
     + call md_abort('Memory deallocation failed for istat',9999)
      if(.not.ma_pop_stack(l_ihop))
     + call md_abort('Memory deallocation failed for ihop',9999)
      if(.not.ma_pop_stack(l_lseq))
     + call md_abort('Memory deallocation failed for lseq',9999)
      if(.not.ma_pop_stack(l_csb))
     + call md_abort('Memory deallocation failed for cdsb',9999)
      if(.not.ma_pop_stack(l_isb))
     + call md_abort('Memory deallocation failed for idsb',9999)
      if(.not.ma_pop_stack(l_qsa))
     + call md_abort('Memory deallocation failed for qsa',9999)
      if(.not.ma_pop_stack(l_sml))
     + call md_abort('Memory deallocation failed for qsa',9999)
      if(.not.ma_pop_stack(l_sfr))
     + call md_abort('Memory deallocation failed for sfr',9999)
      if(.not.ma_pop_stack(l_sgm))
     + call md_abort('Memory deallocation failed for sgm',9999)
      if(.not.ma_pop_stack(l_sat))
     + call md_abort('Memory deallocation failed for sat',9999)
      if(.not.ma_pop_stack(l_wat))
     + call md_abort('Memory deallocation failed for wat',9999)
      if(.not.ma_pop_stack(l_csa))
     + call md_abort('Memory deallocation failed for csa',9999)
      if(.not.ma_pop_stack(l_qwa))
     + call md_abort('Memory deallocation failed for qwa',9999)
      if(.not.ma_pop_stack(l_cwa))
     + call md_abort('Memory deallocation failed for cwa',9999)
      if(.not.ma_pop_stack(l_mas))
     + call md_abort('Memory deallocation failed for amass',9999)
      if(.not.ma_pop_stack(l_num))
     + call md_abort('Memory deallocation failed for anum',9999)
      if(.not.ma_pop_stack(l_dnoe))
     + call md_abort('Memory deallocation failed for dnoe',9999)
      if(.not.ma_pop_stack(l_inoe))
     + call md_abort('Memory deallocation failed for inoe',9999)
      if(.not.ma_pop_stack(l_qu))
     + call md_abort('Memory deallocation failed for qu',9999)
c
      argos_prep_mknoe=.true.
      return
c
 9997 continue
      write(lfnout,1017) filcmd(1:index(filcmd,' ')-1)
 1017 format(' Error reading commands',t40,a,/)
      argos_prep_mknoe=.false.
      return
c
 9999 continue
      argos_prep_mknoe=.false.
      return
      end

      logical function argos_prep_mkrst(irtdb,title,mdold,source,lfnout,
     + lfntop,filtop,lfnpdb,filpdb,lfnrst,filrst,lfntmp,filtmp,
     + lfncmd,filcmd,lfnslv,slvmdl,slvnam,mgrid,nrgrid,iogrid,rogrid,
     + gdist,mnoe,sysnam,altloc,
     + chain,model,lfnxyz,filxyz,lfnqqq,filqqq,lfnpov,filpov,scale,cpk,
     + lfnmat,lfnmrg,nmerge,xmerge,filmrg,fcount,lfnmod,filmod,itopol)
c
c $Id$
c
      implicit none
c
#include "mafdecls.fh"
#include "util.fh"
c
      logical argos_prep_topsiz,argos_prep_toprd,argos_prep_rdpdb,
     + argos_prep_wrtrst,argos_prep_center
      external argos_prep_topsiz,argos_prep_toprd,argos_prep_rdpdb,
     + argos_prep_wrtrst,argos_prep_center
      logical argos_prep_rdgeom,argos_prep_orient,argos_prep_slvsiz,
     + argos_prep_slvrd,argos_prep_slvnum
      external argos_prep_rdgeom,argos_prep_orient,argos_prep_slvsiz,
     + argos_prep_slvrd,argos_prep_slvnum
      logical argos_prep_pdbsiz,argos_prep_solvat,argos_prep_boxsiz,
     + argos_prep_misfit,argos_prep_wrtpdb
      external argos_prep_pdbsiz,argos_prep_solvat,argos_prep_boxsiz,
     + argos_prep_misfit,argos_prep_wrtpdb
      logical argos_prep_rdrst,argos_prep_disres,argos_prep_fix,
     + argos_prep_rstsiz,argos_prep_wtrst
      external argos_prep_rdrst,argos_prep_disres,argos_prep_fix,
     + argos_prep_rstsiz,argos_prep_wtrst
      logical argos_prep_rrst,argos_prep_rstsz,
     + argos_prep_wrtxyz
      external argos_prep_rrst,argos_prep_rstsz,
     + argos_prep_wrtxyz
      logical argos_prep_align,argos_prep_repeat,argos_prep_transl,
     + argos_prep_removw,argos_prep_wrtpov
      external argos_prep_align,argos_prep_repeat,argos_prep_transl,
     + argos_prep_removw,argos_prep_wrtpov
      logical argos_prep_rot,argos_prep_crop,argos_prep_his
      external argos_prep_rot,argos_prep_crop,argos_prep_his
c
      integer irtdb,lfnout,lfntop,lfnpdb,lfnrst,lfntmp,lfncmd,lfnmat
      integer lfnslv,model,mdold,lfnxyz,lfnqqq,lfnpov,lfnmrg,lfnmod
      character*255 filtop,filpdb,filrst,filtmp,filcmd,filnam,filslv
      character*255 filxyz,filqqq,filpov,filmrg(100),filmod
      character*10 slvmdl
      character*3 slvnam
      character*80 sysnam,source,card,card2,cfix,title(2,3)
      character*1 altloc,chain
      real*8 fcount
      real*8 scale,cpk
c
      integer nsa,nwa,nwm,nmis,npbtyp,nbxtyp,noe,nwms,nwas,mwmc,nwmc
      integer mat,msa,mwa,mwm,mnoe,msm,nsm,mwms,msb,nsb,mgrid,ngrid,iopt
      integer nnoe,ifix,jfix,nnwm,nnwa,nnsa,nnwmc,nqu,iropt,irrand,irrep
      integer invert,nmerge,mseq,nseq,itopol
      integer nrgrid,iogrid(5),rogrid(2,5)
      real*8 box(3),boxs(3),touch,xmax(3),edge,rshell,gdist,rfix,drep(3)
      real*8 boxset(3),xmerge(3,100)
      real*8 edgex,edgey,edgez
c
      integer l_cwa,i_cwa,l_mas,i_mas,l_num,i_num,l_sat,i_sat
      integer l_csa,i_csa,l_sar,i_sar,l_sgm,i_sgm,l_qsa,i_qsa
      integer l_sml,i_sml,l_sfnd,i_sfnd,l_xs,i_xs,l_vs,i_vs
      integer l_xw,i_xw,l_vw,i_vw,l_iwmr,i_iwmr,l_sfr,i_sfr
      integer l_inoe,i_inoe,l_dnoe,i_dnoe,l_qu,i_qu
      integer l_xslv,i_xslv,l_vslv,i_vslv,l_xwc,i_xwc,l_vwc,i_vwc
      integer l_wfnd,i_wfnd,l_wcmr,i_wcmr,l_qwa,i_qwa
      integer l_isb,i_isb,l_csb,i_csb,l_grid,i_grid
      integer l_wat,i_wat,l_lang,i_lang,l_ndx,i_ndx,l_vec,i_vec
      integer l_lseq,i_lseq,l_ips,i_ips,l_ihop,i_ihop,l_istat,i_istat
c
      logical lslvnt,ldistr,lsolva,lcryst,lnames,lpdbbx
      integer length,len,islv,idum1,idum2,idum3,idum4,iwater,irenum
      integer nxrep,nyrep,nzrep,itran,jtran,i,j,k,irot,jrot,mlang,lrgpdb
      real*8 xtran(3),xrot(3),arot
      real*8 zdist,xr(3,2),boxmax
c
      integer iskip(3,100),nskip,nmoves
      integer mcolgr,ncolgr
      parameter(mcolgr=20)
      integer icolgr(2,mcolgr)
c
      ncolgr=0
      irrep=0
      lpdbbx=.false.
c
      if(util_print('restart',print_debug)) then
      write(lfnout,1000)
 1000 format(' RESTART FILE GENERATION')
      endif
c
      if(.not.argos_prep_topsiz(lfntop,filtop,lfnout,mat,msa,
     + mwa,msb,nqu,nseq))
     + call md_abort('argos_prep_topsiz failed',9999)
c
      if(util_print('restart',print_high)) then
      write(lfnout,1001) filtop(1:index(filtop,' ')-1),mat,msa,mwa,msb
 1001 format(' Topology',t40,a,//,
     + ' Number of atom types',t40,i8,/,
     + ' Number of solute atoms',t40,i8,/,
     + ' Number of solvent atoms',t40,i8,/,
     + ' Number of solute bonds',t40,i8,/)
      endif
      nwm=0
      nwmc=0
      nwms=0
      nwa=0
      nsa=0
      nsb=0
      nnoe=0
      xmax(1)=0.0d0
      xmax(2)=0.0d0
      xmax(3)=0.0d0
      rshell=1.2d0
      nbxtyp=0
      nxrep=0
      nyrep=0
      nzrep=0
      iropt=0
      irrand=0
      invert=0
      nskip=0
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
      if(.not.ma_push_get(mt_int,msa,'sar',l_sar,i_sar))
     + call md_abort('Memory allocation failed for sar',9999)
      if(.not.ma_push_get(mt_int,mwa,'wat',l_wat,i_wat))
     + call md_abort('Memory allocation failed for wat',9999)
      if(.not.ma_push_get(mt_int,msa,'sat',l_sat,i_sat))
     + call md_abort('Memory allocation failed for sat',9999)
      if(.not.ma_push_get(mt_int,msa,'sgm',l_sgm,i_sgm))
     + call md_abort('Memory allocation failed for sgm',9999)
      if(.not.ma_push_get(mt_int,msa,'sml',l_sml,i_sml))
     + call md_abort('Memory allocation failed for sml',9999)
      if(.not.ma_push_get(mt_int,msa,'sfr',l_sfr,i_sfr))
     + call md_abort('Memory allocation failed for sfr',9999)
      if(.not.ma_push_get(mt_int,msa,'sfnd',l_sfnd,i_sfnd))
     + call md_abort('Memory allocation failed for sfnd',9999)
      if(.not.ma_push_get(mt_dbl,msa,'qsa',l_qsa,i_qsa))
     + call md_abort('Memory allocation failed for qsa',9999)
      if(.not.ma_push_get(mt_dbl,3*msa,'xs',l_xs,i_xs))
     + call md_abort('Memory allocation failed for xs',9999)
      if(.not.ma_push_get(mt_dbl,3*msa,'vs',l_vs,i_vs))
     + call md_abort('Memory allocation failed for xs',9999)
      if(.not.ma_push_get(mt_int,2*msb,'idsb',l_isb,i_isb))
     + call md_abort('Memory allocation failed for idsb',9999)
      if(.not.ma_push_get(mt_dbl,msb,'cdsb',l_csb,i_csb))
     + call md_abort('Memory allocation failed for cdsb',9999)
      if(.not.ma_push_get(mt_int,mseq,'lseq',l_lseq,i_lseq))
     + call md_abort('Memory allocation failed for lseq',9999)
      if(.not.ma_push_get(mt_int,msa,'ips',l_ips,i_ips))
     + call md_abort('Memory allocation failed for ips',9999)
      if(.not.ma_push_get(mt_int,msa,'ihop',l_ihop,i_ihop))
     + call md_abort('Memory allocation failed for ihop',9999)
      if(.not.ma_push_get(mt_int,msa,'istat',l_istat,i_istat))
     + call md_abort('Memory allocation failed for istat',9999)
c
      lslvnt=.false.
      ldistr=.false.
      lsolva=.false.
      lcryst=.false.
c
      nwm=0
      nwmc=0
      nwms=0
      nwa=0
      nsa=0
      nsm=0
      noe=0
      npbtyp=0
c
c     allocate memory for solvent
c
      if(.not.lcryst) then
      nwmc=0
      mwmc=1
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'xwc',l_xwc,i_xwc))
     + call md_abort('ma_alloc_get failed for xwc',9999)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'vwc',l_vwc,i_vwc))
     + call md_abort('ma_alloc_get failed for vwc',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wfnd',l_wfnd,i_wfnd))
     + call md_abort('ma_alloc_get failed for wfnd',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wcmr',l_wcmr,i_wcmr))
     + call md_abort('ma_alloc_get failed for wcmr',9999)
      lcryst=.true.
      endif
      if(.not.lsolva) then
      nwm=0
      mwm=1
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwm,'xw',l_xw,i_xw))
     + call md_abort('ma_alloc_get failed for xw',9999)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwm,'vw',l_vw,i_vw))
     + call md_abort('ma_alloc_get failed for vw',9999)
      if(.not.ma_alloc_get(mt_int,mwm,'iwmr',l_iwmr,i_iwmr))
     + call md_abort('ma_alloc_get failed for iwmr',9999)
      lsolva=.true.
      endif
c
      if(util_print('restart',print_debug)) then
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
c     set protonation index
c     ---------------------
c
c     returned in ips: number of sets to skip to get to next
c     protonation state if this proton hops off.
c
      if(itopol.eq.3) then
      call qhop_pindex(msa,nsa,int_mb(i_sgm),int_mb(i_ihop),
     + int_mb(i_ips))
      endif
c
c     open the command file
c     ---------------------
c
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old',err=99)
      goto 11
   99 continue
      if(util_print('restart',print_default)) then
      write(lfnout,1025)
 1025 format(' No command file found: Default restart directives',/)
      endif
c
c     write restart when no command file is found
c     ===========================================
c
      if(util_print('restart',print_debug)) then
      write(lfnout,1004) source(1:index(source,' ')-1)
 1004 format(' No command file for',t40,a,/)
      endif
      if(source(1:8).eq.'geometry') then
      if(.not.argos_prep_rdgeom(irtdb,byte_mb(i_csa),int_mb(i_sar),
     + int_mb(i_sgm),int_mb(i_sfnd),dbl_mb(i_xs),dbl_mb(i_vs),
     + msa,nsa,nmis,msm,nsm))
     + call md_abort('argos_prep_rdgeom failed',9999)
      if(util_print('restart',print_debug)) then
      write(lfnout,1005)
 1005 format(' Geometry read from rtdb')
      endif
      else
      if(.not.argos_prep_pdbsiz(lfnpdb,filpdb,lfnout,model,
     + idum1,idum2,idum3,idum4,slvnam,nwmc,altloc,chain))
     + call md_abort('argos_prep_pdbsiz failed',9999)
c
c     allocate memory for crystal solvent
c     -----------------------------------
c
      if(lcryst) then
      if(.not.ma_free_heap(l_xwc))
     + call md_abort('ma_free_heap failed on l_xwc',9999)
      if(.not.ma_free_heap(l_vwc))
     + call md_abort('ma_free_heap failed on l_vwc',9999)
      if(.not.ma_free_heap(l_wfnd))
     + call md_abort('ma_free_heap failed on l_wfnd',9999)
      if(.not.ma_free_heap(l_wcmr))
     + call md_abort('ma_free_heap failed on l_wcmr',9999)
      lcryst=.false.
      endif
      mwmc=nwmc
      if(nwmc.eq.0) mwmc=1
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'xwc',l_xwc,i_xwc))
     + call md_abort('ma_alloc_get failed for xwc',9999)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'vwc',l_vwc,i_vwc))
     + call md_abort('ma_alloc_get failed for vwc',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wfnd',l_wfnd,i_wfnd))
     + call md_abort('ma_alloc_get failed for wfnd',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wcmr',l_wcmr,i_wcmr))
     + call md_abort('ma_alloc_get failed for wcmr',9999)
      lcryst=.true.
c
c     read the pdb file
c     -----------------
c
      if(.not.argos_prep_rdpdb(lfnout,filpdb,lfnpdb,filtmp,
     + lfntmp,.true.,
     + byte_mb(i_csa),int_mb(i_sar),int_mb(i_sgm),int_mb(i_sfnd),
     + dbl_mb(i_xs),dbl_mb(i_vs),msa,nsa,nmis,msm,nsm,
     + int_mb(i_wcmr),int_mb(i_wfnd),dbl_mb(i_xwc),dbl_mb(i_vwc),
     + byte_mb(i_cwa),mwmc,mwa,nwmc,nwa,slvnam,iwater,altloc,
     + chain,model,box,lpdbbx))
     + call md_abort('argos_prep_rdpdb failed',9999)
      if(util_print('restart',print_debug)) then
      write(lfnout,1006) filpdb(1:index(filpdb,' ')-1),nmis,iwater
 1006 format(' PDB file',t40,a,2i7/)
      endif
      if(nmis.gt.0.or.iwater.gt.0) then
      if(.not.ma_push_get(mt_dbl,4*mgrid*mgrid*mgrid,'grid',l_grid,
     + i_grid)) call md_abort('Memory allocation failed for grid',9999)
       if(.not.argos_prep_misfit(lfnout,lfnmat,int_mb(i_sfnd),
     + dbl_mb(i_xs),
     + byte_mb(i_csa),dbl_mb(i_qsa),int_mb(i_sgm),msa,nsa,
     + int_mb(i_isb),dbl_mb(i_csb),msb,nsb,
     + dbl_mb(i_grid),mgrid,ngrid,gdist,nmis,iwater,int_mb(i_wfnd),
     + dbl_mb(i_xwc),dbl_mb(i_qwa),mwmc,mwa,nwmc,nwa,npbtyp,box,
     + fcount,nrgrid,iogrid,rogrid))
     + call md_abort('argos_prep_misfit failed',9999)
      if(.not.ma_pop_stack(l_grid))
     + call md_abort('Memory deallocation failed for grid',9999)
      endif
      endif
c
      if(.not.lpdbbx) then
      if(nsa.gt.0) then
c
c     center the solute
c     -----------------
c
      if(.not.argos_prep_center(lfnout,dbl_mb(i_xw),mwm,mwa,nwm,nwa,
     + dbl_mb(i_xwc),mwmc,nwmc,dbl_mb(i_xs),msa,nsa,0))
     + call md_abort('argos_prep_center failed',9999)
c
c     orient the solute
c     -----------------
c
cx      if(.not.argos_prep_orient(lfnout,dbl_mb(i_xw),mwm,mwa,nwm,nwa,
cx     + dbl_mb(i_xwc),mwmc,nwmc,dbl_mb(i_xs),msa,nsa))
cx     + call md_abort('argos_prep_orient failed',9999)
c
      endif
c
c     determine the boxsize
c     ---------------------
c
      if(.not.argos_prep_boxsiz(dbl_mb(i_xw),dbl_mb(i_xwc),mwm,mwmc,mwa,
     + nwm,nwmc,nwa,dbl_mb(i_xs),msa,nsa,box,0))
     +  call md_abort('argos_prep_boxsiz failed',9999)
      if(util_print('restart',print_medium)) then
      write(lfnout,1007) box
 1007 format(/,' Boxsize determined to ',t40,3f12.6,/)
      endif
      endif
c
c     write the restart file
c     ----------------------
c
      if(.not.argos_prep_his(int_mb(i_num),mat,
     + byte_mb(i_csa),int_mb(i_sat),int_mb(i_sgm),int_mb(i_sml),
     + int_mb(i_sfr),dbl_mb(i_xs),msa,nsa))
     +  call md_abort('Error in argos_prep_his',0)
      if(mdold.gt.0) then
      if(.not.argos_prep_wtrst(lfnout,lfnrst,filrst,npbtyp,nbxtyp,box,
     + int_mb(i_iwmr),dbl_mb(i_xw),dbl_mb(i_vw),mwm,mwa,nwm,nwa,
     + int_mb(i_wcmr),dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc,
     + int_mb(i_sar),dbl_mb(i_xs),dbl_mb(i_vs),msa,nsa,
     + int_mb(i_inoe),dbl_mb(i_dnoe),mnoe,noe,msm,nsm))
     + call md_abort('argos_prep_wtrst failed',9999)
      else
      if(.not.argos_prep_wrtrst(lfnout,lfnrst,filrst,title,
     + npbtyp,nbxtyp,box,
     + int_mb(i_iwmr),dbl_mb(i_xw),dbl_mb(i_vw),mwm,mwa,nwm,nwa,
     + int_mb(i_wcmr),dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc,
     + int_mb(i_sar),int_mb(i_sgm),dbl_mb(i_xs),dbl_mb(i_vs),
     + msa,nsa,msm,nsm,irrep,nxrep,nyrep,nzrep,
     + mseq,nseq,int_mb(i_lseq),int_mb(i_ips),int_mb(i_istat),
     + lfnmod,filmod,itopol))
     + call md_abort('argos_prep_wrtrst failed',9999)
      endif
      goto 9
c
c     loop through commands from command file
c     =======================================
c
   11 continue
c
      if(util_print('restart',print_default)) then
      write(lfnout,1026)
 1026 format(' Restart generation using command file directives',/)
      endif
c
    1 continue
c
      read(lfncmd,1101,end=9,err=9997) card
 1101 format(a)
c
      if(util_print('restart',print_high)) then
      write(lfnout,1099) card
 1099 format(' Command option executing: ',a,/)
      endif
c
c     read distance restraints
c     ========================
c
      if(card(1:6).eq.'disres') then
      read(lfncmd,1101,end=9,err=9997) card2
      if(mdold.gt.0) then
      if(.not.argos_prep_disres(card,card2,int_mb(i_sgm),
     + byte_mb(i_csa),msa,
     + nsa,int_mb(i_inoe),dbl_mb(i_dnoe),mnoe,nnoe))
     + call md_abort('argos_prep_disres failed',9999)
      endif
      endif
c
c     read solvent coordinates
c     ========================
c
      if(card(1:6).eq.'solvnt') then
c
c     free solvent memory if previously allocated
c     -------------------------------------------
c
      if(lslvnt) then
      if(.not.ma_free_heap(l_xslv))
     + call md_abort('ma_free_heap failed on l_xslv',9999)
      if(.not.ma_free_heap(l_vslv))
     + call md_abort('ma_free_heap failed on l_vslv',9999)
      lslvnt=.false.
      endif
c
c     check if solvent restart file exists locally
c     --------------------------------------------
c
      filslv=slvmdl(1:index(slvmdl,' ')-1)//'.rst '
      open(unit=lfnslv,file=filslv(1:index(filslv,' ')-1),
     + status='old',err=1020)
      close(unit=lfnslv)
      write(lfnout,1008) filslv(1:index(filslv,' ')-1)
      goto 1021
 1020 continue
c
c     get restart file name for the solvent
c     -------------------------------------
c
      if(.not.util_nwchemrc_get(slvmdl,filslv))
     + call md_abort('util_nwchemrc_get failed on slvmdl',9999)
      if(util_print('topology',print_debug)) then
      write(lfnout,1008) filslv(1:index(filslv,' ')-1)
 1008 format(' Solvent coordinates from ',t40,a,/)
      endif
c
c     check if solvent restart file exists
c     ------------------------------------
c
      open(unit=lfnslv,file=filslv(1:index(filslv,' ')-1),
     + status='old',err=9998)
      close(unit=lfnslv)
c
 1021 continue
c
c     get dimensions from solvent restart file
c     ----------------------------------------
c
      if(.not.argos_prep_slvsiz(lfnslv,filslv,nwms,nwas))
     + call md_abort('argos_prep_slvsiz failed',9999)
      if(nwas.ne.nwa) call md_abort('incompatible solvent',9999)
      mwms=nwms
c
c     allocate memory for solvent coordinates and velocities
c     ------------------------------------------------------
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwms,'xslv',l_xslv,i_xslv))
     + call md_abort('ma_alloc_get failed for xslv',9999)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwms,'vslv',l_vslv,i_vslv))
     + call md_abort('ma_alloc_get failed for vslv',9999)
c
c     read solvent coordinates and velocities
c     ---------------------------------------
c
      if(.not.argos_prep_slvrd(lfnslv,filslv,dbl_mb(i_xslv),
     + dbl_mb(i_vslv),
     + mwms,mwa,nwms,nwa,boxs))
     + call md_abort('argos_prep_slvrd failed',9999)
      lslvnt=.true.
      if(util_print('restart',print_medium)) then
      write(lfnout,1009) filslv(1:index(filslv,' ')-1),nwms,boxs
 1009 format(/,' Solvent file',t40,a,/,
     + ' Number of solvent molecules',t40,i8,/,
     + ' Solvent box size',t40,3f12.6)
      endif
      endif
c
c     read pdb file
c     =============
c
      if(card(1:6).eq.'rd_pdb') then
      filnam=card(8:80)
      length=index(filnam,' ')-1
      if(length.eq.0) filnam=filpdb
      len=index(filpdb,' ')-1
      lnames=length.eq.len.and.filnam(1:len).eq.filpdb(1:len)
c
c     read the pdb file
c     -----------------
c
      if(.not.argos_prep_pdbsiz(lfnpdb,filnam,lfnout,model,
     + idum1,idum2,idum3,idum4,slvnam,nwmc,altloc,chain))
     + call md_abort('argos_prep_pdbsiz failed',9999)
c
c     allocate memory for crystal solvent
c     -----------------------------------
c
      if(lcryst) then
      if(.not.ma_free_heap(l_xwc))
     + call md_abort('ma_free_heap failed on l_xwc',9999)
      if(.not.ma_free_heap(l_vwc))
     + call md_abort('ma_free_heap failed on l_vwc',9999)
      if(.not.ma_free_heap(l_wfnd))
     + call md_abort('ma_free_heap failed on l_wfnd',9999)
      if(.not.ma_free_heap(l_wcmr))
     + call md_abort('ma_free_heap failed on l_wcmr',9999)
      lcryst=.false.
      endif
      mwmc=nwmc
      if(nwmc.eq.0) mwmc=1
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'xwc',l_xwc,i_xwc))
     + call md_abort('ma_alloc_get failed for xwc',9999)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'vwc',l_vwc,i_vwc))
     + call md_abort('ma_alloc_get failed for vwc',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wfnd',l_wfnd,i_wfnd))
     + call md_abort('ma_alloc_get failed for wfnd',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wcmr',l_wcmr,i_wcmr))
     + call md_abort('ma_alloc_get failed for wcmr',9999)
      lcryst=.true.
c
      if(.not.argos_prep_rdpdb(lfnout,filnam,lfnpdb,
     + filtmp,lfntmp,lnames,
     + byte_mb(i_csa),int_mb(i_sar),int_mb(i_sgm),int_mb(i_sfnd),
     + dbl_mb(i_xs),dbl_mb(i_vs),msa,nsa,nmis,msm,nsm,
     + int_mb(i_wcmr),int_mb(i_wfnd),dbl_mb(i_xwc),dbl_mb(i_vwc),
     + byte_mb(i_cwa),mwmc,mwa,nwmc,nwa,slvnam,iwater,altloc,
     + chain,model,box,lpdbbx))
     + call md_abort('argos_prep_rdpdb failed',9999)
      if(util_print('restart',print_default)) then
      write(lfnout,1010) filnam(1:index(filnam,' ')-1),nmis,nsa
 1010 format(' PDB file',t40,a,/,
     + ' Number of missing coordinates',t40,i8,/,
     + ' Total number of atoms',t40,i8,/)
      endif
      if(nmis.gt.0.or.iwater.gt.0) then
      if(.not.ma_push_get(mt_dbl,4*mgrid*mgrid*mgrid,'grid',l_grid,
     + i_grid)) call md_abort('Memory allocation failed for grid',9999)
      if(.not.argos_prep_misfit(lfnout,lfnmat,int_mb(i_sfnd),
     + dbl_mb(i_xs),
     + byte_mb(i_csa),dbl_mb(i_qsa),int_mb(i_sgm),msa,nsa,
     + int_mb(i_isb),dbl_mb(i_csb),msb,nsb,
     + dbl_mb(i_grid),mgrid,ngrid,gdist,nmis,iwater,int_mb(i_wfnd),
     + dbl_mb(i_xwc),dbl_mb(i_qwa),mwmc,mwa,nwmc,nwa,npbtyp,box,
     + fcount,nrgrid,iogrid,rogrid))
     + call md_abort('argos_prep_misfit failed',9999)
      if(util_print('restart',print_default)) then
      write(lfnout,1023)
 1023 format(' Coordinates generated for missing atoms',/)
      endif
      if(.not.ma_pop_stack(l_grid))
     + call md_abort('Memory deallocation failed for grid',9999)
      endif
      endif
c
c     set box
c
      if(card(1:6).eq.'setbox') then
      read(card(7:42),'(3f12.6)') boxset
      if(boxset(1).gt.0.0d0) box(1)=boxset(1)
      if(boxset(2).gt.0.0d0) box(2)=boxset(2)
      if(boxset(3).gt.0.0d0) box(3)=boxset(3)
      if(util_print('restart',print_medium)) then
      write(lfnout,1050) box
 1050 format(' Box dimensions set to ',3f12.6)
      endif
      endif
c
c     read restart file
c     ------------------
c
      if(card(1:6).eq.'rd_rst') then
      filnam=card(8:80)
      length=index(filnam,' ')-1
      if(length.eq.0) filnam=filrst
      if(.not.argos_prep_rstsiz(lfnrst,filnam,nnwm,nnwa,nnsa,nnwmc))
     + call md_abort('argos_prep_rstsiz failed',9999)
      if(nnwm.gt.mwm) then
      if(lsolva) then
      if(.not.ma_free_heap(l_xw))
     + call md_abort('ma_free_heap failed on l_xw',9999)
      if(.not.ma_free_heap(l_vw))
     + call md_abort('ma_free_heap failed on l_vw',9999)
      if(.not.ma_free_heap(l_iwmr))
     + call md_abort('ma_free_heap failed on l_iwmr',9999)
      lsolva=.false.
      endif
      mwm=nnwm
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwm,'xw',l_xw,i_xw))
     + call md_abort('ma_alloc_get failed for xw',9999)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwm,'vw',l_vw,i_vw))
     + call md_abort('ma_alloc_get failed for vw',9999)
      if(.not.ma_alloc_get(mt_int,mwm,'iwmr',l_iwmr,i_iwmr))
     + call md_abort('ma_alloc_get failed for iwmr',9999)
      lsolva=.true.
      endif
      if(lcryst) then
      if(.not.ma_free_heap(l_xwc))
     + call md_abort('ma_free_heap failed on l_xwc',9999)
      if(.not.ma_free_heap(l_vwc))
     + call md_abort('ma_free_heap failed on l_vwc',9999)
      if(.not.ma_free_heap(l_wfnd))
     + call md_abort('ma_free_heap failed on l_wfnd',9999)
      if(.not.ma_free_heap(l_wcmr))
     + call md_abort('ma_free_heap failed on l_wcmr',9999)
      lcryst=.false.
      endif
      nwmc=0
      mwmc=nnwmc
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'xwc',l_xwc,i_xwc))
     + call md_abort('ma_alloc_get failed for xwc',9999)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'vwc',l_vwc,i_vwc))
     + call md_abort('ma_alloc_get failed for vwc',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wfnd',l_wfnd,i_wfnd))
     + call md_abort('ma_alloc_get failed for wfnd',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wcmr',l_wcmr,i_wcmr))
     + call md_abort('ma_alloc_get failed for wcmr',9999)
      lcryst=.true.
      if(.not.argos_prep_rdrst(lfnout,lfnrst,filnam,npbtyp,nbxtyp,box,
     + int_mb(i_iwmr),dbl_mb(i_xw),dbl_mb(i_vw),mwm,mwa,nwm,nwa,
     + int_mb(i_wcmr),dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc,
     + int_mb(i_sar),dbl_mb(i_xs),dbl_mb(i_vs),msa,nsa,
     + int_mb(i_inoe),dbl_mb(i_dnoe),mnoe,noe,msm,nsm,
     + mseq,nseq,int_mb(i_lseq),int_mb(i_ips)))
     + call md_abort('argos_prep_rdrst failed',9999)
      lpdbbx=.true.
      endif
c
c     read restart file
c     ------------------
c
      if(card(1:6).eq.'rd_old') then
      filnam=card(8:80)
      length=index(filnam,' ')-1
      if(length.eq.0) filnam=filrst
      if(.not.argos_prep_rstsz(lfnrst,filnam,nnwm,nnwa,nnsa,nnwmc))
     + call md_abort('argos_prep_rstsiz failed',9999)
      if(nnwm.gt.mwm) then
      if(lsolva) then
      if(.not.ma_free_heap(l_xw))
     + call md_abort('ma_free_heap failed on l_xw',9999)
      if(.not.ma_free_heap(l_vw))
     + call md_abort('ma_free_heap failed on l_vw',9999)
      if(.not.ma_free_heap(l_iwmr))
     + call md_abort('ma_free_heap failed on l_iwmr',9999)
      lsolva=.false.
      endif
      mwm=nnwm
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwm,'xw',l_xw,i_xw))
     + call md_abort('ma_alloc_get failed for xw',9999)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwm,'vw',l_vw,i_vw))
     + call md_abort('ma_alloc_get failed for vw',9999)
      if(.not.ma_alloc_get(mt_int,mwm,'iwmr',l_iwmr,i_iwmr))
     + call md_abort('ma_alloc_get failed for iwmr',9999)
      lsolva=.true.
      endif
      if(lcryst) then
      if(.not.ma_free_heap(l_xwc))
     + call md_abort('ma_free_heap failed on l_xwc',9999)
      if(.not.ma_free_heap(l_vwc))
     + call md_abort('ma_free_heap failed on l_vwc',9999)
      if(.not.ma_free_heap(l_wfnd))
     + call md_abort('ma_free_heap failed on l_wfnd',9999)
      if(.not.ma_free_heap(l_wcmr))
     + call md_abort('ma_free_heap failed on l_wcmr',9999)
      lcryst=.false.
      endif
      nwmc=0
      mwmc=nnwmc
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'xwc',l_xwc,i_xwc))
     + call md_abort('ma_alloc_get failed for xwc',9999)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'vwc',l_vwc,i_vwc))
     + call md_abort('ma_alloc_get failed for vwc',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wfnd',l_wfnd,i_wfnd))
     + call md_abort('ma_alloc_get failed for wfnd',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wcmr',l_wcmr,i_wcmr))
     + call md_abort('ma_alloc_get failed for wcmr',9999)
      lcryst=.true.
      if(.not.argos_prep_rrst(lfnout,lfnrst,filnam,npbtyp,nbxtyp,box,
     + int_mb(i_iwmr),dbl_mb(i_xw),dbl_mb(i_vw),mwm,mwa,nwm,nwa,
     + int_mb(i_wcmr),dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc,
     + int_mb(i_sar),dbl_mb(i_xs),dbl_mb(i_vs),msa,nsa,
     + int_mb(i_inoe),dbl_mb(i_dnoe),mnoe,noe,msm,nsm))
     + call md_abort('argos_prep_rdrst failed',9999)
      lpdbbx=.true.
      endif
c
c     write pdb file
c     --------------
c
      if(card(1:6).eq.'wr_pdb') then
      read(card(7:8),'(2i1)') irenum,lrgpdb
      if(card(9:14).eq.'solute') then
      read(card(16:25),'(i10)') iopt
      if(iopt.lt.0) iopt=nwmc
      filnam=card(26:80)
      else
      filnam=card(9:80)
      iopt=nwm+nwmc
      endif
      length=index(filnam,' ')-1
      if(length.eq.0) filnam=filpdb
      len=index(filpdb,' ')-1
      if(irrand.gt.0) then
      if(.not.argos_prep_boxsiz(dbl_mb(i_xw),dbl_mb(i_xwc),mwm,mwmc,mwa,
     + nwm,nwmc,nwa,dbl_mb(i_xs),msa,nsa,box,irrand))
     +  call md_abort('argos_prep_boxsiz failed',9999)
      if((irrand.eq.1.or.irrand.eq.4).and.drep(1).lt.box(1))
     + drep(1)=box(1)
      if((irrand.eq.2.or.irrand.eq.4).and.drep(2).lt.box(2))
     + drep(2)=box(2)
      if((irrand.eq.3.or.irrand.eq.4).and.drep(3).lt.box(3))
     + drep(3)=box(3)
      endif
      mlang=nxrep*nyrep*iabs(nzrep)
      if(.not.ma_push_get(mt_int,mlang,'lang',l_lang,i_lang))
     + call md_abort('ma_push_get failed for lang',9999)
      if(.not.argos_prep_wrtpdb(lfnout,lfnpdb,lrgpdb,filnam,iopt,box,
     + int_mb(i_num),dbl_mb(i_mas),
     + mat,byte_mb(i_csa),int_mb(i_sat),int_mb(i_sgm),int_mb(i_sml),
     + int_mb(i_sfr),
     + dbl_mb(i_xs),dbl_mb(i_vs),msa,nsa,byte_mb(i_cwa),
     + int_mb(i_wat),dbl_mb(i_xw),dbl_mb(i_vw),
     + mwm,mwa,nwm,nwa,dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc,slvnam,
     + iropt,irrand,nxrep,nyrep,nzrep,drep,msb,nsb,int_mb(i_isb),zdist,
     + nskip,iskip,int_mb(i_lang),lfnmrg,nmerge,xmerge,filmrg,
     + irenum,invert,int_mb(i_ihop),int_mb(i_ips)))
     +  call md_abort('argos_prep_wrtpdb failed',9999)
      if(.not.ma_pop_stack(l_lang))
     + call md_abort('ma_pop_stack failed for lang',9999)
      nmerge=0
      endif
c
      if(card(1:6).eq.'merge ') then
      nmerge=nmerge+1
      read(card(7:30),'(3f8.3)') (xmerge(i,nmerge),i=1,3)
      filmrg(nmerge)=card(31:80)
      endif
c
c     write pov file
c     --------------
c
      if(card(1:6).eq.'wr_pov') then
      filnam=card(8:80)
      length=index(filnam,' ')-1
      if(length.eq.0) filnam=filpov
      len=index(filpov,' ')-1
      if(.not.argos_prep_wrtpov(lfnout,lfnpov,filnam,iopt,box,
     + int_mb(i_num),dbl_mb(i_mas),
     + mat,byte_mb(i_csa),int_mb(i_sat),int_mb(i_sgm),dbl_mb(i_xs),
     + dbl_mb(i_vs),msa,nsa,byte_mb(i_cwa),
     + int_mb(i_wat),dbl_mb(i_xw),dbl_mb(i_vw),
     + mwm,mwa,nwm,nwa,dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc,slvnam,
     + nxrep,nyrep,nzrep,drep,msb,nsb,int_mb(i_isb),zdist,scale,cpk))
     + call md_abort('argos_prep_wrtpov failed',9999)
      endif
c
c     read coordinates
c     ================
c
      if(card(1:6).eq.'rdcoor') then
c
c     determine the coordinate source
c     -------------------------------
c
      if(source(1:8).eq.'geometry') then
c
c     read coordinates from rtdb
c     --------------------------
c
      if(.not.argos_prep_rdgeom(irtdb,byte_mb(i_csa),int_mb(i_sar),
     + int_mb(i_sgm),int_mb(i_sfnd),dbl_mb(i_xs),dbl_mb(i_vs),
     + msa,nsa,nmis,msm,nsm))
     + call md_abort('argos_prep_rdgeom failed',9999)
      else
c
c     free memory for crystal solvent if already allocated
c     ----------------------------------------------------
c
      if(lcryst) then
      if(.not.ma_free_heap(l_xwc))
     + call md_abort('ma_free_heap failed on l_xwc',9999)
      if(.not.ma_free_heap(l_vwc))
     + call md_abort('ma_free_heap failed on l_vwc',9999)
      if(.not.ma_free_heap(l_wfnd))
     + call md_abort('ma_free_heap failed on l_wfnd',9999)
      if(.not.ma_free_heap(l_wcmr))
     + call md_abort('ma_free_heap failed on l_wcmr',9999)
      lcryst=.false.
      endif
c
c     determine the dimensions on the pdb file
c     ----------------------------------------
c
      if(.not.argos_prep_pdbsiz(lfnpdb,filpdb,lfnout,0,
     + idum1,idum2,idum3,idum4,slvnam,nwmc,altloc,chain))
     + call md_abort('argos_prep_pdbsiz failed',9999)
c
c     allocate memory for crystal solvent
c     -----------------------------------
c
      if(nwmc.gt.0) then
      if(lcryst) then
      if(.not.ma_free_heap(l_xwc))
     + call md_abort('ma_free_heap failed on l_xwc',9999)
      if(.not.ma_free_heap(l_vwc))
     + call md_abort('ma_free_heap failed on l_vwc',9999)
      if(.not.ma_free_heap(l_wfnd))
     + call md_abort('ma_free_heap failed on l_wfnd',9999)
      if(.not.ma_free_heap(l_wcmr))
     + call md_abort('ma_free_heap failed on l_wcmr',9999)
      lcryst=.false.
      endif
      mwmc=nwmc
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'xwc',l_xwc,i_xwc))
     + call md_abort('ma_alloc_get failed for xwc',9999)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwmc,'vwc',l_vwc,i_vwc))
     + call md_abort('ma_alloc_get failed for vwc',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wfnd',l_wfnd,i_wfnd))
     + call md_abort('ma_alloc_get failed for wfnd',9999)
      if(.not.ma_alloc_get(mt_int,mwmc,'wcmr',l_wcmr,i_wcmr))
     + call md_abort('ma_alloc_get failed for wcmr',9999)
      lcryst=.true.
      endif
c
c     read the pdb file
c     -----------------
c
      if(.not.argos_prep_rdpdb(lfnout,filpdb,lfnpdb,filtmp,
     + lfntmp,.true.,
     + byte_mb(i_csa),int_mb(i_sar),int_mb(i_sgm),int_mb(i_sfnd),
     + dbl_mb(i_xs),dbl_mb(i_vs),msa,nsa,nmis,msm,nsm,
     + int_mb(i_wcmr),int_mb(i_wfnd),dbl_mb(i_xwc),dbl_mb(i_vwc),
     + byte_mb(i_cwa),mwmc,mwa,nwmc,nwa,slvnam,iwater,altloc,
     + chain,model,box,lpdbbx))
     + call md_abort('argos_prep_rdpdb failed',9999)
      if(util_print('restart',print_default)) then
      write(lfnout,1011) filpdb(1:index(filpdb,' ')-1),nmis,nsa
 1011 format(' PDB file',t40,a,/,
     + ' Number of missing coordinates',t40,i8,/,
     + ' Total number of atoms',t40,i8,/)
      endif
      if(nmis.gt.0.or.iwater.gt.0) then
      if(util_print('restart',print_default)) then
      write(lfnout,1027) nmis,iwater
 1027 format(' Generating missing atom coordinates: ',2i5,/)
      endif
      if(.not.ma_push_get(mt_dbl,4*mgrid*mgrid*mgrid,'grid',l_grid,
     + i_grid)) call md_abort('Memory allocation failed for grid',9999)
      if(.not.argos_prep_misfit(lfnout,lfnmat,int_mb(i_sfnd),
     + dbl_mb(i_xs),
     + byte_mb(i_csa),dbl_mb(i_qsa),int_mb(i_sgm),msa,nsa,
     + int_mb(i_isb),dbl_mb(i_csb),msb,nsb,
     + dbl_mb(i_grid),mgrid,ngrid,gdist,nmis,iwater,int_mb(i_wfnd),
     + dbl_mb(i_xwc),dbl_mb(i_qwa),mwmc,mwa,nwmc,nwa,npbtyp,box,
     + fcount,nrgrid,iogrid,rogrid))
     + call md_abort('argos_prep_misfit failed',9999)
      if(util_print('restart',print_default)) then
      write(lfnout,1024)
 1024 format(' Coordinates generated for missing atoms',/)
      endif
      if(.not.ma_pop_stack(l_grid))
     + call md_abort('Memory deallocation failed for grid',9999)
      endif
      endif
c
c     determine the boxsize
c     ---------------------
c
      if(.not.lpdbbx) then
      if(.not.argos_prep_boxsiz(dbl_mb(i_xw),dbl_mb(i_xwc),mwm,mwmc,mwa,
     + nwm,nwmc,nwa,dbl_mb(i_xs),msa,nsa,box,0))
     + call md_abort('argos_prep_boxsiz failed',9999)
      if(util_print('restart',print_medium)) then
      write(lfnout,1012) box
 1012 format(/,' Boxsize determined to ',t40,3f12.6,/)
      endif
      endif
      endif
c
c     write restart file
c     ------------------
c
      if(card(1:6).eq.'wr_rst') then
      filnam=card(8:80)
      length=index(filnam,' ')-1
      if(length.eq.0) filnam=filrst
      if(.not.argos_prep_wrtrst(lfnout,lfnrst,filnam,title,
     + npbtyp,nbxtyp,box,
     + int_mb(i_iwmr),dbl_mb(i_xw),dbl_mb(i_vw),mwm,mwa,nwm,nwa,
     + int_mb(i_wcmr),dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc,
     + int_mb(i_sar),int_mb(i_sgm),dbl_mb(i_xs),dbl_mb(i_vs),
     + msa,nsa,msm,nsm,irrep,nxrep,nyrep,nzrep,
     + mseq,nseq,int_mb(i_lseq),int_mb(i_ips),int_mb(i_istat),
     + lfnmod,filmod,itopol))
     + call md_abort('argos_prep_wrtrst failed',9999)
      irrep=0
      endif
c
c     write restart file
c     ------------------
c
      if(card(1:6).eq.'wr_old') then
      filnam=card(8:80)
      length=index(filnam,' ')-1
      if(length.eq.0) filnam=filrst
      if(.not.argos_prep_wtrst(lfnout,lfnrst,filnam,npbtyp,nbxtyp,box,
     + int_mb(i_iwmr),dbl_mb(i_xw),dbl_mb(i_vw),mwm,mwa,nwm,nwa,
     + int_mb(i_wcmr),dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc,
     + int_mb(i_sar),dbl_mb(i_xs),dbl_mb(i_vs),msa,nsa,
     + int_mb(i_inoe),dbl_mb(i_dnoe),mnoe,noe,msm,nsm))
     + call md_abort('argos_prep_wtrst failed',9999)
      endif
c
c
c     write xyz file
c     --------------
c
      if(card(1:6).eq.'wr_xyz') then
      if(card(8:13).eq.'solute') then
      read(card(15:24),'(i10)') iopt
      if(iopt.lt.0) iopt=nwmc
      if(iopt.gt.nwmc) iopt=nwmc
      filnam=card(25:80)
      else
      filnam=card(8:80)
      endif
      length=index(filnam,' ')-1
      if(length.eq.0) filnam=filxyz
      filxyz=filnam
      if(.not.argos_prep_wrtxyz(lfnout,lfnxyz,filxyz,lfnqqq,filqqq,
     + int_mb(i_num),mat,
     + int_mb(i_sat),dbl_mb(i_xs),dbl_mb(i_qsa),msa,nsa,box))
     + call md_abort('argos_prep_wrtxyz failed',9999)
      endif
c
c     center solute
c     -------------
c
      if(card(1:6).eq.'center') then
      if(.not.argos_prep_center(lfnout,dbl_mb(i_xw),mwm,mwa,nwm,nwa,
     + dbl_mb(i_xwc),mwmc,nwmc,dbl_mb(i_xs),msa,nsa,0))
     + call md_abort('argos_prep_center failed',9999)
      endif
c
      if(card(1:6).eq.'centrx') then
      if(.not.argos_prep_center(lfnout,dbl_mb(i_xw),mwm,mwa,nwm,nwa,
     + dbl_mb(i_xwc),mwmc,nwmc,dbl_mb(i_xs),msa,nsa,1))
     + call md_abort('argos_prep_center failed',9999)
      endif
c
      if(card(1:6).eq.'centry') then
      if(.not.argos_prep_center(lfnout,dbl_mb(i_xw),mwm,mwa,nwm,nwa,
     + dbl_mb(i_xwc),mwmc,nwmc,dbl_mb(i_xs),msa,nsa,2))
     + call md_abort('argos_prep_center failed',9999)
      endif
c
      if(card(1:6).eq.'centrz') then
      if(.not.argos_prep_center(lfnout,dbl_mb(i_xw),mwm,mwa,nwm,nwa,
     + dbl_mb(i_xwc),mwmc,nwmc,dbl_mb(i_xs),msa,nsa,3))
     + call md_abort('argos_prep_center failed',9999)
      endif
c
c     orient solute
c     -------------
c
      if(card(1:6).eq.'orient') then
      if(.not.argos_prep_orient(lfnout,dbl_mb(i_xw),mwm,mwa,nwm,nwa,
     + dbl_mb(i_xwc),mwmc,nwmc,dbl_mb(i_xs),msa,nsa))
     + call md_abort('argos_prep_orient failed',9999)
      endif
c
c     translate molecule
c     ------------------
c
      if(card(1:6).eq.'tramol') then
      read(card(8:53),'(2i5,3f12.6)',err=9999) itran,jtran,xtran
      if(.not.argos_prep_transl(itran,jtran,xtran,
     + int_mb(i_sml),dbl_mb(i_xs),msa,nsa))
     + call md_abort('argos_prep_transl failed',9999)
      endif
c
c     translate segment
c     -----------------
c
      if(card(1:6).eq.'trasgm') then
      read(card(8:53),'(2i5,3f12.6)',err=9999) itran,jtran,xtran
      if(.not.argos_prep_transl(itran,jtran,xtran,
     + int_mb(i_sgm),dbl_mb(i_xs),msa,nsa))
     + call md_abort('argos_prep_transl failed',9999)
      endif
c
c     translate atom
c     --------------
c
      if(card(1:6).eq.'traatm') then
      read(card(8:53),'(2i5,3f12.6)',err=9999) itran,jtran,xtran
      if(.not.argos_prep_transl(itran,jtran,xtran,
     + int_mb(i_sat),dbl_mb(i_xs),msa,nsa))
     + call md_abort('argos_prep_transl failed',9999)
      endif
c
c     rotate molecule
c     ------------------
c
      if(card(1:6).eq.'rotmol') then
      read(card(8:65),'(2i5,4f12.6)',err=9999) irot,jrot,arot,xrot
      if(.not.argos_prep_rot(irot,jrot,arot,xrot,
     + int_mb(i_sml),dbl_mb(i_xs),msa,nsa))
     + call md_abort('argos_prep_rot failed',9999)
      endif
c
c     rotate segment
c     -----------------
c
      if(card(1:6).eq.'rotsgm') then
      read(card(8:65),'(2i5,4f12.6)',err=9999) irot,jrot,arot,xrot
      if(.not.argos_prep_rot(irot,jrot,arot,xrot,
     + int_mb(i_sgm),dbl_mb(i_xs),msa,nsa))
     + call md_abort('argos_prep_rot failed',9999)
      endif
c
c     rotate atom
c     --------------
c
      if(card(1:6).eq.'rotatm') then
      read(card(8:65),'(2i5,4f12.6)',err=9999) irot,jrot,arot,xrot
      if(.not.argos_prep_rot(irot,jrot,arot,xrot,
     + int_mb(i_sat),dbl_mb(i_xs),msa,nsa))
     + call md_abort('argos_prep_rot failed',9999)
      endif
c
c     crop atom
c     ---------
c
      if(card(1:4).eq.'crop') then
      if(.not.argos_prep_crop(dbl_mb(i_xs),msa,nsa,lfnout))
     + call md_abort('argos_prep_rot failed',9999)
      endif
c
c     collapse solute
c     ---------------
c
      if(card(1:10).eq.'collapsexy') then
      read(card(11:15),'(i5)') nmoves
      if(.not.ma_push_get(mt_int,3*msa,'ndx',l_ndx,i_ndx))
     + call md_abort('Memory allocation failed for ndx',9999)
      if(.not.ma_push_get(mt_dbl,6*msa,'vec',l_vec,i_vec))
     + call md_abort('Memory allocation failed for ndx',9999)
      call argos_prep_collaps(dbl_mb(i_xs),int_mb(i_sml),int_mb(i_ndx),
     + dbl_mb(i_vec),msa,nsa,touch,3,nmoves,ncolgr,icolgr,
     + lfnout,lfnpdb,lrgpdb,sysnam,filnam,iopt,box,
     + int_mb(i_num),dbl_mb(i_mas),
     + mat,byte_mb(i_csa),int_mb(i_sat),int_mb(i_sgm),int_mb(i_sml),
     + int_mb(i_sfr),
     + dbl_mb(i_vs),byte_mb(i_cwa),
     + int_mb(i_wat),dbl_mb(i_xw),dbl_mb(i_vw),
     + mwm,mwa,nwm,nwa,dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc,slvnam,
     + iropt,irrand,nxrep,nyrep,nzrep,drep,msb,nsb,int_mb(i_isb),zdist,
     + nskip,iskip,int_mb(i_lang),lfnmrg,nmerge,xmerge,filmrg,
     + irenum,invert,int_mb(i_ihop),int_mb(i_ips))
      if(.not.ma_pop_stack(l_vec))
     + call md_abort('Memory deallocation failed for vec',9999)
      if(.not.ma_pop_stack(l_ndx))
     + call md_abort('Memory deallocation failed for ndx',9999)
      endif
      if(card(1:9).eq.'collapsez') then
      read(card(11:15),'(i5)') nmoves
      if(.not.ma_push_get(mt_int,3*msa,'ndx',l_ndx,i_ndx))
     + call md_abort('Memory allocation failed for ndx',9999)
      if(.not.ma_push_get(mt_dbl,6*msa,'vec',l_vec,i_vec))
     + call md_abort('Memory allocation failed for ndx',9999)
      call argos_prep_collaps(dbl_mb(i_xs),int_mb(i_sml),int_mb(i_ndx),
     + dbl_mb(i_vec),msa,nsa,touch,1,nmoves,ncolgr,icolgr,
     + lfnout,lfnpdb,lrgpdb,sysnam,filnam,iopt,box,
     + int_mb(i_num),dbl_mb(i_mas),
     + mat,byte_mb(i_csa),int_mb(i_sat),int_mb(i_sgm),int_mb(i_sml),
     + int_mb(i_sfr),
     + dbl_mb(i_vs),byte_mb(i_cwa),
     + int_mb(i_wat),dbl_mb(i_xw),dbl_mb(i_vw),
     + mwm,mwa,nwm,nwa,dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc,slvnam,
     + iropt,irrand,nxrep,nyrep,nzrep,drep,msb,nsb,int_mb(i_isb),zdist,
     + nskip,iskip,int_mb(i_lang),lfnmrg,nmerge,xmerge,filmrg,
     + irenum,invert,int_mb(i_ihop),int_mb(i_ips))
      if(.not.ma_pop_stack(l_vec))
     + call md_abort('Memory deallocation failed for vec',9999)
      if(.not.ma_pop_stack(l_ndx))
     + call md_abort('Memory deallocation failed for ndx',9999)
      endif
      if(card(1:9).eq.'collapsgr') then
      ncolgr=ncolgr+1
      if(ncolgr.gt.10) call md_abort('Increase dimension mcolgr',ncolgr)
      read(card(11:20),'(2i5)') icolgr(1,ncolgr),icolgr(2,ncolgr)
      endif
c
c     remove solvent
c     --------------
c
      if(card(1:6).eq.'removw') then
      read(card(8:80),'(6f12.6,i1)') ((xr(i,j),j=1,2),i=1,3),k
      if(.not.argos_prep_removw(k,xr,dbl_mb(i_xw),dbl_mb(i_vw),
     + mwm,mwa,nwm,nwa,
     + dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc,dbl_mb(i_xs),msa,nsa))
     + call md_abort('Error in argos_prep_removw',0)
      endif
c
c     solvate solute
c     --------------
c
      if(card(1:6).eq.'solvat') then
c
c     determine type of solvation box
c     -------------------------------
c
      islv=0
      nbxtyp=0
      if(card(8:11).eq.'cube') then
      if(card(17:17).eq.' ') then
      box(1)=box(1)+rshell
      box(2)=box(2)+rshell
      box(3)=box(3)+rshell
      edge=max(box(1),box(2),box(3))
      else
      read(card(12:23),'(f12.6)',err=9999) edge
      endif
      box(1)=edge
      box(2)=edge
      box(3)=edge
      islv=1
      endif
      if(card(8:10).eq.'box') then
      if(card(16:16).eq.' ') then
      box(1)=box(1)+rshell
      box(2)=box(2)+rshell
      box(3)=box(3)+rshell
      else
      read(card(11:22),'(f12.6)',err=9999) box(1)
      read(card(23:34),'(f12.6)',err=9999) box(2)
      read(card(35:46),'(f12.6)',err=9999) box(3)
      endif
      islv=1
      endif
      if(card(8:13).eq.'sphere') then
      if(card(19:19).eq.' ') then
      box(1)=box(1)+rshell
      box(2)=box(2)+rshell
      box(3)=box(3)+rshell
      edge=max(box(1),box(2),box(3))
      else
      read(card(14:25),'(f12.6)',err=9999) edge
      endif
      box(1)=edge
      box(2)=edge
      box(3)=edge
      islv=2
      endif
      if(card(8:12).eq.'troct') then
      if(card(18:18).eq.' ') then
      box(1)=box(1)+rshell
      box(2)=box(2)+rshell
      box(3)=box(3)+rshell
      edge=max(box(1),box(2),box(3))
      else
      read(card(13:24),'(f12.6)',err=9999) edge
      endif
      box(1)=edge
      box(2)=edge
      box(3)=edge
      islv=3
      nbxtyp=2
      endif
      if(islv.eq.0) call md_abort('Illegal solvation type',9999)
      if(util_print('restart',print_debug)) then
      write(lfnout,1013) islv,box
 1013 format(' Solvent type',t40,i8,/,
     + ' Box size',t40,3f12.6,/)
      endif
c
c     free memory for system solvent if allocated previously
c     ------------------------------------------------------
c
      if(lsolva) then
      if(.not.ma_free_heap(l_xw))
     + call md_abort('ma_free_heap failed on l_xw',9999)
      if(.not.ma_free_heap(l_vw))
     + call md_abort('ma_free_heap failed on l_vw',9999)
      if(.not.ma_free_heap(l_iwmr))
     + call md_abort('ma_free_heap failed on l_iwmr',9999)
      lsolva=.false.
      endif
c
c     estimate the maximum number of system solvent molecules
c     -------------------------------------------------------
c
      if(.not.argos_prep_slvnum(islv,boxs,nwms,box,mwm))
     + call md_abort('argos_prep_slvnum failed',9999)
      nwm=0
      if(util_print('restart',print_debug)) then
      write(lfnout,1014) mwm
 1014 format(' Maximum number of solvent molecules',t40,i8,/)
      endif
c
c     allocate memory for system solvent
c     ----------------------------------
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwm,'xw',l_xw,i_xw))
     + call md_abort('ma_alloc_get failed for xw',9999)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwm,'vw',l_vw,i_vw))
     + call md_abort('ma_alloc_get failed for vw',9999)
      if(.not.ma_alloc_get(mt_int,mwm,'iwmr',l_iwmr,i_iwmr))
     + call md_abort('ma_alloc_get failed for iwmr',9999)
      lsolva=.true.
c
c     solvate the system
c     ------------------
c
      if(.not.argos_prep_solvat(islv,lfnout,boxs,dbl_mb(i_xslv),
     + dbl_mb(i_vslv),
     + mwms,mwa,nwms,nwa,touch,box,
     + int_mb(i_iwmr),dbl_mb(i_xw),dbl_mb(i_vw),mwm,nwm,
     + dbl_mb(i_xs),msa,nsa,
     + dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc))
     + call md_abort('argos_prep_solvat failed',9999)
      if(util_print('restart',print_low)) then
      write(lfnout,1015) nwm
 1015 format(/,' Number of solvent molecules',t40,i8,/)
      endif
      endif
c
c     set crystal solvent
c     -------------------
c
      if(card(1:6).eq.'crysts') then
      endif
c
c     expand box
c     ----------
c
      if(card(1:6).eq.'expand') then
      if(card(8:9).eq.'to') then
      read(card(11:22),'(f12.6)',end=9999,err=9999) edgex
      box(1)=edgex
      read(card(23:34),'(f12.6)',end=9999,err=9999) edgey
      box(2)=edgey
      read(card(35:46),'(f12.6)',end=9999,err=9999) edgez
      box(3)=edgez
      if(util_print('restart',print_medium)) then
      write(lfnout,1060) box
 1060 format(' Box expanded with to ',t40,3f12.6)
      endif
      elseif(card(8:11).eq.'with') then
      read(card(13:24),'(f12.6)',end=9999,err=9999) edgex
      box(1)=box(1)+edgex
      read(card(25:36),'(f12.6)',end=9999,err=9999) edgey
      box(2)=box(2)+edgey
      read(card(37:48),'(f12.6)',end=9999,err=9999) edgez
      box(3)=box(3)+edgez
      endif
      if(util_print('restart',print_medium)) then
      write(lfnout,1061) edgex,edgey,edgez,box
 1061 format(' Box expanded with ',t40,3f12.6,/,
     +       '                to ',t40,3f12.6,/)
      endif
      endif
c
c     set minimum box size
c     --------------------
c
c      if(.not.lpdbbx) then
      if(card(1:6).eq.'boxsiz') then
      if(.not.argos_prep_boxsiz(dbl_mb(i_xw),dbl_mb(i_xwc),mwm,mwmc,mwa,
     + nwm,nwmc,nwa,dbl_mb(i_xs),msa,nsa,box,0))
     +  call md_abort('argos_prep_boxsiz failed',9999)
      if(util_print('restart',print_medium)) then
      write(lfnout,1033) box
 1033 format(/,' Boxsize determined to ',t40,3f12.6,/)
      endif
      endif
c      endif
c
c     set cubic box
c     -------------
c
      if(card(1:5).eq.'cubic') then
      boxmax=max(box(1),box(2),box(3))
      box(1)=boxmax
      box(2)=boxmax
      box(3)=boxmax
      if(util_print('restart',print_medium)) then
      write(lfnout,1034) box
 1034 format(/,' Boxsize set to ',t40,3f12.6,/)
      endif
      endif
c
c     set periodic boundary conditions
c     --------------------------------
c
      if(card(1:6).eq.'period') then
      npbtyp=1
      if(util_print('restart',print_medium)) then
      write(lfnout,1031)
 1031 format(/,' Periodic boundary conditions on')
      endif
      endif
c
c     set vacuum conditions
c     ---------------------
c
      if(card(1:6).eq.'vacuum') then
      npbtyp=0
      if(util_print('restart',print_medium)) then
      write(lfnout,1032)
 1032 format(/,' Periodic boundary conditions off')
      endif
      endif
c
c     set velocities to zero
c     ----------------------
c
      if(card(1:6).eq.'zerovl') then
      endif
c
c     set coordinates to zero
c     -----------------------
c
      if(card(1:6).eq.'zerocd') then
      endif
c
c     set touch distance
c     ------------------
c
      if(card(1:5).eq.'touch') then
      read(card(7:19),'(f12.6)') touch
      if(util_print('restart',print_medium)) then
      write(lfnout,1016) touch
 1016 format(' Touch distance set to ',t40,f12.6,' nm')
      endif
      endif
c
c     rshell
c     ------
c
      if(card(1:6).eq.'rshell') then
      read(card(8:19),'(f12.6)') rshell
      if(util_print('restart',print_medium)) then
      write(lfnout,1022) rshell
 1022 format(' Rshell set to ',t40,f12.6,' nm',/)
      endif
      endif
c
c     fix
c     ---
c
      if(card(1:3).eq.'fix') then
      if(card(4:7).eq.'out ') then
      ifix=1
      else if(card(4:7).eq.'win ') then
      ifix=-1
      else
      call md_abort('Error in fix command',9999)
      endif
      if(card(8:11).eq.'atm ') then
      jfix=1
      else if(card(8:11).eq.'sgm ') then
      jfix=-1
      else
      call md_abort('Error in fix command',9999)
      endif
      read(card(12:23),'(f12.6)') rfix
      cfix=card(25:80)
      if(.not.argos_prep_fix(lfnout,ifix,jfix,rfix,cfix,
     + int_mb(i_iwmr),dbl_mb(i_xw),mwm,mwa,nwm,nwa,
     + int_mb(i_wcmr),dbl_mb(i_xwc),mwmc,nwmc,
     + int_mb(i_sar),int_mb(i_sgm),byte_mb(i_csa),dbl_mb(i_xs),msa,nsa))
     + call md_abort('argos_prep_fix failed',9999)
      endif
c
c     align
c     =====
c
      if(card(1:6).eq.'align ') then
      if(.not.argos_prep_align(card,lfnout,
     + byte_mb(i_csa),int_mb(i_sat),int_mb(i_sgm),dbl_mb(i_xs),
     + dbl_mb(i_vs),msa,nsa,dbl_mb(i_xw),dbl_mb(i_vw),
     + mwm,mwa,nwm,nwa,dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc))
     + call md_abort('argos_prep_align failed',9999)
      endif
c
c     repeat
c     ======
c
      if(card(1:6).eq.'repeat') then
      read(card(8:12),'(i5)') nxrep
      read(card(13:17),'(i5)') nyrep
      read(card(18:22),'(i5)') nzrep
      read(card(35:46),'(f12.6)') zdist
      read(card(47:51),'(i5)') iropt
      read(card(52:56),'(i5)') irrand
      read(card(57:61),'(i5)') invert
      read(card(62:66),'(i5)') irrep
      if(invert.gt.0) nzrep=-iabs(nzrep)
      if(irrep.eq.0) then
      write(*,1019) nxrep,nyrep,nzrep
 1019 format(/,' Repeating ',3i5,/)
      if(.not.argos_prep_repeat(card,lfnout,drep,box,
     + byte_mb(i_csa),int_mb(i_sat),int_mb(i_sgm),dbl_mb(i_xs),
     + dbl_mb(i_vs),msa,nsa,dbl_mb(i_xw),dbl_mb(i_vw),
     + mwm,mwa,nwm,nwa,dbl_mb(i_xwc),dbl_mb(i_vwc),mwmc,nwmc))
     + call md_abort('argos_prep_repeat failed',9999)
      nskip=0
      endif
      endif
c
c     skip
c     ====
c
      if(card(1:6).eq.'skip  ') then
      nskip=nskip+1
      if(nskip.gt.100) call md_abort('Increase iskip dimension',nskip)
      read(card(8:22),'(3i5)') (iskip(j,nskip),j=1,3)
      endif
c
c     end
c     ---
c
      if(card(1:3).ne.'end') goto 1
c
    9 continue
      close(unit=lfncmd,err=999)
  999 continue
c
c     deallocate memory
c     -----------------
c
      if(.not.ma_pop_stack(l_istat))
     + call md_abort('Memory deallocation failed for istat',9999)
      if(.not.ma_pop_stack(l_ihop))
     + call md_abort('Memory deallocation failed for ihop',9999)
      if(.not.ma_pop_stack(l_ips))
     + call md_abort('Memory deallocation failed for ips',9999)
      if(.not.ma_pop_stack(l_lseq))
     + call md_abort('Memory deallocation failed for lseq',9999)
      if(.not.ma_pop_stack(l_csb))
     + call md_abort('Memory deallocation failed for cdsb',9999)
      if(.not.ma_pop_stack(l_isb))
     + call md_abort('Memory deallocation failed for idsb',9999)
      if(.not.ma_pop_stack(l_vs))
     + call md_abort('Memory deallocation failed for vs',9999)
      if(.not.ma_pop_stack(l_xs))
     + call md_abort('Memory deallocation failed for xs',9999)
      if(.not.ma_pop_stack(l_qsa))
     + call md_abort('Memory deallocation failed for qsa',9999)
      if(.not.ma_pop_stack(l_sfnd))
     + call md_abort('Memory deallocation failed for sfnd',9999)
      if(.not.ma_pop_stack(l_sfr))
     + call md_abort('Memory deallocation failed for sfr',9999)
      if(.not.ma_pop_stack(l_sml))
     + call md_abort('Memory deallocation failed for sml',9999)
      if(.not.ma_pop_stack(l_sgm))
     + call md_abort('Memory deallocation failed for sgm',9999)
      if(.not.ma_pop_stack(l_sat))
     + call md_abort('Memory deallocation failed for sat',9999)
      if(.not.ma_pop_stack(l_wat))
     + call md_abort('Memory deallocation failed for wat',9999)
      if(.not.ma_pop_stack(l_sar))
     + call md_abort('Memory deallocation failed for sar',9999)
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
      if(lslvnt) then
      if(.not.ma_free_heap(l_xslv))
     + call md_abort('ma_free_heap failed on l_xslv',9999)
      if(.not.ma_free_heap(l_vslv))
     + call md_abort('ma_free_heap failed on l_vslv',9999)
      lslvnt=.false.
      endif
      if(lsolva) then
      if(.not.ma_free_heap(l_xw))
     + call md_abort('ma_free_heap failed on l_xw',9999)
      if(.not.ma_free_heap(l_vw))
     + call md_abort('ma_free_heap failed on l_vw',9999)
      if(.not.ma_free_heap(l_iwmr))
     + call md_abort('ma_free_heap failed on l_iwmr',9999)
      lsolva=.false.
      endif
      if(lcryst) then
      if(.not.ma_free_heap(l_xwc))
     + call md_abort('ma_free_heap failed on l_xwc',9999)
      if(.not.ma_free_heap(l_vwc))
     + call md_abort('ma_free_heap failed on l_vwc',9999)
      if(.not.ma_free_heap(l_wfnd))
     + call md_abort('ma_free_heap failed on l_wfnd',9999)
      if(.not.ma_free_heap(l_wcmr))
     + call md_abort('ma_free_heap failed on l_wcmr',9999)
      lcryst=.false.
      endif
c
      argos_prep_mkrst=.true.
      return
c
 9997 continue
      write(lfnout,1017) filcmd(1:index(filcmd,' ')-1)
 1017 format(' Error reading commands',t40,a,/)
      argos_prep_mkrst=.false.
      return
c
 9998 continue
      write(lfnout,1018) filslv(1:index(filslv,' ')-1)
 1018 format(' Error reading solvent',t40,a,/)
      argos_prep_mkrst=.false.
      return
c
 9999 continue
      argos_prep_mkrst=.false.
      return
      end

      logical function argos_prep_mktop(lfnout,title,nparms,mparms,
     + lfnseq,filseq,lfntop,filtop,lfnsgm,lfnpar,
     + lfnmod,filmod,lfnmat,lfnhdb,lfnhop,filhop,
     + lfnbsg,filbsg,ignore,slvnam,itopol,icomb)
c
c $Id$
c
c     in : integer lfnout    = logical file number output file
c          char*80 ffield    = force field from [amber]
c          integer lfnseq    = logical file number for seq file
c          char*80 filseq    = file name of seq file
c          integer lfntop    = logical file number for top file
c          char*80 filtop    = file name of top file
c          integer lfnsgm    = logical file number for sgm file
c
      implicit none
c
#include "mafdecls.fh"
#include "util.fh"
#include "argos_prep_common.fh"
c
      logical argos_prep_seqsiz,argos_prep_rdseq,argos_prep_term,
     + argos_prep_dimens,argos_prep_mklist
      external argos_prep_seqsiz,argos_prep_rdseq,argos_prep_term,
     + argos_prep_dimens,argos_prep_mklist
      logical argos_prep_params,argos_prep_natyps,argos_prep_wrttop,
     + argos_prep_third,argos_prep_excl
      external argos_prep_params,argos_prep_natyps,argos_prep_wrttop,
     + argos_prep_third,argos_prep_excl
      logical argos_prep_nonbon,argos_prep_modify
      external argos_prep_nonbon,argos_prep_modify
c
      integer lfnout,lfnseq,lfntop,lfnsgm,lfnpar,lfnmod,ignore,lfnmat
      integer lfnhdb,lfnhop,lfnbsg
      integer nparms,mparms,itopol,icomb
      character*80 title(2,3)
      character*255 filseq,filtop,filmod,filhop,filbsg
      character*3 slvnam
c
      integer mseq,nseq
      integer l_lseq,i_lseq,l_cseq,i_cseq
c
      integer mlnk,nlnk
      integer l_llnk,i_llnk,l_clnk,i_clnk
c
      integer matm,natm
      integer l_latm,i_latm,l_catm,i_catm,l_qatm,i_qatm
c
      integer mbnd,nbnd
      integer l_lbnd,i_lbnd,l_rbnd,i_rbnd
c
      integer mang,nang
      integer l_lang,i_lang,l_rang,i_rang
c
      integer mdih,ndih
      integer l_ldih,i_ldih,l_kdih,i_kdih,l_rdih,i_rdih,l_ndih,i_ndih
c
      integer mimp,nimp
      integer l_limp,i_limp,l_rimp,i_rimp,l_kimp,i_kimp
c
      integer matt,natt,mats,nats
      integer l_latt,i_latt,l_catt,i_catt,l_patt,i_patt,l_ratt,i_ratt
      integer l_lats,i_lats
c
      integer nval
      integer l_ival,i_ival,l_rval,i_rval,l_ndx,i_ndx
c
      integer m3rd,n3rd
      integer i_l3rd,l_l3rd
c
      integer mexc,nexc
      integer i_lexc,l_lexc
c
      integer mnon,nnon
      integer i_lnon,l_lnon
c
      integer natmt,nbndt,nangt,ndiht,nimpt,n3rdt,mqu
c
      real*8 releps,q14fac,wcorr(10)
c
      logical lupdat
c
      real*8 timer_wall_total
      external timer_wall_total
c
      lupdat=.false.
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1000)
 1000 format('TOPOLOGY GENERATION')
      endif
c
      call timer_init()
      if(.not.argos_prep_seqsiz(lfnout,lfnseq,filseq,nseq,nlnk))
     + call md_abort('argos_prep_seqsize failed',9999)
c
      mseq=nseq+2
      mlnk=nlnk+1
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1001) nseq
 1001 format('sequence length is',i10)
      endif
c
c     allocate memory for sequence
c     ----------------------------
c
c     integer lseq(1,mseq) : number of sequence entry
c                  2       : link type
c                  3       : index of first atom of sequence entry
c                  4       : molecule number
c
      if(.not.ma_push_get(mt_int,6*mseq,'lseq',l_lseq,i_lseq))
     + call md_abort('Memory allocation failed for lseq',9999)
c
c     char*5  cseq(1,mseq) : name of sequence entry on pdb
c                  2       : name of sequence entry on top
c
      if(.not.ma_push_get(mt_byte,10*mseq,'cseq',l_cseq,i_cseq))
     + call md_abort('Memory allocation failed for cseq',9999)
c
c     allocate memory for link list
c     -----------------------------
c
c     integer llnk(1,mlnk) : link segment 1
c                  2       : link segment 2
c                  3       : link type : 0: normal
c                                        1: forced
c     char*4  clnk(1,mlnk) : link segment 1 atom name
c                  2       : link segment 2 atom name
c
      if(.not.ma_push_get(mt_int,3*mlnk,'llnk',l_llnk,i_llnk))
     + call md_abort('Memory allocation failed for llnk',9999)
      if(.not.ma_push_get(mt_byte,8*mlnk,'clnk',l_clnk,i_clnk))
     + call md_abort('Memory allocation failed for clnk',9999)
c
c     read the sequence file
c     ----------------------
c
      if(.not.argos_prep_rdseq(lfnout,lfnseq,filseq,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_llnk),byte_mb(i_clnk),mlnk,nlnk))
     + call md_abort('argos_prep_rdseq failed',9999)
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1002)
 1002 format('sequence is read')
      endif
c
c     determine termini
c     -----------------
c
      if(.not.argos_prep_term(lfnout,lfnsgm,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_llnk),byte_mb(i_clnk),mlnk,nlnk))
     + call md_abort('argos_prep_term failed',9999)
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1003)
 1003 format('termini are determined')
      endif
c
c     determine array dimensions
c     --------------------------
c
      if(.not.argos_prep_dimens(lfnout,lfnsgm,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + natm,nbnd,nang,ndih,nimp,nparms))
     + call md_abort('argos_prep_dimens failed',9999)
      mparms=nparms+1
c
c     adjust dimension to account for links
c     -------------------------------------
c
      matm=natm+1
      mbnd=nbnd+3*nseq+1
      mang=nang+27*nseq+1
      mdih=ndih+45*nseq+1
      mimp=nimp+6*nseq+1
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1004) mseq,matm,mbnd,mang,mdih,mimp
 1004 format('list dimensions are determined:',/,
     + ' number of segments  ',i10,/,
     + ' number of atoms     ',i10,/,
     + ' number of bonds     ',i10,/,
     + ' number of angles    ',i10,/,
     + ' number of dihedrals ',i10,/,
     + ' number of impropers ',i10)
      endif
c
c     allocate memory for atomic lists
c     --------------------------------
c
c     integer latm(1,matm)   : charge group
c                  2         : polarization group
c                  3         : link number
c                  4         : center type
c                  5         : segment number
c                  6         : molecule number
c                  7         : atomic number
c                  8         : fraction
c                  9         : ???
c
c
c     char*6  catm(  1,matm) : atom name
c                  i+1       : atom type set i, i=1,nparms
c
      if(.not.ma_push_get(mt_int,11*matm,'latm',l_latm,i_latm))
     + call md_abort('Memory allocation failed for latm',9999)
      if(.not.ma_push_get(mt_byte,6*mparms*matm,'catm',l_catm,i_catm))
     + call md_abort('Memory allocation failed for catm',9999)
c
c     real*8 qatm(nparms,1,matm) : partial atomic charges sets 1,..,nparms
c                        2       : polarizabilities sets 1,..,nparms
c
      mqu=3*matm
      if(.not.ma_push_get(mt_dbl,2*nparms*mqu,'qatm',l_qatm,i_qatm))
     + call md_abort('Memory allocation failed for qatm',9999)
c
c     allocate memory for bonded lists
c     --------------------------------
c
c     integer lbnd(1:2,mbnd)        : bond indices
c                  3                : parameter source
c                  4                : constraint type
c     real*8  rbnd(nparms,1,mbnd)   : bond length set 1
c                         2         : force constant set 1
c
      if(.not.ma_push_get(mt_int,4*mbnd,'lbnd',l_lbnd,i_lbnd))
     + call md_abort('Memory allocation failed for lbnd',9999)
      if(.not.ma_push_get(mt_dbl,2*nparms*mbnd,'rbnd',l_rbnd,i_rbnd))
     + call md_abort('Memory allocation failed for rbnd',9999)
c
c     integer lang(1:3,mang)        : angle indices
c                  4                : parameter source
c                  5                : constraint type
c     real*8  rang(nparms,1,mang)   : angle set 1
c                         2         : force constant set 1
c
      if(.not.ma_push_get(mt_int,5*mang,'lang',l_lang,i_lang))
     + call md_abort('Memory allocation failed for lang',9999)
      if(.not.ma_push_get(mt_dbl,4*nparms*mang,'rang',l_rang,i_rang))
     + call md_abort('Memory allocation failed for rang',9999)
c
c     integer ldih(1:4,mdih) : torsion indices
c                  5         : parameter source
c                  6         : constraint type
c
c     integer nfdih(nparms,mdih) : number of functions
c
c     integer kdih(6,nparms,mdih) : multiplicity set 1
c
c     real*8  rdih(6,nparms,1,mdih) : phase angle set 1
c                           2       : force constant set 1
c
      if(.not.ma_push_get(mt_int,6*mdih,'ldih',l_ldih,i_ldih))
     + call md_abort('Memory allocation failed for ldih',9999)
      if(.not.ma_push_get(mt_int,nparms*mdih,'nfdih',l_ndih,i_ndih))
     + call md_abort('Memory allocation failed for ldih',9999)
      if(.not.ma_push_get(mt_int,6*nparms*mdih,'kdih',l_kdih,i_kdih))
     + call md_abort('Memory allocation failed for ldih',9999)
      if(.not.ma_push_get(mt_dbl,12*nparms*mdih,'rdih',l_rdih,i_rdih))
     + call md_abort('Memory allocation failed for rdih',9999)
c
c     integer limp(1:4,mimp) : improper torsion indices
c                  5         : parameter source
c                  6         : constraint type
c
c     integer kimp(nparms,mimp) : multiplicity
c
c     real*8  rimp(nparms,1,mimp)   : phase angle set 1
c                         2         : force constant set 1
c
      if(.not.ma_push_get(mt_int,6*mimp,'limp',l_limp,i_limp))
     + call md_abort('Memory allocation failed for limp',9999)
      if(.not.ma_push_get(mt_int,nparms*mimp,'kimp',l_kimp,i_kimp))
     + call md_abort('Memory allocation failed for limp',9999)
      if(.not.ma_push_get(mt_dbl,2*nparms*mimp,'rimp',l_rimp,i_rimp))
     + call md_abort('Memory allocation failed for rimp',9999)
c
c     generate the atomic and bonded lists
c     ------------------------------------
c
      natm=0
      nbnd=0
      nang=0
      ndih=0
      nimp=0
c
      nval=max(mbnd,mang,mdih,mimp,1)
      if(.not.ma_push_get(mt_int,nval,'ival',l_ival,i_ival))
     + call md_abort('Memory allocation failed for ival',9999)
      if(.not.ma_push_get(mt_int,nval,'ndx',l_ndx,i_ndx))
     + call md_abort('Memory allocation failed for ndx',9999)
      if(.not.ma_push_get(mt_dbl,nval,'rval',l_rval,i_rval))
     + call md_abort('Memory allocation failed for rval',9999)
      if(.not.argos_prep_mklist(lfnout,lfnsgm,lfnmat,nparms,mparms,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_llnk),byte_mb(i_clnk),mlnk,nlnk,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_lbnd),dbl_mb(i_rbnd),mbnd,nbnd,
     + int_mb(i_lang),dbl_mb(i_rang),mang,nang,
     + int_mb(i_ldih),int_mb(i_ndih),int_mb(i_kdih),
     + dbl_mb(i_rdih),mdih,ndih,
     + int_mb(i_limp),int_mb(i_kimp),dbl_mb(i_rimp),mimp,nimp,
     + natmt,nbndt,nangt,ndiht,nimpt,wcorr,
     + nval,int_mb(i_ival),dbl_mb(i_rval),int_mb(i_ndx),itopol))
     + call md_abort('argos_prep_mklist failed',9999)
      if(.not.ma_pop_stack(l_rval))
     + call md_abort('Memory deallocation failed for rval',9999)
      if(.not.ma_pop_stack(l_ndx))
     + call md_abort('Memory deallocation failed for ndx',9999)
      if(.not.ma_pop_stack(l_ival))
     + call md_abort('Memory deallocation failed for ival',9999)
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1005)
 1005 format('lists are generated')
      endif
c
c     apply topology modifications to atom types
c     ------------------------------------------
      if(.not.argos_prep_modify(1,lfnout,lfnmod,filmod,nparms,mparms,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_lbnd),dbl_mb(i_rbnd),mbnd,nbnd,nbndt,
     + int_mb(i_lang),dbl_mb(i_rang),mang,nang,nangt,
     + int_mb(i_ldih),int_mb(i_ndih),int_mb(i_kdih),
     + dbl_mb(i_rdih),mdih,ndih,ndiht,
     + int_mb(i_limp),int_mb(i_kimp),dbl_mb(i_rimp),mimp,nimp,nimpt,
     + lupdat))
     + call md_abort('argos_prep_modify failed',9999)
c
c     get number of atom types
c     ------------------------
c
      if(.not.argos_prep_natyps(lfnout,nparms,mparms,
     + byte_mb(i_catm),matm,natmt,natt,nats))
     + call md_abort('argos_prep_natyps failed',9999)
      matt=natt
      mats=nats
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1006) natt
 1006 format(' Number of atom types is ',i10)
      endif
c
c     allocate memory for atom type lists
c     -----------------------------------
c
      if(.not.ma_push_get(mt_int,3*matt,'latt',l_latt,i_latt))
     + call md_abort('Memory allocation failed for latt',9999)
      if(.not.ma_push_get(mt_int,nparms*mats,'lats',l_lats,i_lats))
     + call md_abort('Memory allocation failed for lats',9999)
      if(.not.ma_push_get(mt_byte,36*matt,'catt',l_catt,i_catt))
     + call md_abort('Memory allocation failed for catt',9999)
      if(.not.ma_push_get(mt_dbl,12*matt*matt,'patt',l_patt,i_patt))
     + call md_abort('Memory allocation failed for patt',9999)
      if(.not.ma_push_get(mt_dbl,3*matt,'ratt',l_ratt,i_ratt))
     + call md_abort('Memory allocation failed for ratt',9999)
c
c     substitute force field parameters
c     ---------------------------------
c
      if(.not.argos_prep_params(lfnpar,lfnout,nparms,mparms,
     + releps,q14fac,ignore,
     + int_mb(i_latm),byte_mb(i_catm),matm,natmt,
     + int_mb(i_lbnd),dbl_mb(i_rbnd),mbnd,nbndt,
     + int_mb(i_lang),dbl_mb(i_rang),mang,nangt,
     + int_mb(i_ldih),int_mb(i_ndih),int_mb(i_kdih),
     + dbl_mb(i_rdih),mdih,ndiht,
     + int_mb(i_limp),int_mb(i_kimp),dbl_mb(i_rimp),mimp,nimpt,
     + int_mb(i_latt),int_mb(i_lats),byte_mb(i_catt),dbl_mb(i_patt),
     + dbl_mb(i_ratt),matt,natt,mats,nats,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,icomb))
     + call md_abort('argos_prep_params failed',9999)
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1007)
 1007 format('force field parameters are substituted')
      endif
c
c     allocate memory for third neighbor list
c     ---------------------------------------
c
      m3rd=6*mang
      n3rd=0
      if(.not.ma_push_get(mt_int,2*m3rd,'l3rd',l_l3rd,i_l3rd))
     + call md_abort('Memory allocation failed for l3rd',9999)
      nval=max(m3rd,matm,1)
      if(.not.ma_push_get(mt_int,nval,'ival',l_ival,i_ival))
     + call md_abort('Memory allocation failed for ival',9999)
      if(.not.ma_push_get(mt_int,nval,'ndx',l_ndx,i_ndx))
     + call md_abort('Memory allocation failed for ndx',9999)
c
c     construct third neighbor list
c     -----------------------------
c
      if(.not.argos_prep_third(int_mb(i_lbnd),mbnd,nbnd,int_mb(i_lang),
     + mang,nang,int_mb(i_l3rd),m3rd,n3rd,1,1,1,
     + nval,int_mb(i_ival),int_mb(i_ndx),matm))
     + call md_abort('argos_prep_third failed',9999)
c
      n3rdt=n3rd
      if(nbndt.gt.nbnd.and.nangt.gt.nang) then
      if(.not.argos_prep_third(int_mb(i_lbnd),mbnd,nbndt,int_mb(i_lang),
     + mang,nangt,int_mb(i_l3rd),m3rd,n3rdt,nbnd+1,nang+1,n3rd+1,
     + nval,int_mb(i_ival),int_mb(i_ndx),matm))
     + call md_abort('argos_prep_third failed',9999)
      endif
      if(.not.ma_pop_stack(l_ndx))
     + call md_abort('Memory deallocation failed for ndx',9999)
      if(.not.ma_pop_stack(l_ival))
     + call md_abort('Memory deallocation failed for ival',9999)
c
      if(util_print('topology',print_high)) then
      write(lfnout,1008) n3rd
 1008 format(' Solute third neighbor list length is',i10)
      endif
c
c     allocate memory for excluded pair list
c     --------------------------------------
c
      mexc=n3rd+nbnd+nang+300
      nexc=0
      if(.not.ma_push_get(mt_int,2*mexc,'lexc',l_lexc,i_lexc))
     + call md_abort('Memory allocation failed for lexc',9999)
c
c     construct excluded pair list
c     ----------------------------
c
      nval=max(mexc,1)
      if(.not.ma_push_get(mt_int,nval,'ival',l_ival,i_ival))
     + call md_abort('Memory allocation failed for ival',9999)
      if(.not.ma_push_get(mt_int,nval,'ndx',l_ndx,i_ndx))
     + call md_abort('Memory allocation failed for ndx',9999)
      if(.not.argos_prep_excl(int_mb(i_l3rd),m3rd,n3rd,int_mb(i_lbnd),
     + mbnd,nbnd,int_mb(i_lang),mang,nang,int_mb(i_lexc),mexc,nexc,
     + nval,int_mb(i_ival),int_mb(i_ndx),byte_mb(i_catt),matt,
     + int_mb(i_lats),nparms,mats,int_mb(i_latm),matm,natm,lupdat))
     + call md_abort('argos_prep_excl failed',9999)
      if(.not.ma_pop_stack(l_ndx))
     + call md_abort('Memory deallocation failed for ndx',9999)
      if(.not.ma_pop_stack(l_ival))
     + call md_abort('Memory deallocation failed for ival',9999)
c
      if(util_print('topology',print_high)) then
      write(lfnout,1009) nexc
 1009 format(' Solute excluded pair list length is',i10)
      endif
c
c     allocate memory for solvent non-bonded list
c     -------------------------------------------
c
      mnon=max(1,(natmt-natm)*(natmt-natm))
      nnon=0
      if(.not.ma_push_get(mt_int,2*mnon,'lnon',l_lnon,i_lnon))
     + call md_abort('Memory allocation failed for lnon',9999)
c
c     construct solvent non-bonded list
c     ---------------------------------
c
      if(.not.argos_prep_nonbon(natm+1,natmt,int_mb(i_lbnd),
     + mbnd,nbnd+1,nbndt,
     + int_mb(i_lang),mang,nang+1,nangt,
     + int_mb(i_l3rd),m3rd,n3rd+1,n3rdt,
     + int_mb(i_lnon),mnon,nnon))
     + call md_abort('argos_prep_nonbon failed',9999)
c
      if(util_print('topology',print_high)) then
      write(lfnout,1010) nnon
 1010 format(' Solvent non-bonded list length is',i10)
      endif
c
c     apply topology modifications to bonded parameters
c     -------------------------------------------------
c
      if(.not.argos_prep_modify(2,lfnout,lfnmod,filmod,nparms,mparms,
     + int_mb(i_latm),byte_mb(i_catm),dbl_mb(i_qatm),matm,natm,
     + int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,
     + int_mb(i_lbnd),dbl_mb(i_rbnd),mbnd,nbnd,nbndt,
     + int_mb(i_lang),dbl_mb(i_rang),mang,nang,nangt,
     + int_mb(i_ldih),int_mb(i_ndih),int_mb(i_kdih),
     + dbl_mb(i_rdih),mdih,ndih,ndiht,
     + int_mb(i_limp),int_mb(i_kimp),dbl_mb(i_rimp),mimp,nimp,nimpt,
     + lupdat))
     + call md_abort('argos_prep_modify failed',9999)
c
c     write the topology file
c     -----------------------
c
      if(.not.argos_prep_wrttop(lfnout,title,lfntop,filtop,lfnhdb,
     + lfnhop,filhop,lfnbsg,filbsg,releps,q14fac,
     + nparms,mparms,int_mb(i_lseq),byte_mb(i_cseq),mseq,nseq,mqu,
     + int_mb(i_latt),int_mb(i_lats),byte_mb(i_catt),dbl_mb(i_patt),
     + dbl_mb(i_ratt),matt,natt,mats,nats,int_mb(i_latm),
     + byte_mb(i_catm),
     + dbl_mb(i_qatm),matm,natm,int_mb(i_lbnd),dbl_mb(i_rbnd),mbnd,
     + nbnd,int_mb(i_lang),
     + dbl_mb(i_rang),mang,nang,int_mb(i_ldih),int_mb(i_ndih),
     + int_mb(i_kdih),dbl_mb(i_rdih),mdih,ndih,
     + int_mb(i_limp),int_mb(i_kimp),dbl_mb(i_rimp),mimp,nimp,
     + int_mb(i_l3rd),
     + m3rd,n3rd,int_mb(i_lexc),mexc,nexc,int_mb(i_lnon),mnon,nnon,
     + natmt,nbndt,nangt,ndiht,nimpt,n3rdt,wcorr,slvnam,itopol))
     + call md_abort('argos_prep_wrttop failed',9999)
c
      if(util_print('topology',print_debug)) then
      write(lfnout,1011)
 1011 format(' Topology file is written')
      endif
c
c     deallocate memory
c     -----------------
c
      if(.not.ma_pop_stack(l_lnon))
     + call md_abort('Memory deallocation failed for lnon',9999)
      if(.not.ma_pop_stack(l_lexc))
     + call md_abort('Memory deallocation failed for lexc',9999)
      if(.not.ma_pop_stack(l_l3rd))
     + call md_abort('Memory deallocation failed for l3rd',9999)
      if(.not.ma_pop_stack(l_ratt))
     + call md_abort('Memory deallocation failed for ratt',9999)
      if(.not.ma_pop_stack(l_patt))
     + call md_abort('Memory deallocation failed for patt',9999)
      if(.not.ma_pop_stack(l_catt))
     + call md_abort('Memory deallocation failed for catt',9999)
      if(.not.ma_pop_stack(l_lats))
     + call md_abort('Memory deallocation failed for lats',9999)
      if(.not.ma_pop_stack(l_latt))
     + call md_abort('Memory deallocation failed for latt',9999)
      if(.not.ma_pop_stack(l_rimp))
     + call md_abort('Memory deallocation failed for rimp',9999)
      if(.not.ma_pop_stack(l_kimp))
     + call md_abort('Memory deallocation failed for limp',9999)
      if(.not.ma_pop_stack(l_limp))
     + call md_abort('Memory deallocation failed for limp',9999)
      if(.not.ma_pop_stack(l_rdih))
     + call md_abort('Memory deallocation failed for rdih',9999)
      if(.not.ma_pop_stack(l_kdih))
     + call md_abort('Memory deallocation failed for kdih',9999)
      if(.not.ma_pop_stack(l_ndih))
     + call md_abort('Memory deallocation failed for ldih',9999)
      if(.not.ma_pop_stack(l_ldih))
     + call md_abort('Memory deallocation failed for ldih',9999)
      if(.not.ma_pop_stack(l_rang))
     + call md_abort('Memory deallocation failed for rang',9999)
      if(.not.ma_pop_stack(l_lang))
     + call md_abort('Memory deallocation failed for lang',9999)
      if(.not.ma_pop_stack(l_rbnd))
     + call md_abort('Memory deallocation failed for rbnd',9999)
      if(.not.ma_pop_stack(l_lbnd))
     + call md_abort('Memory deallocation failed for lbnd',9999)
      if(.not.ma_pop_stack(l_qatm))
     + call md_abort('Memory deallocation failed for qatm',9999)
      if(.not.ma_pop_stack(l_catm))
     + call md_abort('Memory deallocation failed for catm',9999)
      if(.not.ma_pop_stack(l_latm))
     + call md_abort('Memory deallocation failed for latm',9999)
      if(.not.ma_pop_stack(l_clnk))
     + call md_abort('Memory deallocation failed for clnk',9999)
      if(.not.ma_pop_stack(l_llnk))
     + call md_abort('Memory deallocation failed for llnk',9999)
      if(.not.ma_pop_stack(l_cseq))
     + call md_abort('Memory deallocation failed for cseq',9999)
      if(.not.ma_pop_stack(l_lseq))
     + call md_abort('Memory deallocation failed for lseq',9999)
c
      argos_prep_mktop=.true.
      return
      end

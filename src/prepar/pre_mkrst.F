      logical function pre_mkrst(source,lfnout,lprint,lfntop,filtop,
     + lfnpdb,filpdb,lfnrst,filrst,lfntmp,filtmp,lfncmd,filcmd)
c
      implicit none
c
#include "mafdecls.fh"
c
      logical pre_topsiz,pre_toprd,pre_rdpdb
      external pre_topsiz,pre_toprd,pre_rdpdb
c
      integer lfnout,lprint,lfntop,lfnpdb,lfnrst,lfntmp,lfncmd
      character*255 filtop,filpdb,filrst,filtmp,filcmd,filnam
      character*80 source,card
c
      integer nsa,nwa,nwm,nmis
      integer msa,mwa
c
      integer l_cwa,i_cwa
      integer l_csa,i_csa,l_sat,i_sat,l_sgm,i_sgm,l_qsa,i_qsa
      integer l_sfnd,i_sfnd,l_xs,i_xs
c
      logical lslvnt,ldistr,lnames
      integer length,len
c
      if(lprint.ge.5) write(lfnout,1000)
 1000 format('RESTART FILE GENERATION')
c
      if(.not.pre_topsiz(lfntop,filtop,lfnout,lprint,msa,mwa))
     + call errquit('pre_topsiz failed',9999)
c
      nwm=0
      nwa=0
      nsa=0
c
c     allocate memory
c     ---------------
c
c     character*16 cwa(mwa)  : solvent atom names
c     character*16 csa(msa)  : solute atom names
c     integer isat(msa)      : solute atom types
c     integer isgm(msa)      : solute segment numbers
c     integer isfnd(msa)     : solute atom found flags
c     real*8 qsa(msa)        : solute atom charges
c     real*8 xs(3,msa)       : solute atom coordinates
c
      if(.not.ma_push_get(mt_byte,16*mwa,"cwa",l_cwa,i_cwa))
     + call errquit('Memory allocation failed for cwa',9999)
      if(.not.ma_push_get(mt_byte,16*msa,"csa",l_csa,i_csa))
     + call errquit('Memory allocation failed for csa',9999)
      if(.not.ma_push_get(mt_int,msa,"sat",l_sat,i_sat))
     + call errquit('Memory allocation failed for sat',9999)
      if(.not.ma_push_get(mt_int,msa,"sgm",l_sgm,i_sgm))
     + call errquit('Memory allocation failed for sgm',9999)
      if(.not.ma_push_get(mt_int,msa,"sfnd",l_sfnd,i_sfnd))
     + call errquit('Memory allocation failed for sfnd',9999)
      if(.not.ma_push_get(mt_dbl,msa,"qsa",l_qsa,i_qsa))
     + call errquit('Memory allocation failed for qsa',9999)
      if(.not.ma_push_get(mt_dbl,3*msa,"xs",l_xs,i_xs))
     + call errquit('Memory allocation failed for xs',9999)
c
      lslvnt=.false.
      ldistr=.false.
c
c     read topology file
c     ------------------
c
      if(.not.pre_toprd(lfntop,filtop,lfnout,lprint,
     + byte_mb(i_cwa),mwa,nwa,int_mb(i_sat),int_mb(i_sgm),
     + byte_mb(i_csa),dbl_mb(i_qsa),msa,nsa))
     + call errquit('pre_toprd failed',9999)
      if(lprint.ge.5) write(lfnout,2000) filtop(1:index(filtop,' ')-1)
 2000 format('Topology successfully read from ',a)
c
c     open the command file
c     ---------------------
c
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old',err=99)
      goto 1
   99 continue
      if(.not.pre_rdpdb(filpdb,lfnpdb,filtmp,lfntmp,.true.,
     + byte_mb(i_csa),int_mb(i_sgm),int_mb(i_sfnd),dbl_mb(i_xs),
     + msa,nsa,nmis))
     + call errquit('pre_rdpdb failed',9999)
      if(lprint.ge.5) write(lfnout,2003) filpdb(1:index(filpdb,' ')-1)
      write(lfnout,2103) nmis,nsa
      goto 9
c
    1 continue
c
      read(lfncmd,1001,end=9,err=9999) card
 1001 format(a)
c
c     read distance restraints
c     ------------------------
c
      if(card(1:6).eq.'disres') then
      endif
c
c     read solvent coordinates
c     ------------------------
c
      if(card(1:6).eq.'solvnt') then
      endif
c
c     read pdb file
c     -------------
c
      if(card(1:6).eq.'rd_pdb') then
      filnam=card(8:80)
      length=index(filnam,' ')-1
      if(length.eq.0) filnam=filpdb
      len=index(filpdb,' ')-1
      lnames=length.eq.len.and.filnam(1:len).eq.filpdb(1:len)
      if(.not.pre_rdpdb(filnam,lfnpdb,filtmp,lfntmp,lnames,
     + byte_mb(i_csa),int_mb(i_sgm),int_mb(i_sfnd),dbl_mb(i_xs),
     + msa,nsa,nmis))
     + call errquit('pre_rdpdb failed',9999)
      if(lprint.ge.5) write(lfnout,2003) filnam(1:index(filnam,' ')-1)
 2003 format('PDB file successfully read from ',a)
      write(lfnout,2103) nmis,nsa
 2103 format('Missing',i5,' atoms out of',i5)
      endif
c
c     read restart file
c     ------------------
c
      if(card(1:6).eq.'rd_rst') then
      endif
c
c     write pdb file
c     --------------
c
      if(card(1:6).eq.'wr_pdb') then
      endif
c
c     write restart file
c     ------------------
c
      if(card(1:6).eq.'wr_rst') then
      endif
c
c     center solute
c     -------------
c
      if(card(1:6).eq.'center') then
      endif
c
c     orient solute
c     -------------
c
      if(card(1:6).eq.'orient') then
      endif
c
c     solvate solute
c     --------------
c
      if(card(1:6).eq.'solvat') then
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
      endif
c
c     set minimum box size
c     --------------------
c
      if(card(1:6).eq.'boxsiz') then
      endif
c
c     set periodic boundary conditions
c     --------------------------------
c
      if(card(1:6).eq.'period') then
      endif
c
c     set vacuum conditions
c     ---------------------
c
      if(card(1:6).eq.'vacuum') then
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
c     end
c     ---
c
      if(card(1:3).ne.'end') goto 1
c
c      if(.not.pre_wrtrst(lfnrst
c
    9 continue
      close(unit=lfncmd,status='delete')
  999 continue
c
c     deallocate memory
c     -----------------
c
      if(.not.ma_pop_stack(l_xs))
     + call errquit('Memory deallocation failed for xs',9999)
      if(.not.ma_pop_stack(l_qsa))
     + call errquit('Memory deallocation failed for qsa',9999)
      if(.not.ma_pop_stack(l_sfnd))
     + call errquit('Memory deallocation failed for sfnd',9999)
      if(.not.ma_pop_stack(l_sgm))
     + call errquit('Memory deallocation failed for sgm',9999)
      if(.not.ma_pop_stack(l_sat))
     + call errquit('Memory deallocation failed for sat',9999)
      if(.not.ma_pop_stack(l_csa))
     + call errquit('Memory deallocation failed for csa',9999)
      if(.not.ma_pop_stack(l_cwa))
     + call errquit('Memory deallocation failed for cwa',9999)
c
      pre_mkrst=.true.
      return
c
 9999 continue
      pre_mkrst=.false.
      return
      end



      subroutine argos_rdrest(lfn,fil)
c
      implicit none
c
#include "argos.fh"
c
      integer lfn
      character*255 fil
c
      character*13 string
c
      if(me.eq.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfn)
c
    1 continue
      read(lfn,1000,end=9997) string
 1000 format(a13)
      if(string.ne.'restart input') goto 1
      read(lfn,1001,end=9998,err=9998) ntype,mdtype
 1001 format(11i7)
      read(lfn,1001) nfpair,nflong
      read(lfn,1001) lwtype,lstype,nfrest,keepr
      if(nserie.eq.1) then
      read(lfn,1002) krun,kequi,kdacq,mrun,mequi,mdacq,ldacq
 1002 format(7i7)
      else
      read(lfn,1002) krun,kequi,kdacq
      endif
      read(lfn,1003) stime,tstep
 1003 format(2f12.6)
      read(lfn,1004) rshort,rlong
 1004 format(2f12.6)
      read(lfn,1005) mshitw,tlwsha
      read(lfn,1005) mshits,tlssha
 1005 format(i7,f12.6)
      read(lfn,1006) ipscal,prsext,prsrlx,compr,ipopt
 1006 format(i5,e12.5,f12.6,e12.5,i5)
      read(lfn,1007) itscal,tmpext1,tmprlx,tmsrlx,tmpext2,
     + tann1,tann2
 1007 format(i5,6f12.6)
      read(lfn,1008) nfgaus,ivopt,tgauss,iseed
 1008 format(2i7,f12.6,i12)
      read(lfn,1009) nfoutp,nfstat,nfprop,nfnewf,ibatch
 1009 format(11i7)
      read(lfn,1009) ibinar,iformt
      read(lfn,1009) nfcoor,nfscoo,nfvelo,nfsvel,nfforc,nfsfor
      read(lfn,1009) nffree
      read(lfn,1009) nfcntr,nfslow
      read(lfn,1010) nfrdf,numrdf,ngc,ngr,ngl,ngrww,ngrsw,ngrss
 1010 format(8i7)
      read(lfn,1011) rrdf,drdf
 1011 format(2f12.6)
      read(lfn,1012) numdis,lendis
 1012 format(2i7)
      read(lfn,1013) numhis,lenhis
 1013 format(11i7)
      read(lfn,1014) nfdip,ndip,rdip
 1014 format(2i7,f12.6)
      read(lfn,1015) nfkirk,nkirk,rkirk
 1015 format(2i7,f12.6)
      endif
c
 9997 continue
      if(me.eq.0) close(unit=lfn)
      return
c
 9998 continue
      call md_abort('Unable to read restart file in argos_rest ',me)
      return
 9999 continue
      call md_abort('Unable to open restart file in argos_rest ',me)
      return
      end
c $Id$

      subroutine qhop_setup(mparms,nparms,matm,natm,mseq,nseq,
     + catm,latm,cseq,lseq)
c
c $Id$
c
      implicit none
c
c     input variables
c     ---------------
c
c     matm            : dimension of atom arrays
c     natm            : number of atoms
c     mseq            : dimension of residue arrays
c     nseq            : number of residues
c     nparms          : number of parameter sets
c     mparms          : number of parameter sets plus 1
c     cseq            : residue names
c     catm(1,1:natm)  : atom names
c     catm(k,1:natm)  : atom types in parameter set k-1
c     latm(5,1:natm)  : residue number
c     latm(10,1:natm) : + heavy atom attached to protonatable hydrogen
c                       - protonatable hydrogen, relative index to heavy atom
c     lseq(5,1:nseq)  : number of protonation states of residue
c     lseq(6,1:nseq)  : current protonation state of residue
c
      integer matm,natm,mseq,nseq,nparms,mparms
      character*10 cseq(mseq)
      character*6 catm(mparms,matm)
      integer latm(10,matm),lseq(6,mseq)
c
      integer i,j
c
      do 1 i=1,natm
c
      write(*,1000) i,catm(1,i),latm(5,i),cseq(latm(5,i)),latm(10,i),
     + lseq(5,latm(5,i)),lseq(6,latm(5,i))
c
 1000 format(i5,1x,a6,1x,i5,1xa10,3i5)
c
    1 continue
c
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c Program gensym:
c
c   Given an input set of generators for any group in 3 dimensional
c  space or less, this program will produce the matrix representations
c  of the operators in the group. The input operators must be chosen in 
c  a particular order and according to the conventions set forth
c  on page 729 of the International Tables of Crystallography vol. 3.
c
c Algorithm:
c    The algorithm works in the following manner: given a set of
c   generators, {G}n choose and initial generator, G(i). Multiply
c   this generator by some operator from the set {O}m. In the first
c   pass of the seqence O(j)=G(i) therefore O(k)=G(i)*O(j)=G(i)*G(i).
c   Which is added to the current symmetry operation list (the matrix
c   symops) if G(i)^2.ne.E. A new generator is then selected from the 
c   generator list and added to the bottom of symops. All previous 
c   operators in the list are then pre-multiplied by this new generator.
c    For example, if the final generator is G3, then sympos may contain:
c
c   G1, G1^2, G2, G2G1, G2G1^2, G2^2, G2^2G1, G2^2G1^2, G3, G3G1, 
c   G3G1^2, G3G2, G3G2G1, G3G2G1^2, G3G2^2, G3G2^2G1, G3G2^2G1^2,
c   G3^2, G3^2G1, G3^2G1^2, G3^2G2, G3^2G2G1, G3^2G2G1^2, G3^2G2^2,
c   G3^2G2^2G1 and G3^2G2^2G1^2.
c
c    Note: Identities and repeated operators generated along the way
c   are trapped and not added to the list, so the actual number of 
c   operators in symops may be considerably less than shown above.
c
c Matrix Mults:
c    The appropriate matrix multiplications are defined in Seitz 
c   notation as:
c
c            {R|t}{S|u}r={R|t}(Sr+u)=RSr+Ru+t={RS}|Ru+t}r
c     
c   R and S are normal point operations t and u are the associated 
c   translation components for whatever dimensional system you are
c   working with. Expainations of this notation can be found in "Space
c   Groups for Solid State Scientists", G. Burns and A.M. Glazer.
c
c
c    The matrix SYMOPS contains the matrix reps. of all group operators
c   except the identity. The variable NOPS holds the number of operators
c   in SYMOPS.
c
c   
c
c
c                                        A.C. Hess 
c                                        D.G. Clerc
c                                        Solid State Theory Group
c                                        MSRC/PNL
c                                        9/13/93
c***********************************************************************
      subroutine gensym(itype,numgrp,numset,symops,nops,oprint,
     $     group_name,geom,rtdb)
C$Id$
      implicit real*8 (a-h,o-z) 

      parameter(maxops=192,tol=1.0d-07,max_gen=6)
      character*2 kpos(-1:3),kneg(-3:1),rotoop(maxops)
      dimension capr(3,4),caps(3,4),symops(maxops*3,4)
*      integer indx(3)
      double precision deter3
      external deter3
      dimension resop(3,4),gens(18,4),detres(3,3),cntvec(3,3)
      character*(*) group_name
      integer geom,rtdb
      logical oprint,use_primitive
      double precision s_vec(max_gen,3)
      data kpos/' 2',' 3',' 4',' 6',' 1'/,kneg/'-1','-6','-4','-3',' m'/

      logical  geom_use_primitive
      external geom_use_primitive

c
c-->call spgen with correct system type flag to make generators
c
         call spgen(itype,numgrp,numset,gens,cntvec,ngen,numvec,
     $              group_name,max_gen,s_vec,oprint)

c
c-----------------------------------------------------------------------
c
c--> outer loop picks up input generators (begins at 5000) inner implied
c    loop mults. current generator by all previous operators generated
c    and stored in symops (begins at 6000).
c
c--> pointers:
c             ipos indexes the current working generator in symops
c             igpos points at the next element to be used in the 
c                  input generator list
c             jpos points to the current operator in the list symops
c             nops holds the current number of operators in symops
c
c----------------------------------------------------------------------
      igen=1
      icnt1=1
      isquare=0
      nops=0
      if (oprint) then
         write(*,27)
         write(*,29) ngen
 27      format(/,16x,'---------------',' GROUP GENERATORS ','---------'
     $        ,'------')
 28      format(/,23x,'GROUP NUMBER AND NAME: ',a12)
 29      format(/,22x,i1,' GENERATORS USED TO FORM THE GROUP')
      endif
5000  ipos=(icnt1-1)*3+1
      igpos=(igen-1)*3+1
c
c--> pickup an input generator and store it in the symop list
c--> set first operator equal to this generator, G(i).
c
      do 200 i=1,3
         do 220 j=1,4
            symops(ipos+(i-1),j)=gens(igpos+(i-1),j)
            capr(i,j)=gens(igpos+(i-1),j)
220      continue
200   continue
      nops=nops+1
c
c--> compute trace and determinant of the generator
c
      do 560 i=1,3
         do 570 j=1,3
            detres(i,j)=capr(i,j)
570      continue
560   continue
      trace=0.0d+00
      do 580 i=1,3
         trace=trace+capr(i,i)
580   continue
*
      det = deter3(detres)
*      call ludcmp(detres,3,3,indx,det)
*      do 590 i=1,3
*         det=det*detres(i,i)
*590   continue
c
c--> store type information
c
c**********************************************************************
c  Screen out all point groups from the next calc   DGC 3/10/94
c**********************************************************************
      if(itype.ge.1.and.itype.le.3) then
         itrace=idint(trace)
c
         if(det.lt.0.0d+00) then
            rotoop(nops)=kneg(itrace)
            if (oprint) write(*,30) kneg(itrace),(s_vec(igen,j),j=1,3)
         else
            rotoop(nops)=kpos(itrace)
            if (oprint) write(*,30) kpos(itrace),(s_vec(igen,j),j=1,3)
         endif
c
      elseif(itype.eq.0) then
         if (oprint) write(*,31) 'PT'
      endif
c
30    format(/,8x,a2,' fold Rotoinversion Operator at (',2(f9.6,','),f9.
     &6,')')
31    format(/,25x,a2,' fold Rotoinversion Operator')
c
c--> write out input generators
c
      if (oprint) call mprint(capr,3,4)
c
c--> get the current operator to do mults with, O(j)
c
5999  icnt2=1
6000  jpos=(icnt2-1)*3+1      
      do 320 i=1,3
         do 330 j=1,4
            caps(i,j)=symops(jpos+(i-1),j)
330      continue
320   continue
c
c--> determine new rotational portion of operator [resop=capr*caps]
c
      do 130 i=1,3
         do 140 j=1,3
            sum=0.0d+00
            do 150 k=1,3
               sum = sum + capr(i,k)*caps(k,j)
150         continue
            resop(i,j)=sum
140      continue
130   continue
c
c--> compute new translational components [Ru+t]
c
      do 250 i=1,3
         sum=0.0d+00
         do 260 j=1,3
            sum=sum+capr(i,j)*caps(j,4)
260      continue
         resop(i,4)=sum
250   continue
      do 270 i=1,3
         resop(i,4)=resop(i,4)+capr(i,4)
270   continue
c
c--> clean up translational components (put in standard convention)
c
      do 390 i=1,3
         if(resop(i,4).ge.-tol.and.resop(i,4).le.tol) then
            resop(i,4)=0.0d+00
         elseif(resop(i,4).ge.1.0d+00-tol.and.resop(i,4).le.
     &          1.0d+00+tol) then
            resop(i,4)=1.0d+00
         endif
390   continue
c
c
      do 400 i=1,3
         if(resop(i,4).lt.0.0d+00) then
            resop(i,4)=resop(i,4)+1.00000000
         elseif(resop(i,4).ge.1.0d+00) then
            resop(i,4)=resop(i,4)-1.00000000
         endif
400   continue
c
c---> Set translation components 1/3 and 2/3 of the generator
c      products to double precision.
c      Otherwise get 1/3=0.33333334 and 2/3=0.66666666.
c      (A similar loop for the generators appears in "spgen.f")
c
      do 247 i=1,3
         if(resop(i,4).le..6667.and.resop(i,4).ge..6666) then
             resop(i,4)=2.0d+00/3.0d+00
         elseif(resop(i,4).le..3334.and.resop(i,4).ge..3333) then
             resop(i,4)=1.0d+00/3.0d+00
         endif
247   continue
c
c--> determine what the newly generated operator is
c
c**********************************************************************
c     Move this after screens? (rotoop(j) is not changed here!)
c**********************************************************************
      do 160 i=1,3
         do 170 j=1,3
            detres(i,j)=resop(i,j)
170      continue
160   continue
c
c--> compute trace and determinant of the new operator
c
      trace=0.0d+00
      do 180 i=1,3
         trace=trace+resop(i,i)
180   continue
*
      det = deter3(detres)
**      call ludcmp(detres,3,3,indx,det)
**      do 190 i=1,3
**         det=det*detres(i,i)
**190   continue
c
c-->  Check whether {RS|Ru+t} is identical to any of the previous
c       operators stored in the "symops" matrix.
c     The computation of Gi^2 (i>1) produces identical operators
c       in groups #75-#230 (e.g. group #75(P4)).  Exclusion of Gi^2
c       (i>1) from the algorithm fixes the problem for groups #75-#194,
c       but it results in the loss of many operators from the cubic
c       groups (e.g. group #195(P23)).
c
      irow=0
      ichecks=1
191   isame=0
      do 192 i=1,3
          do 193 j=1,4
              if(resop(i,j).ge.symops(irow+i,j)-tol.and.resop(i,j).le.sy
     &mops(irow+i,j)+tol) then
                isame=isame+1
              endif
193       continue
192    continue
      if(isame.eq.12.and.jpos.lt.ipos) then
          icnt2=icnt2+1
          goto 6000
      elseif(isame.eq.12.and.jpos.eq.ipos.and.igen.lt.ngen) then
          icnt1=nops+1
          igen=igen+1
          goto 5000
      elseif(isame.eq.12.and.jpos.eq.ipos.and.igen.eq.ngen) then
          goto 223
      elseif(ichecks.lt.nops) then
          irow=irow+3
          ichecks=ichecks+1
          goto 191
      endif
c
c--> add the new operator to the bottom of the symop list if it is not
c    the identity element, increment pointer containing the number of 
c    operations currently in the symops list
c
c**********************************************************************
c  Add the next screen since "idint" merely extracts the integer part
c  of a real number. (E.g. idint(2.999999998912563)=2).
c  Note: The addition is unnecessary for all 230 space groups and
c        point groups 1-13. Point group #14 (D5d), where C5^5=E is
c        computed, is the first case where "idint" causes problems.
c                           DGC 2-25-94
c**********************************************************************
      if(trace.ge.3.00d+00-tol.and.trace.le.3.00d+00+tol) then
        trace=3.00d+00
      endif
      itrace=idint(trace)
      if(itrace.ne.3) then
         isym=nops*3
         do 350 i=1,3
            do 360 j=1,4
               symops(isym+i,j)=resop(i,j)
360         continue
350      continue
         nops=nops+1
c
c--> assign labels to the rotational portion of the generated
c    symmetry operation.
c
c**********************************************************************
c  Screen out all point groups from the next calc   DGC 3/10/94
c**********************************************************************
c
         if(det.lt.0.0d+00.and.itype.ne.0) then
            rotoop(nops)=kneg(itrace)
         elseif(itype.ne.0) then
            rotoop(nops)=kpos(itrace)
         endif
      endif
c
c --> debug prints (print all matrices and new mats determinants)
c
c      write(*,12) 'trace = ', trace
c      write(*,12) 'det   = ', det
c      call mprint(resop,3,4)
c----------------------------------------------------------------------
c--> get ready to do next pass
c
c    If block: 
c       If ipos.lt.jpos then there are more operators in the list symops
c      that need to be multiplied by the current generator.
c
c--->  If an identity was generated before exhausting all possible 
c     multipliers in symops, then get a new multiplier from symops.
c      If an identity was generated after exhausting all possible
c     multipliers in symops, then get a new generator if there is
c     still one present in the generator list, and quit the algorithm
c     if there aren't any generators left.
c----------------------------------------------------------------------
      if(itrace.eq.3.and.jpos.lt.ipos) then
         icnt2=icnt2+1
         goto 6000
      elseif(itrace.eq.3.and.jpos.eq.ipos.and.igen.lt.ngen) then
         icnt1=nops+1
         igen=igen+1
         goto 5000
      elseif(itrace.eq.3.and.jpos.eq.ipos.and.igen.eq.ngen) then
         goto 223
      endif
c
c---------------------------------------------------------------------
c      A new operator has been generated at this point, and if the
c     list of possible multipliers in the list symops has not been
c     exhausted, then get a new multiplier. If the list has been 
c     exhausted (jpos.eq.ipos), then get a new generator unless the
c     space group has only one generator (ngen.eq.1).
c      To compute terms such as G3^2G1, G3^2G1^1,...G3^2G2^2G1^2,
c     the sequence jpos = 1,4,7,...,ipos of multipliers from symops
c     is repeated - but this time they are left-multiplied by G3^2 
c     instead of G3. The indicator "isquare" = 0 or 1 when the left-
c     multiplier is G3 or G3^2, rspt. At the end of the sequence 
c     having G1 as the left-multiplier, the matrix R (capr) is set
c     equal to Gi^2.
c---------------------------------------------------------------------
      if(jpos.lt.ipos) then
         icnt2=icnt2+1
         goto 6000
      elseif(igen.eq.1.and.igen.lt.ngen) then
         igen=igen+1
         icnt1=nops+1
         goto 5000
      elseif(ngen.eq.1) then
         goto 223
      elseif(igen.eq.ngen.and.isquare.eq.1) then
         goto 223
      elseif(isquare.eq.1) then
         isquare=0
         igen=igen+1
         icnt1=nops+1
         goto 5000
      endif
      do 201 i=1,3
         do 221 j=1,4
            capr(i,j)=symops(isym+i,j)
221      continue
201   continue
      isquare=1
      goto 5999
c
12    format(a8,f10.6)
c
c--> Determine if this is a centered space group, if so add the centering
c    vectors to those generated for the point (0,0,0), ie to the current list
c    of symops
c
c--> fill in roto-inversion operators (duplicate from top part of symops
c    to the bottom) then add centering vectors to the 4th col. of symops
c
c
c

223   if(numvec.gt.0) then
c------> fill in rotoinversion ops
         icnt=nops*3
         do 830 ij=1,numvec
c------> Fill in the identity first
            do 827 kiki=1,3
               do 828 ikik=1,3
                   symops(icnt+kiki,ikik)=0
828            continue
827         continue
            do 829 kiki=1,3
               symops(icnt+kiki,kiki)=1
829         continue
c---> label to rotop for identity
            rotoop(nops*ij+ij)=' 1'
            icnt=icnt+3
            do 832 ktop=1,nops*3
               icnt=icnt+1
               do 834 kcol=1,3
                  symops(icnt,kcol)=symops(ktop,kcol)
834            continue
832         continue
830      continue
c------> add centering vectors to col. 4            
         nops1=nops
         do 849 i=1,numvec
           isym=nops*3
c------> Add cntvec to the identity first
             ibot=isym+3*(i-1)
             do 853 j=1,3
                symops(ibot+j,4)=cntvec(i,j)
853          continue
            do 850 iold=1,nops1
               itop=(iold-1)*3
               ibot=isym+itop +3*i
               do 860 j=1,3
                  symops(ibot+j,4)=symops(itop+j,4)+cntvec(i,j)
860            continue
               rotoop(nops+i+iold)=rotoop(iold)
850         continue
            nops=nops1*(i+1)
849      continue
      endif
c
c--> add labels to centered symops
25    format(a14,f10.6)
      nops=nops+numvec
c
c--> print the matrix reps in operator form, with labels
      if(oprint) then
         write(*,1423) nops
         call opprint(symops,rotoop,maxops,nops,itype)
      endif
c
c dgc --> decenter if necessary
      if(numvec.gt.0) then
         !write(*,1420)
         use_primitive = geom_use_primitive(geom)
         if (use_primitive) then
            write(*,1421) group_name
            call dctr(symops,nops1,numgrp,group_name,numvec,cntvec,
     >                numset,geom)
            write(*,1424)
            write(*,1425)
            write(*,1426)
            write(*,1427)
            write(*,1428)
            write(*,1429) nops+1,nops1+1
            nops=nops1
            if (oprint) then
               write(*,1431)
               call opprint(symops,rotoop,maxops,nops,itype)
            end if
            !write(*,1430)
         else
            write(*,1422) group_name
         end if
      endif
c
1420  format(/' Primitive cell exists. ')
1421  format(/' Primitive cell exists. Converting the symmetry',
     & ' operators and the unit cell to primitive cell.'/,
     & ' To not do this conversion',
     & ' add the "conventional" keyword to the symmetry input, e.g.'//,
     & ' geometry'/," ..."/,' symmetry ',A,' conventional'/,
     & ' ...'/,' end'/)
1422  format(/' Primitve cell exists, but not converting to it.',
     & ' Using the conventional cell',
     & ' instead.'/,' To turn this conversion on',
     & ' add the "primitive" keyword to the symmetry input, e.g.'//,
     & ' geometry'/," ..."/,' symmetry ',A,' primitive'/,
     & ' ...'/,' end'/)
1423  format(//'The ',i3,' symmetry operators (excl. E)',
     &' are as follows:'/)
1424  format(/'After de-centering the following are redefined:'/)
1425  format(10X,' (1) lattice parameters a,b,c,alpha,beta, and gamma')
1426  format(10X,' (2) lattice vectors a1,a2,a3')
1427  format(10X,' (3) reciprocal lattice vectors b1,b2,b3')
1428  format(10X,' (4) fractional coordinates')
1429  format(/i3,' operators converted to ',i3,' symmetry operators.')
1430  format(/'DEBUG:primitive cell exists, but dctr was not called.'/)
1431  format(//i3,'The symmetry operators (excl. E) in the de-centered',
     &' basis are as follows:'/)
c
c rjh hack to fix C1
      if (numgrp .eq. 1) nops = 0
c
      end


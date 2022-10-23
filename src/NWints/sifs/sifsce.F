      function sifsce( nenrgy, energy, ietype )
c
c  select and sum the core energies in the energy(*) array.
c
c  usage: real*8   sifsce
c         external sifsce
c         total_core = sifsce(...)
c
c  input: nenrgy = number of energy(*) values.
c         ietype(1:nenrgy) = energy types.
c         energy(1:nenrgy) = energy array.
c
c  output: sifsce = total core energy such that
c                       total_potential = total_electronic + sifsce(...)
c                   is the total clamped-nucleus, born-oppenheimer
c                   potential.
c
      implicit none
c
      integer  nenrgy, ietype(nenrgy)
      real*8   sifsce, energy(nenrgy)
c
      integer  i,      itypea
      real*8   ecore
c
      ecore = (0)
      do 10 i = 1, nenrgy
c
c        # 0 <= ietype elements are 1-e hamiltonian frozen core terms.
c        # ietype < 0  elements with itypea=0 are other core terms to
c                      be added.  e.g. ietype=-1 = nuclear repulsion.
c
         itypea = ietype(i) / 1024
         if ( (ietype(i) .ge. 0 ) .or.
     &    ( (ietype(i) .lt. 0) .and. (itypea .eq. 0) ) ) then
            ecore = ecore + energy(i)
         endif
10    continue
c
      sifsce = ecore
c
      return
      end

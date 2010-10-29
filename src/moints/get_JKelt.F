       double precision function get_JKelt( g_a, nbf, i, j, a, b )
C$Id$
       implicit none
#include "global.fh"
       integer g_a
       integer nbf, i, j, a, b

       integer ij, ab
       double precision xx

c
c Handle permutational symmetry.
c
c
c       ij       ji       ij        ji
c      J    =   J     =  J     =   J
c       ab       ab       ba        ba
c
c
c       ij       ji
c      K    =   K  
c       ab       ba
c
c
       if (i.ge.j) then
         ab = (a-1)*nbf + b
         ij = (i*(i-1))/2 + j
       else
         ab = (b-1)*nbf + a
         ij = (j*(j-1))/2 + i
       endif
       call ga_get(g_a,ab,ab,ij,ij,xx,1)
       get_JKelt = xx
       return
       end



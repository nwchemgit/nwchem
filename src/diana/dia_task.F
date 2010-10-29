      subroutine dia_task
c
c $Id$
c
      implicit none
c
#include "dia_common.fh"
#include "dia_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      call dia_init()
c
    1 continue
c
      call dia_getcmd()
c
c     read reference coordinates
c     --------------------------
c
      if(cmd.eq.'refer ') then
      call dia_getref()
c
c     read topology file
c     ------------------
c
      elseif(cmd.eq.'topol ') then
      call dia_gettop()
c
c     set io mode
c     -----------
c
      elseif(cmd.eq.'iomode') then
      call dia_setiomode()
c
c     read trajectory file header
c     ---------------------------
c
      elseif(cmd.eq.'file  ') then
      call dia_gethdr()
c
c     select which frames to read
c     ---------------------------
c
      elseif(cmd.eq.'frame ') then
      call dia_gettrj()
c
c     atom selection
c     --------------
c
      elseif(cmd.eq.'select') then
      call dia_getsel()
c
c     save atom selection
c     -------------------
c
      elseif(cmd.eq.'savsel') then
      call dia_savsel()
      elseif(cmd.eq.'savsls') then
      call dia_savssl()
c
c     read atom selection
c     -------------------
c
      elseif(card(1:5).eq.'rdsel') then
      call dia_rdsel()
      elseif(cmd.eq.'rdsels') then
      call dia_rdssl()
c
c     geometric selections
c
      elseif(cmd.eq.'distan') then
      call dia_adddistance()
      elseif(cmd.eq.'angle') then
      call dia_addangle()
      elseif(cmd.eq.'torsio') then
      call dia_addtorsion()
      elseif(cmd.eq.'improp') then
      elseif(cmd.eq.'order') then
      call dia_addorder()
c
c     superimpose all stored frames
c     -----------------------------
c
      elseif(card(1:5).eq.'super') then
      call dia_super()
c
c     calculate average coordinates
c     -----------------------------
c
      elseif(cmd.eq.'xaver ') then
      call dia_average()
c
c     copy the trajectory that is currently in memory
c     -----------------------------------------------
c
      elseif(cmd.eq.'copy ') then
      call dia_copy()
c
c     calculate one-dimensional root mean square deviation
c     ----------------------------------------------------
c
      elseif(cmd.eq.'rmsdev') then
      call dia_rmsd()
c
c     build covariance matrix
c     -----------------------
c
      elseif(cmd.eq.'covar ') then
      call dia_covar()
c
c     pca projection
c     --------------
c
      elseif(cmd.eq.'projec') then
      call dia_project()
c
      elseif(cmd.eq.'evalua') then
      call dia_evaluate()
c
      elseif(cmd.eq.'espplt') then
      call dia_espplt()
c
      elseif(cmd.eq.'rama  ') then
      call dia_rama()
      else
      goto 9
      endif
      goto 1
c
    9 continue
      call dia_final()
c
      return
      end

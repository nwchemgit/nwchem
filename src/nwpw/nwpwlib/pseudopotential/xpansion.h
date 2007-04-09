/*
 $Id: xpansion.h,v 1.2 2007-04-09 22:55:52 d3p708 Exp $
*/
#ifndef _GET_CS_H_
#define _GET_CS_H_
/* get_cs.h
*/

extern  void	p_xpansion();
extern  void	dp_xpansion();
extern  void	ddp_xpansion();
extern  void	dddp_xpansion();
extern  void	ddddp_xpansion();
extern	void	psp_xpansion();
extern	void	psi_xpansion();
extern	void	dpsi_xpansion();
extern	void	init_xpansion(int,int,double,double,
                              double*, double*, double*,
                              double*, double*, double*);

#endif

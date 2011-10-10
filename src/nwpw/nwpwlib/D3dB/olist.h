#ifndef _OLIST_H_
#define _OLIST_H_
/* olist.h -
   Author - Eric Bylaska

*/

typedef struct {
	int	max_index;
        int	*list;
} OList_Type;

extern	void	create_olist();
extern	void	insert_olist();
extern	int	index_olist();
extern	void	destroy_olist();
extern	void	print_olist();

#endif
/* $Id$ */

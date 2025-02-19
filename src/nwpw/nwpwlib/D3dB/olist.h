#ifndef _OLIST_H_
#define _OLIST_H_
/* olist.h -
   Author - Eric Bylaska

*/

typedef struct {
	int	max_index;
        int	*list;
} OList_Type;

extern	void	create_olist(OList_Type*, int);
extern	void	insert_olist(OList_Type*, int);
extern	int	index_olist(OList_Type*, int);
extern	void	destroy_olist(OList_Type*);
extern	void	print_olist(OList_Type*);

#endif
/* $Id$ */

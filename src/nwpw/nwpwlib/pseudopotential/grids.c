/* grid.c
   author - Eric Bylaska
   $Id$

*/

#include	<stdlib.h>
#include	<stdio.h>
#include	"grids.h"

/* use a stack for the unused grids */
#define	NIL	((void *) 0)
#define	Push_Stack(s, e)	s = (e->next = s)   ? e : e
#define Pop_Stack(s, e)		s = ( (e=s) == NIL) ? NIL : s->next




/***************************/
/* the grid data structure */
/***************************/

/* use a linked list to keep tract of grids */
typedef	struct tt_struct {
    struct tt_struct	*next;
    double			*grid;
} *Grids_List_Type;

static	int		Ngrid_points;
static	Grids_List_Type	 using_grids_list;
static	Grids_List_Type  unused_grids_stack;


/********************************
 *				*
 *         init_Grids		*
 *				*
 ********************************/

void	init_Grids(Np)
int	Np;
{
    Ngrid_points = Np;
    using_grids_list  = NIL;
    unused_grids_stack = NIL;

} /* init_Grids */

/********************************
 *				*
 *        end_Grids		*
 *				*
 ********************************/
void	end_Grids()
{
    Grids_List_Type node;

    while (using_grids_list != NIL)
    {
        Pop_Stack(using_grids_list,node);
        free(node->grid);
        free(node);
    }

    while (unused_grids_stack != NIL)
    {
        Pop_Stack(unused_grids_stack,node);
        free(node->grid);
        free(node);
    }
}


/********************************
 *				*
 *        alloc_Grid		*
 *				*
 ********************************/

double	*alloc_Grid()
{
    double 	   *grid_return;
    Grids_List_Type node;


    /* get grid form unused grid stack */
    if (unused_grids_stack != NIL)
    {
        Pop_Stack(unused_grids_stack,node);
    }

    /* allocate memory */
    else
    {
        node       = (Grids_List_Type) malloc(sizeof(struct tt_struct));
        node->grid = (double *) malloc(Ngrid_points*sizeof(double));

    }

    /* put node on using grids list */
    node->next       = using_grids_list;
    using_grids_list = node;

    /* access the pointer to the grid array */
    grid_return = node->grid;

    return grid_return;
} /* alloc_Grid */



/********************************
 *				*
 *      dealloc_Grid		*
 *				*
 ********************************/

void	dealloc_Grid(double * grid)
{
    int	    done;
    Grids_List_Type  cur,prev;

    /* find grid on using grids list */
    cur  = using_grids_list;
    prev = using_grids_list;
    done = 0; /*false */
    while ((cur != NIL) && (!done))
    {
        if ((cur->grid) == grid)
            done = 1; /* true */
        else
        {
            prev = cur;
            cur  = cur->next;
        }
    } /*while*/

    /* error: grid not found */
    if (!done)
    {
        printf("Error in dealloc_Grid: grid not found\n");
        printf("              address: %lx\n", grid);
        cur = using_grids_list;
        printf("Address List:");
        while (cur != NIL)
        {
            printf("%lx ", cur->grid);
            cur = cur->next;
        }
        printf("\n");
        exit(93);
    }

    /* remove cur from using grids list */

    /* special case, grid on the head of the list */
    if (prev==cur)
        using_grids_list = using_grids_list->next;

    /* regular case */
    else
        prev->next = cur->next;


    /* put cur on unused grids stack */
    cur->next  = NIL;
    Push_Stack(unused_grids_stack, cur);

} /* dealloc_Grid */


/* 
*******************************************************************************
*
* File:         dist.h
* RCS:          $Header: /tmp/mss/nwchem/src/perfm/dist.h,v 1.1 2006-05-19 18:59:26 edo Exp $
* Description:  Distribution table implemented with an array
* Author:       Fabrizio Petrini
* Created:
* Modified:
* Language:     C
* Package:      N/A
* Status:       Experimental (Do Not Distribute)
*
*******************************************************************************
*/


#ifndef __dist__
#define __dist__


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>


/* the queue event: contains the event value */
typedef struct event
{
  /* the event value */
  unsigned int value;

  /* the event count */
  unsigned int count;

  /* next pointer */
  struct event* next;
} event;


/* the distribution data structure: simply an array of ints and a
   queue to deal with the large elements that overflow the array */
typedef struct distr
{
  /* the distribution name */
  char* name;

  /* the array of counters */
  unsigned int* array;
  
  /* the array size */
  unsigned int size;

  /* the maximum array size: larger elements are stored in the queue */
  unsigned int max_size;

  /* the queue head */
  event* head;

  /* the queue tail */
  event* tail;

  /* the number of events in the queue */
  unsigned int queue_size;
} distribution;


/* create the distribution: take the name and the initial size and the
   threshold of the array size, for which larger elements are stored
   in the queue */
distribution* create_distribution( const char* dname, unsigned int size, 
				   unsigned int max_size );

/* delete the distribution: take the distribution pointer */
void delete_distribution( distribution* d );

/* update the distribution: take the distribution pointer and the
   value */
void update_distribution( distribution* d, unsigned int val );

/* as above plus a counter */
void update_distribution_count( distribution* d, unsigned int val, unsigned int count );

/* print the given distributuion: take the distribution pointer, the
   file prefix, a index to differentiate distributions in the same
   family and a coefficient that is used to scale the values of the
   distribution */
void print_distribution( distribution* d, const char* file_name, unsigned int id, 
			 double f );

/* file open and close */
FILE* open_file( const char* name, const char* mode ); 
void close_file( FILE* fp );
FILE* create_data_file( const char* prefix, const char* graph );


/*
*******************************************************************************
*
*                          user-level C interface
*
*******************************************************************************
*/


/* initialize the monitoring system: we access the distributions
   through integer handles */
void initialize_dist();

/* finalize the monitoring system, flushing all the distribution files
   and writing some usage statistics; take the process id */
void finalize_dist( unsigned int );

/* create a new distribution: take the file name and initial
   distribution size and the threshold to use a queue rather than an
   array; return the distribution handle, a unique index inside the
   distribution array which identifies the distribution */
unsigned int newdist( const char* filename, const unsigned int dist_size, 
		      const unsigned int max_dist_size );

/* return the distribution index, given the the distribution name as a
   string: creates a distribution with a default size, if there is no
   distribution available */
unsigned int getdist( const char* distname );

/* update the distribution */
void updist( const unsigned int handle, const unsigned int value );
void updistcount( const unsigned int handle, const unsigned int value, 
		  const unsigned int count );

/* print the distribution */
void printdist( const unsigned int handle, const char* file_name, 
		unsigned int id, double f );


/*
*******************************************************************************
*
*                         user-level FORTRAN interface
*
*******************************************************************************
*/


/* initialize the monitoring system: we access the distributions
   through integer handles */
void initialize_dist_();

/* finalize the monitoring system, flushing all the distribution files
   and writing some usage statistics */
void finalize_dist_( const unsigned int* );

/* return the distribution index, given the the distribution name as a
   string */
unsigned int getdist_( const char* distname,  
		       const unsigned int distname_length );

/* create a new distribution: take the file name in FORTRAN format,
   the distribution size and the maximum distribution size; return the
   distribution handle, a unique index inside the distribution array
   which identifies the distribution */
unsigned int newdist_( const char* filename, const unsigned int* dist_size, 
		       const unsigned int* max_dist_size, 
		       const unsigned int filename_length );

/* update the distribution */
void updist_( const unsigned int* handle, const unsigned int* value );
void updistcount_( const unsigned int* handle, const unsigned int* value, 
		   const unsigned int* count );

/* print the distribution */ 
void printdist_( const unsigned int* handle, const char* filename, 
		 unsigned int* fd, double* f, 
		 const unsigned int filename_length );

/* Timers */

/* start the timer */
void starttimer_( const unsigned int* handle );

/* stop the timer and stores its value in the distribution: by
   default, it converts it to microseconds */
void endtimer_( const unsigned int* handle );

/* as above, but in milliseconds */
void endtimermilli_( const unsigned int* handle );

/* and, just in case, in seconds */
void endtimersec_( const unsigned int* handle );


#endif /* __dist__ */

/* $Id$ */

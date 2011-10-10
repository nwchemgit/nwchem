/*
*******************************************************************************
*
* File:         dist.c
* RCS:          $Header: /tmp/mss/nwchem/src/perfm/dist.c,v 1.3 2006-05-23 18:33:05 edo Exp $
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


#include "dist.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/*#define MPI_TIME 1*/


/*
*******************************************************************************
*
*                               constants 
*
*******************************************************************************
*/


/* define to control several invariants */
/*#define CHECKING*/
#undef CHECKING

#ifdef CHECKING
#define CHECK(x)  x
#else
#define CHECK(x)
#endif

#define MAX_STRING                      1024
#define MAX_ID                            10
#define TABLE_SIZE                      1024
#define DEFAULT_DISTRIBUTION_SIZE      16384
#define DEFAULT_MAX_DISTRIBUTION_SIZE 131072

#ifdef FDIST
/*
*******************************************************************************
*
*                               diagnostic
*
*******************************************************************************
*/


void dist_error( char* func, char* msg )
{
  printf( "Error in %s: %s\n", func, msg );
  exit(0);
}


/*
*******************************************************************************
*
*                                Timers
*
*******************************************************************************
*/


/* define to use MPI_Wtime, undefine to use gettimeofday */
/*#define MPI_TIME 1*/

#define timersub(a, b, result)                                                \
  do {                                                                        \
    (result)->tv_sec = (a)->tv_sec - (b)->tv_sec;                             \
    (result)->tv_usec = (a)->tv_usec - (b)->tv_usec;                          \
    if ((result)->tv_usec < 0) {                                              \
      --(result)->tv_sec;                                                     \
      (result)->tv_usec += 1000000;                                           \
    }                                                                         \
  } while (0)


#define set_time( time, tv_1, tv_2, tv_diff )  \
        timersub(&tv_2, &tv_1, &tv_diff); \
        *time = ((double)(tv_diff.tv_sec) * 1.0E+6 + \
                (double)(tv_diff.tv_usec)); \
        *time = *time / 1.0E+6

/* get_time set the time based on the timer type, that can be either a
   double, when we use MPI or a struct timeval with gettimeof day 

   diff( t, t2, t1 ) set the first parameter, which must be a double,
   with the difference between the two timer values, t2 and t1.

   the time is returned in seconds
*/


#ifdef MPI_TIME

#include <mpi.h>

#define timer_t double
#define get_noclock(t) t = 0;
#define get_clock(t) t = MPI_Wtime();
#define diff( t, t2, t1 ) { t = (t2 - t1); }
#else
#define timer_t struct timeval
#define get_clock(t) gettimeofday( &t, &tz );
#define diff( t, t2, t1 ) { set_time( &t, t1, t2, t_diff ); }
#endif

timer_t start_time, end_time;

#ifndef MPI_TIME
/* the timezone, used by gettimeofday */
struct timezone tz;
struct timeval t_diff;
#endif


/*
*******************************************************************************
*
*                          distribution queue
*
*******************************************************************************
*/


/* create a new event */
event* create_event( unsigned int value, unsigned int count )
{
  event* e;
  
  /* malloc the data structure */
  if ( !( e = (event*) malloc( sizeof( event) ) ) )
    dist_error( "create_event", "unable to malloc event" );
  
  e->value = value;
  e->count = count;
  e->next = 0;

  /* return the event */
  return e;
}


void check_empty( distribution* d ) 
{
  if ( !d->queue_size )
    dist_error( "check_empty", "the queue is empty" );
}


void check_status( distribution* d ) 
{
  event* e;
  unsigned int i;

  if ( ( d->queue_size && (!d->head || !d->tail) ) || 
       ( !d->queue_size && ( d->head || d->tail ) ) )
    dist_error("check_empty", "wrong internal status");
  
  /* scan the list */
  e = d->head;
  for( i = 0; i < d->queue_size; i++ ) 
    {
    if ( !e )
      dist_error("check_status", "the list is too short");
    else
      e = e->next;
  }
  if ( e ) 
    dist_error("check_status", "the list is too long");
}


/* return the number of events in the queue */
unsigned int get_queue_size( distribution* d ) 
{ 
  return d->queue_size; 
}


/* return a pointer to the first event of the queue */
event* top_event(  distribution* d ) 
{ 
  return d->head; 
}


void put_event( distribution* d, unsigned int value, unsigned int count ) 
{
  event* scan;
  event* e;

  CHECK(check_status( d );)

  if ( d->queue_size ) 
    {
      /* priority queue: scan the list and insert the event according
	 to its value */
      if ( d->head->value == value )
	d->head->count++;
      else if ( d->head->value > value ) 
	{
	  /* insertion at the top */
	  e = create_event( value, count );
	  e->next = d->head;
	  d->head = e;

	  /* increment queue_size */
	  d->queue_size++;
	}
      else 
	{
	  /* scan the list */
	  scan = d->head;
	  while ( scan->next && ( scan->next->value < value ) )
	    scan = scan->next;
	  
	  if ( scan->next && scan->next->value == value )
	    scan->next->count++;
	  else
	    {
	      e = create_event( value, count );
	      e->next = scan->next;
	      scan->next =  e;

	      /* set the tail pointer, if necessary */
	      if ( !e->next )
		d->tail = e;
	      
	      /* increment queue_size */
	      d->queue_size++;
	    }
	  
	}
    }
  else
    {
      /* the queue is empty */
      d->head = d->tail = create_event( value, count );
      
      /* increment queue_size */
      d->queue_size++;
    }	

  CHECK(check_status( d );)
}


event* get_event( distribution* d )
{
  event* tmp;

  /* there must be at least an event */
  CHECK(check_empty( d );)
  CHECK(check_status( d );)

  tmp = d->head;
  d->head = d->head->next;

  /* reset the tail pointer, if this is the last event */
  if ( !d->head )
    d->tail = d->head;
  d->queue_size--;
  CHECK(check_status( d );)
  return tmp;
}


/*
*******************************************************************************
*
*                              distributions
*
*******************************************************************************
*/

static unsigned int ENABLE_MONITORING = 0;
static distribution* dtable[TABLE_SIZE];
static unsigned int distribution_count = 0;
static timer_t table_timer[TABLE_SIZE];


FILE* open_file( const char* name, const char* mode ) 
{
  FILE* fp;

  if ( ( fp = fopen( name, mode ) ) == NULL )
    {
      printf( "open_file: unable to open file %s, %s mode, exiting\n", 
	      name, mode );
      exit(1);
    }
  else
    return fp;
}


void close_file( FILE* fp )
{
  if ( fclose( fp ) )
    {
      printf( "close_file: unable to close file descriptor, exiting\n" );
      exit(1);
    }
}


FILE* create_data_file( const char* prefix, const char* graph )
{
  FILE* fp;
  char* directory_name;
  char file_name[MAX_STRING];

  /* check the env variable */
  directory_name = getenv ( "DIST_DIRECTORY" );

  /* create the file name  */
  if ( directory_name ) 
    sprintf( file_name, "%s/%s%s.dat", directory_name, prefix, graph );
  else
    sprintf( file_name, "%s%s.dat", prefix, graph );
  
  return fp = open_file( file_name, "w" );
}


distribution* create_distribution( const char* dname, unsigned int size, 
				   unsigned int max_size )
{
  unsigned int i;
  distribution* d;

  if ( !dname )
    dist_error( "create_distribution", "distribution name is a null pointer" );

  if ( strlen( dname ) > MAX_STRING )
    dist_error( "create_distribution", "distribution name is too long" );

  if ( size <= 0 )
    dist_error( "create_distribution", 
		"distribution size must be a positive number" );

  if ( max_size <= 0 || max_size < size )
    dist_error( "create_distribution", 
		"distribution threshold must be larger than size" );


  /* allocate the data structure */
  d = (distribution*) malloc ( sizeof( distribution ) );
  d->name = (char*) malloc( strlen( dname ) + 1 );

  /* initialize the fields */
  strcpy( d->name, dname );
  d->size = size;
  d->max_size = max_size;
  d->array = (unsigned int*) malloc( sizeof( unsigned int ) * size );

  for( i = 0; i < size; i++ )
    d->array[i] = 0;

  /* initialize the queue */
  d->head = d->tail = 0;
  d->queue_size = 0;

  return d;
}


void delete_distribution( distribution* d )
{
  free( d->array );
  free( d );
}


void rehash( distribution* d, unsigned int val )
{
  unsigned int new_size, i;
  unsigned int* new_array;

  new_size = d->size;
  
  /* we need to resize the array */
  while ( new_size <= val )
    new_size = 2 * new_size;
  
  /* allocate the new array */
  new_array = (unsigned int*) malloc ( sizeof(unsigned int) * new_size );
  
  /* copy the old values */
  for( i = 0; i < d->size; i++ )
    new_array[i] = d->array[i];
  
  /* and reset the remaining part of the array */
  for( i = d->size; i < new_size; i++ )
    new_array[i] = 0;
  
  /* free the old array */
  free( d->array );
  
  /* update the internal fields */
  d->size = new_size;
  d->array = new_array;
}


void update_distribution( distribution* d, unsigned int val )
{
  if ( val >= d->size && val < d->max_size )
    rehash( d, val );
  else if ( val >= d->max_size )
    {
      put_event( d, val, 1 );
      return;
    }

  /* update the array */
  d->array[val]++;
}


void update_distribution_count( distribution* d, unsigned int val, 
				unsigned int count )
{
  if ( val >= d->size )
    rehash( d, val );
  else if ( val >= d->max_size )
    {
      put_event( d, val, count );
      return;
    }

 
  if ( count <= 0 )
    {
      printf( "update_distribution_count increment value is <= 0 %u\n", 
	      count );
      return; 
    }

  /* update the array */
  d->array[val] += count;
}


void print_distribution( distribution* d, const char* file_name, 
			 unsigned int id, double f )
{
  unsigned int i;
  char string_id[MAX_ID];
  FILE* distribution_fp;
  unsigned int count = 0, event_count = 0, count_bins = 0, 
    min = 0xffffffff, max = 0;
  long long unsigned int sum = 0;

  


  sprintf( string_id, "_%u", id );

  distribution_fp = create_data_file( file_name, string_id );

  /* scan the array */
  for( i = 0; i < d->size; i++ )
    {
      if ( d->array[i] )
	{
	  fprintf( distribution_fp, "%u %f\n", i, 
		   f * ((double) d->array[i] ) );
	  count += d->array[i];
	  count_bins ++;
	  sum += i * d->array[i];
	  
	  if ( i > max ) 
	    max = i;

	  if ( i < min )
	    min = i;
	}
    }

  /* scan the queue of events */
  while( get_queue_size( d ) )
    {
      event* ev = get_event( d );

      fprintf( distribution_fp, "%u %f\n", ev->value, 
	       f * ((double) ev->count) );

      count += ev->count;
      count_bins ++;
      sum += ev->value * ev->count;
      
      if ( ev->value > max ) 
	max = ev->value;
      
      if ( ev->value < min )
	min = ev->value;

      event_count++;

      /* delete the event */
      free( ev );
    }

  /* print the stats at the end of the file */
  fprintf( distribution_fp, "# Count\t%u\n", count );

  if ( count ) 
    {
      fprintf( distribution_fp, "# Bins\t%u\n", count_bins );
      fprintf( distribution_fp, "# Min\t%u\n", min );
      fprintf( distribution_fp, "# Max\t%u\n", max );
      fprintf( distribution_fp, "# Sum\t%lld\n", sum );
      fprintf( distribution_fp, "# Avg\t%f\n", (double) sum/(double) count );
      fprintf( distribution_fp, "# Integral\t%f\n", 
	       (double) sum/((double) max - min + 1) );
    }

  fprintf( distribution_fp, "# Memory Usage\t%u bytes\n", 
	   d->size * sizeof(unsigned int) + event_count * sizeof( event ) );
  
  /* close the data file */
  close_file( distribution_fp );
}


/*
*******************************************************************************
*
*                           C user interface
*
*******************************************************************************
*/


void enable_dist()
{
  ENABLE_MONITORING = 1;
}


void disable_dist()
{
  ENABLE_MONITORING = 0;
}


void initialize_dist()
{
  char *monitoring;

  if ( ( monitoring = getenv( "ENABLE_MONITORING" ) ) )
    {
      if ( !strcmp( monitoring, "ON" ) )
	enable_dist();
      else
	disable_dist();
    }
  else
    enable_dist();
}


void finalize_dist( unsigned int id )
{
  unsigned int i;

  if ( !ENABLE_MONITORING )
    return;

  for( i = 0; i < distribution_count; i++ )
    {
      if ( !dtable[i] )
	{
	  printf( "finalize_dist: skipping null dist %i, table size %u\n",
		  i, distribution_count );
	  continue;
	}
	 
      print_distribution( dtable[i], dtable[i]->name, id, 1.0 );
    }
}


unsigned int newdist( const char* file_name, const unsigned int dist_size, 
	     const unsigned int max_dist_size )
{
  if ( !ENABLE_MONITORING )
    return 0;

  if ( distribution_count >= TABLE_SIZE )
    {
      printf( "newdist: overflowing distribution table %u\n", 
	      distribution_count );
      exit(0); 
    }

  if ( !file_name )
    {
      printf( "newdist: file_name is a null pointer\n" );
      exit(0);
    }
  
  if ( dist_size <= 0 )
    {
      printf( "newdist: *dist_size is not a positive number %u\n", 
	      dist_size );
      exit(0);
    }
  
  dtable[distribution_count++] = create_distribution( file_name, dist_size, 
						      max_dist_size );

  return distribution_count -1;
}

/* return the distribution index, given the the distribution name as a
   string: creates a distribution with a default size, if there is no
   distribution available */
unsigned int getdist( const char* dist_name )
{
  unsigned int i;

  if ( !ENABLE_MONITORING )
    return 0;

  if ( !dist_name )
   {
     printf( "getdist: dist_name is a null pointer\n" );
     exit(0);
   }

  for( i =0; i < distribution_count; i++ )
      if ( !strcmp( dist_name, dtable[i]->name ) )
	return i;

  /* no matching distribution: create a new distribution */
  return newdist( dist_name, DEFAULT_DISTRIBUTION_SIZE, 
		  DEFAULT_MAX_DISTRIBUTION_SIZE );
}


void updist( const unsigned int handle, const unsigned int value )
{
  if ( !ENABLE_MONITORING )
    return;

  if ( handle >= TABLE_SIZE )
    {
      printf( "updist: handle %u is larger than table size %u\n",
	      handle, TABLE_SIZE );
      exit(0);
    }
  
  update_distribution( dtable[handle], value );
}


void updistcount( const unsigned int handle, const unsigned int value, 
		  const unsigned int count )
{
  if ( !ENABLE_MONITORING )
    return;

  if ( handle >= TABLE_SIZE )
    {
      printf( "updistcount: handle %u is larger than the table size %u\n", 
	      handle, TABLE_SIZE );
      exit(0);
    }

  if ( count < 1 )
    {
      printf( "updistcount: count %u is smaller than 1\n", count );
      exit(0);
    }

  update_distribution_count( dtable[handle], value, count );
}


void printdist( const unsigned int handle, const char* file_name, 
		unsigned int id, double f )
{
  if ( !ENABLE_MONITORING )
    return;

  if ( handle >= TABLE_SIZE )
    {
      printf( "printdist: handle %u is larger than the table size %u\n", 
	      handle, TABLE_SIZE );
      exit(0);
    }

  if ( !file_name )
    {
      printf( "printdist: file_name is a null pointer\n" );
      exit(0);
    }

  print_distribution( dtable[handle], file_name, id, f );
}


/*
*******************************************************************************
*
*                            FORTRAN user interface
*
*******************************************************************************
*/


/* converts a FORTRAN string in a C string */
static int f2c_string( const char *f, int flen, char *buf,
		       const int buflen )
{
  while (flen-- && f[flen] == ' ')
    ;
  
  if ((flen+1) >= buflen)
    /* won't fit */
    return 0;                   
  
  flen++;
  buf[flen] = 0;
  while(flen--)
    buf[flen] = f[flen];
  
  return 1;
}


void enable_dist_()
{
  enable_dist();
}


void disable_dist_()
{
  disable_dist();
}


void initialize_dist_()
{
  initialize_dist();
}



void finalize_dist_( const unsigned int* id )
{
  if ( !id )
    {
      printf("finalize_dist_: id %d is a null pointer\n" );
      exit(0);
    }

  finalize_dist( *id );
}

/* create a new distribution: take the file name in FORTRAN format and
   the distribution size; return the distribution handle, a unique
   index inside the distribution array */
unsigned int newdist_( const char* file_name, const unsigned int* dist_size, 
		       const unsigned int* max_dist_size,
		       const unsigned int file_name_length )
{
  char buf[MAX_STRING];

  if ( !ENABLE_MONITORING )
    return;

  if ( !file_name )
    {
      printf( "newdist_: file_name is a null pointer\n" );
      exit(0);
    }

  if ( file_name_length == 0 )
    {
      printf( "newdist_: file_name_length is not a positive number %d\n", 
	      file_name_length );
      exit(0);
    }
  
  if ( !dist_size )
    {
      printf( "newdist_: dist_size is a null pointer\n" );
      exit(0);
    }

  if ( *dist_size == 0 )
    {
      printf( "newdist_: *dist_size is not a positive number %d\n", 
	      *dist_size );
      exit(0);
    }

  /* converts the string in C format */
  if ( !f2c_string( file_name, file_name_length, buf, MAX_STRING ) )
    {
      printf( "newdist_: string does not fit in buffer\n" );
      exit(0);
    }
  
  dtable[distribution_count++] = create_distribution( buf, *dist_size, 
						      *max_dist_size );

  return distribution_count - 1;
}


/* return the distribution index, given the distribution name as a
   string: creates a distribution with a default size, if there is no
   distribution available */
unsigned int getdist_( const char* dist_name, 
		       const unsigned int dist_name_length )
{
  int i;
  char buf[MAX_STRING];

  if ( !ENABLE_MONITORING )
    return;

  if ( !dist_name )
   {
     printf( "getdist_: dist_name is a null pointer\n" );
     exit(0);
   }

  if ( dist_name_length == 0 )
    {
     printf( "getdist_: dist_name_length %d is not a positive number\n",
	     dist_name_length );
     exit(0);
    }

  /* converts the string in C format */
  if ( !f2c_string( dist_name, dist_name_length, buf, MAX_STRING ) )
    {
      printf( "getdist_: string does not fit in buffer\n" );
      exit(0);
    }
  
  for( i = 0; i < distribution_count; i++ )
      if ( !strcmp( buf, dtable[i]->name ) )
	return i;

  /* no matching distribution: create a new distribution */
  return newdist( buf, DEFAULT_DISTRIBUTION_SIZE, 
		  DEFAULT_MAX_DISTRIBUTION_SIZE );
}


void updist_( const unsigned int* handle, const unsigned int* value )
{
  if ( !ENABLE_MONITORING )
    return;

  if ( !handle )
    {
      printf( "updist_: handle is a null pointer\n" );
      exit(0);
    }

  if ( *handle >= TABLE_SIZE )
    {
      printf( "updist_: *handle %d is larger than the table size %d\n",
	      *handle, TABLE_SIZE );
      exit(0);
    }
  
  if ( !value )
    {
      printf( "updist_: value is a null pointer\n" );
      exit(0);
    }
  
  update_distribution( dtable[*handle], *value );
}


void updistcount_( const unsigned int* handle, const unsigned int* value, 
		   const unsigned int* count )
{
  if ( !ENABLE_MONITORING )
    return;

  if ( !handle )
    {
      printf( "updistcount_: handle is a null pointer\n" );
      exit(0);
    }

  if ( *handle >= TABLE_SIZE )
    {
      printf( "updistcount_: *handle %d is larger than the table size %d\n", 
	      *handle, TABLE_SIZE );
      exit(0);
    }
  
  if ( !value )
    {
      printf( "updistcount_: value is a null pointer\n" );
      exit(0);
    }

  if ( !count )
    {
      printf( "updistcount_: count is a null pointer\n" );
      exit(0);
    }

  update_distribution_count( dtable[*handle], *value, *count );
}


void printdist_( const unsigned int* handle, const char* file_name, 
		 unsigned int* id, double* f,
		 const unsigned int file_name_length )
{
  char buf[MAX_STRING];

  if ( !ENABLE_MONITORING )
    return;

  /* converts the string in C format */
  if ( !f2c_string( file_name, file_name_length, buf, MAX_STRING ) )
    {
      printf( "printdist_: *string does not fit in buffer\n" );
      exit(0);
    }

  if ( !handle )
    {
      printf( "printdist_: handle is a null pointer\n" );
      exit(0);
    }

  if ( *handle < 0 )
    {
      printf( "printdist_: *handle is a negative number %d\n", 
	      *handle );
      exit(0);
    }

  if ( *handle >= TABLE_SIZE )
    {
      printf( "printdist_: *handle %d is larger than the table size %d\n", 
	      *handle, TABLE_SIZE );
      exit(0);
    }

  if ( !file_name )
    {
      printf( "printdist_: file_name is a null pointer\n" );
      exit(0);
    }

  if ( file_name_length <= 0 )
    {
      printf( "printdist_: file_name_length is not a positive number %d\n", 
	      file_name_length );
      exit(0);
    }
  
  if ( !id )
    {
      printf( "printdist_: id is a null pointer" ); 
      exit(0);
    }

  if ( !f )
    {
      printf( "printdist_: f is a null pointer" ); 
      exit(0);
    }

  print_distribution( dtable[*handle], buf, *id, *f );
}


/* Timers */

/* star the timer */
void starttimer_( const unsigned int* handle )
{
  if ( !ENABLE_MONITORING )
    return;

  if ( !handle )
    {
      printf( "starttimer_: handle is a null pointer\n" );
      exit(0);
    }

  if ( *handle < 0 )
    {
      printf( "starttimer_: *handle is a negative number %d\n", 
	      *handle );
      exit(0);
    }

  if ( *handle >= TABLE_SIZE )
    {
      printf( "starttimer_: *handle is out of range %d TABLE_SIZE\n", 
	      *handle, TABLE_SIZE );
      exit(0);
    }

  get_clock( table_timer[*handle] );
}


/* stop the timer and stores its value in the distribution: by
   default, it converts it to microseconds */
void endtimer_( const unsigned int* handle )
{
  double t;

  if ( !ENABLE_MONITORING )
    return;

  get_clock( end_time );

  if ( !handle )
    {
      printf( "endtimer_: handle is a null pointer\n" );
      exit(0);
    }

  if ( *handle < 0 )
    {
      printf( "endtimer_: *handle is a negative number %d\n", 
	      *handle );
      exit(0);
    }
 
  if ( *handle >= TABLE_SIZE )
    {
      printf( "endtimer_: *handle is out of range %d TABLE_SIZE\n", 
	      *handle, TABLE_SIZE );
      exit(0);
    }

  diff( t, end_time, table_timer[*handle] );

  update_distribution( dtable[*handle], (unsigned int) (t*1000000) );
}


/* as above, but in milliseconds */
void endtimermilli_( const unsigned int* handle )
{
  double t;

  if ( !ENABLE_MONITORING )
    return;

  get_clock( end_time );

  if ( !handle )
    {
      printf( "endtimermilli_: handle is a null pointer\n" );
      exit(0);
    }

  if ( *handle < 0 )
    {
      printf( "endtimermilli_: *handle is a negative number %d\n", 
	      *handle );
      exit(0);
    }

  if ( *handle >= TABLE_SIZE )
    {
      printf( "endtimermilli_: *handle is out of range %d TABLE_SIZE\n", 
	      *handle, TABLE_SIZE );
      exit(0);
    }
 
  diff( t, end_time, table_timer[*handle] );

  update_distribution( dtable[*handle], (unsigned int) (t*1.0E+3) );
}


/* and, just in case, in seconds */
void endtimersec_( const unsigned int* handle )
{
  double t;

  if ( !ENABLE_MONITORING )
    return;

  get_clock( end_time );

  if ( !handle )
    {
      printf( "endtimersec_: handle is a null pointer\n" );
      exit(0);
    }

  if ( *handle < 0 )
    {
      printf( "endtimersec_: *handle is a negative number %d\n", 
	      *handle );
      exit(0);
    }

   if ( *handle >= TABLE_SIZE )
    {
      printf( "endtimersec_: *handle is out of range %d TABLE_SIZE\n", 
	      *handle, TABLE_SIZE );
      exit(0);
    }

  diff( t, end_time, table_timer[*handle] );

  update_distribution( dtable[*handle], (unsigned int) t );
}


//#define DEBUG_DIST
#ifdef  DEBUG_DIST

main()
{
  unsigned int handle, handle_time, i, j, id = 20;
  double f = 1.0;
  
  initialize_dist_();

  /* C functions */
  handle = getdist( "pluto" ) ;
 
  for ( i = 0; i < DEFAULT_MAX_DISTRIBUTION_SIZE + 100; i++ )
    {
      updist( handle, i );
    }
  
  //  printdist( handle, "pluto", 10, 1.0 );
  
  /* Fortran functions */
  handle = getdist_( "pippo", strlen( "pippo" ) + 1 );
  handle_time = getdist_( "pippo_time", strlen( "pippo_time" ) + 1 ) ;
  
  for ( i = 0; i < DEFAULT_MAX_DISTRIBUTION_SIZE + 100; i++ )
    {
      starttimer_( &handle_time );
      updist_( &handle, &i );
      endtimer_( &handle_time );
    }
  
  //printdist_( &handle, "pippo", &id, &f, strlen( "pippo" ) + 1 );
  //printdist_( &handle_time, "pippo_time", &id, &f, 
  //       strlen( "pippo_time" ) + 1 );
  finalize_dist_( &id );
}


#endif /* DEBUG_DIST */
#else
/* stubs */
void starttimer_( const unsigned int* handle )
{
    return;
}
void endtimer_( const unsigned int* handle )
{
    return;
}
void updist_( const unsigned int* handle, const unsigned int* value )
{
    return;
}
unsigned int getdist_( const char* dist_name, 
		       const unsigned int dist_name_length )
{
    return (unsigned int) 0;
}
void initialize_dist_()
{
  return;
}
void finalize_dist_(const unsigned int* id )
{
  return;
}
#endif
/* $Id$ */

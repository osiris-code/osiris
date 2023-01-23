/*****************************************************************************************

Simple C interface to PAPI routines. Currently only the following hardware counters are
supported:
x86: PAPI_TOT_INS, PAPI_FP_INS and PAPI_L2_DCM
BlueGene/Q: PAPI_TOT_INS, PAPI_FP_OPS and PAPI_L1_DCM

*****************************************************************************************/

#include "fortran.h"

#include <stdlib.h>
#include <stdio.h>

#include <papi.h>


#ifdef __bgq__

/* Events for the BlueGene/Q */
#define NUM_PAPI_EVENTS 3
char *papi_event_name[NUM_PAPI_EVENTS] = { "PAPI_TOT_INS" , "PAPI_FP_OPS" , "PAPI_L1_DCM" }; 

// I couldn't get L2Unit:::PEVT_L2_MISSES to work, we get a PAPI_EINVAL error when we
// call PAPI_start_counters

#else

/* Default events */
#define NUM_PAPI_EVENTS 3
char *papi_event_name[NUM_PAPI_EVENTS] = { "PAPI_TOT_INS" , "PAPI_FP_INS" , "PAPI_L2_DCM" }; 

#endif

static int papi_events[ NUM_PAPI_EVENTS ];

typedef long long t_counter;

static t_counter counters[NUM_PAPI_EVENTS];

void FNAME(init_papi_counters)
( int *ierr )
{
  
  unsigned i;
  PAPI_event_info_t info;
  
  /* Initialize PAPI library */
  if ( PAPI_library_init( PAPI_VER_CURRENT ) != PAPI_VER_CURRENT ) {
    fprintf( stderr, "(*error*) Unable to initialize PAPI library.\n");
	*ierr = -1;
	return;
  }
  
  /* Initialize PAPI events */
  for( i = 0; i < NUM_PAPI_EVENTS; i ++ ) {
    
    /* Get event code */
    if ( PAPI_event_name_to_code(papi_event_name[i], &papi_events[i]) != PAPI_OK ) {
      fprintf( stderr, "(*error*) Unable to get event code for event %s \n", papi_event_name[i] );
      *ierr = -2;
      return;
    }
    
    /* Verify that event is supported */
    if ( PAPI_get_event_info(papi_events[i], &info) != PAPI_OK ) {
      fprintf( stderr, "(*error*) Event %s is not supported\n", papi_event_name[i] );
      *ierr = -2;
      return;
    }
    
    counters[i] = 0;
  }

  /* Start the counters */
  if ( PAPI_start_counters( papi_events , NUM_PAPI_EVENTS ) != PAPI_OK) {
    fprintf( stderr, "(*error*) Unable to start PAPI counters.\n");
	*ierr = -3;
	return;
  }
  
  /* Initialization completed ok */
  *ierr = 0;
}


void FNAME(get_papi_info)
( t_counter* rcy, t_counter* ucy, t_counter* rus, t_counter* uus, t_counter* hwc )
{
  unsigned i;
  
  *rcy = PAPI_get_real_cyc();
  *ucy = PAPI_get_virt_cyc();
  *rus = PAPI_get_real_usec();
  *uus = PAPI_get_virt_usec();
  
  PAPI_accum_counters( counters, NUM_PAPI_EVENTS );
  
  for( i = 0; i < NUM_PAPI_EVENTS; i ++ ) hwc[i] = counters[i];
  
}

void FNAME(papi_start_counters)
()
{
  PAPI_start_counters( papi_events , NUM_PAPI_EVENTS );
}

void FNAME(papi_stop_counters)
()
{
  unsigned i;
  t_counter temp[NUM_PAPI_EVENTS];
  
  PAPI_stop_counters( temp , NUM_PAPI_EVENTS );
  for( i = 0; i < NUM_PAPI_EVENTS; i ++ ) counters[i] += temp[i];
}


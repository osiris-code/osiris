/*****************************************************************************************

Relativistic particle pusher, SSE optimized version (double precision)

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__


#include "os-spec-push-ssed.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "fortran.h"

#include "vector-sse.h"
#include "splines-sse.h"
#include "os-spec-current-ssed.h"

#include "split-vec.h"

/*****************************************************************************************
Single pass particle advance deposit

This version of the code does the full push in a single pass: field interpolation, du/dt and dx/dt. 
This allows for savings on memory accesses, since particle positions need only be read once and it 
avoids storing temporary field and momenta values. 

Each group of 4 particles is immediatly split to avoid storing old position values in a special
buffer optimized for the current deposition routines.

The routine then deposits the current for all virtual (split) particles.
*****************************************************************************************/


/*****************************************************************************************
vmntrim

Returns the required shift (trim) for particles that must remain in [ -0.5, 0.5 [
to the nearest cell point:
          x >= 0.5 : +1
  -0.5 <= x <  0.5 :  0
  -0.5 <  x        : -1

*****************************************************************************************/

static inline __m128d vmntrim( const __m128d vx )
{
  __m128d va, vb;

  va = _mm_cmplt_pd( vx, _mm_set1_pd( -0.5 ) );
  va = _mm_and_pd(   va, _mm_set1_pd( -1.0 ) );
  
  vb = _mm_cmpge_pd( vx, _mm_set1_pd( +0.5 ) );
  vb = _mm_and_pd(   vb, _mm_set1_pd( +1.0 ) );
  
  return _mm_add_pd( va, vb );
}

/*****************************************************************************************
LOADFLD2

Loads 2 field values corresponding to 2 particle positions into a vector
variable.
*****************************************************************************************/


#define LOADFLD2( fp, shift ) \
  _mm_set_pd( (fp[1])[shift], (fp[0])[shift] )

/*****************************************************************************************
vhp_shift
  
Gets the required h shift and cell index for interpolating quantities on a staggered grid
(quantities located at the half cell position).

  h  = ( dx < 0 ) ? 1.0 : 0.0
  ixh = ix - h * delta

*****************************************************************************************/

#define vhp_shift(dx,ix, h, ixh,delta ) {							                 \
  __m128d cmp  = _mm_cmplt_pd( dx, _mm_setzero_pd());                                    \
  h = _mm_and_pd( cmp, _mm_set1_pd( 1.0 ) );                              \
  ixh  = _mm_sub_epi32( ix, _mm_and_si128(                                               \
               _mm_shuffle_epi32( _mm_castpd_si128( cmp ), _MM_SHUFFLE( 0, 0, 2, 0 ) ),  \
               _mm_set1_epi32( delta )) );                                               \
}

/*****************************************************************************************
vdudt_boris

Use a Boris push to advance particle velocities using interpolated
fields.
*****************************************************************************************/

#define VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene )           \
{                                                                                        \
   __m128d const c1 = _mm_set1_pd( 1.0 ); \
   __m128d vut1, vut2, vut3, vutsq;\
   __m128d vgamma, vtem_gamma, votsq;\
   \
   /* Perform first half of electric field acceleration.*/\
   ve1 = _mm_mul_pd( ve1, vtem );\
   ve2 = _mm_mul_pd( ve2, vtem );\
   ve3 = _mm_mul_pd( ve3, vtem );\
   \
   vut1 = _mm_add_pd( vu1, ve1 );\
   vut2 = _mm_add_pd( vu2, ve2 );\
   vut3 = _mm_add_pd( vu3, ve3 );\
   \
   /* Perform first half of the rotation and store in u. */\
   vutsq = _mm_add_pd( _mm_mul_pd( vut3, vut3 ),\
                          _mm_add_pd( _mm_mul_pd( vut2, vut2 ),\
                                         _mm_mul_pd( vut1, vut1 ) ) );\
   \
   vgamma = _mm_sqrt_pd( _mm_add_pd( vutsq, c1 ) );\
   \
   vene = _mm_div_pd( vutsq, _mm_add_pd( vgamma, c1 ) );\
   \
   vtem_gamma = _mm_div_pd( vtem, vgamma );\
   \
   vb1 = _mm_mul_pd( vb1, vtem_gamma );\
   vb2 = _mm_mul_pd( vb2, vtem_gamma );\
   vb3 = _mm_mul_pd( vb3, vtem_gamma );\
   \
   vu1 = _mm_add_pd( _mm_mul_pd(vb3, vut2), vut1 );\
   vu2 = _mm_add_pd( _mm_mul_pd(vb1, vut3), vut2 );\
   vu3 = _mm_add_pd( _mm_mul_pd(vb2, vut1), vut3 );\
   \
   vu1 = _mm_sub_pd( vu1, _mm_mul_pd(vb2, vut3)  );\
   vu2 = _mm_sub_pd( vu2, _mm_mul_pd(vb3, vut1)  );\
   vu3 = _mm_sub_pd( vu3, _mm_mul_pd(vb1, vut2)  );\
   \
   /* Perform second half of the rotation. */\
  votsq = _mm_add_pd( _mm_mul_pd(vb1, vb1),    c1 );\
  votsq = _mm_add_pd( _mm_mul_pd(vb2, vb2), votsq );\
  votsq = _mm_add_pd( _mm_mul_pd(vb3, vb3), votsq );\
  \
  votsq = _mm_div_pd( _mm_set1_pd( 2.0 ), votsq ) ;\
   \
   vb1 = _mm_mul_pd( vb1, votsq );\
   vb2 = _mm_mul_pd( vb2, votsq );\
   vb3 = _mm_mul_pd( vb3, votsq );\
   \
  vut1 = _mm_add_pd( _mm_mul_pd(vb3, vu2), vut1 );\
  vut2 = _mm_add_pd( _mm_mul_pd(vb1, vu3), vut2 );\
  vut3 = _mm_add_pd( _mm_mul_pd(vb2, vu1), vut3 );\
   \
  vut1 = _mm_sub_pd( vut1, _mm_mul_pd(vb2, vu3) );\
  vut2 = _mm_sub_pd( vut2, _mm_mul_pd(vb3, vu1) );\
  vut3 = _mm_sub_pd( vut3, _mm_mul_pd(vb1, vu2) );\
   \
   /* Perform second half of electric field acceleration.*/\
   vu1 = _mm_add_pd( vut1, ve1 );\
   vu2 = _mm_add_pd( vut2, ve2 );\
   vu3 = _mm_add_pd( vut3, ve3 );\
}

/*****************************************************************************************
VRGAMMA

Get 1 / Lorentz gamma. Using a macro is beneficial because this is 
called from a function that is already inlined.
*****************************************************************************************/

#define VRGAMMA( vrg, vu1, vu2, vu3 ) {  \
  const __m128d c1 = _mm_set1_pd( 1.0 );                    \
  (vrg) = _mm_add_pd( _mm_mul_pd((vu3),(vu3)),                                     \
                        _mm_add_pd(_mm_mul_pd((vu2),(vu2)),                       \
                          _mm_add_pd(_mm_mul_pd((vu1),(vu1)),c1 )));  \
  (vrg) = _mm_div_pd( c1, _mm_sqrt_pd(vrg) );                                      \
}

/****************************************************************************************

  Generate specific particle push functions for 2D and 3D, 1st to 4th order

****************************************************************************************/

#define __TEMPLATE__

// These macros will append _s1, _s2, etc to function names.
// These cannot be used to define the fortran function names because of recursive macro
// expansion issues.
#define ONAME(f, o) OJOIN(f, o)
#define OJOIN(f, o) f ## _s ## o


/********************************** Linear interpolation ********************************/

#define ADVANCE_DEPOSIT_1D     FNAME( vadvance_deposit_1d_s1 )
#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s1 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s1 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s1 )

// Interpolation order
#define ORDER 1
// Number of interpolation points ( ORDER + 1 )
#define NP 2
// Grid offset ( 1 + ORDER/2 )
#define OFFSET 1


#include __FILE__

/******************************** Quadratic interpolation *******************************/

#define ADVANCE_DEPOSIT_1D     FNAME( vadvance_deposit_1d_s2 )
#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s2 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s2 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s2 )

#define ORDER 2
#define NP 3
#define OFFSET 2

#include __FILE__

/********************************** Cubic interpolation *********************************/

#define ADVANCE_DEPOSIT_1D     FNAME( vadvance_deposit_1d_s3 )
#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s3 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s3 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s3 )

#define ORDER 3
#define NP 4
#define OFFSET 2

#include __FILE__

/********************************* Quartic interpolation ********************************/

#define ADVANCE_DEPOSIT_1D     FNAME( vadvance_deposit_1d_s4 )
#define ADVANCE_DEPOSIT_2D     FNAME( vadvance_deposit_2d_s4 )
#define ADVANCE_DEPOSIT_2D_CYL FNAME( vadvance_deposit_2d_cyl_s4 )
#define ADVANCE_DEPOSIT_3D     FNAME( vadvance_deposit_3d_s4 )

#define ORDER 4
#define NP 5
#define OFFSET 3

#include __FILE__

/****************************************************************************************/

#else

/****************************************************************************************

  Template function definitions for 2D and 3D particle push

****************************************************************************************/

/***************   Generate Function names based on interpolation order  ****************/

#define SPLINE   ONAME( vsplined, ORDER )
#define SPLINEH  ONAME( vsplinehd, ORDER )

#define DEP_CURRENT_1D ONAME( vdepcurrent_1d, ORDER )
#define DEP_CURRENT_2D ONAME( vdepcurrent_2d, ORDER )
#define DEP_CURRENT_3D ONAME( vdepcurrent_3d, ORDER )

extern void ADVANCE_DEPOSIT_1D
( int *ix, double *x, double *u, double *q, int *i0, int *i1, double *rqm,
  double *efield, double *bfield, int *emf_size, int *emf_offset,  
  double *j, int *j_size, int *j_offset, 
  double *dx, double *dt, double *ene )
{
  unsigned int k, i, np, np_total;
  
  DECLARE_ALIGNED_16( double jnorm[1] );
  DECLARE_ALIGNED_16( int dix[4] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf1D vpbuf;

  // Simulation constants
  __m128d const vtem    = _mm_set1_pd( 0.5 * (*dt) / (*rqm) );
  __m128d const vdt_dx1 = _mm_set1_pd( (*dt) / dx[0] );

  // get pointers to position 0 of each field component
  unsigned int const deltaX = 3; // 3 field components

  double* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX;
  double* const e2 = e1 + 1;
  double* const e3 = e1 + 2;

  double* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX;
  double* const b2 = b1  + 1;
  double* const b3 = b1  + 2;
  
  
  // Normalization for currents
  jnorm[0] = dx[0] / (*dt);
  
  
  // jump to 1st particle
  x  += (*i0-1);
  ix += (*i0-1);
  u  += 3*(*i0-1);
  q  += (*i0-1);
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {

    for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH; i ++ ) {
      ix[ (np_total + i)     ] = 1;
      x[  (np_total + i)     ] = 0.;
      u[ 3 * (np_total + i)     ] = 0.;
      u[ 3 * (np_total + i) + 1 ] = 0.;
      u[ 3 * (np_total + i) + 2 ] = 0.;
      q[ (np_total + i ) ] = 0.;
    }
    
    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
      printf("(*error*) Number of particles to push is still not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
      exit(1);
    }

  }

  if ( p_cache_size_1D % VEC_WIDTH ) {
    printf("(*error*) p_cache_size_1D is not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
    exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_1D ) {

    // Initialize energy
    __m128d vene_group = _mm_setzero_pd();

      // find number of particles to process
    np = ( np_total - k < p_cache_size_1D ) ? np_total - k : p_cache_size_1D;

    // Initialize virtual particle buffer
    vpbuf.np = 0;

    // Push all particles in group
    for( i = 0; i < np; i += VEC_WIDTH ) {
      __m128d vene;

      __m128d vx0, vx1;
      __m128i vix0, vix1;

      __m128d ve1, ve2, ve3;
      __m128d vb1, vb2, vb3;

      __m128d vu1, vu2, vu3;

      // Load particle positions 
      vx0  = _mm_load_pd( (double const*) x ); 
      vix0 = _mm_loadu_si128( (__m128i const* ) ix );
    
      // Interpolate fields
      {
        __m128d  hx;
        __m128i vidxi, vidxih;
        __m128d vwx[NP], vwxh[NP];

        // idxi = 3 * ix0
        vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );

        SPLINE( vx0, vwx );

        vhp_shift( vx0, vidxi, hx, vidxih, 3 ); 

        SPLINEH( vx0, hx, vwxh );

        // Interpolate E field
        {
          unsigned int k1;

          double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
          DECLARE_ALIGNED_16( int idx1[4] ); // We'll need to store 4 integers
          DECLARE_ALIGNED_16( int idx2[4] );
          DECLARE_ALIGNED_16( int idx3[4] );

          // get pointers to fields 
          _mm_store_si128( (__m128i *) idx1, vidxih );
          _mm_store_si128( (__m128i *) idx2,  vidxi );
          _mm_store_si128( (__m128i *) idx3,  vidxi );

          // This cannot be done vectorially since pointers are 64 bit
          for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
            fld1[k1] = e1 + idx1[k1];
            fld2[k1] = e2 + idx2[k1];
            fld3[k1] = e3 + idx3[k1];
          }

          ve1 = _mm_setzero_pd();
          ve2 = _mm_setzero_pd();
          ve3 = _mm_setzero_pd();

          __m128d f1point[NP], f2point[NP], f3point[NP];

          for ( k1 = 0; k1 < NP; k1++ ) {
            unsigned int shift = k1*3;

            f1point[k1] = LOADFLD2( fld1, shift ); 
            f2point[k1] = LOADFLD2( fld2, shift ); 
            f3point[k1] = LOADFLD2( fld3, shift ); 
          }

          for ( k1 = 0; k1 < NP; k1 ++ ) {
            ve1 = _mm_add_pd( ve1, _mm_mul_pd( f1point[k1], vwxh[k1] ) );
            ve2 = _mm_add_pd( ve2, _mm_mul_pd( f2point[k1],  vwx[k1] ) );
            ve3 = _mm_add_pd( ve3, _mm_mul_pd( f3point[k1],  vwx[k1] ) );
          }

        }     

        // Interpolate B field
        {
          unsigned int k1;

          double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
          DECLARE_ALIGNED_16( int idx1[4] ); // We'll need to store 4 integers
          DECLARE_ALIGNED_16( int idx2[4] );
          DECLARE_ALIGNED_16( int idx3[4] );

          // get pointers to fields 
          _mm_store_si128( (__m128i *) idx1,  vidxi );
          _mm_store_si128( (__m128i *) idx2, vidxih );
          _mm_store_si128( (__m128i *) idx3, vidxih );

          // This cannot be done vectorially since pointers are 64 bit
          for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
            fld1[k1] = b1 + idx1[k1];
            fld2[k1] = b2 + idx2[k1];
            fld3[k1] = b3 + idx3[k1];
          }

          vb1 = _mm_setzero_pd();
          vb2 = _mm_setzero_pd();
          vb3 = _mm_setzero_pd();

          __m128d f1point[NP], f2point[NP], f3point[NP];

          for ( k1 = 0; k1 < NP; k1++ ) {
            unsigned int shift = k1*3;

            f1point[k1] = LOADFLD2( fld1, shift ); 
            f2point[k1] = LOADFLD2( fld2, shift ); 
            f3point[k1] = LOADFLD2( fld3, shift ); 
          }

          for ( k1 = 0; k1 < NP; k1 ++ ) {
            vb1 = _mm_add_pd( vb1, _mm_mul_pd( f1point[k1],  vwx[k1] ) );
            vb2 = _mm_add_pd( vb2, _mm_mul_pd( f2point[k1], vwxh[k1] ) );
            vb3 = _mm_add_pd( vb3, _mm_mul_pd( f3point[k1], vwxh[k1] ) );
          }

        }     
      }

      // load momenta
      _MM_LOAD2v3_PD( vu1, vu2, vu3, u );

      // advance momenta
      VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );

      // Store Results
      _MM_STORE2v3_PD(u, vu1, vu2, vu3); 

      // ---------- advance position and get velocities
      {
        __m128d vrg;

      // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        __m128d vq = _mm_load_pd( q );

        // Accumulate energy
        vene_group = _mm_add_pd( _mm_mul_pd( vq, vene ), vene_group );

        // get velocities
        __m128d vv1 = _mm_mul_pd( vu1, vrg );
        __m128d vv2 = _mm_mul_pd( vu2, vrg ); // this is required for current deposition
        __m128d vv3 = _mm_mul_pd( vu3, vrg ); // this is required for current deposition

        // Push positions
        __m128d vx1 = _mm_add_pd( _mm_mul_pd( vv1, vdt_dx1 ), vx0 );

        // Store virtual particles with positions still indexed to the original cell
        STOREU2P1D( vpbuf, vpbuf.np, vx0, vx1, vq, vv2, vv3, vix0 );

        // Trim positions and store results
        __m128d vtr1 = vmntrim(vx1);

        vx1  = _mm_sub_pd( vx1, vtr1 );
        _mm_store_pd( x, vx1 );

        __m128i vitr1 = _mm_cvttpd_epi32( vtr1 ) ;
        _mm_store_si128( (__m128i *) dix, vitr1 );

        __m128i vix1 = _mm_add_epi32( vix0, vitr1 ); 
        _mm_storeu_si128( (__m128i *) ix, vix1 );
      }

      // ---------- split trajectories for current deposition
      vsplit1D( &vpbuf, dix );

      // ---------- advance pointers
      x  += VEC_WIDTH;  
      ix += VEC_WIDTH;
      u  += 3 * VEC_WIDTH; 
      q  += VEC_WIDTH;
    }

    // Deposit current from all virtual particles
    DEP_CURRENT_1D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    // this is done in double precision
    *ene += _mm_reduce_add_pd( vene_group );

  }

}

/********************************** 2D advance deposit **********************************/

extern void ADVANCE_DEPOSIT_2D
( int *ix, double *x, double *u, double *q, int *i0, int *i1, double *rqm,
  double *efield, double *bfield, int *emf_size, int *emf_offset,  
  double *j, int *j_size, int *j_offset, 
  double *dx, double *dt, double *ene )
{
  unsigned int k, i, np, np_total;
  
  DECLARE_ALIGNED_16( double jnorm[2] );
  
  DECLARE_ALIGNED_16( unsigned int cross[4] );
  DECLARE_ALIGNED_16( int dix[4] );
  DECLARE_ALIGNED_16( int diy[4] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf2D vpbuf;

  // Simulation constants
  __m128d const vtem     = _mm_set1_pd( 0.5 * (*dt) / (*rqm) );
  __m128d const vrdx1_dt = _mm_set1_pd( (*dt) / dx[0] );
  __m128d const vrdx2_dt = _mm_set1_pd( (*dt) / dx[1] );
  
  // get pointers to position 0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components
  unsigned int const deltaY = emf_size[0] * deltaX;

  double* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX + 
                              (emf_offset[1] - OFFSET) * deltaY;
  double* const e2 = e1 + 1;
  double* const e3 = e1 + 2;

  double* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX + 
                              (emf_offset[1] - OFFSET) * deltaY;
  double* const b2 = b1  + 1;
  double* const b3 = b1  + 2;
  
  // Normalization for currents
  jnorm[0] = dx[0] / (*dt) / 2;
  jnorm[1] = dx[1] / (*dt) / 2;
    
  // jump to 1st particle
  x  += 2*(*i0-1);
  ix += 2*(*i0-1);
  u  += 3*(*i0-1);
  q  += *i0-1;
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {
	
	for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH; i ++ ) {
	  ix[ 2 * (np_total + i)     ] = 1;
	  ix[ 2 * (np_total + i) + 1 ] = 1;
	  x[ 2 * (np_total + i)     ] = 0.;
	  x[ 2 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i)     ] = 0.;
	  u[ 3 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i) + 2 ] = 0.;
	  q[ (np_total + i ) ] = 0.;
	}
    
    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
       printf("(*error*) Number of particles to push is still not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
       exit(1);
    }
    
  }
  
  if ( p_cache_size_2D % VEC_WIDTH ) {
       printf("(*error*) p_cache_size_2D is not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
       exit(1);
  }
  
  for ( k = 0; k < np_total; k += p_cache_size_2D ) {

	// Initialize energy
	__m128d vene_group = _mm_setzero_pd();
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;

	// Initialize virtual particle buffer
	vpbuf.np = 0;
		
	// Push all particles in group
	for( i = 0; i < np; i+= VEC_WIDTH ) {
  	  __m128d vene;

	  __m128d vx0, vy0;
	  __m128i vix0, viy0;
  
	  __m128d ve1, ve2, ve3;
	  __m128d vb1, vb2, vb3;
  
	  __m128d vu1, vu2, vu3;
	  
  
	  // Load particle positions 
	  _MM_LOAD2v2_PD( vx0, vy0, x );
	  _MM_LOAD2v2_EPI32( vix0, viy0, ix );
  
	  // Interpolate fields
	  {
		 __m128d hx, hy;
		 __m128i vidxi, vidxj, vidxih, vidxjh;
		 __m128d vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];
	 
	     // idxi = 3 * ix0
	     vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
	     // idxj = deltaY * iy0
	     vidxj = vsmul( viy0, deltaY );
	     	     
		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy);
	 
		 vhp_shift( vx0, vidxi, hx, vidxih, 3 ); 
		 vhp_shift( vy0, vidxj, hy, vidxjh, deltaY ); 
	 
		 SPLINEH( vx0, hx, vwxh );
		 SPLINEH( vy0, hy, vwyh );

		 // Interpolate E field
         {
			unsigned int k1, k2;
			
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[4] ); // We'll need to store 4 integers
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32(  vidxi,  vidxj ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
              fld1[k1] = e1 + idx1[k1];
              fld2[k1] = e2 + idx2[k1];
              fld3[k1] = e3 + idx3[k1];
			}
			
			ve1 = _mm_setzero_pd();
			ve2 = _mm_setzero_pd();
			ve3 = _mm_setzero_pd();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m128d f1line, f2line, f3line;
			  f1line = _mm_setzero_pd();
			  f2line = _mm_setzero_pd();
			  f3line = _mm_setzero_pd();

			  __m128d f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 unsigned int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD2( fld1, shift ); 
				 f2point[k1] = LOADFLD2( fld2, shift ); 
				 f3point[k1] = LOADFLD2( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm_add_pd( f1line, _mm_mul_pd( f1point[k1], vwxh[k1] ) );
				 f2line = _mm_add_pd( f2line, _mm_mul_pd( f2point[k1],  vwx[k1] ) );
				 f3line = _mm_add_pd( f3line, _mm_mul_pd( f3point[k1],  vwx[k1] ) );
			  }

			  ve1 = _mm_add_pd( ve1, _mm_mul_pd( f1line,  vwy[k2] ) );
			  ve2 = _mm_add_pd( ve2, _mm_mul_pd( f2line, vwyh[k2] ) );
			  ve3 = _mm_add_pd( ve3, _mm_mul_pd( f3line,  vwy[k2] ) );
			}
			        
         }		 

		 // Interpolate B field
         {
			unsigned int k1, k2;
			
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[4] ); // We'll need to store 4 integers
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32( vidxih, vidxjh ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < VEC_WIDTH; k1++) {
              fld1[k1] = b1 + idx1[k1];
              fld2[k1] = b2 + idx2[k1];
              fld3[k1] = b3 + idx3[k1];
			}
			
			vb1 = _mm_setzero_pd();
			vb2 = _mm_setzero_pd();
			vb3 = _mm_setzero_pd();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m128d f1line, f2line, f3line;
			  f1line = _mm_setzero_pd();
			  f2line = _mm_setzero_pd();
			  f3line = _mm_setzero_pd();

			  __m128d f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 unsigned int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD2( fld1, shift ); 
				 f2point[k1] = LOADFLD2( fld2, shift ); 
				 f3point[k1] = LOADFLD2( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm_add_pd( f1line, _mm_mul_pd( f1point[k1],  vwx[k1] ) );
				 f2line = _mm_add_pd( f2line, _mm_mul_pd( f2point[k1], vwxh[k1] ) );
				 f3line = _mm_add_pd( f3line, _mm_mul_pd( f3point[k1], vwxh[k1] ) );
			  }

			  vb1 = _mm_add_pd( vb1, _mm_mul_pd( f1line, vwyh[k2] ) );
			  vb2 = _mm_add_pd( vb2, _mm_mul_pd( f2line,  vwy[k2] ) );
			  vb3 = _mm_add_pd( vb3, _mm_mul_pd( f3line, vwyh[k2] ) );
			}
			        
         }		 
	  }

	   
      // Load momenta
      _MM_LOAD2v3_PD( vu1, vu2, vu3, u );

	  // ---------- advance momenta
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );
	  
	  // Store Momenta
	  _MM_STORE2v3_PD(u, vu1, vu2, vu3); 
	  
   
	  // ---------- advance position and get velocities
	  {
		__m128d vrg;
  
        // Get 1 / \gamma
		VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        __m128d vq = _mm_load_pd( q );

		// Accumulate energy
		vene_group = _mm_add_pd( _mm_mul_pd( vq, vene ), vene_group );
	  
		// get velocities
		__m128d vv1 = _mm_mul_pd( vu1, vrg );
		__m128d vv2 = _mm_mul_pd( vu2, vrg );
		__m128d vv3 = _mm_mul_pd( vu3, vrg ); // this is required for current deposition
			
		// Push positions
		__m128d vx1 = _mm_add_pd( _mm_mul_pd( vv1, vrdx1_dt ), vx0 );
		__m128d vy1 = _mm_add_pd( _mm_mul_pd( vv2, vrdx2_dt ), vy0 );

        // Store virtual particles with positions still indexed to the original cell
        STOREU2P2D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );
	  
	    // Trim positions and store results
		__m128d vtr1 = vmntrim(vx1);
		__m128d vtr2 = vmntrim(vy1);

		vx1  = _mm_sub_pd( vx1, vtr1 );
		vy1  = _mm_sub_pd( vy1, vtr2 );
		_MM_STORE2v2_PD( x, vx1, vy1 );
		
		__m128i vitr1 = _mm_cvttpd_epi32( vtr1 ) ;
		__m128i vitr2 = _mm_cvttpd_epi32( vtr2 ) ;
		_mm_store_si128( (__m128i *) dix, vitr1 );
		_mm_store_si128( (__m128i *) diy, vitr2 );
		
		__m128i vix1 = _mm_add_epi32( vix0, vitr1 ); 
		__m128i viy1 = _mm_add_epi32( viy0, vitr2 ); 
		
		_MM_STORE2v2_EPI32( ix, vix1, viy1 );

        // calculate crossings
        _mm_store_si128( (__m128i *) cross, 
                    _mm_add_epi32(_mm_andnot_si128(_mm_cmpeq_epi32(vix0, vix1),_mm_set1_epi32(1)), 
                                  _mm_andnot_si128(_mm_cmpeq_epi32(viy0, viy1),_mm_set1_epi32(2))));
	  }
	  
	  
	  // ---------- split trajectories for current deposition
	  
	  // The new split uses positions indexed to the original cell
      vsplit2D( &vpbuf, cross, dix, diy );

	  
	  // ---------- advance pointers
	  x  += 2 * VEC_WIDTH;  
	  ix += 2 * VEC_WIDTH;
	  u  += 3 * VEC_WIDTH; 
	  q  += VEC_WIDTH;
  
	}
  
	// Deposit current from all virtual particles
	DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );

	// Accumulate energy from group
	*ene += _mm_reduce_add_pd( vene_group );
		
  }
    
}

extern void ADVANCE_DEPOSIT_2D_CYL
( int *ix, double *x, double *u, double *q, int *i0, int *i1, double *rqm,
  int *ilb2 , 
  double *efield, double *bfield, int *emf_size, int *emf_offset,  
  double *j, int *j_size, int *j_offset, 
  double *dx, double *dt, double *ene )
{
  unsigned int k, i, np, np_total;
    
  DECLARE_ALIGNED_16( double jnorm[2] );
  
  DECLARE_ALIGNED_16( unsigned int cross[4] );
  DECLARE_ALIGNED_16( int dix[4] );
  DECLARE_ALIGNED_16( int diy[4] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf2D vpbuf;

  // Simulation Constants
  __m128d const vtem     = _mm_set1_pd( 0.5 * (*dt) / (*rqm) );
  __m128d const vrdx1_dt = _mm_set1_pd( (*dt) / dx[0] );
  
  // get pointers to position 0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components
  unsigned int const deltaY = emf_size[0] * deltaX;

  double* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX + 
                              (emf_offset[1] - OFFSET) * deltaY;
  double* const e2 = e1 + 1;
  double* const e3 = e1 + 2;

  double* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX + 
                              (emf_offset[1] - OFFSET) * deltaY;
  double* const b2 = b1  + 1;
  double* const b3 = b1  + 2;
  
  // Normalization for currents
  jnorm[0] = dx[0] / (*dt) / 2;
  jnorm[1] = dx[1] / (*dt) / 2;

  #if ( ORDER == 1 ) || ( ORDER == 3 )
  __m128d const gshift2 = _mm_set1_pd( (*ilb2 - 2 ) );
  #else
  __m128d const gshift2 = _mm_set1_pd( (*ilb2 - 2 ) - 0.5 );
  #endif
  
  __m128d const dr      = _mm_set1_pd( dx[1] );
  __m128d const rdr     = _mm_set1_pd( 1.0/dx[1] );
  __m128d const vdt     = _mm_set1_pd( *dt );
  
  // jump to 1st particle
  x  += 2*(*i0-1);
  ix += 2*(*i0-1);
  u  += 3*(*i0-1);
  q  += *i0-1;
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {
	
	for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH; i ++ ) {
	  ix[ 2 * (np_total + i)     ] = 1;
	  ix[ 2 * (np_total + i) + 1 ] = 1;
	  x[ 2 * (np_total + i)     ] = 0.;
	  x[ 2 * (np_total + i) + 1 ] = 0.25; // Avoid the axis
	  u[ 3 * (np_total + i)     ] = 0.;
	  u[ 3 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i) + 2 ] = 0.;
	  q[ (np_total + i ) ] = 0.;
	}
    
    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
       printf("(*error*) Number of particles to push is still not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
       exit(1);
    }
    
  }
  
  if ( p_cache_size_2D % VEC_WIDTH ) {
       printf("(*error*) p_cache_size_2D is not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {

	// Initialize energy
	__m128d vene_group = _mm_setzero_pd();
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;

	// Initialize virtual particle buffer
	vpbuf.np = 0;
		
	// Push all particles in group
	for( i = 0; i < np; i+=VEC_WIDTH ) {
 	  __m128d vene;

	  __m128d vx0, vy0;
	  __m128i vix0, viy0;
  
	  __m128d ve1, ve2, ve3;
	  __m128d vb1, vb2, vb3;
  
	  __m128d vu1, vu2, vu3;
	    
	  // Load particle positions 
	  _MM_LOAD2v2_PD( vx0, vy0, x );
	  _MM_LOAD2v2_EPI32( vix0, viy0, ix );
  
	  // Interpolate fields
	  {
		 __m128d hx, hy;
		 __m128i vidxi, vidxj, vidxih, vidxjh;
		 __m128d vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];
	 
	     // idxi = 3 * ix0
	     vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
	     // idxj = deltaY * iy0
	     vidxj = vsmul( viy0, deltaY );
	     	     
		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy);
	 
		 vhp_shift( vx0, vidxi, hx, vidxih, 3 ); 
		 vhp_shift( vy0, vidxj, hy, vidxjh, deltaY ); 
	 
		 SPLINEH( vx0, hx, vwxh );
		 SPLINEH( vy0, hy, vwyh );

		 // Interpolate E field
         {
			unsigned int k1, k2;
			
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[4] ); // We'll need to store 4 integers
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32(  vidxi,  vidxj ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < 2; k1++) {
              fld1[k1] = e1 + idx1[k1];
              fld2[k1] = e2 + idx2[k1];
              fld3[k1] = e3 + idx3[k1];
			}
			
			ve1 = _mm_setzero_pd();
			ve2 = _mm_setzero_pd();
			ve3 = _mm_setzero_pd();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m128d f1line, f2line, f3line;
			  f1line = _mm_setzero_pd();
			  f2line = _mm_setzero_pd();
			  f3line = _mm_setzero_pd();

			  __m128d f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 unsigned int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD2( fld1, shift ); 
				 f2point[k1] = LOADFLD2( fld2, shift ); 
				 f3point[k1] = LOADFLD2( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm_add_pd( f1line, _mm_mul_pd( f1point[k1], vwxh[k1] ) );
				 f2line = _mm_add_pd( f2line, _mm_mul_pd( f2point[k1],  vwx[k1] ) );
				 f3line = _mm_add_pd( f3line, _mm_mul_pd( f3point[k1],  vwx[k1] ) );
			  }

			  ve1 = _mm_add_pd( ve1, _mm_mul_pd( f1line,  vwy[k2] ) );
			  ve2 = _mm_add_pd( ve2, _mm_mul_pd( f2line, vwyh[k2] ) );
			  ve3 = _mm_add_pd( ve3, _mm_mul_pd( f3line,  vwy[k2] ) );
			}
			        
         }		 

		 // Interpolate B field
         {
			unsigned int k1, k2;
			
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[4] ); // We'll need to store 4 integers
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, _mm_add_epi32(  vidxi, vidxjh ) );
			_mm_store_si128( (__m128i *) idx2, _mm_add_epi32( vidxih,  vidxj ) );
			_mm_store_si128( (__m128i *) idx3, _mm_add_epi32( vidxih, vidxjh ) );
            
            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < 2; k1++) {
              fld1[k1] = b1 + idx1[k1];
              fld2[k1] = b2 + idx2[k1];
              fld3[k1] = b3 + idx3[k1];
			}
			
			vb1 = _mm_setzero_pd();
			vb2 = _mm_setzero_pd();
			vb3 = _mm_setzero_pd();
		    		    
			for ( k2 = 0; k2 < NP; k2++ ) {
			  
			  register __m128d f1line, f2line, f3line;
			  f1line = _mm_setzero_pd();
			  f2line = _mm_setzero_pd();
			  f3line = _mm_setzero_pd();

			  __m128d f1point[NP], f2point[NP], f3point[NP];

			  for ( k1 = 0; k1 < NP; k1++ ) {
				 unsigned int shift = k1*3 + k2*deltaY;
				 
				 f1point[k1] = LOADFLD2( fld1, shift ); 
				 f2point[k1] = LOADFLD2( fld2, shift ); 
				 f3point[k1] = LOADFLD2( fld3, shift ); 
              }
                       
              for ( k1 = 0; k1 < NP; k1 ++ ) {
				 f1line = _mm_add_pd( f1line, _mm_mul_pd( f1point[k1],  vwx[k1] ) );
				 f2line = _mm_add_pd( f2line, _mm_mul_pd( f2point[k1], vwxh[k1] ) );
				 f3line = _mm_add_pd( f3line, _mm_mul_pd( f3point[k1], vwxh[k1] ) );
			  }

			  vb1 = _mm_add_pd( vb1, _mm_mul_pd( f1line, vwyh[k2] ) );
			  vb2 = _mm_add_pd( vb2, _mm_mul_pd( f2line,  vwy[k2] ) );
			  vb3 = _mm_add_pd( vb3, _mm_mul_pd( f3line, vwyh[k2] ) );
			}
			        
         }		 
	  }

	   
      // Load momenta
      _MM_LOAD2v3_PD( vu1, vu2, vu3, u );

	  // ---------- advance momenta
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );
	     
	  // ---------- advance position in cylindrical geometry
	  {
		__m128d vrg;
		__m128d vx1, vy1;

                  // Get 1 / \gamma
		VRGAMMA( vrg, vu1, vu2, vu3 );

        __m128d vq = _mm_load_pd( q );

		// Accumulate energy
		vene_group = _mm_add_pd( _mm_mul_pd( vq, vene ), vene_group );
	  
		// get velocities
		__m128d vv1 = _mm_mul_pd( vu1, vrg );
		__m128d vv2 = _mm_mul_pd( vu2, vrg );
		__m128d vv3 = _mm_mul_pd( vu3, vrg ); // this is required for current deposition
			
		// Push positions
		vx1 = _mm_add_pd( _mm_mul_pd( vv1, vrdx1_dt ), vx0 );
        
        // x2 - 3D like push
        // This further changes u2 and u3
        {  
		   
		   __m128d gix2, tmp;
		   __m128d r_old, r_new;
		   __m128d x2_new, x3_new; 
		   
		   // gshift2 = shift_ix2 - 0.5d0
   
		   // gix2   = viy0 + gshift2;          // global cell
		   gix2   = _mm_add_pd( _mm_cvtepi32_pd( viy0 ), gshift2 );
		   r_old  = _mm_mul_pd( _mm_add_pd( vy0, gix2 ), dr );

		   x2_new = _mm_add_pd( _mm_mul_pd( vv2, vdt ), r_old );
		   x3_new =                    _mm_mul_pd( vv3, vdt ) ;

		   r_new  = _mm_sqrt_pd( _mm_add_pd( _mm_mul_pd( x2_new, x2_new ),
											 _mm_mul_pd( x3_new, x3_new ) ) ); 
		   
		   // This is a protection against roundoff for cold plasmas
		   // vy1 = ( r_old == r_new ) ? vy0 : r_new * rdr - gix2);
		   vy1 = _mm_blendv_pd( _mm_sub_pd( _mm_mul_pd( r_new, rdr ), gix2 ), vy0, 
		                        _mm_cmpeq_pd( r_old, r_new ) );
		   		   
		   // Correct p_r and p_\theta to conserve angular momentum
		   tmp  = _mm_div_pd( _mm_set1_pd(1.0), r_new );
		   vu2  = _mm_mul_pd( _mm_add_pd( _mm_mul_pd( vu2, x2_new ), 
										  _mm_mul_pd( vu3, x3_new ) ), 
										  tmp );
		   vu3  = _mm_mul_pd( vu3 , _mm_mul_pd( r_old, tmp ) );
   
		}        
        
		// Store Momenta
		_MM_STORE2v3_PD(u, vu1, vu2, vu3) 
        
        // Store virtual particles with positions still indexed to the original cell
        STOREU2P2D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );
	  
	    // Trim positions and store results
		__m128d vtr1 = vmntrim(vx1);
		__m128d vtr2 = vmntrim(vy1);

		vx1  = _mm_sub_pd( vx1, vtr1 );
		vy1  = _mm_sub_pd( vy1, vtr2 );
		_MM_STORE2v2_PD( x, vx1, vy1 );
		
		__m128i vitr1 = _mm_cvttpd_epi32( vtr1 ) ;
		__m128i vitr2 = _mm_cvttpd_epi32( vtr2 ) ;
		_mm_store_si128( (__m128i *) dix, vitr1 );
		_mm_store_si128( (__m128i *) diy, vitr2 );
		
		__m128i vix1 = _mm_add_epi32( vix0, vitr1 ); 
		__m128i viy1 = _mm_add_epi32( viy0, vitr2 ); 
		_MM_STORE2v2_EPI32( ix, vix1, viy1 );

        // calculate crossings
        _mm_store_si128( (__m128i *) cross, 
                    _mm_add_epi32(_mm_andnot_si128(_mm_cmpeq_epi32(vix0, vix1),_mm_set1_epi32(1)), 
                                  _mm_andnot_si128(_mm_cmpeq_epi32(viy0, viy1),_mm_set1_epi32(2))));
	  }
	  
	  
	  // ---------- split trajectories for current deposition
	  
	  // The new split uses positions indexed to the original cell
      vsplit2D( &vpbuf, cross, dix, diy );

	  
	  // ---------- advance pointers
	  x  += 2 * VEC_WIDTH;  
	  ix += 2 * VEC_WIDTH;
	  u  += 3 * VEC_WIDTH; 
	  q  += VEC_WIDTH;
  
	}
  
	// Deposit current from all virtual particles
	DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );

	// Accumulate energy from group
	*ene += _mm_reduce_add_pd( vene_group );
		
  }
    
}

/********************************** 3D advance deposit **********************************/

extern void ADVANCE_DEPOSIT_3D
( int *ix, double *x, double *u, double *q, int *i0, int *i1, double *rqm,
  double *efield, double *bfield, int *emf_size, int *emf_offset,  
  double *j, int *j_size, int *j_offset, 
  double *dx, double *dt, double *ene )
{
  unsigned int k, i, np, np_total;
    
  DECLARE_ALIGNED_16( double jnorm[3] );

  DECLARE_ALIGNED_16( unsigned int cross[4] );
  DECLARE_ALIGNED_16( int dix[4] );
  DECLARE_ALIGNED_16( int diy[4] );
  DECLARE_ALIGNED_16( int diz[4] );
  
  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf3D vpbuf;

  // Simulation constants
  __m128d const vtem     = _mm_set1_pd( 0.5 * (*dt) / (*rqm) );
  __m128d const vrdx1_dt = _mm_set1_pd( (*dt) / dx[0] );
  __m128d const vrdx2_dt = _mm_set1_pd( (*dt) / dx[1] );
  __m128d const vrdx3_dt = _mm_set1_pd( (*dt) / dx[2] );
  
  // get pointers to position 0,0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components
  unsigned int const deltaY = emf_size[0] * deltaX;
  unsigned int const deltaZ = emf_size[1] * deltaY;

  double* const e1 = efield + (emf_offset[0] - OFFSET)*deltaX +
                              (emf_offset[1] - OFFSET)*deltaY + 
                              (emf_offset[2] - OFFSET)*deltaZ;
  double* const e2 = e1 + 1;
  double* const e3 = e1 + 2;

  double* const b1 = bfield + (emf_offset[0] - OFFSET)*deltaX + 
                              (emf_offset[1] - OFFSET)*deltaY + 
                              (emf_offset[2] - OFFSET)*deltaZ;
  double* const b2 = b1  + 1;
  double* const b3 = b1  + 2;
    
  // Normalization for currents
  jnorm[0] = dx[0] / (*dt) / 3;
  jnorm[1] = dx[1] / (*dt) / 3;
  jnorm[2] = dx[2] / (*dt) / 3;
  
  // jump to 1st particle
  x  += 3*(*i0-1);
  ix += 3*(*i0-1);
  u  += 3*(*i0-1);
  q  += *i0-1;
  
  
  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of 4 add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {
	

	for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH; i ++ ) {
	  ix[ 3 * (np_total + i)     ] = 1;
	  ix[ 3 * (np_total + i) + 1 ] = 1;
	  ix[ 3 * (np_total + i) + 2 ] = 1;
	  x[ 3 * (np_total + i)     ] = 0.;
	  x[ 3 * (np_total + i) + 1 ] = 0.;
	  x[ 3 * (np_total + i) + 2 ] = 0.;
	  u[ 3 * (np_total + i)     ] = 0.;
	  u[ 3 * (np_total + i) + 1 ] = 0.;
	  u[ 3 * (np_total + i) + 2 ] = 0.;
	  q[ (np_total + i ) ] = 0.;
	}
    
    np_total += ( VEC_WIDTH - np_total%VEC_WIDTH );
    if ( np_total % VEC_WIDTH ) {
       printf("(*error*) Number of particles to push is still not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
       exit(1);
    }
    
  }
  
  if ( p_cache_size_3D % VEC_WIDTH ) {
       printf("(*error*) p_cache_size_3D is not a multiple of VECWIDTH = %d!\n", VEC_WIDTH );
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_3D ) {

	// Initialize energy
	__m128d vene_group = _mm_setzero_pd();
        
    // find number of particles to process
    np = ( np_total - k < p_cache_size_3D ) ? np_total - k : p_cache_size_3D;
    
	// Initialize virtual particle buffer
	vpbuf.np = 0;
		
	// Push all particles in group
	for( i = 0; i < np; i+=VEC_WIDTH  ) {

      __m128d vene; 

	  __m128d vx0, vy0, vz0;
	  __m128i vix0, viy0, viz0;
  
	  __m128d ve1, ve2, ve3;
	  __m128d vb1, vb2, vb3;
  
	  __m128d vu1, vu2, vu3;
	    
	  // Load particle positions 
	  _MM_LOAD2v3_PD( vx0, vy0, vz0, x );
	  _MM_LOAD2v3_EPI32( vix0, viy0, viz0, ix  );
  
	  // Interpolate fields
	  {
		 __m128d hx, hy, hz;
		 __m128i vidxi, vidxj, vidxk, vidxih, vidxjh, vidxkh;
		 __m128d vwx[NP], vwy[NP], vwz[NP], vwxh[NP], vwyh[NP], vwzh[NP];

	     // idxi = 3 * ix0
	     vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
	     // idxj = deltaY * iy0
	     vidxj = vsmul( viy0, deltaY );
	     // idxk = deltaZ * iz0
	     vidxk = vsmul( viz0, deltaZ );
	 
		 SPLINE( vx0, vwx );
		 SPLINE( vy0, vwy );
		 SPLINE( vz0, vwz );
	 
		 vhp_shift( vx0, vidxi, hx, vidxih, 3 ); 
		 vhp_shift( vy0, vidxj, hy, vidxjh, deltaY ); 
		 vhp_shift( vz0, vidxk, hz, vidxkh, deltaZ ); 
	 
		 SPLINEH( vx0, hx, vwxh );
		 SPLINEH( vy0, hy, vwyh );
		 SPLINEH( vz0, hz, vwzh );
	 
         // Interpolate E field
         {
			unsigned int k1, k2, k3;
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[4] ); // We'll need to store 4 integers
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, 
			                  _mm_add_epi32( _mm_add_epi32( vidxih, vidxj ), vidxk ) );
			_mm_store_si128( (__m128i *) idx2, 
			                  _mm_add_epi32( _mm_add_epi32( vidxi, vidxjh ), vidxk ) );
			_mm_store_si128( (__m128i *) idx3, 
			                  _mm_add_epi32( _mm_add_epi32( vidxi, vidxj ), vidxkh ) );

            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < 2; k1++) {
              fld1[k1] = e1 + idx1[k1];
              fld2[k1] = e2 + idx2[k1];
              fld3[k1] = e3 + idx3[k1];
			}
			
			ve1 = _mm_setzero_pd();
			ve2 = _mm_setzero_pd();
			ve3 = _mm_setzero_pd();
		  
			for ( k3 = 0; k3 < NP; k3++ ) {
			   register __m128d f1plane, f2plane, f3plane;
			   f1plane = _mm_setzero_pd();
			   f2plane = _mm_setzero_pd();
			   f3plane = _mm_setzero_pd();
				 
			   for ( k2 = 0; k2 < NP; k2++ ) {
				 
				 register __m128d f1line, f2line, f3line;
				 f1line = _mm_setzero_pd();
				 f2line = _mm_setzero_pd();
				 f3line = _mm_setzero_pd();

				 __m128d f1point[NP], f2point[NP], f3point[NP];
				 for ( k1 = 0; k1 < NP; k1++ ) {
					unsigned int shift = k1*3 + k2*deltaY + k3*deltaZ;
					
					f1point[k1] = LOADFLD2( fld1, shift ); 
					f2point[k1] = LOADFLD2( fld2, shift ); 
					f3point[k1] = LOADFLD2( fld3, shift ); 
				 }
				 
				 for ( k1 = 0; k1 < NP; k1 ++ ) {
					f1line = _mm_add_pd( f1line, _mm_mul_pd( f1point[k1], vwxh[k1] ) );
					f2line = _mm_add_pd( f2line, _mm_mul_pd( f2point[k1],  vwx[k1] ) );
					f3line = _mm_add_pd( f3line, _mm_mul_pd( f3point[k1],  vwx[k1] ) );
				 }
				 
			     f1plane = _mm_add_pd( f1plane, _mm_mul_pd( f1line,  vwy[k2] ) );
			     f2plane = _mm_add_pd( f2plane, _mm_mul_pd( f2line, vwyh[k2] ) );
			     f3plane = _mm_add_pd( f3plane, _mm_mul_pd( f3line,  vwy[k2] ) );
			   }
			   
			   ve1 = _mm_add_pd( ve1, _mm_mul_pd( f1plane,  vwz[k3] ) );
			   ve2 = _mm_add_pd( ve2, _mm_mul_pd( f2plane,  vwz[k3] ) );
			   ve3 = _mm_add_pd( ve3, _mm_mul_pd( f3plane, vwzh[k3] ) );
			}
         
         }		 

         // Interpolate B field
	     {
			unsigned int k1, k2, k3;
			double *fld1[VEC_WIDTH], *fld2[VEC_WIDTH], *fld3[VEC_WIDTH];
			DECLARE_ALIGNED_16( int idx1[4] ); // We'll need to store 4 integers
			DECLARE_ALIGNED_16( int idx2[4] );
			DECLARE_ALIGNED_16( int idx3[4] );
			
			// get pointers to fields 
			_mm_store_si128( (__m128i *) idx1, 
			                  _mm_add_epi32( _mm_add_epi32( vidxi, vidxjh ), vidxkh ) );
			_mm_store_si128( (__m128i *) idx2, 
			                  _mm_add_epi32( _mm_add_epi32( vidxih, vidxj ), vidxkh ) );
			_mm_store_si128( (__m128i *) idx3, 
			                  _mm_add_epi32( _mm_add_epi32( vidxih, vidxjh ), vidxk ) );

            // This cannot be done vectorially since pointers are 64 bit
            for ( k1 = 0; k1 < 2; k1++) {
              fld1[k1] = b1 + idx1[k1];
              fld2[k1] = b2 + idx2[k1];
              fld3[k1] = b3 + idx3[k1];
			}
			
			vb1 = _mm_setzero_pd();
			vb2 = _mm_setzero_pd();
			vb3 = _mm_setzero_pd();
		  
			for ( k3 = 0; k3 < NP; k3++ ) {
			   register __m128d f1plane, f2plane, f3plane;
			   f1plane = _mm_setzero_pd();
			   f2plane = _mm_setzero_pd();
			   f3plane = _mm_setzero_pd();
				 
			   for ( k2 = 0; k2 < NP; k2++ ) {
				 
				 register __m128d f1line, f2line, f3line;
				 f1line = _mm_setzero_pd();
				 f2line = _mm_setzero_pd();
				 f3line = _mm_setzero_pd();

				 __m128d f1point[NP], f2point[NP], f3point[NP];
				 for ( k1 = 0; k1 < NP; k1++ ) {
					unsigned int shift = k1*3 + k2*deltaY + k3*deltaZ;
					
					f1point[k1] = LOADFLD2( fld1, shift ); 
					f2point[k1] = LOADFLD2( fld2, shift ); 
					f3point[k1] = LOADFLD2( fld3, shift ); 
				 }
				 
				 for ( k1 = 0; k1 < NP; k1 ++ ) {
					f1line = _mm_add_pd( f1line, _mm_mul_pd( f1point[k1],  vwx[k1] ) );
					f2line = _mm_add_pd( f2line, _mm_mul_pd( f2point[k1], vwxh[k1] ) );
					f3line = _mm_add_pd( f3line, _mm_mul_pd( f3point[k1], vwxh[k1] ) );
				 }
				 
			     f1plane = _mm_add_pd( f1plane, _mm_mul_pd( f1line, vwyh[k2] ) );
			     f2plane = _mm_add_pd( f2plane, _mm_mul_pd( f2line,  vwy[k2] ) );
			     f3plane = _mm_add_pd( f3plane, _mm_mul_pd( f3line, vwyh[k2] ) );
			   }
			   
			   vb1 = _mm_add_pd( vb1, _mm_mul_pd( f1plane, vwzh[k3] ) );
			   vb2 = _mm_add_pd( vb2, _mm_mul_pd( f2plane, vwzh[k3] ) );
			   vb3 = _mm_add_pd( vb3, _mm_mul_pd( f3plane,  vwz[k3] ) );
			}
         
         }		 
		 
	  }	  
	  
      // Load momenta
      _MM_LOAD2v3_PD( vu1, vu2, vu3, u );

	  // ---------- advance momenta
	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );
	  
	  // Store Momenta
	  _MM_STORE2v3_PD(u, vu1, vu2, vu3);
   
	  // ---------- advance position and get velocities
	  {
		__m128d vrg;
  
          // Get 1 / \gamma
		VRGAMMA( vrg, vu1, vu2, vu3 );

	  	          // Load charge
		__m128d vq = _mm_load_pd( q );

				// Accumulate energy
		vene_group = _mm_add_pd( _mm_mul_pd( vq, vene ), vene_group );
	  
        // Push positions
		__m128d vx1 = _mm_add_pd( vx0, _mm_mul_pd( _mm_mul_pd( vu1, vrg ), vrdx1_dt ) );
		__m128d vy1 = _mm_add_pd( vy0, _mm_mul_pd( _mm_mul_pd( vu2, vrg ), vrdx2_dt ) );
		__m128d vz1 = _mm_add_pd( vz0, _mm_mul_pd( _mm_mul_pd( vu3, vrg ), vrdx3_dt ) );
	  
	          // Store virtual particles with positions still indexed to the original cell
		STOREU2P3D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix0, viy0, viz0 );

       // Find position trim
		__m128d vtrx = vmntrim(vx1);
		__m128d vtry = vmntrim(vy1);
		__m128d vtrz = vmntrim(vz1);
				
        // Trim positions and store results
		vx1   = _mm_sub_pd( vx1, vtrx );
		vy1   = _mm_sub_pd( vy1, vtry );
		vz1   = _mm_sub_pd( vz1, vtrz );
        _MM_STORE2v3_PD( x, vx1, vy1, vz1 );

       // find cell crossings and store
		__m128i vitrx = _mm_cvttpd_epi32( vtrx );
		__m128i vitry = _mm_cvttpd_epi32( vtry );
		__m128i vitrz = _mm_cvttpd_epi32( vtrz );
		_mm_store_si128( (__m128i *) dix, vitrx );
		_mm_store_si128( (__m128i *) diy, vitry );
		_mm_store_si128( (__m128i *) diz, vitrz );
		
        // Trim cell indexes and store
		__m128i vix1 = _mm_add_epi32( vix0, vitrx ); 
		__m128i viy1 = _mm_add_epi32( viy0, vitry ); 
		__m128i viz1 = _mm_add_epi32( viz0, vitrz ); 
		_MM_STORE2v3_EPI32( ix, vix1, viy1, viz1 );

		// use vitrx to save registers
		vitrx =                       _mm_andnot_si128(_mm_cmpeq_epi32(vix0, vix1),_mm_set1_epi32(1));
		vitrx = _mm_add_epi32( vitrx, _mm_andnot_si128(_mm_cmpeq_epi32(viy0, viy1),_mm_set1_epi32(2)));
		vitrx = _mm_add_epi32( vitrx, _mm_andnot_si128(_mm_cmpeq_epi32(viz0, viz1),_mm_set1_epi32(4)));
		_mm_store_si128( (__m128i *) cross, vitrx );                  

	  }
	  
	  
	  // ---------- split trajectories for current deposition

      vsplit3D( &vpbuf, cross, dix, diy, diz );

	  
	  // ---------- advance pointers
	  x  += 3*VEC_WIDTH;  
	  ix += 3*VEC_WIDTH;
	  u  += 3*VEC_WIDTH; 
	  q  += VEC_WIDTH;
  
	}

	// Deposit current from all virtual particles
	DEP_CURRENT_3D( j, j_size, j_offset, jnorm, &vpbuf );

	// Accumulate energy from group
	*ene += _mm_reduce_add_pd( vene_group );
		
  }
    
}


/****************************************************************************************
   Clear template definitions
****************************************************************************************/

#undef ORDER
#undef ADVANCE_DEPOSIT_1D
#undef ADVANCE_DEPOSIT_2D
#undef ADVANCE_DEPOSIT_2D_CYL
#undef ADVANCE_DEPOSIT_3D
#undef OFFSET
#undef SPLINE
#undef SPLINEH
#undef NP
#undef DEP_CURRENT_1D
#undef DEP_CURRENT_2D
#undef DEP_CURRENT_3D

#endif

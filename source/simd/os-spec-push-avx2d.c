/*****************************************************************************************

Relativistic particle pusher, AVX2 optimized version (double precision)

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__


#include "os-spec-push-avx2d.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "fortran.h"

#include "vector-avx2.h"
#include "splines-avx2.h"
#include "os-spec-current-avx2d.h"

#include "split-vec.h"


/*****************************************************************************************
Single pass particle advance deposit

This version of the code does the full push in a single pass: field interpolation, du/dt and dx/dt.
This allows for savings on memory accesses, since particle positions need only be read once and it
avoids storing temporary field and momenta values.

Each group of VEC_WIDTH particles is immediatly split to avoid storing old position values in a special
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

static inline __m256d vmntrim( const __m256d vx )
{
  __m256d va, vb;

  va = _mm256_cmp_pd( vx, _mm256_set1_pd( -0.5 ), _CMP_LT_OS );
  va = _mm256_and_pd( va, _mm256_set1_pd( +1.0 ) );

  vb = _mm256_cmp_pd( vx, _mm256_set1_pd( +0.5 ), _CMP_GE_OS );
  vb = _mm256_and_pd( vb, _mm256_set1_pd( +1.0 ) );

  return  _mm256_sub_pd( vb, va );
}

/*****************************************************************************************
vhp_shift

Gets the required h shift and cell index for interpolating quantities on a staggered grid
(quantities located at the half cell position).

  h  = ( dx < 0 ) ? 1.0 : 0.0
  ixh = ix - h * delta

*****************************************************************************************/


#define VHP_SHIFT(dx,ix, h, ixh,delta ) {                              \
   __m256d cmp  = _mm256_cmp_pd( (dx), _mm256_setzero_pd(), _CMP_LT_OS ); \
   (h) =  _mm256_and_pd( cmp, _mm256_set1_pd( 1.0 ) ); \
   cmp = _mm256_castps_pd( _mm256_permute_ps( _mm256_castpd_ps( cmp ), 0x08 ) ); \
   __m128i cmpi = _mm_castpd_si128( \
        _mm_unpacklo_pd( _mm256_castpd256_pd128( cmp  ), \
                 _mm256_extractf128_pd( cmp , 1 ) ) ); \
   (ixh)  = _mm_sub_epi32( (ix), _mm_and_si128( cmpi,_mm_set1_epi32( delta )) ); \
}


/*****************************************************************************************
vdudt_boris

Use a Boris push to advance particle velocities using interpolated
fields.
*****************************************************************************************/

// This version uses full precision division and square root

#define VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene ) \
{\
  __m256d const c1 = _mm256_set1_pd( 1.0 );\
  __m256d vut1, vut2, vut3, vutsq;\
  __m256d vgamma, vtem_gamma, votsq;\
\
  ve1 = _mm256_mul_pd( ve1, vtem );\
  ve2 = _mm256_mul_pd( ve2, vtem );\
  ve3 = _mm256_mul_pd( ve3, vtem );\
\
  vut1 = _mm256_add_pd( vu1, ve1 );\
  vut2 = _mm256_add_pd( vu2, ve2 );\
  vut3 = _mm256_add_pd( vu3, ve3 );\
\
  vutsq = _mm256_fmadd_pd( vut3, vut3,\
         _mm256_fmadd_pd( vut2, vut2, _mm256_mul_pd( vut1, vut1 ) ) );\
\
  vgamma = _mm256_sqrt_pd( _mm256_add_pd( vutsq, c1 ) );\
\
  vene = _mm256_div_pd( vutsq, _mm256_add_pd( vgamma, c1) );\
\
  vtem_gamma = _mm256_div_pd( vtem, vgamma );\
\
  vb1 = _mm256_mul_pd( vb1, vtem_gamma );\
  vb2 = _mm256_mul_pd( vb2, vtem_gamma );\
  vb3 = _mm256_mul_pd( vb3, vtem_gamma );\
\
  vu1 = _mm256_fmadd_pd( vb3, vut2, vut1 );\
  vu2 = _mm256_fmadd_pd( vb1, vut3, vut2 );\
  vu3 = _mm256_fmadd_pd( vb2, vut1, vut3 );\
\
  vu1 = _mm256_fnmadd_pd( vb2, vut3, vu1 );\
  vu2 = _mm256_fnmadd_pd( vb3, vut1, vu2 );\
  vu3 = _mm256_fnmadd_pd( vb1, vut2, vu3 );\
\
  votsq = _mm256_fmadd_pd( vb1, vb1,    c1 );\
  votsq = _mm256_fmadd_pd( vb2, vb2, votsq );\
  votsq = _mm256_fmadd_pd( vb3, vb3, votsq );\
\
  votsq = _mm256_div_pd( _mm256_set1_pd( 2.0 ), votsq ) ;\
\
  vb1 = _mm256_mul_pd( vb1, votsq );\
  vb2 = _mm256_mul_pd( vb2, votsq );\
  vb3 = _mm256_mul_pd( vb3, votsq );\
\
  vut1 = _mm256_fmadd_pd( vb3, vu2, vut1 );\
  vut2 = _mm256_fmadd_pd( vb1, vu3, vut2 );\
  vut3 = _mm256_fmadd_pd( vb2, vu1, vut3 );\
\
  vut1 = _mm256_fnmadd_pd( vb2, vu3, vut1 );\
  vut2 = _mm256_fnmadd_pd( vb3, vu1, vut2 );\
  vut3 = _mm256_fnmadd_pd( vb1, vu2, vut3 );\
\
  vu1 = _mm256_add_pd( vut1, ve1 );\
  vu2 = _mm256_add_pd( vut2, ve2 );\
  vu3 = _mm256_add_pd( vut3, ve3 );\
}

/*****************************************************************************************
VRGAMMA

Get 1 / Lorentz gamma.
*****************************************************************************************/


#define VRGAMMA( vrg, vu1, vu2, vu3 ) {                            		                  \
   const __m256d c1 = _mm256_set1_pd( 1.0 );                       		                  \
   (vrg) = _mm256_div_pd( c1, _mm256_sqrt_pd( _mm256_fmadd_pd( (vu3), (vu3),              \
                                              _mm256_fmadd_pd( (vu2), (vu2),              \
                                              _mm256_fmadd_pd( (vu1), (vu1), c1 ) ) ) ) );\
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

#define SPLINE  ONAME( vsplined, ORDER )
#define SPLINEH  ONAME( vsplinehd, ORDER )

#define DEP_CURRENT_1D ONAME( vdepcurrent_1d, ORDER )
#define DEP_CURRENT_2D ONAME( vdepcurrent_2d, ORDER )
#define DEP_CURRENT_3D ONAME( vdepcurrent_3d, ORDER )

/********************************** 1D advance deposit **********************************/


extern void ADVANCE_DEPOSIT_1D
( int *ix, double *x, double *u, double *q, int *i0, int *i1, double *rqm,
  double *efield, double *bfield, int *emf_size, int *emf_offset,
  double *j, int *j_size, int *j_offset,
  double *dx, double *dt, double *ene )
{

  int k, i, np, np_total;

  DECLARE_ALIGNED_32( double jnorm[1] );
  DECLARE_ALIGNED_32( int dix[VEC_WIDTH] );

  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf1D vpbuf;

  // Simulation constants
  __m256d const vtem    = _mm256_set1_pd( 0.5 * (*dt) / (*rqm) );
  __m256d const vdt_dx1 = _mm256_set1_pd( (*dt) / dx[0] );

  // get pointers to position 0 of each field component
  unsigned int const deltaX = 3; // 3 field components

  double* const e1 = efield + (emf_offset[0] - OFFSET ) * deltaX;
  double* const e2 = e1 + 1;
  double* const e3 = e1 + 2;

  double* const b1 = bfield + (emf_offset[0] - OFFSET ) * deltaX;
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

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {

    for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH; i ++ ) {
      ix[ (np_total + i)        ] = 1;
      x[  (np_total + i)        ] = 0.;
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
    __m256d vene_group = _mm256_setzero_pd();

    // find number of particles to process
    np = ( np_total - k < p_cache_size_1D ) ? np_total - k : p_cache_size_1D;

    // Initialize virtual particle buffer
    vpbuf.np = 0;

    // Push all particles in group
    for( i = 0; i < np; i+=VEC_WIDTH ) {
      __m256d vene;
      __m256d vx0;
      __m128i vix0;

      __m256d ve1, ve2, ve3;
      __m256d vb1, vb2, vb3;

      __m256d vu1, vu2, vu3;

      // Load particle positions
      vx0  = _mm256_load_pd( (double const*) x );
      vix0 = _mm_load_si128( (__m128i const* ) ix );

      // Interpolate fields
      {
        __m256d  hx;
        __m128i vidxi, vidxih;
        __m256d vwx[NP], vwxh[NP];

        // idxi = 3 * ix0
        vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );

        SPLINE( vx0, vwx );

        VHP_SHIFT( vx0, vidxi, hx, vidxih, 3 );

        SPLINEH( vx0, hx, vwxh );

        // Interpolate E field
        {
          unsigned int k1;

          ve1 = _mm256_setzero_pd();
          ve2 = _mm256_setzero_pd();
          ve3 = _mm256_setzero_pd();

          for ( k1 = 0; k1 < NP; k1 ++ ) {
            unsigned shift = k1*3;
            ve1 = _mm256_fmadd_pd( _mm256_i32gather_pd( e1 + shift, vidxih, 8 ), vwxh[k1], ve1 );
            ve2 = _mm256_fmadd_pd( _mm256_i32gather_pd( e2 + shift, vidxi , 8 ),  vwx[k1], ve2 );
            ve3 = _mm256_fmadd_pd( _mm256_i32gather_pd( e3 + shift, vidxi , 8 ),  vwx[k1], ve3 );
          }

        }

        // Interpolate B field
        {
          unsigned int k1;

          vb1 = _mm256_setzero_pd();
          vb2 = _mm256_setzero_pd();
          vb3 = _mm256_setzero_pd();

          for ( k1 = 0; k1 < NP; k1 ++ ) {
            unsigned shift = k1*3;
            vb1 = _mm256_fmadd_pd( _mm256_i32gather_pd( b1 + shift,  vidxi, 8 ),  vwx[k1], vb1 );
            vb2 = _mm256_fmadd_pd( _mm256_i32gather_pd( b2 + shift, vidxih, 8 ), vwxh[k1], vb2 );
            vb3 = _mm256_fmadd_pd( _mm256_i32gather_pd( b3 + shift, vidxih, 8 ), vwxh[k1], vb3 );
          }

        }

      }

      // Load momenta
      _MM256_LOAD4v3_PD( vu1, vu2, vu3, u );

      // ---------- advance momenta
      VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );

      // Store Results
      _MM256_STORE4v3_PD(u, vu1, vu2, vu3);

      // ---------- advance position and get velocities
      {
        __m256d vrg;

        // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        __m256d vq = _mm256_load_pd( (double const*) q );

        // Accumulate energy
        vene_group = _mm256_fmadd_pd( vq, vene, vene_group );

        // get velocities
        __m256d vv1 = _mm256_mul_pd( vu1, vrg );
        __m256d vv2 = _mm256_mul_pd( vu2, vrg ); // this is required for current deposition
        __m256d vv3 = _mm256_mul_pd( vu3, vrg ); // this is required for current deposition

        // Push positions
        __m256d vx1 = _mm256_fmadd_pd( vv1, vdt_dx1, vx0 );

        // Store virtual particles with positions still indexed to the
        // original cell
        STOREU4P1D( vpbuf, vpbuf.np, vx0, vx1, vq, vv2, vv3, vix0 );

        // Trim positions and store results
        __m256d vtr1 = vmntrim(vx1);

        vx1  = _mm256_sub_pd( vx1, vtr1 );
        _mm256_store_pd( x, vx1 );

        __m128i vitr1 = _mm256_cvttpd_epi32( vtr1 ) ;
        _mm_store_si128( (__m128i *) dix, vitr1 );

        __m128i vix1 = _mm_add_epi32( vix0, vitr1 );
        _mm_store_si128( (__m128i *) ix, vix1 );

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
    *ene += _mm256_reduce_add_pd( vene_group );

  }

}


/********************************** 2D advance deposit **********************************/


extern void ADVANCE_DEPOSIT_2D
( int *ix, double *x, double *u, double *q, int *i0, int *i1, double *rqm,
  double *efield, double *bfield, int *emf_size, int *emf_offset,
  double *j, int *j_size, int *j_offset,
  double *dx, double *dt, double *ene )
{

  /* debug and sanity checks */

/*
  printf("In ADVANCE_DEPOSIT_2D, vec_width = %d \n", VEC_WIDTH);

  CHECK_ALIGN(ix);
  CHECK_ALIGN(x);
  CHECK_ALIGN(u);
  CHECK_ALIGN(q);

  CHECK_ALIGN(efield);
  CHECK_ALIGN(bfield);
*/

  /* main code start */

  int k, i, np, np_total;

  DECLARE_ALIGNED_32( double jnorm[2] );

  DECLARE_ALIGNED_32( unsigned int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int diy[VEC_WIDTH] );

  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf2D vpbuf;

  // Simulation constants
  __m256d const vtem    = _mm256_set1_pd( 0.5f * (*dt) / (*rqm) );
  __m256d const vdt_dx1 = _mm256_set1_pd( (*dt) / dx[0] );
  __m256d const vdt_dx2 = _mm256_set1_pd( (*dt) / dx[1] );

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

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
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
       printf("(*error*) Number of particles to push is still not a multiple of %d!\n", VEC_WIDTH);
       exit(1);
    }

  }

  if ( p_cache_size_2D % VEC_WIDTH ) {
       printf("(*error*) p_cache_size_2D is not a multiple of VEC_WIDTH = %d!\n", VEC_WIDTH);
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {

    // Initialize energy
    __m256d vene_group = _mm256_setzero_pd();

    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;

  	// Initialize virtual particle buffer
  	vpbuf.np = 0;

  	// Push all particles in group
  	for( i = 0; i < np; i+=VEC_WIDTH ) {
  	  __m256d vene;
      __m256d vx0, vy0;
  	  __m128i vix0, viy0;

  	  __m256d ve1, ve2, ve3;
  	  __m256d vb1, vb2, vb3;

  	  __m256d vu1, vu2, vu3;

  	  // Load particle positions
      _MM256_LOAD4v2_PD( vx0, vy0, x );
      _MM_LOAD4v2_EPI32( vix0, viy0, ix  );

  	  // Interpolate fields
  	  {
        __m256d  hx, hy;
        __m128i vidxi, vidxj, vidxih, vidxjh;
        __m256d vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];

        // idxi = 3 * ix0
        vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
        // idxj = deltaY * iy0
        vidxj = _mm_mullo_epi32( _mm_set1_epi32( deltaY ), viy0 );

        SPLINE( vx0, vwx );
        SPLINE( vy0, vwy);

        VHP_SHIFT( vx0, vidxi, hx, vidxih, 3 );
        VHP_SHIFT( vy0, vidxj, hy, vidxjh, deltaY );

        SPLINEH( vx0, hx, vwxh );
        SPLINEH( vy0, hy, vwyh );

        // Interpolate E field
        {
          unsigned int k1, k2;

          __m128i idx1, idx2, idx3;

          idx1 = _mm_add_epi32( vidxih,  vidxj );
          idx2 = _mm_add_epi32(  vidxi, vidxjh );
          idx3 = _mm_add_epi32(  vidxi,  vidxj );

          ve1 = _mm256_setzero_pd();
          ve2 = _mm256_setzero_pd();
          ve3 = _mm256_setzero_pd();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m256d f1line, f2line, f3line;
            f1line = _mm256_setzero_pd();
            f2line = _mm256_setzero_pd();
            f3line = _mm256_setzero_pd();

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              unsigned shift = k1*3 + k2*deltaY;
              f1line = _mm256_fmadd_pd( _mm256_i32gather_pd( e1 + shift, idx1, 8 ), vwxh[k1], f1line );
              f2line = _mm256_fmadd_pd( _mm256_i32gather_pd( e2 + shift, idx2, 8 ),  vwx[k1], f2line );
              f3line = _mm256_fmadd_pd( _mm256_i32gather_pd( e3 + shift, idx3, 8 ),  vwx[k1], f3line );
            }

            ve1 = _mm256_fmadd_pd( f1line,  vwy[k2], ve1 );
            ve2 = _mm256_fmadd_pd( f2line, vwyh[k2], ve2 );
            ve3 = _mm256_fmadd_pd( f3line,  vwy[k2], ve3 );
          }

        }

        // Interpolate B field
        {
          unsigned int k1, k2;

          __m128i idx1, idx2, idx3;

          idx1 = _mm_add_epi32(  vidxi, vidxjh );
          idx2 = _mm_add_epi32( vidxih,  vidxj );
          idx3 = _mm_add_epi32( vidxih, vidxjh );

          vb1 = _mm256_setzero_pd();
          vb2 = _mm256_setzero_pd();
          vb3 = _mm256_setzero_pd();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m256d f1line, f2line, f3line;
            f1line = _mm256_setzero_pd();
            f2line = _mm256_setzero_pd();
            f3line = _mm256_setzero_pd();

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              unsigned shift = k1*3 + k2*deltaY;
              f1line = _mm256_fmadd_pd( _mm256_i32gather_pd( b1 + shift, idx1, 8 ),  vwx[k1], f1line );
              f2line = _mm256_fmadd_pd( _mm256_i32gather_pd( b2 + shift, idx2, 8 ), vwxh[k1], f2line );
              f3line = _mm256_fmadd_pd( _mm256_i32gather_pd( b3 + shift, idx3, 8 ), vwxh[k1], f3line );
            }

            vb1 = _mm256_fmadd_pd( f1line, vwyh[k2], vb1 );
            vb2 = _mm256_fmadd_pd( f2line,  vwy[k2], vb2 );
            vb3 = _mm256_fmadd_pd( f3line, vwyh[k2], vb3 );
          }
        }

      }

  	  // Load momenta
  	  _MM256_LOAD4v3_PD( vu1, vu2, vu3, u );

  	  // ---------- advance momenta
  	  VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );

      // Store Results
      _MM256_STORE4v3_PD(u, vu1, vu2, vu3);

      // ---------- advance position and get velocities
      {
        __m256d vrg;

        // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        __m256d vq = _mm256_load_pd( q );

        // Accumulate energy
        vene_group = _mm256_fmadd_pd( vq, vene, vene_group );

        // get velocities
        __m256d vv1 = _mm256_mul_pd( vu1, vrg );
        __m256d vv2 = _mm256_mul_pd( vu2, vrg );
        __m256d vv3 = _mm256_mul_pd( vu3, vrg ); // this is required for current deposition

        // Push positions
        __m256d vx1 = _mm256_fmadd_pd( vv1, vdt_dx1, vx0 );
        __m256d vy1 = _mm256_fmadd_pd( vv2, vdt_dx2, vy0 );

        // Store virtual particles with positions still indexed to the
        // original cell
        STOREU4P2D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );

        // Find position trim
        __m256d vtr1 = vmntrim(vx1);
        __m256d vtr2 = vmntrim(vy1);

        // Trim positions and store results
        vx1  = _mm256_sub_pd( vx1, vtr1 );
        vy1  = _mm256_sub_pd( vy1, vtr2 );
        _MM256_STORE4v2_PD( x, vx1, vy1 );

        // find cell crossings and store
        __m128i vitr1 = _mm256_cvtpd_epi32( vtr1 ) ;
        __m128i vitr2 = _mm256_cvtpd_epi32( vtr2 ) ;

        _mm_store_si128( (__m128i *) dix, vitr1 );
        _mm_store_si128( (__m128i *) diy, vitr2 );

        // Trim cell indexes and store
        __m128i vix1 = _mm_add_epi32( vix0, vitr1 );
        __m128i viy1 = _mm_add_epi32( viy0, vitr2 );
        _MM_STORE4v2_EPI32( ix, vix1, viy1 );

        // calculate crossings
        _mm_store_si128( (__m128i *) cross,
                    _mm_or_si128(_mm_andnot_si128(_mm_cmpeq_epi32(vix0, vix1),_mm_set1_epi32(1)),
                                 _mm_andnot_si128(_mm_cmpeq_epi32(viy0, viy1),_mm_set1_epi32(2))));

      }

      // ---------- split trajectories for current deposition
      vsplit2D( &vpbuf, cross, dix, diy );

      // ---------- advance pointers
      x  += 2 * VEC_WIDTH;
      u  += 3 * VEC_WIDTH;
      ix += 2 * VEC_WIDTH;
      q  += VEC_WIDTH;

  	}

    // Deposit current from all virtual particles
    DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    *ene += _mm256_reduce_add_pd( vene_group );

  }

}

extern void ADVANCE_DEPOSIT_2D_CYL
( int *ix, double *x, double *u, double *q, int *i0, int *i1, double *rqm,
  int *ilb2 ,
  double *efield, double *bfield, int *emf_size, int *emf_offset,
  double *j, int *j_size, int *j_offset,
  double *dx, double *dt, double *ene )
{

  /* debug and sanity checks */

/*
  printf("In ADVANCE_DEPOSIT_2D, vec_width = %d \n", VEC_WIDTH);

  CHECK_ALIGN(ix);
  CHECK_ALIGN(x);
  CHECK_ALIGN(u);
  CHECK_ALIGN(q);

  CHECK_ALIGN(efield);
  CHECK_ALIGN(bfield);
*/

  /* main code start */

  int k, i, np, np_total;

  DECLARE_ALIGNED_32( double jnorm[2] );

  DECLARE_ALIGNED_32( unsigned int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int diy[VEC_WIDTH] );

  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf2D vpbuf;

  // Simulation constants
  __m256d const vtem    = _mm256_set1_pd( 0.5 * (*dt) / (*rqm) );
  __m256d const vdt_dx1 = _mm256_set1_pd( (*dt) / dx[0] );

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
  __m256d const gshift2 = _mm256_set1_pd( (*ilb2 - 2 ) );
  #else
  __m256d const gshift2 = _mm256_set1_pd( (*ilb2 - 2 ) - 0.5 );
  #endif

  __m256d const dr      = _mm256_set1_pd( dx[1] );
  __m256d const rdr     = _mm256_set1_pd( 1.0/dx[1] );
  __m256d const vdt     = _mm256_set1_pd( *dt );

  // jump to 1st particle
  x  += 2*(*i0-1);
  ix += 2*(*i0-1);
  u  += 3*(*i0-1);
  q  += *i0-1;

  // Get total number of particles to push
  np_total = *i1 - *i0 + 1;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
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
       printf("(*error*) Number of particles to push is still not a multiple of %d!\n", VEC_WIDTH);
       exit(1);
    }

  }

  if ( p_cache_size_2D % VEC_WIDTH ) {
       printf("(*error*) p_cache_size_2D is not a multiple of VEC_WIDTH = %d!\n", VEC_WIDTH);
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {

    // Initialize energy
    __m256d vene_group = _mm256_setzero_pd();

    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;

    // Initialize virtual particle buffer
    vpbuf.np = 0;

    // Push all particles in group
    for( i = 0; i < np; i+=VEC_WIDTH ) {
      __m256d vene;
      __m256d vx0, vy0;
      __m128i vix0, viy0;

      __m256d ve1, ve2, ve3;
      __m256d vb1, vb2, vb3;

      __m256d vu1, vu2, vu3;

      // Load particle positions
      _MM256_LOAD4v2_PD( vx0, vy0, x );
      _MM_LOAD4v2_EPI32( vix0, viy0, ix  );

      // Interpolate fields
      {
        __m256d  hx, hy;
        __m128i vidxi, vidxj, vidxih, vidxjh;
        __m256d vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];

        // idxi = 3 * ix0
        vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
        // idxj = deltaY * iy0
        vidxj = _mm_mullo_epi32( _mm_set1_epi32( deltaY ), viy0 );

        SPLINE( vx0, vwx );
        SPLINE( vy0, vwy);

        VHP_SHIFT( vx0, vidxi, hx, vidxih, 3 );
        VHP_SHIFT( vy0, vidxj, hy, vidxjh, deltaY );

        SPLINEH( vx0, hx, vwxh );
        SPLINEH( vy0, hy, vwyh );

        // Interpolate E field
        {
          unsigned int k1, k2;

          __m128i idx1, idx2, idx3;

          idx1 = _mm_add_epi32( vidxih,  vidxj );
          idx2 = _mm_add_epi32(  vidxi, vidxjh );
          idx3 = _mm_add_epi32(  vidxi,  vidxj );

          ve1 = _mm256_setzero_pd();
          ve2 = _mm256_setzero_pd();
          ve3 = _mm256_setzero_pd();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m256d f1line, f2line, f3line;
            f1line = _mm256_setzero_pd();
            f2line = _mm256_setzero_pd();
            f3line = _mm256_setzero_pd();

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              unsigned shift = k1*3 + k2*deltaY;
              f1line = _mm256_fmadd_pd( _mm256_i32gather_pd( e1 + shift, idx1, 8 ), vwxh[k1], f1line );
              f2line = _mm256_fmadd_pd( _mm256_i32gather_pd( e2 + shift, idx2, 8 ),  vwx[k1], f2line );
              f3line = _mm256_fmadd_pd( _mm256_i32gather_pd( e3 + shift, idx3, 8 ),  vwx[k1], f3line );
            }

            ve1 = _mm256_fmadd_pd( f1line,  vwy[k2], ve1 );
            ve2 = _mm256_fmadd_pd( f2line, vwyh[k2], ve2 );
            ve3 = _mm256_fmadd_pd( f3line,  vwy[k2], ve3 );
          }

        }

        // Interpolate B field
        {
          unsigned int k1, k2;

          __m128i idx1, idx2, idx3;

          idx1 = _mm_add_epi32(  vidxi, vidxjh );
          idx2 = _mm_add_epi32( vidxih,  vidxj );
          idx3 = _mm_add_epi32( vidxih, vidxjh );

          vb1 = _mm256_setzero_pd();
          vb2 = _mm256_setzero_pd();
          vb3 = _mm256_setzero_pd();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m256d f1line, f2line, f3line;
            f1line = _mm256_setzero_pd();
            f2line = _mm256_setzero_pd();
            f3line = _mm256_setzero_pd();

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              unsigned shift = k1*3 + k2*deltaY;
              f1line = _mm256_fmadd_pd( _mm256_i32gather_pd( b1 + shift, idx1, 8 ),  vwx[k1], f1line );
              f2line = _mm256_fmadd_pd( _mm256_i32gather_pd( b2 + shift, idx2, 8 ), vwxh[k1], f2line );
              f3line = _mm256_fmadd_pd( _mm256_i32gather_pd( b3 + shift, idx3, 8 ), vwxh[k1], f3line );
            }

            vb1 = _mm256_fmadd_pd( f1line, vwyh[k2], vb1 );
            vb2 = _mm256_fmadd_pd( f2line,  vwy[k2], vb2 );
            vb3 = _mm256_fmadd_pd( f3line, vwyh[k2], vb3 );
          }
        }

      }

      // Load momenta
      _MM256_LOAD4v3_PD( vu1, vu2, vu3, u );

      // ---------- advance momenta
      VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );


      // ---------- advance position in cylindrical geometry
      {
        __m256d vrg;
        __m256d vx1, vy1;

        // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        __m256d vq = _mm256_load_pd( q );

        // Accumulate energy
        vene_group = _mm256_fmadd_pd( vq, vene, vene_group );

        // get velocities
        __m256d vv1 = _mm256_mul_pd( vu1, vrg );
        __m256d vv2 = _mm256_mul_pd( vu2, vrg );
        __m256d vv3 = _mm256_mul_pd( vu3, vrg ); // this is required for current deposition

        // Push positions
        vx1 = _mm256_fmadd_pd( vv1, vdt_dx1, vx0 );

        {

          __m256d gix2, tmp;
          __m256d r_old, r_new;
          __m256d x2_new, x3_new;

          // gix2   = viy0 + gshift2;          // global cell
          gix2 = _mm256_add_pd( _mm256_cvtepi32_pd( viy0 ), gshift2 );

          // r_old = ( vy0 + gix2 ) * dr;      // global radial position
          r_old = _mm256_mul_pd( _mm256_add_pd( vy0, gix2 ), dr );

          // x2_new = r_old + vv2 * dt;
          // x3_new =         vv3 * dt;
          x2_new = _mm256_fmadd_pd( vv2, vdt, r_old );
          x3_new = _mm256_mul_pd( vv3, vdt ) ;

          // r_new = sqrt( x2_new*x2_new + x3_new*x3_new );
          r_new = _mm256_sqrt_pd( _mm256_fmadd_pd( x2_new, x2_new,
                                                    _mm256_mul_pd( x3_new, x3_new ) ) );

          // This is a protection against roundoff for cold plasmas
          // vy1 = ( r_old == r_new ) ? vy0 : r_new * rdr - gix2);
          vy1 = _mm256_blendv_pd( _mm256_fmsub_pd( r_new, rdr, gix2 ), vy0,
                                   _mm256_cmp_pd( r_old, r_new, _CMP_EQ_OQ ) );

          // Correct p_r and p_\theta to conserve angular momentum
          // tmp = 1.0 / r_new;
          // vu2 = ( vu2 * x2_new + vu3 * x3_new ) * tmp;
          // vu3 = vu3 * r_old * tmp;

          tmp = _mm256_div_pd( _mm256_set1_pd(1.0), r_new );

          vu2 = _mm256_mul_pd( _mm256_fmadd_pd( vu2, x2_new,
                                                 _mm256_mul_pd( vu3, x3_new ) ),
                                tmp );

          vu3 = _mm256_mul_pd( vu3, _mm256_mul_pd( r_old, tmp ) );

        }

        // Store Momenta
        _MM256_STORE4v3_PD(u, vu1, vu2, vu3);


        // Store virtual particles with positions still indexed to the
        // original cell
        STOREU4P2D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );

        // Find position trim
        __m256d vtr1 = vmntrim(vx1);
        __m256d vtr2 = vmntrim(vy1);

        // Trim positions and store results
        vx1  = _mm256_sub_pd( vx1, vtr1 );
        vy1  = _mm256_sub_pd( vy1, vtr2 );
        _MM256_STORE4v2_PD( x, vx1, vy1 );

        // find cell crossings and store
        __m128i vitr1 = _mm256_cvtpd_epi32( vtr1 ) ;
        __m128i vitr2 = _mm256_cvtpd_epi32( vtr2 ) ;

        _mm_store_si128( (__m128i *) dix, vitr1 );
        _mm_store_si128( (__m128i *) diy, vitr2 );

        // Trim cell indexes and store
        __m128i vix1 = _mm_add_epi32( vix0, vitr1 );
        __m128i viy1 = _mm_add_epi32( viy0, vitr2 );
        _MM_STORE4v2_EPI32( ix, vix1, viy1 );

        // calculate crossings
        _mm_store_si128( (__m128i *) cross,
                    _mm_or_si128(_mm_andnot_si128(_mm_cmpeq_epi32(vix0, vix1),_mm_set1_epi32(1)),
                                 _mm_andnot_si128(_mm_cmpeq_epi32(viy0, viy1),_mm_set1_epi32(2))));

      }

      // ---------- split trajectories for current deposition
      vsplit2D( &vpbuf, cross, dix, diy );

      // ---------- advance pointers
      x  += 2 * VEC_WIDTH;
      u  += 3 * VEC_WIDTH;
      ix += 2 * VEC_WIDTH;
      q  += VEC_WIDTH;

    }

    // Deposit current from all virtual particles
    DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    *ene += _mm256_reduce_add_pd( vene_group );

  }
}

/********************************** 3D advance deposit **********************************/


extern void ADVANCE_DEPOSIT_3D
( int *ix, double *x, double *u, double *q, int *i0, int *i1, double *rqm,
  double *efield, double *bfield, int *emf_size, int *emf_offset,
  double *j, int *j_size, int *j_offset,
  double *dx, double *dt, double *ene )
{

  /* debug and sanity checks */

/*
  printf("In ADVANCE_DEPOSIT_2D, vec_width = %d \n", VEC_WIDTH);

  CHECK_ALIGN(ix);
  CHECK_ALIGN(x);
  CHECK_ALIGN(u);
  CHECK_ALIGN(q);

  CHECK_ALIGN(efield);
  CHECK_ALIGN(bfield);
*/

  /* main code start */

  int k, i, np, np_total;

  DECLARE_ALIGNED_32( double jnorm[3] );

  DECLARE_ALIGNED_32( unsigned int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int diy[VEC_WIDTH] );
  DECLARE_ALIGNED_32( int diz[VEC_WIDTH] );

  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf3D vpbuf;

  // Simulation constants
  __m256d const vtem    = _mm256_set1_pd( 0.5 * (*dt) / (*rqm) );
  __m256d const vdt_dx1 = _mm256_set1_pd( (*dt) / dx[0] );
  __m256d const vdt_dx2 = _mm256_set1_pd( (*dt) / dx[1] );
  __m256d const vdt_dx3 = _mm256_set1_pd( (*dt) / dx[2] );

  // get pointers to position 0,0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components
  unsigned int const deltaY = emf_size[0] * deltaX;
  unsigned int const deltaZ = emf_size[1] * deltaY;

  double* const e1 = efield + (emf_offset[0]-OFFSET)*deltaX +
                             (emf_offset[1]-OFFSET)*deltaY +
                             (emf_offset[2]-OFFSET)*deltaZ;
  double* const e2 = e1 + 1;
  double* const e3 = e1 + 2;

  double* const b1 = bfield + (emf_offset[0]-OFFSET)*deltaX +
                             (emf_offset[1]-OFFSET)*deltaY +
                             (emf_offset[2]-OFFSET)*deltaZ;
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

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np_total % VEC_WIDTH ) {

	for( i = 0; i < VEC_WIDTH - np_total%VEC_WIDTH; i ++ ) {
	  // printf(" idx = %d \n ", (np_total + i) );

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
       printf("(*error*) p_cache_size_3D is not a multiple of VEC_WIDTH = %d!\n", VEC_WIDTH );
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_3D ) {

    // Initialize energy
    __m256d vene_group = _mm256_setzero_pd();

    // find number of particles to process
    np = ( np_total - k < p_cache_size_3D ) ? np_total - k : p_cache_size_3D;

  	// Initialize virtual particle buffer
  	vpbuf.np = 0;

  	// Push all particles in group
  	for( i = 0; i < np; i+=VEC_WIDTH ) {
  	  __m256d vene;
      __m256d vx0, vy0, vz0;
  	  __m128i vix0, viy0, viz0;

  	  __m256d ve1, ve2, ve3;
  	  __m256d vb1, vb2, vb3;

  	  __m256d vu1, vu2, vu3;

  	  // Load particle positions
      _MM256_LOAD4v3_PD( vx0, vy0, vz0, x );
      _MM_LOAD4v3_EPI32( vix0, viy0, viz0, ix  );

      // Interpolate fields
      {
        __m256d  hx, hy, hz;
        __m128i vidxi, vidxj, vidxk, vidxih, vidxjh, vidxkh;
        __m256d vwx[NP], vwy[NP], vwz[NP], vwxh[NP], vwyh[NP], vwzh[NP];

        // idxi = 3 * ix0
        vidxi = _mm_add_epi32( _mm_add_epi32( vix0, vix0 ), vix0 );
        // idxj = deltaY * iy0
        vidxj = _mm_mullo_epi32( _mm_set1_epi32( deltaY ), viy0 );
        // idxk = deltaZ * iz0
        vidxk = _mm_mullo_epi32( _mm_set1_epi32( deltaZ ), viz0 );

        SPLINE( vx0, vwx );
        SPLINE( vy0, vwy );
        SPLINE( vz0, vwz );

        VHP_SHIFT( vx0, vidxi, hx, vidxih, 3 );
        VHP_SHIFT( vy0, vidxj, hy, vidxjh, deltaY );
        VHP_SHIFT( vz0, vidxk, hz, vidxkh, deltaZ );

        SPLINEH( vx0, hx, vwxh );
        SPLINEH( vy0, hy, vwyh );
        SPLINEH( vz0, hz, vwzh );


        // Interpolate E field
        {
          unsigned int k1, k2, k3;

          __m128i idx1, idx2, idx3;

          idx1 = _mm_add_epi32(_mm_add_epi32( vidxih, vidxj  ), vidxk  );
          idx2 = _mm_add_epi32(_mm_add_epi32( vidxi,  vidxjh ), vidxk  );
          idx3 = _mm_add_epi32(_mm_add_epi32( vidxi,  vidxj  ), vidxkh );

          ve1 = _mm256_setzero_pd();
          ve2 = _mm256_setzero_pd();
          ve3 = _mm256_setzero_pd();

          for ( k3 = 0; k3 < NP; k3++ ) {
            register __m256d f1plane, f2plane, f3plane;
            f1plane = _mm256_setzero_pd();
            f2plane = _mm256_setzero_pd();
            f3plane = _mm256_setzero_pd();

            for ( k2 = 0; k2 < NP; k2++ ) {

              register __m256d f1line, f2line, f3line;
              f1line = _mm256_setzero_pd();
              f2line = _mm256_setzero_pd();
              f3line = _mm256_setzero_pd();

              for ( k1 = 0; k1 < NP; k1++ ) {
                unsigned shift =  k1*3 + k2*deltaY + k3*deltaZ ;
                f1line = _mm256_fmadd_pd( _mm256_i32gather_pd( e1 + shift, idx1, 8 ), vwxh[k1], f1line );
                f2line = _mm256_fmadd_pd( _mm256_i32gather_pd( e2 + shift, idx2, 8 ),  vwx[k1], f2line );
                f3line = _mm256_fmadd_pd( _mm256_i32gather_pd( e3 + shift, idx3, 8 ),  vwx[k1], f3line );
              }

              f1plane = _mm256_fmadd_pd( f1line,  vwy[k2], f1plane );
              f2plane = _mm256_fmadd_pd( f2line, vwyh[k2], f2plane );
              f3plane = _mm256_fmadd_pd( f3line,  vwy[k2], f3plane );

            }

            ve1 = _mm256_fmadd_pd( f1plane,  vwz[k3], ve1 );
            ve2 = _mm256_fmadd_pd( f2plane,  vwz[k3], ve2 );
            ve3 = _mm256_fmadd_pd( f3plane, vwzh[k3], ve3 );
          }

        }

        // Interpolate B field
        {
          unsigned k1, k2, k3;

          __m128i idx1, idx2, idx3;

          idx1 = _mm_add_epi32(_mm_add_epi32( vidxi,  vidxjh ), vidxkh );
          idx2 = _mm_add_epi32(_mm_add_epi32( vidxih, vidxj  ), vidxkh );
          idx3 = _mm_add_epi32(_mm_add_epi32( vidxih, vidxjh ), vidxk  );

          vb1 = _mm256_setzero_pd();
          vb2 = _mm256_setzero_pd();
          vb3 = _mm256_setzero_pd();

          for ( k3 = 0; k3 < NP; k3++ ) {
            register __m256d f1plane, f2plane, f3plane;
            f1plane = _mm256_setzero_pd();
            f2plane = _mm256_setzero_pd();
            f3plane = _mm256_setzero_pd();

            for ( k2 = 0; k2 < NP; k2++ ) {

              register __m256d f1line, f2line, f3line;
              f1line = _mm256_setzero_pd();
              f2line = _mm256_setzero_pd();
              f3line = _mm256_setzero_pd();

              for ( k1 = 0; k1 < NP; k1++ ) {
                unsigned shift =  k1*3 + k2*deltaY + k3*deltaZ ;
                f1line = _mm256_fmadd_pd( _mm256_i32gather_pd( b1 + shift, idx1, 8 ),  vwx[k1], f1line );
                f2line = _mm256_fmadd_pd( _mm256_i32gather_pd( b2 + shift, idx2, 8 ), vwxh[k1], f2line );
                f3line = _mm256_fmadd_pd( _mm256_i32gather_pd( b3 + shift, idx3, 8 ), vwxh[k1], f3line );
              }

              f1plane = _mm256_fmadd_pd( f1line, vwyh[k2], f1plane );
              f2plane = _mm256_fmadd_pd( f2line,  vwy[k2], f2plane );
              f3plane = _mm256_fmadd_pd( f3line, vwyh[k2], f3plane );
            }

            vb1 = _mm256_fmadd_pd( f1plane, vwzh[k3], vb1 );
            vb2 = _mm256_fmadd_pd( f2plane, vwzh[k3], vb2 );
            vb3 = _mm256_fmadd_pd( f3plane,  vwz[k3], vb3 );
          }

        }

      }

      // Load momenta
      _MM256_LOAD4v3_PD( vu1, vu2, vu3, u );

      // ---------- advance momenta
      VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );

      // Store Results
      _MM256_STORE4v3_PD(u, vu1, vu2, vu3);

      // ---------- advance position and get velocities
      {
        __m256d vrg;

        // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        __m256d vq = _mm256_load_pd( q );

        // Accumulate energy
        vene_group = _mm256_fmadd_pd( vq, vene, vene_group );

        // Push positions
        __m256d vx1 = _mm256_fmadd_pd( _mm256_mul_pd( vu1, vrg ), vdt_dx1, vx0 );
        __m256d vy1 = _mm256_fmadd_pd( _mm256_mul_pd( vu2, vrg ), vdt_dx2, vy0 );
        __m256d vz1 = _mm256_fmadd_pd( _mm256_mul_pd( vu3, vrg ), vdt_dx3, vz0 );

        // Store virtual particles with positions still indexed to the
        // original cell
        STOREU4P3D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix0, viy0, viz0 );

        // Find position trim
        __m256d vtrx = vmntrim(vx1);
        __m256d vtry = vmntrim(vy1);
        __m256d vtrz = vmntrim(vz1);

        // Trim positions and store results
        vx1  = _mm256_sub_pd( vx1, vtrx );
        vy1  = _mm256_sub_pd( vy1, vtry );
        vz1  = _mm256_sub_pd( vz1, vtrz );
        _MM256_STORE4v3_PD( x, vx1, vy1, vz1 );

        // find cell crossings and store
        __m128i vitrx = _mm256_cvtpd_epi32( vtrx ) ;
        __m128i vitry = _mm256_cvtpd_epi32( vtry ) ;
        __m128i vitrz = _mm256_cvtpd_epi32( vtrz ) ;

        _mm_store_si128( (__m128i *) dix, vitrx );
        _mm_store_si128( (__m128i *) diy, vitry );
        _mm_store_si128( (__m128i *) diz, vitrz );

        // Trim cell indexes and store
        __m128i vix1 = _mm_add_epi32( vix0, vitrx );
        __m128i viy1 = _mm_add_epi32( viy0, vitry );
        __m128i viz1 = _mm_add_epi32( viz0, vitrz );
        _MM_STORE4v3_EPI32( ix, vix1, viy1, viz1 );

        __m128i vcross;
        vcross =                       _mm_andnot_si128(_mm_cmpeq_epi32(vix0, vix1),_mm_set1_epi32(1));
        vcross = _mm_or_si128( vcross, _mm_andnot_si128(_mm_cmpeq_epi32(viy0, viy1),_mm_set1_epi32(2)));
        vcross = _mm_or_si128( vcross, _mm_andnot_si128(_mm_cmpeq_epi32(viz0, viz1),_mm_set1_epi32(4)));
        _mm_store_si128( (__m128i *) cross, vcross );

      }

      // ---------- split trajectories for current deposition
      vsplit3D( &vpbuf, cross, dix, diy, diz );

      // ---------- advance pointers
      x  += 3 * VEC_WIDTH;
      u  += 3 * VEC_WIDTH;
      ix += 3 * VEC_WIDTH;
      q  += VEC_WIDTH;

    }

    // Deposit current from all virtual particles
    DEP_CURRENT_3D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    *ene += _mm256_reduce_add_pd( vene_group );
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

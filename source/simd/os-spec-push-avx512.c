/*****************************************************************************************

Relativistic particle pusher, AVX-512 optimized version (single precision)

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__


#include "os-spec-push-avx512.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "fortran.h"

#include "vector-avx512.h"
#include "splines-avx512.h"
#include "os-spec-current-avx512.h"

#include "split-vec.h"


#ifdef __AVX512_FAST_EST__
#warning Using AVX512 fast estimates for reciprocal and reciprocal square root
#endif


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

Also returns the "cross" mask
*****************************************************************************************/

// Reference implementation
#define VMNTRIM( vx, trim, cross ) { \
                                                                                         \
  trim = _mm512_setzero_ps();                                                            \
                                                                                         \
  __mmask16 mu = _mm512_cmp_ps_mask( vx, _mm512_set1_ps(  0.5f ), _MM_CMPINT_GE );       \
  __mmask16 ml = _mm512_cmp_ps_mask( vx, _mm512_set1_ps( -0.5f ), _MM_CMPINT_LT );       \
                                                                                         \
  trim = _mm512_mask_mov_ps( trim, mu, _mm512_set1_ps( +1.0f ) );                        \
  trim = _mm512_mask_mov_ps( trim, ml, _mm512_set1_ps( -1.0f ) );                        \
                                                                                         \
  cross = _mm512_kor( mu, ml );                                                          \
}


/*****************************************************************************************
vhp_shift

Gets the required h shift and cell index for interpolating quantities on a staggered grid
(quantities located at the half cell position).

  h  = ( dx < 0 ) ? 1.0 : 0.0
  ixh = ix - h * delta

*****************************************************************************************/

#define VHP_SHIFT(dx,ix, h,ixh,delta ) {                                                 \
   __mmask16 cmp = _mm512_cmp_ps_mask( dx, _mm512_setzero_ps(), _MM_CMPINT_LT );         \
                                                                                         \
   h  =  _mm512_mask_mov_ps( _mm512_setzero_ps(), cmp, _mm512_set1_ps( 1.0f ));          \
   ixh = _mm512_mask_sub_epi32( ix, cmp, ix, _mm512_set1_epi32( delta ) );               \
}

/*****************************************************************************************
vdudt_boris

Use a Boris push to advance particle velocities using interpolated
fields.
*****************************************************************************************/

#ifdef __AVX512_FAST_EST__

// This version uses approximations for 1/x and 1/sqrt(x) (correct to 2^-28) instead of the
// full precision division and square root

#define VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene ) \
{\
  __m512 const c1 = _mm512_set1_ps( 1.0f );\
  register __m512 vut1, vut2, vut3, vutsq;\
  register __m512 vrgamma, vtem_gamma, votsq;\
\
  ve1 = _mm512_mul_ps( ve1, vtem );\
  ve2 = _mm512_mul_ps( ve2, vtem );\
  ve3 = _mm512_mul_ps( ve3, vtem );\
\
  vut1 = _mm512_add_ps( vu1, ve1 );\
  vut2 = _mm512_add_ps( vu2, ve2 );\
  vut3 = _mm512_add_ps( vu3, ve3 );\
\
  vutsq = _mm512_fmadd_ps( vut3, vut3,\
         _mm512_fmadd_ps( vut2, vut2, _mm512_mul_ps( vut1, vut1 ) ) );\
\
  vrgamma = _mm512_rsqrt28_ps( _mm512_add_ps( vutsq, c1 ) );\
\
  vene = _mm512_mul_ps( vutsq, \
         _mm512_rcp28_ps( _mm512_add_ps( _mm512_rcp28_ps(vrgamma), c1) )) ;\
\
  vtem_gamma = _mm512_mul_ps( vtem, vrgamma );\
\
  vb1 = _mm512_mul_ps( vb1, vtem_gamma );\
  vb2 = _mm512_mul_ps( vb2, vtem_gamma );\
  vb3 = _mm512_mul_ps( vb3, vtem_gamma );\
\
  vu1 = _mm512_fmadd_ps( vb3, vut2, vut1 );\
  vu2 = _mm512_fmadd_ps( vb1, vut3, vut2 );\
  vu3 = _mm512_fmadd_ps( vb2, vut1, vut3 );\
\
  vu1 = _mm512_fnmadd_ps( vb2, vut3, vu1 );\
  vu2 = _mm512_fnmadd_ps( vb3, vut1, vu2 );\
  vu3 = _mm512_fnmadd_ps( vb1, vut2, vu3 );\
\
  votsq = _mm512_fmadd_ps( vb1, vb1,    c1 );\
  votsq = _mm512_fmadd_ps( vb2, vb2, votsq );\
  votsq = _mm512_fmadd_ps( vb3, vb3, votsq );\
\
  votsq = _mm512_rcp28_ps(votsq) ;\
  votsq = _mm512_add_ps( votsq, votsq );\
\
  vb1 = _mm512_mul_ps( vb1, votsq );\
  vb2 = _mm512_mul_ps( vb2, votsq );\
  vb3 = _mm512_mul_ps( vb3, votsq );\
\
  vut1 = _mm512_fmadd_ps( vb3, vu2, vut1 );\
  vut2 = _mm512_fmadd_ps( vb1, vu3, vut2 );\
  vut3 = _mm512_fmadd_ps( vb2, vu1, vut3 );\
\
  vut1 = _mm512_fnmadd_ps( vb2, vu3, vut1 );\
  vut2 = _mm512_fnmadd_ps( vb3, vu1, vut2 );\
  vut3 = _mm512_fnmadd_ps( vb1, vu2, vut3 );\
\
  vu1 = _mm512_add_ps( vut1, ve1 );\
  vu2 = _mm512_add_ps( vut2, ve2 );\
  vu3 = _mm512_add_ps( vut3, ve3 );\
}

#else

// This version uses full precision division and square root

#define VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene ) \
{\
  __m512 const c1 = _mm512_set1_ps( 1.0f );\
  __m512 vut1, vut2, vut3, vutsq;\
  __m512 vgamma, vtem_gamma, votsq;\
\
  ve1 = _mm512_mul_ps( ve1, vtem );\
  ve2 = _mm512_mul_ps( ve2, vtem );\
  ve3 = _mm512_mul_ps( ve3, vtem );\
\
  vut1 = _mm512_add_ps( vu1, ve1 );\
  vut2 = _mm512_add_ps( vu2, ve2 );\
  vut3 = _mm512_add_ps( vu3, ve3 );\
\
  vutsq = _mm512_fmadd_ps( vut3, vut3,\
         _mm512_fmadd_ps( vut2, vut2, _mm512_mul_ps( vut1, vut1 ) ) );\
\
  vgamma = _mm512_sqrt_ps( _mm512_add_ps( vutsq, c1 ) );\
\
  vene = _mm512_div_ps( vutsq, _mm512_add_ps( vgamma, c1) );\
\
  vtem_gamma = _mm512_div_ps( vtem, vgamma );\
\
  vb1 = _mm512_mul_ps( vb1, vtem_gamma );\
  vb2 = _mm512_mul_ps( vb2, vtem_gamma );\
  vb3 = _mm512_mul_ps( vb3, vtem_gamma );\
\
  vu1 = _mm512_fmadd_ps( vb3, vut2, vut1 );\
  vu2 = _mm512_fmadd_ps( vb1, vut3, vut2 );\
  vu3 = _mm512_fmadd_ps( vb2, vut1, vut3 );\
\
  vu1 = _mm512_fnmadd_ps( vb2, vut3, vu1 );\
  vu2 = _mm512_fnmadd_ps( vb3, vut1, vu2 );\
  vu3 = _mm512_fnmadd_ps( vb1, vut2, vu3 );\
\
  votsq = _mm512_fmadd_ps( vb1, vb1,    c1 );\
  votsq = _mm512_fmadd_ps( vb2, vb2, votsq );\
  votsq = _mm512_fmadd_ps( vb3, vb3, votsq );\
\
  votsq = _mm512_div_ps( _mm512_set1_ps( 2.0f ), votsq ) ;\
\
  vb1 = _mm512_mul_ps( vb1, votsq );\
  vb2 = _mm512_mul_ps( vb2, votsq );\
  vb3 = _mm512_mul_ps( vb3, votsq );\
\
  vut1 = _mm512_fmadd_ps( vb3, vu2, vut1 );\
  vut2 = _mm512_fmadd_ps( vb1, vu3, vut2 );\
  vut3 = _mm512_fmadd_ps( vb2, vu1, vut3 );\
\
  vut1 = _mm512_fnmadd_ps( vb2, vu3, vut1 );\
  vut2 = _mm512_fnmadd_ps( vb3, vu1, vut2 );\
  vut3 = _mm512_fnmadd_ps( vb1, vu2, vut3 );\
\
  vu1 = _mm512_add_ps( vut1, ve1 );\
  vu2 = _mm512_add_ps( vut2, ve2 );\
  vu3 = _mm512_add_ps( vut3, ve3 );\
}

#endif

/*****************************************************************************************
VRGAMMA

Get 1 / Lorentz gamma.
*****************************************************************************************/

#ifdef __AVX512_FAST_EST__

// This version uses an approximations for 1/sqrt(x) (correct to 2^-28) instead of the
// full precision division and square root

#define VRGAMMA( vrg, vu1, vu2, vu3 ) {                                 \
   const __m512 c1 = _mm512_set1_ps( 1.0f );                            \
   (vrg) = _mm512_rsqrt28_ps( _mm512_fmadd_ps( (vu3), (vu3),              \
                              _mm512_fmadd_ps( (vu2), (vu2),              \
                              _mm512_fmadd_ps( (vu1), (vu1), c1 ) ) ) );  \
}

#else

#define VRGAMMA( vrg, vu1, vu2, vu3 ) {                                                 \
   const __m512 c1 = _mm512_set1_ps( 1.0f );                                            \
   (vrg) = _mm512_div_ps( c1, _mm512_sqrt_ps( _mm512_fmadd_ps( (vu3), (vu3),              \
                                              _mm512_fmadd_ps( (vu2), (vu2),              \
                                              _mm512_fmadd_ps( (vu1), (vu1), c1 ) ) ) ) );\
}

#endif

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

#define SPLINE  ONAME( vspline, ORDER )
#define SPLINEH  ONAME( vsplineh, ORDER )

#define DEP_CURRENT_1D ONAME( vdepcurrent_1d, ORDER )
#define DEP_CURRENT_2D ONAME( vdepcurrent_2d, ORDER )
#define DEP_CURRENT_3D ONAME( vdepcurrent_3d, ORDER )

/********************************** 1D advance deposit **********************************/


extern void ADVANCE_DEPOSIT_1D
( int *ix, float *x, float *u, float *q, int *i0, int *i1, float *rqm,
  float *efield, float *bfield, int *emf_size, int *emf_offset,
  float *j, int *j_size, int *j_offset,
  double *dx, double *dt, double *ene )
{

  /* debug and sanity checks */

/*
  printf("In ADVANCE_DEPOSIT_1D, vec_width = %d \n", VEC_WIDTH);

  CHECK_ALIGN(ix);
  CHECK_ALIGN(x);
  CHECK_ALIGN(u);
  CHECK_ALIGN(q);

  CHECK_ALIGN(efield);
  CHECK_ALIGN(bfield);
*/

  /* main code start */

  int k, i, np, np_total;

  DECLARE_ALIGNED_64( float jnorm[1] );
  DECLARE_ALIGNED_64( int dix[VEC_WIDTH] );

  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf1D vpbuf;

  // Simulation constants
  __m512 const vtem    = _mm512_set1_ps( (float) 0.5 * (*dt) / (*rqm) );
  __m512 const vdt_dx1 = _mm512_set1_ps( (float) (*dt) / dx[0] );

  // get pointers to position 0 of each field component
  unsigned int const deltaX = 3; // 3 field components

  float* const e1 = efield + (emf_offset[0] - OFFSET ) * deltaX;
  float* const e2 = e1 + 1;
  float* const e3 = e1 + 2;

  float* const b1 = bfield + (emf_offset[0] - OFFSET ) * deltaX;
  float* const b2 = b1  + 1;
  float* const b3 = b1  + 2;


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
    __m512 vene_group = _mm512_setzero_ps();

    float* bx   =  &x[ k ];
    int*   bix  = &ix[ k ];
    float* bu   =  &u[ 3*k ];
    float* bq   =  &q[ k ];

    // find number of particles to process
    np = ( np_total - k < p_cache_size_1D ) ? np_total - k : p_cache_size_1D;

    // Initialize virtual particle buffer
    vpbuf.np = 0;

    // Push all particles in group
    for( i = 0; i < np; i+=VEC_WIDTH ) {
      __m512 vene;
      __m512 vx0;
      __m512i vix0;

      __m512 ve1, ve2, ve3;
      __m512 vb1, vb2, vb3;

      // Load particle positions
      vx0  = _mm512_load_ps( &bx[ i ] );
      vix0 = _mm512_load_epi32( &bix[ i ] );

      // Interpolate fields
      {
        __m512  hx;
        __m512i vidxi, vidxih;
        __m512 vwx[NP], vwxh[NP];

        // idxi = 3 * ix0
        vidxi = _mm512_add_epi32( _mm512_add_epi32( vix0, vix0 ), vix0 );

        SPLINE( vx0, vwx );

        VHP_SHIFT( vx0, vidxi, hx, vidxih, 3 );

        SPLINEH( vx0, hx, vwxh );

        // Interpolate E field
        {
          unsigned int k1;

          ve1 = _mm512_setzero_ps();
          ve2 = _mm512_setzero_ps();
          ve3 = _mm512_setzero_ps();

          for ( k1 = 0; k1 < NP; k1 ++ ) {
            unsigned shift = k1*3;
            ve1 = _mm512_fmadd_ps( _mm512_i32gather_ps( vidxih, e1 + shift, _MM_SCALE_4 ), vwxh[k1], ve1 );
            ve2 = _mm512_fmadd_ps( _mm512_i32gather_ps( vidxi , e2 + shift, _MM_SCALE_4 ),  vwx[k1], ve2 );
            ve3 = _mm512_fmadd_ps( _mm512_i32gather_ps( vidxi , e3 + shift, _MM_SCALE_4 ),  vwx[k1], ve3 );
          }

        }

        // Interpolate B field
        {
          unsigned int k1;

          vb1 = _mm512_setzero_ps();
          vb2 = _mm512_setzero_ps();
          vb3 = _mm512_setzero_ps();

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              unsigned shift = k1*3;
              vb1 = _mm512_fmadd_ps( _mm512_i32gather_ps(  vidxi, b1 + shift, _MM_SCALE_4 ),  vwx[k1], vb1 );
              vb2 = _mm512_fmadd_ps( _mm512_i32gather_ps( vidxih, b2 + shift, _MM_SCALE_4 ), vwxh[k1], vb2 );
              vb3 = _mm512_fmadd_ps( _mm512_i32gather_ps( vidxih, b3 + shift, _MM_SCALE_4 ), vwxh[k1], vb3 );
            }

        }

      }

      // Load momenta
      __m512 vu1 = _mm512_load_ps( &bu[ 3*i      ] );
      __m512 vu2 = _mm512_load_ps( &bu[ 3*i + 16 ] );
      __m512 vu3 = _mm512_load_ps( &bu[ 3*i + 32 ] );
      _MM512_TRANSPOSE16v3_PS( vu1, vu2, vu3 );

      // ---------- advance momenta
      VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );

      // ---------- advance position and get velocities
      {
        __m512 vrg;

        // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        __m512 vq = _mm512_load_ps( &bq[ i ] );

        // Accumulate energy
        vene_group = _mm512_fmadd_ps( vq, vene, vene_group );

        // get velocities
        __m512 vv1 = _mm512_mul_ps( vu1, vrg );
        __m512 vv2 = _mm512_mul_ps( vu2, vrg ); // this is required for current deposition
        __m512 vv3 = _mm512_mul_ps( vu3, vrg ); // this is required for current deposition

        // Store momenta
        _MM512_TRANSPOSE3v16_PS( vu1, vu2, vu3 );
        _mm512_store_ps( &bu[ 3*i      ], vu1 );
        _mm512_store_ps( &bu[ 3*i + 16 ], vu2 );
        _mm512_store_ps( &bu[ 3*i + 32 ], vu3 );

        // Push positions
        __m512 vx1 = _mm512_fmadd_ps( vv1, vdt_dx1, vx0 );

        // Store virtual particles with positions still indexed to the
        // original cell
        STOREU16P1D( vpbuf, vpbuf.np, vx0, vx1, vq, vv2, vv3, vix0 );

        // Find position trim, calculate crossings and store result
        __m512 vtr1;
        __mmask16 xcross;
        VMNTRIM( vx1, vtr1, xcross );

        // Trim positions and store results
        vx1 = _mm512_sub_ps( vx1, vtr1 );
        _mm512_store_ps( &bx[ i ], vx1 );

        // find cell crossings and store
        __m512i vitr1 = _mm512_cvtps_epi32( vtr1 ) ;
        _mm512_store_epi32( (__m512i *) dix, vitr1 );

        // Trim cell indexes and store
        vix0 = _mm512_add_epi32( vix0, vitr1 );
        _mm512_store_epi32( &bix[ i ], vix0 );

      }

      // ---------- split trajectories for current deposition
      vsplit1D( &vpbuf, dix );

    }

    // Deposit current from all virtual particles
    DEP_CURRENT_1D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    *ene += _mm512_reduce_add_ps( vene_group );

  }

}


/********************************** 2D advance deposit **********************************/


extern void ADVANCE_DEPOSIT_2D
( int *ix, float *x, float *u, float *q, int *i0, int *i1, float *rqm,
  float *efield, float *bfield, int *emf_size, int *emf_offset,
  float *j, int *j_size, int *j_offset,
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

  DECLARE_ALIGNED_64( float jnorm[2] );

  DECLARE_ALIGNED_64( unsigned int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_64( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_64( int diy[VEC_WIDTH] );

  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf2D vpbuf;

  // Simulation constants
  __m512 const vtem    = _mm512_set1_ps( 0.5f * (*dt) / (*rqm) );
  __m512 const vdt_dx1 = _mm512_set1_ps( (*dt) / dx[0] );
  __m512 const vdt_dx2 = _mm512_set1_ps( (*dt) / dx[1] );

  // get pointers to position 0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components
  unsigned int const deltaY = emf_size[0] * deltaX;

  float* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX +
                             (emf_offset[1] - OFFSET) * deltaY;
  float* const e2 = e1 + 1;
  float* const e3 = e1 + 2;

  float* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX +
                             (emf_offset[1] - OFFSET) * deltaY;
  float* const b2 = b1  + 1;
  float* const b3 = b1  + 2;


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

  }

  if ( p_cache_size_2D % VEC_WIDTH ) {
       printf("(*error*) p_cache_size_2D is not a multiple of VEC_WIDTH = %d!\n", VEC_WIDTH);
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {

    // Initialize energy
    __m512 vene_group = _mm512_setzero_ps();

    float* bx   =  &x[ 2*k ];
    int*   bix  = &ix[ 2*k ];
    float* bu   =  &u[ 3*k ];
    float* bq   =  &q[ k ];

    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;

    // Initialize virtual particle buffer
    vpbuf.np = 0;

    // Push all particles in group
    for( i = 0; i < np; i+=VEC_WIDTH ) {
      __m512 vene;
      __m512 vx0, vx1, vy0, vy1;
      __m512i vix0, viy0;

      register __m512 ve1, ve2, ve3;
      register __m512 vb1, vb2, vb3;

      __m512 vu1, vu2, vu3;

      // Load particle positions
      // _MM512_LOAD16v2_PS( vx0, vy0, x );
      vx0 = _mm512_load_ps( &bx[ 2*i      ] );
      vy0 = _mm512_load_ps( &bx[ 2*i + 16 ] );
      _MM512_TRANSPOSE16v2_PS( vx0, vy0 );

      // _MM512_LOAD16v2_EPI32( vix0, viy0, ix  );
      vix0 = _mm512_load_epi32( &bix[ 2*i      ] );
      viy0 = _mm512_load_epi32( &bix[ 2*i + 16 ] );
      _MM512_TRANSPOSE16v2_EPI32( vix0, viy0 );

      // Interpolate fields
      {
        __m512  hx, hy;
        __m512i vidxi, vidxj, vidxih, vidxjh;
        __m512 vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];

        // idxi = 3 * ix0
        vidxi = _mm512_add_epi32( _mm512_add_epi32( vix0, vix0 ), vix0 );
        // idxj = deltaY * iy0
        vidxj = _mm512_mullo_epi32( _mm512_set1_epi32( deltaY ), viy0 );

        SPLINE( vx0, vwx );
        SPLINE( vy0, vwy);

        VHP_SHIFT( vx0, vidxi, hx, vidxih, 3 );
        VHP_SHIFT( vy0, vidxj, hy, vidxjh, deltaY );

        SPLINEH( vx0, hx, vwxh );
        SPLINEH( vy0, hy, vwyh );

        // Interpolate E field
        {
          unsigned int k1, k2;

          __m512i idx1, idx2, idx3;

          idx1 = _mm512_add_epi32( vidxih,  vidxj );
          idx2 = _mm512_add_epi32(  vidxi, vidxjh );
          idx3 = _mm512_add_epi32(  vidxi,  vidxj );

          ve1 = _mm512_setzero_ps();
          ve2 = _mm512_setzero_ps();
          ve3 = _mm512_setzero_ps();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m512 f1line, f2line, f3line;
            f1line = _mm512_setzero_ps();
            f2line = _mm512_setzero_ps();
            f3line = _mm512_setzero_ps();

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              unsigned shift = k1*3 + k2*deltaY;
              f1line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx1 , e1 + shift, _MM_SCALE_4 ), vwxh[k1], f1line );
              f2line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx2 , e2 + shift, _MM_SCALE_4 ),  vwx[k1], f2line );
              f3line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx3 , e3 + shift, _MM_SCALE_4 ),  vwx[k1], f3line );
            }

            ve1 = _mm512_fmadd_ps( f1line,  vwy[k2], ve1 );
            ve2 = _mm512_fmadd_ps( f2line, vwyh[k2], ve2 );
            ve3 = _mm512_fmadd_ps( f3line,  vwy[k2], ve3 );
          }

        }

        // Interpolate B field
        {
          unsigned int k1, k2;

          __m512i idx1, idx2, idx3;

          idx1 = _mm512_add_epi32(  vidxi, vidxjh );
          idx2 = _mm512_add_epi32( vidxih,  vidxj );
          idx3 = _mm512_add_epi32( vidxih, vidxjh );

          vb1 = _mm512_setzero_ps();
          vb2 = _mm512_setzero_ps();
          vb3 = _mm512_setzero_ps();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m512 f1line, f2line, f3line;
            f1line = _mm512_setzero_ps();
            f2line = _mm512_setzero_ps();
            f3line = _mm512_setzero_ps();

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              unsigned shift = k1*3 + k2*deltaY;
              f1line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx1, b1 + shift, _MM_SCALE_4 ),  vwx[k1], f1line );
              f2line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx2, b2 + shift, _MM_SCALE_4 ), vwxh[k1], f2line );
              f3line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx3, b3 + shift, _MM_SCALE_4 ), vwxh[k1], f3line );
            }

            vb1 = _mm512_fmadd_ps( f1line, vwyh[k2], vb1 );
            vb2 = _mm512_fmadd_ps( f2line,  vwy[k2], vb2 );
            vb3 = _mm512_fmadd_ps( f3line, vwyh[k2], vb3 );
          }
        }

      }

      // Load momenta
      vu1 = _mm512_load_ps( &bu[ 3*i      ] );
      vu2 = _mm512_load_ps( &bu[ 3*i + 16 ] );
      vu3 = _mm512_load_ps( &bu[ 3*i + 32 ] );
      _MM512_TRANSPOSE16v3_PS( vu1, vu2, vu3 );

      // ---------- advance momenta
      VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );

      // ---------- advance position and get velocities
      {
        register __m512 vrg;
        register __m512 vtr1, vtr2;
        register __m512i vitr1, vitr2;
        register __m512 vv1, vv2, vv3;
        register __m512 vq;

        __mmask16 xcross, ycross;

        // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        vq = _mm512_load_ps( &bq[ i ] );

        // Accumulate energy
        vene_group = _mm512_fmadd_ps( vq, vene, vene_group );

        // get velocities
        vv1 = _mm512_mul_ps( vu1, vrg );
        vv2 = _mm512_mul_ps( vu2, vrg );
        vv3 = _mm512_mul_ps( vu3, vrg ); // this is required for current deposition

        // Store momenta
        _MM512_TRANSPOSE3v16_PS( vu1, vu2, vu3 );
        _mm512_store_ps( &bu[ 3*i      ], vu1 );
        _mm512_store_ps( &bu[ 3*i + 16 ], vu2 );
        _mm512_store_ps( &bu[ 3*i + 32 ], vu3 );

        // Push positions
        vx1 = _mm512_fmadd_ps( vv1, vdt_dx1, vx0 );
        vy1 = _mm512_fmadd_ps( vv2, vdt_dx2, vy0 );


        // Store virtual particles with positions still indexed to the
        // original cell
        STOREU16P2D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );

        // Find position trim, calculate crossings and store result
        VMNTRIM( vx1, vtr1, xcross );
        VMNTRIM( vy1, vtr2, ycross );

        __m512i vcross = _mm512_setzero_epi32();

        // x crossings
        vcross = _mm512_mask_or_epi32( vcross, xcross, vcross, _mm512_set1_epi32(1) );
        // y crossings
        vcross = _mm512_mask_or_epi32( vcross, ycross, vcross, _mm512_set1_epi32(2) );
        // Store vcross
        _mm512_store_epi32( (__m512i *) cross, vcross );

        // Trim positions and store results
        vx1  = _mm512_sub_ps( vx1, vtr1 );
        vy1  = _mm512_sub_ps( vy1, vtr2 );

        //    _MM512_STORE16v2_PS( x, vx1, vy1 );
        _MM512_TRANSPOSE2v16_PS( vx1, vy1 );
        _mm512_store_ps( &bx[ 2*i      ], vx1 );
        _mm512_store_ps( &bx[ 2*i + 16 ], vy1 );

        // find cell crossings and store
        vitr1 = _mm512_cvtps_epi32( vtr1 ) ;
        vitr2 = _mm512_cvtps_epi32( vtr2 ) ;

        _mm512_store_epi32( (__m512i *) dix, vitr1 );
        _mm512_store_epi32( (__m512i *) diy, vitr2 );

        // Trim cell indexes and store
        vix0 = _mm512_add_epi32( vix0, vitr1 );
        viy0 = _mm512_add_epi32( viy0, vitr2 );

        //    _MM512_STORE16v2_EPI32( ix, vix1, viy1 );
        _MM512_TRANSPOSE2v16_EPI32( vix0, viy0 );
        _mm512_store_epi32( &bix[ 2*i      ], vix0 );
        _mm512_store_epi32( &bix[ 2*i + 16 ], viy0 );

      }

      // ---------- split trajectories for current deposition
      vsplit2D( &vpbuf, cross, dix, diy );

    }

    // Deposit current from all virtual particles
    DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    *ene += _mm512_reduce_add_ps( vene_group );

  }

}

extern void ADVANCE_DEPOSIT_2D_CYL
( int *ix, float *x, float *u, float *q, int *i0, int *i1, float *rqm,
  int *ilb2 ,
  float *efield, float *bfield, int *emf_size, int *emf_offset,
  float *j, int *j_size, int *j_offset,
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

  DECLARE_ALIGNED_64( float jnorm[2] );

  DECLARE_ALIGNED_64( unsigned int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_64( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_64( int diy[VEC_WIDTH] );

  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf2D vpbuf;

  // Simulation constants
  __m512 const vtem    = _mm512_set1_ps( (float) 0.5f * (*dt) / (*rqm) );
  __m512 const vdt_dx1 = _mm512_set1_ps( (float) (*dt) / dx[0] );

  // get pointers to position 0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components
  unsigned int const deltaY = emf_size[0] * deltaX;

  float* const e1 = efield + (emf_offset[0] - OFFSET) * deltaX +
                             (emf_offset[1] - OFFSET) * deltaY;
  float* const e2 = e1 + 1;
  float* const e3 = e1 + 2;

  float* const b1 = bfield + (emf_offset[0] - OFFSET) * deltaX +
                             (emf_offset[1] - OFFSET) * deltaY;
  float* const b2 = b1  + 1;
  float* const b3 = b1  + 2;


  // Normalization for currents
  jnorm[0] = dx[0] / (*dt) / 2;
  jnorm[1] = dx[1] / (*dt) / 2;

  #if ( ORDER == 1 ) || ( ORDER == 3 )
  __m512d const gshift2 = _mm512_set1_pd( (*ilb2 - 2 ) );
  #else
  __m512d const gshift2 = _mm512_set1_pd( (*ilb2 - 2 ) - 0.5 );
  #endif

  __m512d const dr      = _mm512_set1_pd( dx[1] );
  __m512d const rdr     = _mm512_set1_pd( 1.0/dx[1] );
  __m512d const vdt     = _mm512_set1_pd( *dt );

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
    x[ 2 * (np_total + i) + 1 ] = 0.25f; // Avoid the axis
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
       printf("(*error*) p_cache_size_2D is not a multiple of VEC_WIDTH = %d!\n", VEC_WIDTH);
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_2D ) {

    // Initialize energy
    __m512 vene_group = _mm512_setzero_ps();

    float* bx   =  &x[ 2*k ];
    int*   bix  = &ix[ 2*k ];
    float* bu   =  &u[ 3*k ];
    float* bq   =  &q[ k ];

    // find number of particles to process
    np = ( np_total - k < p_cache_size_2D ) ? np_total - k : p_cache_size_2D;

    // Initialize virtual particle buffer
    vpbuf.np = 0;

    // Push all particles in group
    for( i = 0; i < np; i+=VEC_WIDTH ) {
      __m512 vene;
      __m512 vx0, vy0;
      __m512i vix0, viy0;

      __m512 ve1, ve2, ve3;
      __m512 vb1, vb2, vb3;

      __m512 vu1, vu2, vu3;

      // Load particle positions
      // _MM512_LOAD16v2_PS( vx0, vy0, x );
      vx0 = _mm512_load_ps( &bx[ 2*i      ] );
      vy0 = _mm512_load_ps( &bx[ 2*i + 16 ] );
      _MM512_TRANSPOSE16v2_PS( vx0, vy0 );

      // _MM512_LOAD16v2_EPI32( vix0, viy0, ix  );
      vix0 = _mm512_load_epi32( &bix[ 2*i      ] );
      viy0 = _mm512_load_epi32( &bix[ 2*i + 16 ] );
      _MM512_TRANSPOSE16v2_EPI32( vix0, viy0 );

      // Interpolate fields
      {
        __m512  hx, hy;
        __m512i vidxi, vidxj, vidxih, vidxjh;
        __m512 vwx[NP], vwy[NP], vwxh[NP], vwyh[NP];

        // idxi = 3 * ix0
        vidxi = _mm512_add_epi32( _mm512_add_epi32( vix0, vix0 ), vix0 );
        // idxj = deltaY * iy0
        vidxj = _mm512_mullo_epi32( _mm512_set1_epi32( deltaY ), viy0 );

        SPLINE( vx0, vwx );
        SPLINE( vy0, vwy);

        VHP_SHIFT( vx0, vidxi, hx, vidxih, 3 );
        VHP_SHIFT( vy0, vidxj, hy, vidxjh, deltaY );

        SPLINEH( vx0, hx, vwxh );
        SPLINEH( vy0, hy, vwyh );

        // Interpolate E field
        {
          unsigned int k1, k2;

          __m512i idx1, idx2, idx3;

          idx1 = _mm512_add_epi32( vidxih,  vidxj );
          idx2 = _mm512_add_epi32(  vidxi, vidxjh );
          idx3 = _mm512_add_epi32(  vidxi,  vidxj );

          ve1 = _mm512_setzero_ps();
          ve2 = _mm512_setzero_ps();
          ve3 = _mm512_setzero_ps();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m512 f1line, f2line, f3line;
            f1line = _mm512_setzero_ps();
            f2line = _mm512_setzero_ps();
            f3line = _mm512_setzero_ps();

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              unsigned shift = k1*3 + k2*deltaY;
              f1line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx1 , e1 + shift, _MM_SCALE_4 ), vwxh[k1], f1line );
              f2line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx2 , e2 + shift, _MM_SCALE_4 ),  vwx[k1], f2line );
              f3line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx3 , e3 + shift, _MM_SCALE_4 ),  vwx[k1], f3line );
            }

            ve1 = _mm512_fmadd_ps( f1line,  vwy[k2], ve1 );
            ve2 = _mm512_fmadd_ps( f2line, vwyh[k2], ve2 );
            ve3 = _mm512_fmadd_ps( f3line,  vwy[k2], ve3 );
          }

        }

        // Interpolate B field
        {
          unsigned int k1, k2;

          __m512i idx1, idx2, idx3;

          idx1 = _mm512_add_epi32(  vidxi, vidxjh );
          idx2 = _mm512_add_epi32( vidxih,  vidxj );
          idx3 = _mm512_add_epi32( vidxih, vidxjh );

          vb1 = _mm512_setzero_ps();
          vb2 = _mm512_setzero_ps();
          vb3 = _mm512_setzero_ps();

          for ( k2 = 0; k2 < NP; k2++ ) {

            register __m512 f1line, f2line, f3line;
            f1line = _mm512_setzero_ps();
            f2line = _mm512_setzero_ps();
            f3line = _mm512_setzero_ps();

            for ( k1 = 0; k1 < NP; k1 ++ ) {
              unsigned shift = k1*3 + k2*deltaY;
              f1line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx1, b1 + shift, _MM_SCALE_4 ),  vwx[k1], f1line );
              f2line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx2, b2 + shift, _MM_SCALE_4 ), vwxh[k1], f2line );
              f3line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx3, b3 + shift, _MM_SCALE_4 ), vwxh[k1], f3line );
            }

            vb1 = _mm512_fmadd_ps( f1line, vwyh[k2], vb1 );
            vb2 = _mm512_fmadd_ps( f2line,  vwy[k2], vb2 );
            vb3 = _mm512_fmadd_ps( f3line, vwyh[k2], vb3 );
          }
        }

      }

      // Load momenta
      vu1 = _mm512_load_ps( &bu[ 3*i      ] );
      vu2 = _mm512_load_ps( &bu[ 3*i + 16 ] );
      vu3 = _mm512_load_ps( &bu[ 3*i + 32 ] );
      _MM512_TRANSPOSE16v3_PS( vu1, vu2, vu3 );

      // ---------- advance momenta
      VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );

      // ---------- advance position in cylindrical geometry
      {
        __m512 vrg;
        __m512 vx1, vy1;

        // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        __m512 vq = _mm512_load_ps( &bq[ i ] );

        // Accumulate energy
        vene_group = _mm512_fmadd_ps( vq, vene, vene_group );

        // get velocities
        __m512 vv1 = _mm512_mul_ps( vu1, vrg );
        __m512 vv2 = _mm512_mul_ps( vu2, vrg );
        __m512 vv3 = _mm512_mul_ps( vu3, vrg ); // this is required for current deposition

        // Push positions
        vx1 = _mm512_fmadd_ps( vv1, vdt_dx1, vx0 );

        {  // This section needs to be performed in double precision

          __m512d gix2a, gix2b;
          __m512d vy0a, vy0b, vy1a, vy1b;
          __m512d r_olda, r_oldb, r_newa, r_newb;
          __m512d x2_newa, x2_newb, x3_newa, x3_newb;

          __m512d vv2a, vv2b, vv3a, vv3b;
          __m512d vu2a, vu2b, vu3a, vu3b;

          __m512d tmpa, tmpb;

          // gix2   = viy0 + gshift2;          // global cell
          gix2a = _mm512_add_pd( _mm512_cvtepi32_pd( _mm512_castsi512_si256(viy0) ), gshift2 );
          gix2b = _mm512_add_pd( _mm512_cvtepi32_pd( _mm512_extracti64x4_epi64(viy0, 1) ), gshift2 );

          // r_old = ( vy0 + gix2 ) * dr;      // global radial position
          _MM512_CVTPS_PD( vy0a, vy0b, vy0 );
          r_olda = _mm512_mul_pd( _mm512_add_pd( vy0a, gix2a ), dr );
          r_oldb = _mm512_mul_pd( _mm512_add_pd( vy0b, gix2b ), dr );

          // x2_new = r_old + vv2 * dt;
          // x3_new =         vv3 * dt;

          _MM512_CVTPS_PD( vv2a, vv2b, vv2 );
          _MM512_CVTPS_PD( vv3a, vv3b, vv3 );

          x2_newa = _mm512_fmadd_pd( vv2a, vdt, r_olda );
          x3_newa = _mm512_mul_pd( vv3a, vdt ) ;

          x2_newb = _mm512_fmadd_pd( vv2b, vdt, r_oldb );
          x3_newb = _mm512_mul_pd( vv3b, vdt );

          // r_new = sqrt( x2_new*x2_new + x3_new*x3_new );
          r_newa = _mm512_sqrt_pd( _mm512_fmadd_pd( x2_newa, x2_newa,
                                                    _mm512_mul_pd( x3_newa, x3_newa ) ) );

          r_newb = _mm512_sqrt_pd( _mm512_fmadd_pd( x2_newb, x2_newb,
                                                    _mm512_mul_pd( x3_newb, x3_newb ) ) );

          // This is a protection against roundoff for cold plasmas
          // vy1 = ( r_old == r_new ) ? vy0 : r_new * rdr - gix2);
          vy1a = _mm512_mask_mov_pd( _mm512_fmsub_pd( r_newa, rdr, gix2a ),
                                   _mm512_cmp_pd_mask( r_olda, r_newa, _CMP_EQ_OQ ), vy0a );
          vy1b = _mm512_mask_mov_pd( _mm512_fmsub_pd( r_newb, rdr, gix2b ),
                                   _mm512_cmp_pd_mask( r_oldb, r_newb, _CMP_EQ_OQ ), vy0b );

          // Convert vy1 to single precision
          _MM512_CVTPD_PS( vy1, vy1a, vy1b );

          // Correct p_r and p_\theta to conserve angular momentum
          // tmp = 1.0 / r_new;
          // vu2 = ( vu2 * x2_new + vu3 * x3_new ) * tmp;
          // vu3 = vu3 * r_old * tmp;

          // Convert vu2, vu3 to double precision
          _MM512_CVTPS_PD( vu2a, vu2b, vu2 );
          _MM512_CVTPS_PD( vu3a, vu3b, vu3 );

          tmpa = _mm512_div_pd( _mm512_set1_pd(1.0), r_newa );
          tmpb = _mm512_div_pd( _mm512_set1_pd(1.0), r_newb );

          vu2a = _mm512_mul_pd( _mm512_fmadd_pd( vu2a, x2_newa,
                                                 _mm512_mul_pd( vu3a, x3_newa ) ),
                                tmpa );
          vu2b = _mm512_mul_pd( _mm512_fmadd_pd( vu2b, x2_newb,
                                                 _mm512_mul_pd( vu3b, x3_newb ) ),
                                tmpb );

          vu3a = _mm512_mul_pd( vu3a, _mm512_mul_pd( r_olda, tmpa ) );
          vu3b = _mm512_mul_pd( vu3b, _mm512_mul_pd( r_oldb, tmpb ) );

          // Convert vu2, vu3 to single precision
          _MM512_CVTPD_PS( vu2, vu2a, vu2b );
          _MM512_CVTPD_PS( vu3, vu3a, vu3b );
        }

        // Store momenta
        _MM512_TRANSPOSE3v16_PS( vu1, vu2, vu3 );
        _mm512_store_ps( &bu[ 3*i      ], vu1 );
        _mm512_store_ps( &bu[ 3*i + 16 ], vu2 );
        _mm512_store_ps( &bu[ 3*i + 32 ], vu3 );

        // Store virtual particles with positions still indexed to the
        // original cell
        STOREU16P2D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vq, vv3, vix0, viy0 );

        // Find position trim, calculate crossings and store result
        __mmask16 xcross, ycross;
        __m512 vtr1, vtr2;
        VMNTRIM( vx1, vtr1, xcross );
        VMNTRIM( vy1, vtr2, ycross );

        __m512i vcross = _mm512_setzero_epi32();

        // x crossings
        vcross = _mm512_mask_or_epi32( vcross, xcross, vcross, _mm512_set1_epi32(1) );
        // y crossings
        vcross = _mm512_mask_or_epi32( vcross, ycross, vcross, _mm512_set1_epi32(2) );
        // Store vcross
        _mm512_store_epi32( (__m512i *) cross, vcross );

        // Trim positions and store results
        vx1  = _mm512_sub_ps( vx1, vtr1 );
        vy1  = _mm512_sub_ps( vy1, vtr2 );

        //    _MM512_STORE16v2_PS( x, vx1, vy1 );
        _MM512_TRANSPOSE2v16_PS( vx1, vy1 );
        _mm512_store_ps( &bx[ 2*i      ], vx1 );
        _mm512_store_ps( &bx[ 2*i + 16 ], vy1 );

        // find cell crossings and store
        __m512i vitr1 = _mm512_cvtps_epi32( vtr1 ) ;
        __m512i vitr2 = _mm512_cvtps_epi32( vtr2 ) ;

        _mm512_store_epi32( (__m512i *) dix, vitr1 );
        _mm512_store_epi32( (__m512i *) diy, vitr2 );

        // Trim cell indexes and store
        vix0 = _mm512_add_epi32( vix0, vitr1 );
        viy0 = _mm512_add_epi32( viy0, vitr2 );

        //    _MM512_STORE16v2_EPI32( ix, vix1, viy1 );
        _MM512_TRANSPOSE2v16_EPI32( vix0, viy0 );
        _mm512_store_epi32( &bix[ 2*i      ], vix0 );
        _mm512_store_epi32( &bix[ 2*i + 16 ], viy0 );

      }

      // ---------- split trajectories for current deposition
      vsplit2D( &vpbuf, cross, dix, diy );

    }

    // Deposit current from all virtual particles
    DEP_CURRENT_2D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    *ene += _mm512_reduce_add_ps( vene_group );

  }
}

/********************************** 3D advance deposit **********************************/


extern void ADVANCE_DEPOSIT_3D
( int *ix, float *x, float *u, float *q, int *i0, int *i1, float *rqm,
  float *efield, float *bfield, int *emf_size, int *emf_offset,
  float *j, int *j_size, int *j_offset,
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

  DECLARE_ALIGNED_64( float jnorm[3] );

  DECLARE_ALIGNED_64( unsigned int cross[VEC_WIDTH] );
  DECLARE_ALIGNED_64( int dix[VEC_WIDTH] );
  DECLARE_ALIGNED_64( int diy[VEC_WIDTH] );
  DECLARE_ALIGNED_64( int diz[VEC_WIDTH] );

  // This cannot be made static because of multi-threaded (OpenMP) code
  t_split_buf3D vpbuf;

  // Simulation constants
  __m512 const vtem    = _mm512_set1_ps( 0.5 * (*dt) / (*rqm) );
  __m512 const vdt_dx1 = _mm512_set1_ps( (*dt) / dx[0] );
  __m512 const vdt_dx2 = _mm512_set1_ps( (*dt) / dx[1] );
  __m512 const vdt_dx3 = _mm512_set1_ps( (*dt) / dx[2] );

  // get pointers to position 0,0,0 of each field component
  unsigned int const deltaX = 3; // 3 field components
  unsigned int const deltaY = emf_size[0] * deltaX;
  unsigned int const deltaZ = emf_size[1] * deltaY;

  float* const e1 = efield + (emf_offset[0]-OFFSET)*deltaX +
                             (emf_offset[1]-OFFSET)*deltaY +
                             (emf_offset[2]-OFFSET)*deltaZ;
  float* const e2 = e1 + 1;
  float* const e3 = e1 + 2;

  float* const b1 = bfield + (emf_offset[0]-OFFSET)*deltaX +
                             (emf_offset[1]-OFFSET)*deltaY +
                             (emf_offset[2]-OFFSET)*deltaZ;
  float* const b2 = b1  + 1;
  float* const b3 = b1  + 2;

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

   }

  if ( p_cache_size_3D % VEC_WIDTH ) {
       printf("(*error*) p_cache_size_3D is not a multiple of VEC_WIDTH = %d!\n", VEC_WIDTH );
       exit(1);
  }

  for ( k = 0; k < np_total; k += p_cache_size_3D ) {

    // Initialize energy
    __m512 vene_group = _mm512_setzero_ps();

    float* bx   =  &x[ 3*k ];
    int*   bix  = &ix[ 3*k ];
    float* bu   =  &u[ 3*k ];
    float* bq   =  &q[ k ];

    // find number of particles to process
    np = ( np_total - k < p_cache_size_3D ) ? np_total - k : p_cache_size_3D;

    // Initialize virtual particle buffer
    vpbuf.np = 0;

    // Push all particles in group
    for( i = 0; i < np; i+=VEC_WIDTH ) {
      __m512 vene;
      __m512 vx0, vx1, vy0, vy1, vz0, vz1;
      __m512i vix0, viy0, viz0;

      register __m512 ve1, ve2, ve3;
      register __m512 vb1, vb2, vb3;

      __m512 vu1, vu2, vu3;

      // Load particle positions
      vx0 = _mm512_load_ps( &bx[ 3*i      ] );
      vy0 = _mm512_load_ps( &bx[ 3*i + 16 ] );
      vz0 = _mm512_load_ps( &bx[ 3*i + 32 ] );
      _MM512_TRANSPOSE16v3_PS(vx0, vy0, vz0 );

      vix0 = _mm512_load_epi32( &bix[ 3*i      ] );
      viy0 = _mm512_load_epi32( &bix[ 3*i + 16 ] );
      viz0 = _mm512_load_epi32( &bix[ 3*i + 32 ] );
      _MM512_TRANSPOSE16v3_EPI32(vix0, viy0, viz0 );

      // Interpolate fields
      {
        __m512  hx, hy, hz;
        __m512i vidxi, vidxj, vidxk, vidxih, vidxjh, vidxkh;
        __m512 vwx[NP], vwy[NP], vwz[NP], vwxh[NP], vwyh[NP], vwzh[NP];

        // idxi = 3 * ix0
        vidxi = _mm512_add_epi32( _mm512_add_epi32( vix0, vix0 ), vix0 );
        // idxj = deltaY * iy0
        vidxj = _mm512_mullo_epi32( _mm512_set1_epi32( deltaY ), viy0 );
        // idxk = deltaZ * iz0
        vidxk = _mm512_mullo_epi32( _mm512_set1_epi32( deltaZ ), viz0 );

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

          __m512i idx1, idx2, idx3;

          idx1 = _mm512_add_epi32(_mm512_add_epi32( vidxih, vidxj  ), vidxk  );
          idx2 = _mm512_add_epi32(_mm512_add_epi32( vidxi,  vidxjh ), vidxk  );
          idx3 = _mm512_add_epi32(_mm512_add_epi32( vidxi,  vidxj  ), vidxkh );

          ve1 = _mm512_setzero_ps();
          ve2 = _mm512_setzero_ps();
          ve3 = _mm512_setzero_ps();

          for ( k3 = 0; k3 < NP; k3++ ) {
            register __m512 f1plane, f2plane, f3plane;
            f1plane = _mm512_setzero_ps();
            f2plane = _mm512_setzero_ps();
            f3plane = _mm512_setzero_ps();

            for ( k2 = 0; k2 < NP; k2++ ) {

              register __m512 f1line, f2line, f3line;
              f1line = _mm512_setzero_ps();
              f2line = _mm512_setzero_ps();
              f3line = _mm512_setzero_ps();

              for ( k1 = 0; k1 < NP; k1++ ) {
                unsigned shift =  k1*3 + k2*deltaY + k3*deltaZ ;
                f1line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx1 , e1 + shift, _MM_SCALE_4 ), vwxh[k1], f1line );
                f2line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx2 , e2 + shift, _MM_SCALE_4 ),  vwx[k1], f2line );
                f3line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx3 , e3 + shift, _MM_SCALE_4 ),  vwx[k1], f3line );
              }

              f1plane = _mm512_fmadd_ps( f1line,  vwy[k2], f1plane );
              f2plane = _mm512_fmadd_ps( f2line, vwyh[k2], f2plane );
              f3plane = _mm512_fmadd_ps( f3line,  vwy[k2], f3plane );

            }

            ve1 = _mm512_fmadd_ps( f1plane,  vwz[k3], ve1 );
            ve2 = _mm512_fmadd_ps( f2plane,  vwz[k3], ve2 );
            ve3 = _mm512_fmadd_ps( f3plane, vwzh[k3], ve3 );
          }

        }

        // Interpolate B field
        {
          unsigned k1, k2, k3;

          __m512i idx1, idx2, idx3;

          idx1 = _mm512_add_epi32(_mm512_add_epi32( vidxi,  vidxjh ), vidxkh );
          idx2 = _mm512_add_epi32(_mm512_add_epi32( vidxih, vidxj  ), vidxkh );
          idx3 = _mm512_add_epi32(_mm512_add_epi32( vidxih, vidxjh ), vidxk  );

          vb1 = _mm512_setzero_ps();
          vb2 = _mm512_setzero_ps();
          vb3 = _mm512_setzero_ps();

          for ( k3 = 0; k3 < NP; k3++ ) {
            register __m512 f1plane, f2plane, f3plane;
            f1plane = _mm512_setzero_ps();
            f2plane = _mm512_setzero_ps();
            f3plane = _mm512_setzero_ps();

            for ( k2 = 0; k2 < NP; k2++ ) {

              register __m512 f1line, f2line, f3line;
              f1line = _mm512_setzero_ps();
              f2line = _mm512_setzero_ps();
              f3line = _mm512_setzero_ps();

              for ( k1 = 0; k1 < NP; k1++ ) {
                unsigned shift =  k1*3 + k2*deltaY + k3*deltaZ ;
                f1line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx1 , b1 + shift, _MM_SCALE_4 ),  vwx[k1], f1line );
                f2line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx2 , b2 + shift, _MM_SCALE_4 ), vwxh[k1], f2line );
                f3line = _mm512_fmadd_ps( _mm512_i32gather_ps( idx3 , b3 + shift, _MM_SCALE_4 ), vwxh[k1], f3line );
              }

              f1plane = _mm512_fmadd_ps( f1line, vwyh[k2], f1plane );
              f2plane = _mm512_fmadd_ps( f2line,  vwy[k2], f2plane );
              f3plane = _mm512_fmadd_ps( f3line, vwyh[k2], f3plane );
            }

            vb1 = _mm512_fmadd_ps( f1plane, vwzh[k3], vb1 );
            vb2 = _mm512_fmadd_ps( f2plane, vwzh[k3], vb2 );
            vb3 = _mm512_fmadd_ps( f3plane,  vwz[k3], vb3 );
          }

        }

      }

      // Load momenta
      vu1 = _mm512_load_ps( &bu[ 3*i      ] );
      vu2 = _mm512_load_ps( &bu[ 3*i + 16 ] );
      vu3 = _mm512_load_ps( &bu[ 3*i + 32 ] );
      _MM512_TRANSPOSE16v3_PS( vu1, vu2, vu3 );

      // ---------- advance momenta
      VDUDT_BORIS( vtem, ve1, ve2, ve3, vb1, vb2, vb3, vu1, vu2, vu3, vene );


      // ---------- advance position and get velocities
      {
        register __m512 vrg;
        register __m512 vtrx, vtry, vtrz;
        register __m512i vitrx, vitry, vitrz;
        register __m512 vq;

        __mmask16 xcross, ycross, zcross;

        // Get 1 / \gamma
        VRGAMMA( vrg, vu1, vu2, vu3 );

        // Load charge
        vq = _mm512_load_ps( &bq[ i ] );

        // Accumulate energy
        vene_group = _mm512_fmadd_ps( vq, vene, vene_group );

        // Push positions
        vx1 = _mm512_fmadd_ps( _mm512_mul_ps( vu1, vrg ), vdt_dx1, vx0 );
        vy1 = _mm512_fmadd_ps( _mm512_mul_ps( vu2, vrg ), vdt_dx2, vy0 );
        vz1 = _mm512_fmadd_ps( _mm512_mul_ps( vu3, vrg ), vdt_dx3, vz0 );

        // Store momenta
        _MM512_TRANSPOSE3v16_PS( vu1, vu2, vu3 );
        _mm512_store_ps( &bu[ 3*i      ], vu1 );
        _mm512_store_ps( &bu[ 3*i + 16 ], vu2 );
        _mm512_store_ps( &bu[ 3*i + 32 ], vu3 );

        // Store virtual particles with positions still indexed to the
        // original cell
        STOREU16P3D( vpbuf, vpbuf.np, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix0, viy0, viz0 );

        // Find position trim, calculate crossings and store result
        VMNTRIM( vx1, vtrx, xcross );
        VMNTRIM( vy1, vtry, ycross );
        VMNTRIM( vz1, vtrz, zcross );

        __m512i vcross = _mm512_setzero_epi32();

        // x crossings
        vcross = _mm512_mask_or_epi32( vcross, xcross, vcross, _mm512_set1_epi32(1) );

        // y crossings
        vcross = _mm512_mask_or_epi32( vcross, ycross, vcross, _mm512_set1_epi32(2) );

        // z crossings
        vcross = _mm512_mask_or_epi32( vcross, zcross, vcross, _mm512_set1_epi32(4) );

        // Store vcross
        _mm512_store_epi32( (__m512i *) cross, vcross );

        // Trim positions and store results
        vx1  = _mm512_sub_ps( vx1, vtrx );
        vy1  = _mm512_sub_ps( vy1, vtry );
        vz1  = _mm512_sub_ps( vz1, vtrz );

        _MM512_TRANSPOSE3v16_PS( vx1, vy1, vz1 );
        _mm512_store_ps( &bx[ 3*i      ], vx1 );
        _mm512_store_ps( &bx[ 3*i + 16 ], vy1 );
        _mm512_store_ps( &bx[ 3*i + 32 ], vz1 );

        // find cell crossings and store
        vitrx = _mm512_cvtps_epi32( vtrx ) ;
        vitry = _mm512_cvtps_epi32( vtry ) ;
        vitrz = _mm512_cvtps_epi32( vtrz ) ;

        _mm512_store_epi32( (__m512i *) dix, vitrx );
        _mm512_store_epi32( (__m512i *) diy, vitry );
        _mm512_store_epi32( (__m512i *) diz, vitrz );

        // Trim cell indexes and store
        vix0 = _mm512_add_epi32( vix0, vitrx );
        viy0 = _mm512_add_epi32( viy0, vitry );
        viz0 = _mm512_add_epi32( viz0, vitrz );

        _MM512_TRANSPOSE3v16_EPI32( vix0, viy0, viz0 );
        _mm512_store_epi32( &bix[ 3*i      ], vix0 );
        _mm512_store_epi32( &bix[ 3*i + 16 ], viy0 );
        _mm512_store_epi32( &bix[ 3*i + 32 ], viz0 );

      }

      // ---------- split trajectories for current deposition

      vsplit3D( &vpbuf, cross, dix, diy, diz );

    }

    // Deposit current from all virtual particles
    DEP_CURRENT_3D( j, j_size, j_offset, jnorm, &vpbuf );

    // Accumulate energy from group
    *ene += _mm512_reduce_add_ps( vene_group );
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

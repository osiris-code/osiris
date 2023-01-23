/*****************************************************************************************

Charge conserving current deposition, Intel AVX2 optimized version

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-avx2.h"

#include "vector-avx2.h"
#include "splines-avx2.h"
#include "os-spec-push-avx2.h"

/*****************************************************************************************
vwl_s1
*****************************************************************************************/

inline void vwl_s1( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] )
{
  vwl[0] = _mm256_mul_ps(vqn, _mm256_sub_ps(vx1, vx0));
}

/*****************************************************************************************
vwl_s2
*****************************************************************************************/

inline void vwl_s2( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] )
{

  __m256 const c1_2 = _mm256_set1_ps( 0.5f );

  __m256 d    = _mm256_sub_ps( vx1, vx0 );
  __m256 s1_2 = _mm256_sub_ps( c1_2, vx0 );
  __m256 p1_2 = _mm256_add_ps( c1_2, vx0 );

  __m256 n    = _mm256_mul_ps( vqn, d );

  vwl[0] = _mm256_mul_ps( n, _mm256_fnmadd_ps( c1_2, d, s1_2 ));
  vwl[1] = _mm256_mul_ps( n, _mm256_fmadd_ps(  c1_2, d, p1_2 ));

}

/*****************************************************************************************
vwl_s3
*****************************************************************************************/

inline void vwl_s3( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] )
{
  __m256 const c1_2 = _mm256_set1_ps( 0.5f );
  __m256 const c3_4 = _mm256_set1_ps( 0.75f );
  __m256 const c1_3 = _mm256_set1_ps( 1.0f/3.0f );

  __m256   d = _mm256_sub_ps(  vx1, vx0 );
  __m256   s = _mm256_sub_ps( c1_2, vx0 );
  __m256   p = _mm256_add_ps( c1_2, vx0 );
  __m256 d_3 = _mm256_mul_ps( c1_3,   d );

  vwl[0] = _mm256_mul_ps( c1_2, _mm256_fmsub_ps( s, s, _mm256_mul_ps( d, _mm256_sub_ps( s, d_3 )) ) );
  vwl[1] = _mm256_fnmadd_ps( d, _mm256_add_ps( vx0, d_3 ), _mm256_fnmadd_ps( vx0, vx0, c3_4 ));
  vwl[2] = _mm256_mul_ps( c1_2, _mm256_fmadd_ps( p, p, _mm256_mul_ps( d, _mm256_add_ps( p, d_3 )) ) );

  __m256 n = _mm256_mul_ps( vqn, d );
  vwl[0] = _mm256_mul_ps( n, vwl[0] );
  vwl[1] = _mm256_mul_ps( n, vwl[1] );
  vwl[2] = _mm256_mul_ps( n, vwl[2] );
}

/*****************************************************************************************
vwl_s4
*****************************************************************************************/

inline void vwl_s4( __m256 const vqn, __m256 const vx0, __m256 const vx1, __m256 vwl[] )
{
  __m256 const c1_2 = _mm256_set1_ps( 0.5f );
  __m256 const c1_4 = _mm256_set1_ps( 0.25f );
  __m256 const c1_6 = _mm256_set1_ps( 1.0f/6.0f );
  __m256 const c3_2 = _mm256_set1_ps( 1.5f );
  __m256 const c3_4 = _mm256_set1_ps( 0.75f );
  __m256 const c1_3 = _mm256_set1_ps( 1.0f/3.0f );
  __m256 const c2_3 = _mm256_set1_ps( 2.0f/3.0f );

  __m256 d = _mm256_sub_ps( vx1, vx0 );
  __m256 s = _mm256_sub_ps( c1_2, vx0 );
  __m256 p = _mm256_add_ps( c1_2, vx0 );

  __m256 t = _mm256_sub_ps( vx0, c1_6 );
  __m256 u = _mm256_add_ps( vx0, c1_6 );

  __m256 s2 = _mm256_mul_ps( s, s );
  __m256 s3 = _mm256_mul_ps( s2, s );

  __m256 p2 = _mm256_mul_ps( p, p );
  __m256 p3 = _mm256_mul_ps( p2, p );

  __m256 d_2 = _mm256_mul_ps( d, c1_2 );
  __m256 d_4 = _mm256_mul_ps( d, c1_4 );

  vwl[0] = _mm256_mul_ps( c1_6, _mm256_fmadd_ps( d, _mm256_fmsub_ps( d, _mm256_sub_ps( s, d_4 ), _mm256_mul_ps( c3_2, s2 )  ), s3 ) );
  vwl[1] = _mm256_fmadd_ps( d, _mm256_fmadd_ps(d_2, _mm256_add_ps(t, d_4), _mm256_fmsub_ps(c3_4, _mm256_mul_ps(t, t), c1_3)),

                          _mm256_fmadd_ps(c1_2, p3, _mm256_sub_ps(c2_3, p2)));
  vwl[2] = _mm256_fnmadd_ps(d, _mm256_fmadd_ps(d_2, _mm256_add_ps(u, d_4), _mm256_fmsub_ps(c3_4, _mm256_mul_ps(u, u), c1_3)),
                          _mm256_fmadd_ps(c1_2, s3, _mm256_sub_ps(c2_3, s2)));
  vwl[3] = _mm256_mul_ps( c1_6, _mm256_fmadd_ps( d, _mm256_fmadd_ps( d, _mm256_add_ps( p, d_4 ), _mm256_mul_ps( c3_2, p2 )  ), p3 ) );

  __m256 n = _mm256_mul_ps( vqn, d );
  vwl[0] = _mm256_mul_ps( n, vwl[0] );
  vwl[1] = _mm256_mul_ps( n, vwl[1] );
  vwl[2] = _mm256_mul_ps( n, vwl[2] );
  vwl[3] = _mm256_mul_ps( n, vwl[3] );

}

/****************************************************************************************

  Generate specific current deposition functions for 1D, 2D and 3D, 1st to 4th order

****************************************************************************************/

#define __TEMPLATE__

// These macros will append _s1, _s2, etc to function names.
#define ONAME(f, o) OJOIN(f, o)
#define OJOIN(f, o) f ## _s ## o

/********************************** Linear interpolation ********************************/

// Interpolation order
#define ORDER 1
// Number of interpolation points ( ORDER + 1 )
#define NP 2
// Current grid offset (1 + ORDER/2)
#define OFFSET 1

#include __FILE__

/******************************** Quadratic interpolation *******************************/

#define ORDER 2
#define NP 3
#define OFFSET 2

#include __FILE__

/********************************** Cubic interpolation *********************************/

#define ORDER 3
#define NP 4
#define OFFSET 2

#include __FILE__

/********************************* Quartic interpolation ********************************/

#define ORDER 4
#define NP 5
#define OFFSET 3

#include __FILE__

#else

/****************************************************************************************

  Template function definitions for 2D and 3D current deposition

****************************************************************************************/

/*

Note: In 2D there was a division by 2 removed from the calculation of the perpendicular
      weights, so jnorm must be dx / dt / (2*NORM). In 3D it was a division by 3, so

      jnorm must be dx / dt / (3*NORM). To avoid repeating this every time the

      DEP_CURRENT functions are called, this is done in the ADVANCE_DEPOSIT functions.

*/

/***************   Generate Function names based on interpolation order  ****************/

#define DEP_CURRENT_1D ONAME( vdepcurrent_1d, ORDER )
#define DEP_CURRENT_2D ONAME( vdepcurrent_2d, ORDER )
#define DEP_CURRENT_3D ONAME( vdepcurrent_3d, ORDER )

#define SPLINE  ONAME( vspline, ORDER )
#define WL      ONAME( vwl, ORDER )

/****************************************************************************************
  1D current deposition
*****************************************************************************************/

/*
  Reference implementation, serial deposition
*/

#if 0

inline void DEP_CURRENT_1D

(float * const current, int const * const size, int const * const offset,

 float * const norm, t_split_buf1D * const part)
{

  typedef struct Current { float j1, j2, j3; } t_current;
  t_current *p0;

  fvec j2[NP],j3[NP];
  fvec vwl1[ORDER];

  __m256 vx0, vx1, vq, vvy, vvz;

  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

  __m256 vqnx, vqvy, vqvz;
  __m256 vs0x[NP], vs1x[NP];

  __m256 const c1_2  = _mm256_set1_ps( 0.5f );

  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 );

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);

  int const np = part -> np;

  // If number of particles in buffer is not multiple of 8 add dummy particles to the end

  if ( np % VEC_WIDTH ) {
    for( int i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       int k = np + i;

       part->x0[k] = part->x1[k] = 0.0f;
       part-> q[k] = part->vy[k] = part->vz[k] = 0.0f;

       part->ix[k] = 1;
    }
  }

  for( int i = 0; i < np; i += VEC_WIDTH ) {

    int k, k1;
    __m256i vix;

    // load 8 particles
    LOAD8P1D( part, i, vx0, vx1, vq, vvy, vvz, vix )

    // Store cell index
    _mm256_store_si256( (__m256i *) idx, vix );

    // Get splines

    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );

    vqnx = _mm256_mul_ps( vq, vnorm1 );
    vqvy = _mm256_mul_ps( vq, vvy );
    vqvz = _mm256_mul_ps( vq, vvz );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256 *) vwl1 );

    // get j2,j3 current
    for ( int k1 = 0; k1 < NP; k1++ ) {
      __m256 vwp1 = _mm256_mul_ps( c1_2, _mm256_add_ps(vs0x[k1], vs1x[k1]) );
      j2[k1].v8   = _mm256_mul_ps( vqvy, vwp1 );
      j3[k1].v8   = _mm256_mul_ps( vqvz, vwp1 );
    }

    // Loop by particle on the outside loop
    for ( int k = 0; k < VEC_WIDTH; k ++ ) {

      p0 = pj + idx[k];

      // accumulate j1
      for ( int k1 = 0; k1 < ORDER; k1++ ) {
        p0[k1].j1 += vwl1[k1].v[k] ;
      }

      // accumulate j2, j3
      for ( int k1 = 0; k1 < NP; k1++ ) {
        p0[k1].j2 += j2[k1].v[k];
        p0[k1].j3 += j3[k1].v[k];
      }

    }

  }

}

#endif

/*
  mk. I - Transpose currents, add 1 particle at a time
*/

#if 1

#warning Using mk. I 1D - algorithm

#if ( ORDER == 1 )

#define ACC_CURRENT1D( k, v, p ) {                                      \
  float *p0 = pj + idx[k];                                              \
  __m256 v0, j0;                                                      \
  v0 = _mm256_loadu_ps( (void *) p0 );                                 \
  j0 = _mm256_set_m128( _mm256_extractf128_ps( v[1], p ),              \
                        _mm256_extractf128_ps( v[0], p ) );            \
  v0 = _mm256_add_ps( v0, _mm256_permutevar8x32_ps( j0, comp_mask)); \
  _mm256_storeu_ps( (void *) p0, v0);                                  \
}

#elif ( ORDER == 2 )

#define ACC_CURRENT1D( k, v, p ) {                                      \
  float * const p0 = pj + idx[k];                                       \
  __m256 v0a, j0a;                                                      \
  __m128 v0b, j0b;                                                      \
  v0a = _mm256_loadu_ps( (void *) p0 );                                 \
  v0b = _mm_loadu_ps((void *) (p0+6));                                  \
  j0a = _mm256_set_m128( _mm256_extractf128_ps( v[1], p ),              \
                         _mm256_extractf128_ps( v[0], p ) );            \
  j0b = _mm256_extractf128_ps( v[2], p );                               \
  v0a = _mm256_add_ps( v0a, _mm256_permutevar8x32_ps( j0a, comp_mask)); \
  v0b = _mm_add_ps( v0b, j0b );                                         \
  _mm256_storeu_ps( (void *) p0, v0a);            \
  _mm_storeu_ps( (void *) (p0+6), v0b ); \
}

#elif ( ORDER == 3 )

#define ACC_CURRENT1D( k, v, p ) {                                      \
  float * const p0 = pj + idx[k];                                       \
  __m256 v0a, j0a;                                                      \
  __m256 v0b, j0b;                                                      \
  v0a = _mm256_loadu_ps( (void *) p0 );                                 \
  v0b = _mm256_loadu_ps( (void *) (p0 + 6) );                           \
  j0a = _mm256_set_m128( _mm256_extractf128_ps( v[1], p ),              \
                         _mm256_extractf128_ps( v[0], p ) );            \
  j0b = _mm256_set_m128( _mm256_extractf128_ps( v[3], p ),              \
                         _mm256_extractf128_ps( v[2], p ) );            \
  v0a = _mm256_add_ps( v0a, _mm256_permutevar8x32_ps( j0a, comp_mask)); \
  v0b = _mm256_add_ps( v0b, _mm256_permutevar8x32_ps( j0b, comp_mask)); \
  _mm256_storeu_ps( (void *) p0, v0a);            \
  _mm256_storeu_ps( (void *) (p0+6), v0b);            \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT1D( k, v, p ) {                                      \
  float * const p0 = pj + idx[k];                                       \
  __m256 v0a, j0a;                                                      \
  __m256 v0b, j0b;                                                      \
  __m128 v0c, j0c;                                                      \
  v0a = _mm256_loadu_ps( (void *) p0 );                                 \
  v0b = _mm256_loadu_ps( (void *) (p0 + 6) );                           \
  v0c = _mm_loadu_ps( (void *) (p0 + 12) );                             \
  j0a = _mm256_set_m128( _mm256_extractf128_ps( v[1], p ),              \
                         _mm256_extractf128_ps( v[0], p ) );            \
  j0b = _mm256_set_m128( _mm256_extractf128_ps( v[3], p ),              \
                         _mm256_extractf128_ps( v[2], p ) );            \
  j0c = _mm256_extractf128_ps(                  v[4], p );              \
  v0a = _mm256_add_ps( v0a, _mm256_permutevar8x32_ps( j0a, comp_mask)); \
  v0b = _mm256_add_ps( v0b, _mm256_permutevar8x32_ps( j0b, comp_mask)); \
  v0c = _mm_add_ps( v0c, j0c );                                         \
  _mm256_storeu_ps( (void *) p0, v0a);            \
  _mm256_storeu_ps( (void *) (p0+6), v0b);        \
  _mm_storeu_ps( (void *) (p0+12), v0c );         \
}

#else
#error Unsupported interpolation order
#endif

void DEP_CURRENT_1D

(float * const current, int const * const size, int const * const offset,

                       float * const norm, t_split_buf1D * const part)
{

// Mask used for compressing current values in form (j1a, j2a, j3a, 0, j1b, j2b, j3b, 0) into
// (j1a, j2a, j3a, j1b, j2b, j3b, 0, 0 )
__m256i const comp_mask  = _mm256_set_epi32( 7, 3, 6, 5, 4, 2, 1, 0 );

  __m256 vx0, vx1, vq, vvy, vvz;

  __m256 vqnx, vqvy, vqvz;
  __m256 vs0x[NP], vs1x[NP];

  __m256 const c1_2  = _mm256_set1_ps( 0.5f );

  float* const pj =  current + ( offset[0] - OFFSET ) * 3 ;

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);

  int const np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH ) {
    for( int i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       int k = np + i;

       part->x0[k] = part->x1[k] = 0.0f;
       part-> q[k] = part->vy[k] = part->vz[k] = 0.0f;

       part->ix[k] = 1;
    }
  }

  for( int i = 0; i < np; i += VEC_WIDTH ) {

    __m256 vwl1[NP];

    __m256i vix;
    DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );
    __m256 a[NP], b[NP], c[NP], d[NP];

    // load VEC_WIDTH particles
    LOAD8P1D( part, i, vx0, vx1, vq, vvy, vvz, vix )

    // Store cell index
    vix = _mm256_add_epi32( vix, _mm256_add_epi32( vix, vix ) );

    _mm256_store_si256( (__m256i *) idx, vix );

    // Get splines

    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );

    vqnx = _mm256_mul_ps( vq, vnorm1 );
    vqvy = _mm256_mul_ps( vq, vvy );
    vqvz = _mm256_mul_ps( vq, vvz );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 ); vwl1[NP-1] = _mm256_setzero_ps();

    // get currents and transpose vectors
    for ( int k1 = 0; k1 < NP; k1++ ) {
      __m256 vwp1, j1, j2, j3;

      vwp1 = _mm256_mul_ps( c1_2, _mm256_add_ps(vs0x[k1], vs1x[k1]) );

      j1  = vwl1[k1];
      j2  = _mm256_mul_ps( vqvy, vwp1 );
      j3  = _mm256_mul_ps( vqvz, vwp1 );

      // Do a 8x4 transpose, setting the last term to 0

      __m256 const zero = _mm256_setzero_ps();

      __m256 t0 = _mm256_shuffle_ps( j1, j2, 0x44 );
      __m256 t2 = _mm256_shuffle_ps( j1, j2, 0xEE );
      __m256 t1 = _mm256_shuffle_ps( j3, zero, 0x44 );
      __m256 t3 = _mm256_shuffle_ps( j3, zero, 0xEE );

      a[k1] = _mm256_shuffle_ps( t0, t1, 0x88 );
      b[k1] = _mm256_shuffle_ps( t0, t1, 0xDD );
      c[k1] = _mm256_shuffle_ps( t2, t3, 0x88 );
      d[k1] = _mm256_shuffle_ps( t2, t3, 0xDD );

    }

    ACC_CURRENT1D( 0, a, 0 )
    ACC_CURRENT1D( 1, b, 0 )
    ACC_CURRENT1D( 2, c, 0 )
    ACC_CURRENT1D( 3, d, 0 )
    ACC_CURRENT1D( 4, a, 1 )
    ACC_CURRENT1D( 5, b, 1 )
    ACC_CURRENT1D( 6, c, 1 )
    ACC_CURRENT1D( 7, d, 1 )

  }

}

#undef ACC_CURRENT1D

#endif

/****************************************************************************************
  2D current deposition
*****************************************************************************************/

#if 0

/*
  Reference implementation, serial deposition
*/

#warning Using Reference 2D algorithm

void DEP_CURRENT_2D

(float * const current, int const * const size, int const * const offset,

                       float * const norm, t_split_buf2D * const part)
{

  typedef struct Current { float j1, j2, j3; } t_current;
  t_current *pjpart, *p0;

  fvec j3[NP][NP];
  fvec vwp1[NP], vwp2[NP];

  fvec vwl1[ORDER], vwl2[ORDER];

  __m256 vx0, vx1, vy0, vy1, vq, vvz;

  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

  __m256 vqnx, vqny, vqvz;
  __m256 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP];

  __m256 const c1_3 = _mm256_set1_ps( 1.0f/3.0f );
  __m256 const c1_2 = _mm256_set1_ps( 0.5f );

  int const Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 +

                                                 ( offset[1] - OFFSET ) * 3 * Dy );

  // Setting this to 1/4 removes 1 multiplication from wl2
  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  __m256 const vnorm2 = _mm256_set1_ps(norm[1]);

  int const np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH != 0 ) {
    for( int i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       int k = np + i;

       part->x0[k] = part->x1[k] = 0.0f;
       part->y0[k] = part->y1[k] = 0.0f;
       part-> q[k] = part->vz[k] = 0.0f;

       part->ix[k] = part->iy[k] = 1;
    }
  }

  // Each virtual particle uses 8 doubles, so 4 vp = 32 doubles
  for( int i = 0; i < np; i += VEC_WIDTH ) {

    __m256i vix, viy;

    // load VEC_WIDTH particles
    LOAD8P2D( part, i, vx0, vx1, vy0, vy1, vq, vvz, vix, viy )

    // Store cell index
    _mm256_store_si256( (__m256i *) idx,

            _mm256_add_epi32( vix, _mm256_mullo_epi32( viy, _mm256_set1_epi32( Dy ) ) ) );

    // Get splines

    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );

    vqnx = _mm256_mul_ps( vq, vnorm1 );
    vqny = _mm256_mul_ps( vq, vnorm2 );
    vqvz = _mm256_mul_ps( vq, vvz);
    vqvz = _mm256_mul_ps( vqvz, c1_3 );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256 *) vwl1 );
    WL( vqny, vy0, vy1, (__m256 *) vwl2 );

    // get perpendicular weights
    for( int k = 0; k < NP; k++ ) {
      vwp1[k].v8 = _mm256_add_ps(vs0y[k], vs1y[k]);
      vwp2[k].v8 = _mm256_add_ps(vs0x[k], vs1x[k]);
    }

    // get j3 current
    for( int k2 = 0; k2 < NP; k2++ ) {
      for ( int k1 = 0; k1 < NP; k1++ ) {
        __m256 s00, s10, tmp1, tmp2;

    		s00  = _mm256_mul_ps( vs0x[k1], vs0y[k2] );
    		tmp1 = _mm256_fmadd_ps( vs1y[k2], vs1x[k1], s00 );  // tmp1 = s0x*s0y + s1x*s1y

    		s10  = _mm256_mul_ps( vs0x[k1], vs1y[k2] );
    		tmp2 = _mm256_fmadd_ps( vs0y[k2], vs1x[k1], s10 );  // tmp2 = s0x*s1y + s1x*s0y

    		tmp1 = _mm256_fmadd_ps( tmp2, c1_2, tmp1 );       // tmp1 = tmp1 + 0.5*tmp2

    		j3[k1][k2].v8 = _mm256_mul_ps( vqvz, tmp1 );          // j3 = vqvz * tmp1
      }
    }

    // Loop by particle on the outside loop
    for ( int k = 0; k < VEC_WIDTH; k ++ ) {

      pjpart = pj + idx[k];

      // accumulate j1
      for( int k2 = 0; k2 < NP; k2++ ) {
        p0 = pjpart + k2*Dy;
        for ( int k1 = 0; k1 < ORDER; k1++ ) {
          p0[k1].j1 += vwl1[k1].v[k] * vwp1[k2].v[k];
        }
      }

      // accumulate j2 - making k2 the outside loop gives marginal perf. gain
      for( int k2 = 0; k2 < ORDER; k2++ ) {
        p0 = pjpart + k2*Dy;
        for ( int k1 = 0; k1 < NP; k1++ ) {
          p0[k1].j2 += vwl2[k2].v[k] * vwp2[k1].v[k];
        }
      }

      // accumulate j3
      for( int k2 = 0; k2 < NP; k2++ ) {
        p0 = pjpart + k2*Dy;
        for ( int k1 = 0; k1 < NP; k1++ ) {
          p0[k1].j3 += (j3[k1][k2]).v[k];
        }
      }

    }

  }

}
#endif

#if 1
/*

mk. I - Transpose currents, save 1 particle at a time

process each line sequentially
*/

#warning Using mk. I 2D algorithm

#if ( ORDER == 1 )

#define ACC_CURRENT2D( k, v, p ) { \
  float * const line = pj + idx[k]; \
  for( int k2 = 0; k2 < 2; k2++ ){ \
    __m256 va = _mm256_loadu_ps( &line[ k2*Dy ] ); \
    __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k2][1], p ), \
                                 _mm256_extractf128_ps( v[k2][0], p )); \
    \
    va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask )); \
    _mm256_storeu_ps( &line[ k2*Dy ] , va ); \
  } \
}

#elif ( ORDER == 2 )

#define ACC_CURRENT2D( k, v, p ) { \
  float * const line = pj + idx[k]; \
  for( int k2 = 0; k2 < 3; k2++ ){ \
    __m256 va = _mm256_loadu_ps( &line[ k2*Dy     ] ); \
    __m128 vb = _mm_loadu_ps(    &line[ k2*Dy + 6 ] ); \
    \
    __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k2][1], p ), \
                                 _mm256_extractf128_ps( v[k2][0], p )); \
    __m128 jb = _mm256_extractf128_ps(                  v[k2][2], p ); \
    \
    va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask)); \
    vb = _mm_add_ps(    vb, jb ); \
    _mm256_storeu_ps( &line[ k2*Dy     ] , va ); \
    _mm_storeu_ps(    &line[ k2*Dy + 6 ], vb ); \
  } \
}

#elif ( ORDER == 3 )

#define ACC_CURRENT2D( k, v, p ) { \
  float * const line = pj + idx[k]; \
  for( int k2 = 0; k2 < 4; k2++ ){ \
    __m256 va = _mm256_loadu_ps( &line[ k2*Dy     ] ); \
    __m256 vb = _mm256_loadu_ps( &line[ k2*Dy + 6 ] ); \
    \
    __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k2][1], p ), \
                                 _mm256_extractf128_ps( v[k2][0], p )); \
    __m256 jb = _mm256_set_m128( _mm256_extractf128_ps( v[k2][3], p ), \
                                 _mm256_extractf128_ps( v[k2][2], p )); \
    \
    va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask)); \
    vb = _mm256_add_ps( vb, _mm256_permutevar8x32_ps( jb, comp_mask)); \
    \
    _mm256_storeu_ps( &line[ k2*Dy     ], va ); \
    _mm256_storeu_ps( &line[ k2*Dy + 6 ], vb ); \
  } \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT2D( k, v, p ) { \
  float * const line = pj + idx[k]; \
  \
  for( int k2 = 0; k2 < 5; k2++ ){ \
    __m256 va = _mm256_loadu_ps( &line[ k2*Dy      ] ); \
    __m256 vb = _mm256_loadu_ps( &line[ k2*Dy +  6 ] ); \
    __m128 vc =    _mm_loadu_ps( &line[ k2*Dy + 12 ] ); \
    \
    __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k2][1], p ), \
                                 _mm256_extractf128_ps( v[k2][0], p )); \
    __m256 jb = _mm256_set_m128( _mm256_extractf128_ps( v[k2][3], p ), \
                                 _mm256_extractf128_ps( v[k2][2], p )); \
    __m128 jc = _mm256_extractf128_ps(                  v[k2][4], p ); \
    \
    va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask)); \
    vb = _mm256_add_ps( vb, _mm256_permutevar8x32_ps( jb, comp_mask)); \
    vc = _mm_add_ps(    vc, jc ); \
    \
    _mm256_storeu_ps( &line[ k2*Dy      ], va ); \
    _mm256_storeu_ps( &line[ k2*Dy +  6 ], vb ); \
       _mm_storeu_ps( &line[ k2*Dy + 12 ], vc ); \
  } \
}

#else
#error Unsupported interpolation order
#endif

void DEP_CURRENT_2D

(float * restrict const current, int const * restrict const size, int const * restrict const offset,

                       float * restrict const norm, t_split_buf2D * restrict const part)
{

// Mask used for compressing current values in form (j1a, j2a, j3a, 0, j1b, j2b, j3b, 0) into
// (j1a, j2a, j3a, j1b, j2b, j3b, 0, 0 )
__m256i const comp_mask  = _mm256_set_epi32( 7, 3, 6, 5, 4, 2, 1, 0 );

  __m256 vwp1[NP], vwp2[NP];

  __m256 vwl1[NP], vwl2[NP];

  __m256 vx0, vx1, vy0, vy1, vq, vvz;

  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

  __m256 vqnx, vqny, vqvz;
  __m256 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP];

  __m256 const c1_3 = _mm256_set1_ps( 1.0f/3.0f );
  __m256 const c1_2  = _mm256_set1_ps( 0.5f );

  // Current array has 3 field components
  float* const pj = current + 3 * ( ( offset[0] - OFFSET ) +

                                    ( offset[1] - OFFSET ) * size[0] );

  // Normalization for current
  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  __m256 const vnorm2 = _mm256_set1_ps(norm[1]);

  int const Dy = 3*size[0];
  __m256i const vDy = _mm256_set1_epi32( Dy );

  int const np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH != 0 ) {
    for( int i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       int k = np + i;

       part->x0[k] = part->x1[k] = 0.0f;
       part->y0[k] = part->y1[k] = 0.0f;
       part-> q[k] = part->vz[k] = 0.0f;

       part->ix[k] = part->iy[k] = 1;
    }
  }

  for( int i = 0; i < np; i += VEC_WIDTH ) {

    __m256i vix, viy;

    // load VEC_WIDTH particles
    LOAD8P2D( part, i, vx0, vx1, vy0, vy1, vq, vvz, vix, viy )

    // Store cell index
    // idx = 3 * vix + 3 * size[0] * viy
    vix = _mm256_add_epi32( vix, _mm256_add_epi32( vix, vix ) );
    viy = _mm256_mullo_epi32( viy, vDy );

    _mm256_store_si256( (__m256i *) idx, _mm256_add_epi32( vix, viy ) );

    // Get splines

    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );

    vqnx = _mm256_mul_ps( vq, vnorm1 );
    vqny = _mm256_mul_ps( vq, vnorm2 );
    vqvz = _mm256_mul_ps( vq, vvz);
    vqvz = _mm256_mul_ps( vqvz, c1_3 );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 ); vwl1[NP-1] = _mm256_setzero_ps();

    WL( vqny, vy0, vy1, vwl2 ); vwl2[NP-1] = _mm256_setzero_ps();

    // get perpendicular weights
    for( int k = 0; k < NP; k++ ) {
      vwp1[k] = _mm256_add_ps(vs0y[k], vs1y[k]);
      vwp2[k] = _mm256_add_ps(vs0x[k], vs1x[k]);
    }

    __m256 a[NP][NP], b[NP][NP], c[NP][NP], d[NP][NP];

    for( int k2 = 0; k2 < NP; k2++ ) {
      for ( int k1 = 0; k1 < NP; k1++ ) {
        __m256 j1, j2, j3;

        j1 = _mm256_mul_ps( vwl1[k1] , vwp1[k2] );
        j2 = _mm256_mul_ps( vwl2[k2] , vwp2[k1] );

        __m256 s00, s10, tmp1, tmp2;

        s00  = _mm256_mul_ps( vs0x[k1], vs0y[k2] );
        tmp1 = _mm256_fmadd_ps( vs1y[k2], vs1x[k1], s00 );

        s10  = _mm256_mul_ps( vs0x[k1], vs1y[k2] );
        tmp2 = _mm256_fmadd_ps( vs0y[k2], vs1x[k1], s10 );

        tmp1 = _mm256_fmadd_ps( tmp2, c1_2, tmp1 );

        j3 = _mm256_mul_ps( vqvz, tmp1 );

        // Do a 4x4 transpose inside each 128 bit lane, setting the 4th component to 0
        __m256 const zero = _mm256_setzero_ps();

        __m256 t0 = _mm256_shuffle_ps( j1, j2, 0x44 );
        __m256 t2 = _mm256_shuffle_ps( j1, j2, 0xEE );
        __m256 t1 = _mm256_shuffle_ps( j3, zero, 0x44 );
        __m256 t3 = _mm256_shuffle_ps( j3, zero, 0xEE );

        a[k2][k1] = _mm256_shuffle_ps( t0, t1, 0x88 );
        b[k2][k1] = _mm256_shuffle_ps( t0, t1, 0xDD );
        c[k2][k1] = _mm256_shuffle_ps( t2, t3, 0x88 );
        d[k2][k1] = _mm256_shuffle_ps( t2, t3, 0xDD );
      }
    }

    ACC_CURRENT2D( 0, a, 0 )
    ACC_CURRENT2D( 1, b, 0 )
    ACC_CURRENT2D( 2, c, 0 )
    ACC_CURRENT2D( 3, d, 0 )
    ACC_CURRENT2D( 4, a, 1 )
    ACC_CURRENT2D( 5, b, 1 )
    ACC_CURRENT2D( 6, c, 1 )
    ACC_CURRENT2D( 7, d, 1 )

  }

}

#undef ACC_CURRENT2D

#endif

/****************************************************************************************
  3D current deposition
****************************************************************************************/

/*

  Reference implementation, serial deposition

*/

#if 0

#warning Using reference 3D algorithm

void DEP_CURRENT_3D
(float * const  current, int const * const  size, int const * const  offset,
                       float * const norm, t_split_buf3D * const part)
{

  typedef struct Current { float j1, j2, j3; } t_current;

  __m256 vx0, vx1, vy0, vy1,  vz0, vz1, vq;

  __m256 vqnx, vqny, vqnz;
  __m256 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP];

  fvec vwl1[ORDER], vwl2[ORDER], vwl3[ORDER];
  fvec vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];

  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

  int const Dy = size[0];

  int const Dz = Dy * size[1];

  t_current*  pj = (t_current *) ( current + ( offset[0] - OFFSET  ) * 3 +

                                             ( offset[1] - OFFSET  ) * 3 * Dy +

                                             ( offset[2] - OFFSET  ) * 3 * Dz );

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  __m256 const vnorm2 = _mm256_set1_ps(norm[1]);
  __m256 const vnorm3 = _mm256_set1_ps(norm[2]);

  int np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np % VEC_WIDTH ) {
    for( int i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
      int k = np + i;

      part -> x0[k] = part -> x1[k] = 0.;
      part -> y0[k] = part -> y1[k] = 0.;
      part -> z0[k] = part -> z1[k] = 0.;
      part -> q[k] = 0.;
      part -> ix[k] = part -> iy[k] = part -> iz[k] = 1;
    }
  }

  for( int i = 0; i < np; i+=VEC_WIDTH ) {

    __m256i vix, viy, viz, vidx;

    // load VEC_WIDTH particles
    LOAD8P3D( part, i, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz )

    // Store cell index
    vidx = _mm256_add_epi32( vix,  _mm256_mullo_epi32( viy, _mm256_set1_epi32( Dy ) ) );
    vidx = _mm256_add_epi32( vidx, _mm256_mullo_epi32( viz, _mm256_set1_epi32( Dz ) ) );

    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );

    vqnx = _mm256_mul_ps( vq, vnorm1 );
    vqny = _mm256_mul_ps( vq, vnorm2 );
    vqnz = _mm256_mul_ps( vq, vnorm3 );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256 *) vwl1 );
    WL( vqny, vy0, vy1, (__m256 *) vwl2 );
    WL( vqnz, vz0, vz1, (__m256 *) vwl3 );

    // get perpendicular weights
    for( int k2 = 0; k2 < NP; k2++ ) {
      for( int k1 = 0; k1 < NP; k1++ ) {

        const __m256 c1_2 = _mm256_set1_ps(0.5f);
        __m256 tmp1, tmp2;

        // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] +

        //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )

        tmp1 = _mm256_mul_ps( vs0y[k1], vs0z[k2] );
        tmp2 = _mm256_mul_ps( vs0y[k1], vs1z[k2] );
        tmp1 = _mm256_fmadd_ps( vs1y[k1], vs1z[k2], tmp1 );
        tmp2 = _mm256_fmadd_ps( vs1y[k1], vs0z[k2], tmp2 );

        vwp1[k2][k1].v8 = _mm256_fmadd_ps( c1_2, tmp2, tmp1 );

        // wp2[k2][k1] = s0x[k1]*s0z[k2] + s1x[k1]*s1z[k2] +

        //               0.5*( s0x[k1]*s1z[k2] + s1x[k1]*s0z[k2] )

        tmp1 = _mm256_mul_ps( vs0x[k1], vs0z[k2] );
        tmp2 = _mm256_mul_ps( vs0x[k1], vs1z[k2] );
        tmp1 = _mm256_fmadd_ps( vs1x[k1], vs1z[k2], tmp1 );
        tmp2 = _mm256_fmadd_ps( vs1x[k1], vs0z[k2], tmp2 );

        vwp2[k2][k1].v8 = _mm256_fmadd_ps( c1_2, tmp2, tmp1 );

        // wp3[k2][k1] = s0x[k1]*s0y[k2] + s1x[k1]*s1y[k2] +

        //               0.5*( s0x[k1]*s1y[k2] + s1x[k1]*s0y[k2] )

        tmp1 = _mm256_mul_ps( vs0x[k1], vs0y[k2] );
        tmp2 = _mm256_mul_ps( vs0x[k1], vs1y[k2] );
        tmp1 = _mm256_fmadd_ps( vs1x[k1], vs1y[k2], tmp1 );
        tmp2 = _mm256_fmadd_ps( vs1x[k1], vs0y[k2], tmp2 );

        vwp3[k2][k1].v8 = _mm256_fmadd_ps( c1_2, tmp2, tmp1 );

      }
    }

    // The following code will need to access the individual vidx values
    _mm256_store_si256( (__m256i *) idx, vidx );

    // Loop by particle on the outside loop
    for ( int k = 0; k < VEC_WIDTH; k ++ ) {

      t_current * pjpart = pj + idx[k];

      // accumulate j1
      for( int k3 = 0; k3 < NP; k3++ ) {
        for( int k2 = 0; k2 < NP; k2++ ) {
          t_current *p0 = pjpart + k2*Dy + k3*Dz;
          for ( int k1 = 0; k1 < ORDER; k1++ ) {
            p0[k1].j1 += vwl1[k1].v[k] * vwp1[k3][k2].v[k];
          }
        }
      }

      // accumulate j2
      for( int k3 = 0; k3 < NP; k3++ ) {
        for( int k2 = 0; k2 < ORDER; k2++ ) {
          t_current *p0 = pjpart + k2*Dy + k3*Dz;
          for ( int k1 = 0; k1 < NP; k1++ ) {
            p0[k1].j2 += vwl2[k2].v[k] * vwp2[k3][k1].v[k];
          }
        }
      }

      // accumulate j3
      for( int k3 = 0; k3 < ORDER; k3++ ) {
        for( int k2 = 0; k2 < NP; k2++ ) {
          t_current *p0 = pjpart + k2*Dy + k3*Dz;
          for ( int k1 = 0; k1 < NP; k1++ ) {
            p0[k1].j3 += vwl3[k3].v[k] * vwp3[k2][k1].v[k];
          }
        }
      }

    }

  }
}

#endif

#if 1

/* mk. I

For each vector of particles loop over z coordinate and deposit one current plane at a
time using the same algorithm as 2D mk. I

*/

#warning Using mk. I 3D algorithm

#if ( ORDER == 1 )

#define ACC_CURRENT3D( k, v, p ) { \
  float * const line = pj + idx[k] + k3*Dz; \
  \
  for( int k2 = 0; k2 < 2; k2++ ){ \
    __m256 va = _mm256_loadu_ps( &line[ k2*Dy ] ); \
    __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k2][1], p ), \
                                 _mm256_extractf128_ps( v[k2][0], p )); \
    \
    va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask )); \
    _mm256_storeu_ps( &line[ k2*Dy ] , va ); \
  } \
}

#elif ( ORDER == 2 )

#define ACC_CURRENT3D( k, v, p ) { \
  float * const line = pj + idx[k] + k3*Dz; \
  \
  for( int k2 = 0; k2 < 3; k2++ ){ \
    __m256 va = _mm256_loadu_ps( &line[ k2*Dy     ] ); \
    __m128 vb = _mm_loadu_ps(    &line[ k2*Dy + 6 ] ); \
    \
    __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k2][1], p ), \
                                 _mm256_extractf128_ps( v[k2][0], p )); \
    __m128 jb = _mm256_extractf128_ps(                  v[k2][2], p ); \
    \
    va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask)); \
    vb = _mm_add_ps(    vb, jb ); \
    _mm256_storeu_ps( &line[ k2*Dy     ], va ); \
    _mm_storeu_ps(    &line[ k2*Dy + 6 ], vb ); \
  } \
}

#elif ( ORDER == 3 )

#define ACC_CURRENT3D( k, v, p ) { \
  float * const line = pj + idx[k] + k3*Dy; \
  \
  for( int k2 = 0; k2 < 4; k2++ ){ \
    __m256 va = _mm256_loadu_ps( &line[ k2*Dy     ] ); \
    __m256 vb = _mm256_loadu_ps( &line[ k2*Dy + 6 ] ); \
    \
    __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k2][1], p ),  \
                                 _mm256_extractf128_ps( v[k2][0], p )); \
    __m256 jb = _mm256_set_m128( _mm256_extractf128_ps( v[k2][3], p ),  \
                                 _mm256_extractf128_ps( v[k2][2], p )); \
    \
    va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask)); \
    vb = _mm256_add_ps( vb, _mm256_permutevar8x32_ps( jb, comp_mask)); \
    \
    _mm256_storeu_ps( &line[ k2*Dy     ], va ); \
    _mm256_storeu_ps( &line[ k2*Dy + 6 ], vb ); \
  } \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT3D( k, v, p ) { \
  float * const line = pj + idx[k] + k3*Dy; \
  \
  for( int k2 = 0; k2 < 5; k2++ ){ \
    __m256 va = _mm256_loadu_ps( &line[ k2*Dy      ] ); \
    __m256 vb = _mm256_loadu_ps( &line[ k2*Dy +  6 ] ); \
    __m128 vc =    _mm_loadu_ps( &line[ k2*Dy + 12 ] ); \
    \
    __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k2][1], p ), \
                                 _mm256_extractf128_ps( v[k2][0], p )); \
    __m256 jb = _mm256_set_m128( _mm256_extractf128_ps( v[k2][3], p ), \
                                 _mm256_extractf128_ps( v[k2][2], p )); \
    __m128 jc = _mm256_extractf128_ps(                  v[k2][4], p ); \
    \
    va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask)); \
    vb = _mm256_add_ps( vb, _mm256_permutevar8x32_ps( jb, comp_mask)); \
    vc = _mm_add_ps(    vc, jc ); \
    \
    _mm256_storeu_ps( &line[ k2*Dy      ], va ); \
    _mm256_storeu_ps( &line[ k2*Dy +  6 ], vb ); \
       _mm_storeu_ps( &line[ k2*Dy + 12 ], vc ); \
  } \
}

#else
#error Unsupported interpolation order
#endif

void DEP_CURRENT_3D
(float * const  current, int const * const  size, int const * const  offset,
                       float * const norm, t_split_buf3D * const part)
{

// Mask used for compressing current values in form (j1a, j2a, j3a, 0, j1b, j2b, j3b, 0) into
// (j1a, j2a, j3a, j1b, j2b, j3b, 0, 0 )
__m256i const comp_mask  = _mm256_set_epi32( 7, 3, 6, 5, 4, 2, 1, 0 );

  __m256 vx0, vx1, vy0, vy1,  vz0, vz1, vq;

  __m256 vqnx, vqny, vqnz;
  __m256 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP];

  __m256 vwl1[NP], vwl2[NP], vwl3[NP];
  __m256 vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];

  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

  int const Dy = 3 * size[0];

  int const Dz = Dy * size[1];

  __m256i const vDy = _mm256_set1_epi32( Dy );

  __m256i const vDz = _mm256_set1_epi32( Dz );

  float* const pj = current + ( offset[0] - OFFSET  ) * 3 +

  ( offset[1] - OFFSET  ) * Dy +

  ( offset[2] - OFFSET  ) * Dz;

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  __m256 const vnorm2 = _mm256_set1_ps(norm[1]);
  __m256 const vnorm3 = _mm256_set1_ps(norm[2]);

  int np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np % VEC_WIDTH ) {
    for( int i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
      int k = np + i;

      part -> x0[k] = part -> x1[k] = 0.;
      part -> y0[k] = part -> y1[k] = 0.;
      part -> z0[k] = part -> z1[k] = 0.;
      part -> q[k] = 0.;
      part -> ix[k] = part -> iy[k] = part -> iz[k] = 1;
    }
  }

  for( int i = 0; i < np; i+=VEC_WIDTH ) {

    __m256i vix, viy, viz;

    // load VEC_WIDTH particles
    LOAD8P3D( part, i, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz )

    // Store cell index
    vix = _mm256_add_epi32( vix, _mm256_add_epi32( vix, vix ) );
    viy = _mm256_mullo_epi32( viy, vDy );

    viz = _mm256_mullo_epi32( viz, vDz );

    _mm256_store_si256( (__m256i *) idx, _mm256_add_epi32( _mm256_add_epi32( vix, viy ), viz ) );

    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );

    vqnx = _mm256_mul_ps( vq, vnorm1);
    vqny = _mm256_mul_ps( vq, vnorm2);
    vqnz = _mm256_mul_ps( vq, vnorm3);

    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 ); vwl1[NP-1] = _mm256_setzero_ps();
    WL( vqny, vy0, vy1, vwl2 ); vwl2[NP-1] = _mm256_setzero_ps();
    WL( vqnz, vz0, vz1, vwl3 ); vwl3[NP-1] = _mm256_setzero_ps();

    // get perpendicular weights
    for( int k2 = 0; k2 < NP; k2++ ) {
      for( int k1 = 0; k1 < NP; k1++ ) {

        const __m256 c1_2 = _mm256_set1_ps(0.5f);
        __m256 tmp1, tmp2;

        // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] +

        //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )

        tmp1 = _mm256_mul_ps( vs0y[k1], vs0z[k2] );
        tmp2 = _mm256_mul_ps( vs0y[k1], vs1z[k2] );
        tmp1 = _mm256_fmadd_ps( vs1y[k1], vs1z[k2], tmp1 );
        tmp2 = _mm256_fmadd_ps( vs1y[k1], vs0z[k2], tmp2 );

        vwp1[k2][k1] = _mm256_fmadd_ps( c1_2, tmp2, tmp1 );

        // wp2[k2][k1] = s0x[k1]*s0z[k2] + s1x[k1]*s1z[k2] +

        //               0.5*( s0x[k1]*s1z[k2] + s1x[k1]*s0z[k2] )

        tmp1 = _mm256_mul_ps( vs0x[k1], vs0z[k2] );
        tmp2 = _mm256_mul_ps( vs0x[k1], vs1z[k2] );
        tmp1 = _mm256_fmadd_ps( vs1x[k1], vs1z[k2], tmp1 );
        tmp2 = _mm256_fmadd_ps( vs1x[k1], vs0z[k2], tmp2 );

        vwp2[k2][k1] = _mm256_fmadd_ps( c1_2, tmp2, tmp1 );

        // wp3[k2][k1] = s0x[k1]*s0y[k2] + s1x[k1]*s1y[k2] +

        //               0.5*( s0x[k1]*s1y[k2] + s1x[k1]*s0y[k2] )

        tmp1 = _mm256_mul_ps( vs0x[k1], vs0y[k2] );
        tmp2 = _mm256_mul_ps( vs0x[k1], vs1y[k2] );
        tmp1 = _mm256_fmadd_ps( vs1x[k1], vs1y[k2], tmp1 );
        tmp2 = _mm256_fmadd_ps( vs1x[k1], vs0y[k2], tmp2 );

        vwp3[k2][k1] = _mm256_fmadd_ps( c1_2, tmp2, tmp1 );

      }
    }

    // Accumulate current 1 plane at a time using the 2D algorithm
    for( int k3 = 0; k3 < NP; k3++ ) {
      __m256 a[NP][NP], b[NP][NP], c[NP][NP], d[NP][NP];

      for( int k2 = 0; k2 < NP; k2++ ) {

        for ( int k1 = 0; k1 < NP; k1++ ) {
          __m256 j1, j2, j3;

          // Calculate current components
          j1 = _mm256_mul_ps( vwl1[k1] , vwp1[k3][k2] );
          j2 = _mm256_mul_ps( vwl2[k2] , vwp2[k3][k1] );
          j3 = _mm256_mul_ps( vwl3[k3] , vwp3[k2][k1] );

          // Do a 8x4 transpose, setting the last term to 0

          __m256 const zero = _mm256_setzero_ps();

          __m256 t0 = _mm256_shuffle_ps( j1, j2, 0x44 );
          __m256 t2 = _mm256_shuffle_ps( j1, j2, 0xEE );
          __m256 t1 = _mm256_shuffle_ps( j3, zero, 0x44 );
          __m256 t3 = _mm256_shuffle_ps( j3, zero, 0xEE );

          a[k2][k1] = _mm256_shuffle_ps( t0, t1, 0x88 );
          b[k2][k1] = _mm256_shuffle_ps( t0, t1, 0xDD );
          c[k2][k1] = _mm256_shuffle_ps( t2, t3, 0x88 );
          d[k2][k1] = _mm256_shuffle_ps( t2, t3, 0xDD );
        }
      }

      // Accumulate electric current in k3 plane
      ACC_CURRENT3D( 0, a, 0 )
      ACC_CURRENT3D( 1, b, 0 )
      ACC_CURRENT3D( 2, c, 0 )
      ACC_CURRENT3D( 3, d, 0 )
      ACC_CURRENT3D( 4, a, 1 )
      ACC_CURRENT3D( 5, b, 1 )
      ACC_CURRENT3D( 6, c, 1 )
      ACC_CURRENT3D( 7, d, 1 )

    }

  }

}

#undef ACC_CURRENT3D

#endif

#if 0

/* mk. II

Pre calculate all currents, deposit 1 particle at a time

*/

#warning Using mk. II 3D algorithm

#if ( ORDER == 1 )

#define ACC_CURRENT3D( k, v, p ) { \
  float * const line = pj + idx[k]; \
  \
  for( int k3 = 0; k3 < 2; k3++ ){ \
    for( int k2 = 0; k2 < 2; k2++ ){ \
      __m256 va = _mm256_loadu_ps( &line[ k2*Dy + k3*Dz ] ); \
      __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k3][k2][0], p ), \
                                   _mm256_extractf128_ps( v[k3][k2][1], p )); \
      \
      va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask )); \
      _mm256_storeu_ps( &line[ k2*Dy + k3*Dz ] , va ); \
    } \
  } \
}

#elif ( ORDER == 2 )

#define ACC_CURRENT3D( k, v, p ) { \
  float * const line = pj + idx[k]; \
  \
  for( int k3 = 0; k3 < 3; k3++ ){ \
    for( int k2 = 0; k2 < 3; k2++ ){ \
      __m256 va = _mm256_loadu_ps( &line[ k2*Dy + k3*Dz     ] ); \
      __m128 vb = _mm_loadu_ps(    &line[ k2*Dy + k3*Dz + 6 ] ); \
      \
      __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k3][k2][0], p ), \
                                   _mm256_extractf128_ps( v[k3][k2][1], p )); \
      __m128 jb = _mm256_extractf128_ps(                  v[k3][k2][2], p ); \
      \
      va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask)); \
      vb = _mm_add_ps(    vb, jb ); \
      _mm256_storeu_ps( &line[ k2*Dy + k3*Dz     ] , va ); \
      _mm_storeu_ps(    &line[ k2*Dy + k3*Dz + 6 ], vb ); \
    } \
  } \
}

#elif ( ORDER == 3 )

#define ACC_CURRENT3D( k, v, p ) { \
  float * const line = pj + idx[k]; \
  \
  for( int k3 = 0; k3 < 4; k3++ ){ \
    for( int k2 = 0; k2 < 4; k2++ ){ \
      __m256 va = _mm256_loadu_ps( &line[ k2*Dy + k3*Dz     ] ); \
      __m256 vb = _mm256_loadu_ps( &line[ k2*Dy + k3*Dz + 6 ] ); \
      \
      __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k3][k2][0], p ),  \
                                   _mm256_extractf128_ps( v[k3][k2][1], p )); \
      __m256 jb = _mm256_set_m128( _mm256_extractf128_ps( v[k3][k2][2], p ),  \
                                   _mm256_extractf128_ps( v[k3][k2][3], p )); \
      \
      va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask)); \
      vb = _mm256_add_ps( vb, _mm256_permutevar8x32_ps( jb, comp_mask)); \
      \
      _mm256_storeu_ps( &line[ k2*Dy + k3*Dz     ], va ); \
      _mm256_storeu_ps( &line[ k2*Dy + k3*Dz + 6 ], vb ); \
    } \
  } \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT3D( k, v, p ) { \
  float * const line = pj + idx[k]; \
  \
  for( int k3 = 0; k3 < 5; k3++ ){ \
    for( int k2 = 0; k2 < 5; k2++ ){ \
      __m256 va = _mm256_loadu_ps( &line[ k2*Dy + k3*Dz     ] ); \
      __m256 vb = _mm256_loadu_ps( &line[ k2*Dy + k3*Dz+  6 ] ); \
      __m128 vc =    _mm_loadu_ps( &line[ k2*Dy + k3*Dz+ 12 ] ); \
      \
      __m256 ja = _mm256_set_m128( _mm256_extractf128_ps( v[k3][k2][0], p ), \
                                   _mm256_extractf128_ps( v[k3][k2][1], p )); \
      __m256 jb = _mm256_set_m128( _mm256_extractf128_ps( v[k3][k2][2], p ), \
                                   _mm256_extractf128_ps( v[k3][k2][3], p )); \
      __m128 jc = _mm256_extractf128_ps(                  v[k3][k2][4], p ); \
      \
      va = _mm256_add_ps( va, _mm256_permutevar8x32_ps( ja, comp_mask)); \
      vb = _mm256_add_ps( vb, _mm256_permutevar8x32_ps( jb, comp_mask)); \
      vc = _mm_add_ps(    vc, jc ); \
      \
      _mm256_storeu_ps( &line[ k2*Dy + k3*Dz      ], va ); \
      _mm256_storeu_ps( &line[ k2*Dy + k3*Dz +  6 ], vb ); \
         _mm_storeu_ps( &line[ k2*Dy + k3*Dz + 12 ], vc ); \
    } \
  } \
}

#else
#error Unsupported interpolation order
#endif

void DEP_CURRENT_3D
(float * const  current, int const * const  size, int const * const  offset,
                       float * const norm, t_split_buf3D * const part)
{

// Mask used for compressing current values in form (j1a, j2a, j3a, 0, j1b, j2b, j3b, 0) into
// (j1a, j2a, j3a, j1b, j2b, j3b, 0, 0 )
__m256i const comp_mask  = _mm256_set_epi32( 7, 3, 6, 5, 4, 2, 1, 0 );

  __m256 vx0, vx1, vy0, vy1,  vz0, vz1, vq;

  __m256 vqnx, vqny, vqnz;
  __m256 vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP];

  __m256 vwl1[NP], vwl2[NP], vwl3[NP];
  __m256 vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];

  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

  int const Dy = 3 * size[0];

  int const Dz = Dy * size[1];

  __m256i const vDy = _mm256_set1_epi32( Dy );

  __m256i const vDz = _mm256_set1_epi32( Dz );

  float* const pj = current + ( offset[0] - OFFSET  ) * 3 +

  ( offset[1] - OFFSET  ) * Dy +

  ( offset[2] - OFFSET  ) * Dz;

  __m256 const vnorm1 = _mm256_set1_ps(norm[0]);
  __m256 const vnorm2 = _mm256_set1_ps(norm[1]);
  __m256 const vnorm3 = _mm256_set1_ps(norm[2]);

  int np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np % VEC_WIDTH ) {
    for( int i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
      int k = np + i;

      part -> x0[k] = part -> x1[k] = 0.;
      part -> y0[k] = part -> y1[k] = 0.;
      part -> z0[k] = part -> z1[k] = 0.;
      part -> q[k] = 0.;
      part -> ix[k] = part -> iy[k] = part -> iz[k] = 1;
    }
  }

  for( int i = 0; i < np; i+=VEC_WIDTH ) {

    __m256i vix, viy, viz;

    // load VEC_WIDTH particles
    LOAD8P3D( part, i, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz )

    // Store cell index
    vix = _mm256_add_epi32( vix, _mm256_add_epi32( vix, vix ) );
    viy = _mm256_mullo_epi32( viy, vDy );

    viz = _mm256_mullo_epi32( viz, vDz );

    _mm256_store_si256( (__m256i *) idx, _mm256_add_epi32( _mm256_add_epi32( vix, viy ), viz ) );

    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );

    vqnx = _mm256_mul_ps( vq, vnorm1);
    vqny = _mm256_mul_ps( vq, vnorm2);
    vqnz = _mm256_mul_ps( vq, vnorm3);

    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 ); vwl1[NP-1] = _mm256_setzero_ps();
    WL( vqny, vy0, vy1, vwl2 ); vwl2[NP-1] = _mm256_setzero_ps();
    WL( vqnz, vz0, vz1, vwl3 ); vwl3[NP-1] = _mm256_setzero_ps();

    // get perpendicular weights
    for( int k2 = 0; k2 < NP; k2++ ) {
      for( int k1 = 0; k1 < NP; k1++ ) {

        const __m256 c1_2 = _mm256_set1_ps(0.5f);
        __m256 tmp1, tmp2;

        // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] +

        //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )

        tmp1 = _mm256_mul_ps( vs0y[k1], vs0z[k2] );
        tmp2 = _mm256_mul_ps( vs0y[k1], vs1z[k2] );
        tmp1 = _mm256_fmadd_ps( vs1y[k1], vs1z[k2], tmp1 );
        tmp2 = _mm256_fmadd_ps( vs1y[k1], vs0z[k2], tmp2 );

        vwp1[k2][k1] = _mm256_fmadd_ps( c1_2, tmp2, tmp1 );

        // wp2[k2][k1] = s0x[k1]*s0z[k2] + s1x[k1]*s1z[k2] +

        //               0.5*( s0x[k1]*s1z[k2] + s1x[k1]*s0z[k2] )

        tmp1 = _mm256_mul_ps( vs0x[k1], vs0z[k2] );
        tmp2 = _mm256_mul_ps( vs0x[k1], vs1z[k2] );
        tmp1 = _mm256_fmadd_ps( vs1x[k1], vs1z[k2], tmp1 );
        tmp2 = _mm256_fmadd_ps( vs1x[k1], vs0z[k2], tmp2 );

        vwp2[k2][k1] = _mm256_fmadd_ps( c1_2, tmp2, tmp1 );

        // wp3[k2][k1] = s0x[k1]*s0y[k2] + s1x[k1]*s1y[k2] +

        //               0.5*( s0x[k1]*s1y[k2] + s1x[k1]*s0y[k2] )

        tmp1 = _mm256_mul_ps( vs0x[k1], vs0y[k2] );
        tmp2 = _mm256_mul_ps( vs0x[k1], vs1y[k2] );
        tmp1 = _mm256_fmadd_ps( vs1x[k1], vs1y[k2], tmp1 );
        tmp2 = _mm256_fmadd_ps( vs1x[k1], vs0y[k2], tmp2 );

        vwp3[k2][k1] = _mm256_fmadd_ps( c1_2, tmp2, tmp1 );

      }
    }

    __m256 a[NP][NP][NP], b[NP][NP][NP], c[NP][NP][NP], d[NP][NP][NP];
    for( int k3 = 0; k3 < NP; k3++ ) {
      for( int k2 = 0; k2 < NP; k2++ ) {
        for ( int k1 = 0; k1 < NP; k1++ ) {
          __m256 j1, j2, j3;

          // Calculate current components
          j1 = _mm256_mul_ps( vwl1[k1] , vwp1[k3][k2] );
          j2 = _mm256_mul_ps( vwl2[k2] , vwp2[k3][k1] );
          j3 = _mm256_mul_ps( vwl3[k3] , vwp3[k2][k1] );

          // Do a 8x4 transpose, setting the last term to 0

          __m256 const zero = _mm256_setzero_ps();

          __m256 t0 = _mm256_shuffle_ps( j1, j2, 0x44 );
          __m256 t2 = _mm256_shuffle_ps( j1, j2, 0xEE );
          __m256 t1 = _mm256_shuffle_ps( j3, zero, 0x44 );
          __m256 t3 = _mm256_shuffle_ps( j3, zero, 0xEE );

          a[k3][k2][k1] = _mm256_shuffle_ps( t0, t1, 0x88 );
          b[k3][k2][k1] = _mm256_shuffle_ps( t0, t1, 0xDD );
          c[k3][k2][k1] = _mm256_shuffle_ps( t2, t3, 0x88 );
          d[k3][k2][k1] = _mm256_shuffle_ps( t2, t3, 0xDD );
        }
      }
    }

    ACC_CURRENT3D( 0, a, 0 )
    ACC_CURRENT3D( 1, b, 0 )
    ACC_CURRENT3D( 2, c, 0 )
    ACC_CURRENT3D( 3, d, 0 )
    ACC_CURRENT3D( 4, a, 1 )
    ACC_CURRENT3D( 5, b, 1 )
    ACC_CURRENT3D( 6, c, 1 )
    ACC_CURRENT3D( 7, d, 1 )

  }

}

#undef ACC_CURRENT3D

#endif

#undef DEP_CURRENT_1D
#undef DEP_CURRENT_2D
#undef DEP_CURRENT_3D
#undef SPLINE
#undef WL

#undef OFFSET
#undef ORDER
#undef NP

#endif

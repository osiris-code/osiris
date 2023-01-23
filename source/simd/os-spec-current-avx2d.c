/*****************************************************************************************

Charge conserving current deposition, Intel AVX2 optimized version

*****************************************************************************************/

/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-avx2d.h"
#include "splines-avx2.h"

/*****************************************************************************************
vwl_s1
*****************************************************************************************/

inline void vwl_s1( __m256d const vqn, __m256d const vx0, __m256d const vx1, __m256d vwl[] )
{
  vwl[0] = _mm256_mul_pd(vqn, _mm256_sub_pd(vx1, vx0));
}


/*****************************************************************************************
vwl_s2
*****************************************************************************************/

inline void vwl_s2( __m256d const vqn, __m256d const vx0, __m256d const vx1, __m256d vwl[] )
{
  const __m256d c1_2 = _mm256_set1_pd( 0.5 );

  __m256d d    = _mm256_sub_pd( vx1, vx0 );
  __m256d s1_2 = _mm256_sub_pd( c1_2, vx0 );
  __m256d p1_2 = _mm256_add_pd( c1_2, vx0 );

  __m256d n    = _mm256_mul_pd( vqn, d );

  vwl[0] = _mm256_mul_pd( n, _mm256_fnmadd_pd( c1_2, d, s1_2 ));
  vwl[1] = _mm256_mul_pd( n, _mm256_fmadd_pd(  c1_2, d, p1_2 ));
}


/*****************************************************************************************
vwl_s3
*****************************************************************************************/

inline void vwl_s3( __m256d const vqn, __m256d const vx0, __m256d const vx1, __m256d vwl[] )
{
  __m256d const c1_2 = _mm256_set1_pd( 0.5 );
  __m256d const c3_4 = _mm256_set1_pd( 0.75 );
  __m256d const c1_3 = _mm256_set1_pd( 1.0/3.0 );

  __m256d   d = _mm256_sub_pd(  vx1, vx0 );
  __m256d   s = _mm256_sub_pd( c1_2, vx0 );
  __m256d   p = _mm256_add_pd( c1_2, vx0 );
  __m256d d_3 = _mm256_mul_pd( c1_3,   d );

  vwl[0] = _mm256_mul_pd( c1_2, _mm256_fmsub_pd( s, s, _mm256_mul_pd( d, _mm256_sub_pd( s, d_3 )) ) );
  vwl[1] = _mm256_fnmadd_pd( d, _mm256_add_pd( vx0, d_3 ), _mm256_fnmadd_pd( vx0, vx0, c3_4 ));
  vwl[2] = _mm256_mul_pd( c1_2, _mm256_fmadd_pd( p, p, _mm256_mul_pd( d, _mm256_add_pd( p, d_3 )) ) );

  __m256d n = _mm256_mul_pd( vqn, d );
  vwl[0] = _mm256_mul_pd( n, vwl[0] );
  vwl[1] = _mm256_mul_pd( n, vwl[1] );
  vwl[2] = _mm256_mul_pd( n, vwl[2] );
}

/*****************************************************************************************
vwl_s4
*****************************************************************************************/

inline void vwl_s4( __m256d const vqn, __m256d const vx0, __m256d const vx1, __m256d vwl[] )
{
  __m256d const c1_2 = _mm256_set1_pd( 0.5 );
  __m256d const c1_4 = _mm256_set1_pd( 0.25 );
  __m256d const c1_6 = _mm256_set1_pd( 1.0/6.0 );
  __m256d const c3_2 = _mm256_set1_pd( 1.5 );
  __m256d const c3_4 = _mm256_set1_pd( 0.75 );
  __m256d const c1_3 = _mm256_set1_pd( 1.0/3.0 );
  __m256d const c2_3 = _mm256_set1_pd( 2.0/3.0 );

  __m256d d = _mm256_sub_pd( vx1, vx0 );
  __m256d s = _mm256_sub_pd( c1_2, vx0 );
  __m256d p = _mm256_add_pd( c1_2, vx0 );

  __m256d t = _mm256_sub_pd( vx0, c1_6 );
  __m256d u = _mm256_add_pd( vx0, c1_6 );

  __m256d s2 = _mm256_mul_pd( s, s );
  __m256d s3 = _mm256_mul_pd( s2, s );

  __m256d p2 = _mm256_mul_pd( p, p );
  __m256d p3 = _mm256_mul_pd( p2, p );

  __m256d d_2 = _mm256_mul_pd( d, c1_2 );
  __m256d d_4 = _mm256_mul_pd( d, c1_4 );

  vwl[0] = _mm256_mul_pd( c1_6, _mm256_fmadd_pd( d, _mm256_fmsub_pd( d, _mm256_sub_pd( s, d_4 ), _mm256_mul_pd( c3_2, s2 )  ), s3 ) );
  vwl[1] = _mm256_fmadd_pd( d, _mm256_fmadd_pd(d_2, _mm256_add_pd(t, d_4), _mm256_fmsub_pd(c3_4, _mm256_mul_pd(t, t), c1_3)),
                          _mm256_fmadd_pd(c1_2, p3, _mm256_sub_pd(c2_3, p2)));
  vwl[2] = _mm256_fnmadd_pd(d, _mm256_fmadd_pd(d_2, _mm256_add_pd(u, d_4), _mm256_fmsub_pd(c3_4, _mm256_mul_pd(u, u), c1_3)),
                          _mm256_fmadd_pd(c1_2, s3, _mm256_sub_pd(c2_3, s2)));
  vwl[3] = _mm256_mul_pd( c1_6, _mm256_fmadd_pd( d, _mm256_fmadd_pd( d, _mm256_add_pd( p, d_4 ), _mm256_mul_pd( c3_2, p2 )  ), p3 ) );

  __m256d n = _mm256_mul_pd( vqn, d );
  vwl[0] = _mm256_mul_pd( n, vwl[0] );
  vwl[1] = _mm256_mul_pd( n, vwl[1] );
  vwl[2] = _mm256_mul_pd( n, vwl[2] );
  vwl[3] = _mm256_mul_pd( n, vwl[3] );
}

/****************************************************************************************

  Generate specific current deposition functions for 2D and 3D, 1st to 4th order

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

#define SPLINE  ONAME( vsplined, ORDER )
#define WL      ONAME( vwl, ORDER )

/*
  Reference implementation, serial deposition
*/

#if 1

void DEP_CURRENT_1D
(double * const current, int const * const size, int const * const offset,
 double * const norm, t_split_buf1D * const part)
{

  typedef struct Current { double j1, j2, j3; } t_current;
  t_current *p0;

  dvec j2[NP],j3[NP];
  dvec vwl1[ORDER];

  __m256d vx0, vx1, vq, vvy, vvz;

  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

  __m256d vqnx, vqvy, vqvz;
  __m256d vs0x[NP], vs1x[NP];

  unsigned int i, np;

  __m256d const c1_2  = _mm256_set1_pd( 0.5 );
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 );
  __m256d const vnorm1 = _mm256_set1_pd(norm[0]);
  np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       unsigned k = np + i;

       part->x0[k] = part->x1[k] = 0.0;
       part-> q[k] = part->vy[k] = part->vz[k] = 0.0;

       part->ix[k] = 1;
    }
  }

  for( i = 0; i < np; i += VEC_WIDTH ) {

    unsigned k, k1;
    __m128i vix;

    // load 4 particles
    LOAD4P1D( part, i, vx0, vx1, vq, vvy, vvz, vix )

    // Store cell index
    _mm_store_si128( (__m128i *) idx, vix );

    // Get splines
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );

    vqnx = _mm256_mul_pd( vq, vnorm1 );
    vqvy = _mm256_mul_pd( vq, vvy );
    vqvz = _mm256_mul_pd( vq, vvz );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256d *) vwl1 );

    // get j2,j3 current
    for ( k1 = 0; k1 < NP; k1++ ) {
      __m256d vwp1 = _mm256_mul_pd( c1_2, _mm256_add_pd(vs0x[k1], vs1x[k1]) );
      j2[k1].v4  = _mm256_mul_pd( vqvy, vwp1 );
      j3[k1].v4  = _mm256_mul_pd( vqvz, vwp1 );
    }

    // Loop by particle on the outside loop
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

      p0 = pj + idx[k];

      // accumulate j1
      for ( k1 = 0; k1 < ORDER; k1++ ) {
        p0[k1].j1 += vwl1[k1].v[k] ;
      }

      // accumulate j2, j3
      for ( k1 = 0; k1 < NP; k1++ ) {
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

#if 0

#warning Using mk. I 1D - algorithm

#if ( ORDER == 1 )


#define ACC_CURRENT1D( k, v ) { \
  double *p0 = pj + idx[k]; \
  __m256d j0, j1; \
  j0 = _mm256_loadu_pd( p0 ); \
  j1 = _mm256_loadu_pd( p0 + 3 ); \
  j0 = _mm256_add_pd( j0, v[0] ); \
  j1 = _mm256_add_pd( j1, v[1] ); \
  _mm256_storeu_pd( p0,     j0); \
  _mm256_storeu_pd( p0 + 3, j1); \
}


#elif ( ORDER == 2 )

#define ACC_CURRENT1D( k, v ) { \
  double *p0 = pj + idx[k]; \
  __m256d j0, j1, j2; \
  j0 = _mm256_loadu_pd( p0 ); \
  j1 = _mm256_loadu_pd( p0 + 3 ); \
  j2 = _mm256_loadu_pd( p0 + 6 ); \
  j0 = _mm256_add_pd( j0, v[0] ); \
  j1 = _mm256_add_pd( j1, v[1] ); \
  j2 = _mm256_add_pd( j2, v[2] ); \
  _mm256_storeu_pd( p0,     j0); \
  _mm256_storeu_pd( p0 + 3, j1); \
  _mm256_storeu_pd( p0 + 6, j2); \
}

#elif ( ORDER == 3 )

#define ACC_CURRENT1D( k, v ) { \
  double *p0 = pj + idx[k]; \
  __m256d j0, j1, j2, j3; \
  j0 = _mm256_loadu_pd( p0 ); \
  j1 = _mm256_loadu_pd( p0 + 3 ); \
  j2 = _mm256_loadu_pd( p0 + 6 ); \
  j3 = _mm256_loadu_pd( p0 + 9 ); \
  j0 = _mm256_add_pd( j0, v[0] ); \
  j1 = _mm256_add_pd( j1, v[1] ); \
  j2 = _mm256_add_pd( j2, v[2] ); \
  j3 = _mm256_add_pd( j3, v[3] ); \
  _mm256_storeu_pd( p0,     j0); \
  _mm256_storeu_pd( p0 + 3, j1); \
  _mm256_storeu_pd( p0 + 6, j2); \
  _mm256_storeu_pd( p0 + 9, j3); \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT1D( k, v ) {  \
  double *p0 = pj + idx[k]; \
  __m256d j0, j1, j2, j3, j4; \
  j0 = _mm256_loadu_pd( p0 ); \
  j1 = _mm256_loadu_pd( p0 + 3 ); \
  j2 = _mm256_loadu_pd( p0 + 6 ); \
  j3 = _mm256_loadu_pd( p0 + 9 ); \
  j4 = _mm256_loadu_pd( p0 + 12 ); \
  j0 = _mm256_add_pd( j0, v[0] ); \
  j1 = _mm256_add_pd( j1, v[1] ); \
  j2 = _mm256_add_pd( j2, v[2] ); \
  j3 = _mm256_add_pd( j3, v[3] ); \
  j4 = _mm256_add_pd( j4, v[4] ); \
  _mm256_storeu_pd( p0,     j0); \
  _mm256_storeu_pd( p0 + 3, j1); \
  _mm256_storeu_pd( p0 + 6, j2); \
  _mm256_storeu_pd( p0 + 9, j3); \
  _mm256_storeu_pd( p0 + 12, j4); \
}

#else
#error Unsupported interpolation order
#endif


void DEP_CURRENT_1D
(double * const current, int const * const size, int const * const offset,
                       double * const norm, t_split_buf1D * const part)
{

  __m256d vx0, vx1, vq, vvy, vvz;

  __m256d vqnx, vqvy, vqvz;
  __m256d vs0x[NP], vs1x[NP];

  unsigned i, np;

  __m256d const c1_2  = _mm256_set1_pd( 0.5 );

  double* const pj =  current + ( offset[0] - OFFSET ) * 3 ;

  __m256d const vnorm1 = _mm256_set1_pd(norm[0]);

  np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       unsigned k = np + i;

       part->x0[k] = part->x1[k] = 0.0;
       part-> q[k] = part->vy[k] = part->vz[k] = 0.0;

       part->ix[k] = 1;
    }
  }

  for( i = 0; i < np; i += VEC_WIDTH ) {

    unsigned k1;

    __m256d vwl1[NP];

    __m128i vix;
    DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );
    __m256d a[NP], b[NP], c[NP], d[NP];

    // load VEC_WIDTH particles
    LOAD4P1D( part, i, vx0, vx1, vq, vvy, vvz, vix )

    // Store cell index
    vix = _mm_add_epi32( vix, _mm_add_epi32( vix, vix ) );
    _mm_store_si128( (__m128i *) idx, vix );

    // Get splines
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );

    vqnx = _mm256_mul_pd( vq, vnorm1 );
    vqvy = _mm256_mul_pd( vq, vvy );
    vqvz = _mm256_mul_pd( vq, vvz );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 ); vwl1[NP-1] = _mm256_setzero_pd();

    // get currents and transpose vectors
    for ( k1 = 0; k1 < NP; k1++ ) {
      __m256d vwp1, j1, j2, j3;

      vwp1 = _mm256_mul_pd( c1_2, _mm256_add_pd(vs0x[k1], vs1x[k1]) );

      j1  = vwl1[k1];
      j2  = _mm256_mul_pd( vqvy, vwp1 );
      j3  = _mm256_mul_pd( vqvz, vwp1 );

      // Do a 4x4 transpose, setting the last term to 0

      __m256d const zero = _mm256_setzero_pd();

      __m256d t0 = _mm256_shuffle_pd( j1, j2, 0x0 );
      __m256d t2 = _mm256_shuffle_pd( j1, j2, 0xF );
      __m256d t1 = _mm256_shuffle_pd( j3, zero, 0x0 );
      __m256d t3 = _mm256_shuffle_pd( j3, zero, 0xEE );

      a[k1] = _mm256_permute2f128_pd( t0, t2, 0x20 );
      b[k1] = _mm256_permute2f128_pd( t0, t2, 0x31 );
      c[k1] = _mm256_permute2f128_pd( t1, t3, 0x20 );
      d[k1] = _mm256_permute2f128_pd( t1, t3, 0x31 );

    }

    ACC_CURRENT1D( 0, a )
    ACC_CURRENT1D( 1, b )
    ACC_CURRENT1D( 2, c )
    ACC_CURRENT1D( 3, d )

  }

}

#undef ACC_CURRENT1D

#endif




#if 1

/*
  Reference implementation, serial deposition
*/

#warning Using Reference 2D algorithm

void DEP_CURRENT_2D
(double * const current, int const * const size, int const * const offset,
                       double * const norm, t_split_buf2D * const part)
{

  typedef struct Current { double j1, j2, j3; } t_current;
  t_current *pjpart, *p0;


  dvec j3[NP][NP];
  dvec vwp1[NP], vwp2[NP];
  dvec vwl1[ORDER], vwl2[ORDER];

  __m256d vx0, vx1, vy0, vy1, vq, vvz;

  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

  __m256d vqnx, vqny, vqvz;
  __m256d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP];

  unsigned np;
  unsigned i, k, k1, k2;

  __m256d const c1_3 = _mm256_set1_pd( 1.0/3.0 );
  __m256d const c1_2  = _mm256_set1_pd( 0.5 );

  unsigned const Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 +
                                                 ( offset[1] - OFFSET ) * 3 * Dy );

  // Setting this to 1/4 removes 1 multiplication from wl2
  __m256d const vnorm1 = _mm256_set1_pd(norm[0]);
  __m256d const vnorm2 = _mm256_set1_pd(norm[1]);

  np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH != 0 ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       unsigned k = np + i;

       part->x0[k] = part->x1[k] = 0.0;
       part->y0[k] = part->y1[k] = 0.0;
       part-> q[k] = part->vz[k] = 0.0;

       part->ix[k] = part->iy[k] = 1;
    }
  }

  // Each virtual particle uses 8 doubles, so 4 vp = 32 doubles
  for( i = 0; i < np; i += VEC_WIDTH ) {

    __m128i vix, viy;

    // load VEC_WIDTH particles
    LOAD4P2D( part, i, vx0, vx1, vy0, vy1, vq, vvz, vix, viy )

    // Store cell index
    _mm_store_si128( (__m128i *) idx,
            _mm_add_epi32( vix, _mm_mullo_epi32( viy, _mm_set1_epi32( Dy ) ) ) );

    // Get splines
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );

    vqnx = _mm256_mul_pd( vq, vnorm1 );
    vqny = _mm256_mul_pd( vq, vnorm2 );
    vqvz = _mm256_mul_pd( vq, vvz);
    vqvz = _mm256_mul_pd( vqvz, c1_3 );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256d *) vwl1 );
    WL( vqny, vy0, vy1, (__m256d *) vwl2 );

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k].v4 = _mm256_add_pd(vs0y[k], vs1y[k]);
      vwp2[k].v4 = _mm256_add_pd(vs0x[k], vs1x[k]);
    }

    // get j3 current
    for( k2 = 0; k2 < NP; k2++ ) {
      for ( k1 = 0; k1 < NP; k1++ ) {
        __m256d s00, s10, tmp1, tmp2;

    		s00  = _mm256_mul_pd( vs0x[k1], vs0y[k2] );
    		tmp1 = _mm256_fmadd_pd( vs1y[k2], vs1x[k1], s00 );  // tmp1 = s0x*s0y + s1x*s1y

    		s10  = _mm256_mul_pd( vs0x[k1], vs1y[k2] );
    		tmp2 = _mm256_fmadd_pd( vs0y[k2], vs1x[k1], s10 );  // tmp2 = s0x*s1y + s1x*s0y

    		tmp1 = _mm256_fmadd_pd( tmp2, c1_2, tmp1 );       // tmp1 = tmp1 + 0.5*tmp2

    		j3[k1][k2].v4 = _mm256_mul_pd( vqvz, tmp1 );          // j3 = vqvz * tmp1
      }
    }

    // Loop by particle on the outside loop
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

      pjpart = pj + idx[k];

      // accumulate j1
      for( k2 = 0; k2 < NP; k2++ ) {
        p0 = pjpart + k2*Dy;
        #pragma unroll(ORDER)
        for ( k1 = 0; k1 < ORDER; k1++ ) {
          p0[k1].j1 += vwl1[k1].v[k] * vwp1[k2].v[k];
        }
      }

      // accumulate j2 - making k2 the outside loop gives marginal perf. gain
      for( k2 = 0; k2 < ORDER; k2++ ) {
        p0 = pjpart + k2*Dy;
        #pragma unroll(NP)
        for ( k1 = 0; k1 < NP; k1++ ) {
          p0[k1].j2 += vwl2[k2].v[k] * vwp2[k1].v[k];
        }
      }

      // accumulate j3
      for( k2 = 0; k2 < NP; k2++ ) {
        p0 = pjpart + k2*Dy;
        #pragma unroll(NP)
        for ( k1 = 0; k1 < NP; k1++ ) {
          p0[k1].j3 += (j3[k1][k2]).v[k];
        }
      }

    }

  }

}
#endif

#if 0
/*
mk. I - Transpose currents, save 1 particle at a time

process each line sequentially
*/


#warning Using mk. I 2D algorithm

#if ( ORDER == 1 )

#define ACC_CURRENT2D( k, v ) { \
  for( int j = 0; j < 2; j++ ){ \
    double * const p0 = pj + idx[k] + j*Dy; \
    __m256d j0, j1; \
    j0 = _mm256_loadu_pd( p0 ); \
    j1 = _mm256_loadu_pd( p0 + 3 ); \
    j0 = _mm256_add_pd( j0, v[j][0] ); \
    j1 = _mm256_add_pd( j1, v[j][1] ); \
    _mm256_storeu_pd( p0,     j0); \
    _mm256_storeu_pd( p0 + 3, j1); \
  } \
}

#elif ( ORDER == 2 )

#define ACC_CURRENT2D( k, v ) { \
  for( int j = 0; j < 3; j++ ){ \
    double * const p0 = pj + idx[k] + j*Dy; \
    __m256d j0, j1, j2; \
    j0 = _mm256_loadu_pd( p0 ); \
    j1 = _mm256_loadu_pd( p0 + 3 ); \
    j2 = _mm256_loadu_pd( p0 + 6 ); \
    j0 = _mm256_add_pd( j0, v[j][0] ); \
    j1 = _mm256_add_pd( j1, v[j][1] ); \
    j2 = _mm256_add_pd( j2, v[j][2] ); \
    _mm256_storeu_pd( p0,     j0); \
    _mm256_storeu_pd( p0 + 3, j1); \
    _mm256_storeu_pd( p0 + 6, j2); \
  } \
}


#elif ( ORDER == 3 )

#define ACC_CURRENT2D( k, v ) { \
  for( int j = 0; j < 4; j++ ){ \
    double * const p0 = pj + idx[k] + j*Dy; \
    __m256d j0, j1, j2, j3; \
    j0 = _mm256_loadu_pd( p0 ); \
    j1 = _mm256_loadu_pd( p0 + 3 ); \
    j2 = _mm256_loadu_pd( p0 + 6 ); \
    j3 = _mm256_loadu_pd( p0 + 9 ); \
    j0 = _mm256_add_pd( j0, v[j][0] ); \
    j1 = _mm256_add_pd( j1, v[j][1] ); \
    j2 = _mm256_add_pd( j2, v[j][2] ); \
    j3 = _mm256_add_pd( j3, v[j][3] ); \
    _mm256_storeu_pd( p0,     j0); \
    _mm256_storeu_pd( p0 + 3, j1); \
    _mm256_storeu_pd( p0 + 6, j2); \
    _mm256_storeu_pd( p0 + 9, j3); \
  } \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT2D( k, v ) { \
  for( int j = 0; j < 5; j++ ){ \
    double * const p0 = pj + idx[k] + j*Dy; \
    __m256d j0, j1, j2, j3, j4; \
    j0 = _mm256_loadu_pd( p0 ); \
    j1 = _mm256_loadu_pd( p0 + 3 ); \
    j2 = _mm256_loadu_pd( p0 + 6 ); \
    j3 = _mm256_loadu_pd( p0 + 9 ); \
    j4 = _mm256_loadu_pd( p0 + 12 ); \
    j0 = _mm256_add_pd( j0, v[j][0] ); \
    j1 = _mm256_add_pd( j1, v[j][1] ); \
    j2 = _mm256_add_pd( j2, v[j][2] ); \
    j3 = _mm256_add_pd( j3, v[j][3] ); \
    j4 = _mm256_add_pd( j4, v[j][4] ); \
    _mm256_storeu_pd( p0,     j0); \
    _mm256_storeu_pd( p0 + 3, j1); \
    _mm256_storeu_pd( p0 + 6, j2); \
    _mm256_storeu_pd( p0 + 9, j3); \
    _mm256_storeu_pd( p0 + 12, j4); \
  } \
}

#else
#error Unsupported interpolation order
#endif


void DEP_CURRENT_2D
(double * restrict const current, int const * restrict const size, int const * restrict const offset,
                       double * restrict const norm, t_split_buf2D * restrict const part)
{

  __m256d vwp1[NP], vwp2[NP];
  __m256d vwl1[NP], vwl2[NP];

  __m256d vx0, vx1, vy0, vy1, vq, vvz;

  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

  __m256d vqnx, vqny, vqvz;
  __m256d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP];

  __m256d const c1_3 = _mm256_set1_pd( 1.0f/3.0f );
  __m256d const c1_2  = _mm256_set1_pd( 0.5f );

  // Current array has 3 field components
  double* const pj = current + 3 * ( ( offset[0] - OFFSET ) +
                                    ( offset[1] - OFFSET ) * size[0] );

  // Normalization for current
  __m256d const vnorm1 = _mm256_set1_pd(norm[0]);
  __m256d const vnorm2 = _mm256_set1_pd(norm[1]);

  unsigned const Dy = 3*size[0];
  __m128i const vDy = _mm_set1_epi32( Dy );

  unsigned i;
  unsigned np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH != 0 ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       unsigned k = np + i;

       part->x0[k] = part->x1[k] = 0.0;
       part->y0[k] = part->y1[k] = 0.0;
       part-> q[k] = part->vz[k] = 0.0;

       part->ix[k] = part->iy[k] = 1;
    }
  }

  for( i = 0; i < np; i += VEC_WIDTH ) {

    unsigned k, k1, k2;
    __m128i vix, viy;

    // load VEC_WIDTH particles
    LOAD4P2D( part, i, vx0, vx1, vy0, vy1, vq, vvz, vix, viy )

    // Store cell index
    // idx = 3 * vix + 3 * size[0] * viy
    vix = _mm_add_epi32( vix, _mm_add_epi32( vix, vix ) );
    viy = _mm_mullo_epi32( viy, vDy );
    _mm_store_si128( (__m128i *) idx, _mm_add_epi32( vix, viy ) );

    // Get splines
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );

    vqnx = _mm256_mul_pd( vq, vnorm1 );
    vqny = _mm256_mul_pd( vq, vnorm2 );
    vqvz = _mm256_mul_pd( vq, vvz);
    vqvz = _mm256_mul_pd( vqvz, c1_3 );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 ); vwl1[NP-1] = _mm256_setzero_pd();
    WL( vqny, vy0, vy1, vwl2 ); vwl2[NP-1] = _mm256_setzero_pd();

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k] = _mm256_add_pd(vs0y[k], vs1y[k]);
      vwp2[k] = _mm256_add_pd(vs0x[k], vs1x[k]);
    }


    __m256d a[NP][NP], b[NP][NP], c[NP][NP], d[NP][NP];

    for( k2 = 0; k2 < NP; k2++ ) {
      for ( k1 = 0; k1 < NP; k1++ ) {
        __m256d j1, j2, j3;

        j1 = _mm256_mul_pd( vwl1[k1] , vwp1[k2] );
        j2 = _mm256_mul_pd( vwl2[k2] , vwp2[k1] );

        __m256d s00, s10, tmp1, tmp2;

        s00  = _mm256_mul_pd( vs0x[k1], vs0y[k2] );
        tmp1 = _mm256_fmadd_pd( vs1y[k2], vs1x[k1], s00 );

        s10  = _mm256_mul_pd( vs0x[k1], vs1y[k2] );
        tmp2 = _mm256_fmadd_pd( vs0y[k2], vs1x[k1], s10 );

        tmp1 = _mm256_fmadd_pd( tmp2, c1_2, tmp1 );

        j3 = _mm256_mul_pd( vqvz, tmp1 );

        // Do a 4x4 transpose, setting the last term to 0
        __m256d const zero = _mm256_setzero_pd();

        __m256d t0 = _mm256_shuffle_pd( j1, j2, 0x0 );
        __m256d t2 = _mm256_shuffle_pd( j1, j2, 0xF );
        __m256d t1 = _mm256_shuffle_pd( j3, zero, 0x0 );
        __m256d t3 = _mm256_shuffle_pd( j3, zero, 0xEE );

        a[k2][k1] = _mm256_permute2f128_pd( t0, t2, 0x20 );
        b[k2][k1] = _mm256_permute2f128_pd( t0, t2, 0x31 );
        c[k2][k1] = _mm256_permute2f128_pd( t1, t3, 0x20 );
        d[k2][k1] = _mm256_permute2f128_pd( t1, t3, 0x31 );
      }
    }


    ACC_CURRENT2D( 0, a )
    ACC_CURRENT2D( 1, b )
    ACC_CURRENT2D( 2, c )
    ACC_CURRENT2D( 3, d )

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

#if 1

#warning Using reference 3D algorithm

void DEP_CURRENT_3D
(double * const  current, int const * const  size, int const * const  offset,
                       double * const norm, t_split_buf3D * const part)
{

  typedef struct Current { double j1, j2, j3; } t_current;

  __m256d vx0, vx1, vy0, vy1,  vz0, vz1, vq;

  __m256d vqnx, vqny, vqnz;
  __m256d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP];

  dvec vwl1[ORDER], vwl2[ORDER], vwl3[ORDER];
  dvec vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];

  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

  int const Dy = size[0];
  int const Dz = Dy * size[1];

  t_current*  pj = (t_current *) ( current + ( offset[0] - OFFSET  ) * 3 +
                                             ( offset[1] - OFFSET  ) * 3 * Dy +
                                             ( offset[2] - OFFSET  ) * 3 * Dz );

  __m256d const vnorm1 = _mm256_set1_pd(norm[0]);
  __m256d const vnorm2 = _mm256_set1_pd(norm[1]);
  __m256d const vnorm3 = _mm256_set1_pd(norm[2]);

  unsigned i;
  unsigned np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
      unsigned k = np + i;

      part -> x0[k] = part -> x1[k] = 0.;
      part -> y0[k] = part -> y1[k] = 0.;
      part -> z0[k] = part -> z1[k] = 0.;
      part -> q[k] = 0.;
      part -> ix[k] = part -> iy[k] = part -> iz[k] = 1;
    }
  }

  for( i = 0; i < np; i+=VEC_WIDTH ) {

    unsigned k, k1, k2, k3;
    __m128i vix, viy, viz, vidx;

    // load VEC_WIDTH particles
    LOAD4P3D( part, i, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz )

    // Store cell index
    vidx = _mm_add_epi32( vix,  _mm_mullo_epi32( viy, _mm_set1_epi32( Dy ) ) );
    vidx = _mm_add_epi32( vidx, _mm_mullo_epi32( viz, _mm_set1_epi32( Dz ) ) );

    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );

    vqnx = _mm256_mul_pd( vq, vnorm1 );
    vqny = _mm256_mul_pd( vq, vnorm2 );
    vqnz = _mm256_mul_pd( vq, vnorm3 );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m256d *) vwl1 );
    WL( vqny, vy0, vy1, (__m256d *) vwl2 );
    WL( vqnz, vz0, vz1, (__m256d *) vwl3 );

    // get perpendicular weights
    for( k2 = 0; k2 < NP; k2++ ) {
      for( k1 = 0; k1 < NP; k1++ ) {

        const __m256d c1_2 = _mm256_set1_pd(0.5);
        __m256d tmp1, tmp2;

        // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] +
        //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )

        tmp1 = _mm256_mul_pd( vs0y[k1], vs0z[k2] );
        tmp2 = _mm256_mul_pd( vs0y[k1], vs1z[k2] );
        tmp1 = _mm256_fmadd_pd( vs1y[k1], vs1z[k2], tmp1 );
        tmp2 = _mm256_fmadd_pd( vs1y[k1], vs0z[k2], tmp2 );

        vwp1[k2][k1].v4 = _mm256_fmadd_pd( c1_2, tmp2, tmp1 );

        // wp2[k2][k1] = s0x[k1]*s0z[k2] + s1x[k1]*s1z[k2] +
        //               0.5*( s0x[k1]*s1z[k2] + s1x[k1]*s0z[k2] )

        tmp1 = _mm256_mul_pd( vs0x[k1], vs0z[k2] );
        tmp2 = _mm256_mul_pd( vs0x[k1], vs1z[k2] );
        tmp1 = _mm256_fmadd_pd( vs1x[k1], vs1z[k2], tmp1 );
        tmp2 = _mm256_fmadd_pd( vs1x[k1], vs0z[k2], tmp2 );

        vwp2[k2][k1].v4 = _mm256_fmadd_pd( c1_2, tmp2, tmp1 );

        // wp3[k2][k1] = s0x[k1]*s0y[k2] + s1x[k1]*s1y[k2] +
        //               0.5*( s0x[k1]*s1y[k2] + s1x[k1]*s0y[k2] )

        tmp1 = _mm256_mul_pd( vs0x[k1], vs0y[k2] );
        tmp2 = _mm256_mul_pd( vs0x[k1], vs1y[k2] );
        tmp1 = _mm256_fmadd_pd( vs1x[k1], vs1y[k2], tmp1 );
        tmp2 = _mm256_fmadd_pd( vs1x[k1], vs0y[k2], tmp2 );

        vwp3[k2][k1].v4 = _mm256_fmadd_pd( c1_2, tmp2, tmp1 );

      }
    }

    // The following code will need to access the individual vidx values
    _mm_store_si128( (__m128i *) idx, vidx );

    // Loop by particle on the outside loop
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

      t_current * pjpart = pj + idx[k];

      // accumulate j1
      for( k3 = 0; k3 < NP; k3++ ) {
        for( k2 = 0; k2 < NP; k2++ ) {
          t_current *p0 = pjpart + k2*Dy + k3*Dz;
          #pragma unroll(ORDER)
          for ( k1 = 0; k1 < ORDER; k1++ ) {
            p0[k1].j1 += vwl1[k1].v[k] * vwp1[k3][k2].v[k];
          }
        }
      }

      // accumulate j2
      for( k3 = 0; k3 < NP; k3++ ) {
        for( k2 = 0; k2 < ORDER; k2++ ) {
          t_current *p0 = pjpart + k2*Dy + k3*Dz;
          #pragma unroll(NP)
          for ( k1 = 0; k1 < NP; k1++ ) {
            p0[k1].j2 += vwl2[k2].v[k] * vwp2[k3][k1].v[k];
          }
        }
      }

      // accumulate j3
      for( k3 = 0; k3 < ORDER; k3++ ) {
        for( k2 = 0; k2 < NP; k2++ ) {
          t_current *p0 = pjpart + k2*Dy + k3*Dz;
          #pragma unroll(NP)
          for ( k1 = 0; k1 < NP; k1++ ) {
            p0[k1].j3 += vwl3[k3].v[k] * vwp3[k2][k1].v[k];
          }
        }
      }

    }

  }
}

#endif




#if 0

/* mk. I

For each vector of particles loop over z coordinate and deposit one current plane at a
time using the same algorithm as 2D mk. VI

*/

#warning Using mk. I 3D algorithm

#if ( ORDER == 1 )

#define ACC_CURRENT3D( k, v ) { \
  for( int j = 0; j < 2; j++ ){ \
    double * const p0 = pj + k3*Dz + idx[k] + j*Dy; \
    __m256d j0, j1; \
    j0 = _mm256_loadu_pd( p0 ); \
    j1 = _mm256_loadu_pd( p0 + 3 ); \
    j0 = _mm256_add_pd( j0, v[j][0] ); \
    j1 = _mm256_add_pd( j1, v[j][1] ); \
    _mm256_storeu_pd( p0,     j0); \
    _mm256_storeu_pd( p0 + 3, j1); \
  } \
}


#elif ( ORDER == 2 )

#define ACC_CURRENT3D( k, v ) { \
  for( int j = 0; j < 3; j++ ){ \
    double * const p0 = pj + k3 * Dz + idx[k] + j*Dy; \
    __m256d j0, j1, j2; \
    j0 = _mm256_loadu_pd( p0 ); \
    j1 = _mm256_loadu_pd( p0 + 3 ); \
    j2 = _mm256_loadu_pd( p0 + 6 ); \
    j0 = _mm256_add_pd( j0, v[j][0] ); \
    j1 = _mm256_add_pd( j1, v[j][1] ); \
    j2 = _mm256_add_pd( j2, v[j][2] ); \
    _mm256_storeu_pd( p0,     j0); \
    _mm256_storeu_pd( p0 + 3, j1); \
    _mm256_storeu_pd( p0 + 6, j2); \
  } \
}

#elif ( ORDER == 3 )

#define ACC_CURRENT3D( k, v ) { \
  for( int j = 0; j < 4; j++ ){ \
    double * const p0 = pj + k3 * Dz + idx[k] + j*Dy; \
    __m256d j0, j1, j2, j3; \
    j0 = _mm256_loadu_pd( p0 ); \
    j1 = _mm256_loadu_pd( p0 + 3 ); \
    j2 = _mm256_loadu_pd( p0 + 6 ); \
    j3 = _mm256_loadu_pd( p0 + 9 ); \
    j0 = _mm256_add_pd( j0, v[j][0] ); \
    j1 = _mm256_add_pd( j1, v[j][1] ); \
    j2 = _mm256_add_pd( j2, v[j][2] ); \
    j3 = _mm256_add_pd( j3, v[j][3] ); \
    _mm256_storeu_pd( p0,     j0); \
    _mm256_storeu_pd( p0 + 3, j1); \
    _mm256_storeu_pd( p0 + 6, j2); \
    _mm256_storeu_pd( p0 + 9, j3); \
  } \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT3D( k, v ) { \
  for( int j = 0; j < 5; j++ ){ \
    double * const p0 = pj + k3 * Dz + idx[k] + j*Dy; \
    __m256d j0, j1, j2, j3, j4; \
    j0 = _mm256_loadu_pd( p0 ); \
    j1 = _mm256_loadu_pd( p0 + 3 ); \
    j2 = _mm256_loadu_pd( p0 + 6 ); \
    j3 = _mm256_loadu_pd( p0 + 9 ); \
    j4 = _mm256_loadu_pd( p0 + 12 ); \
    j0 = _mm256_add_pd( j0, v[j][0] ); \
    j1 = _mm256_add_pd( j1, v[j][1] ); \
    j2 = _mm256_add_pd( j2, v[j][2] ); \
    j3 = _mm256_add_pd( j3, v[j][3] ); \
    j4 = _mm256_add_pd( j4, v[j][4] ); \
    _mm256_storeu_pd( p0,     j0); \
    _mm256_storeu_pd( p0 + 3, j1); \
    _mm256_storeu_pd( p0 + 6, j2); \
    _mm256_storeu_pd( p0 + 9, j3); \
    _mm256_storeu_pd( p0 + 12, j4); \
  } \
}


#else
#error Unsupported interpolation order
#endif


void DEP_CURRENT_3D
(double * const  current, int const * const  size, int const * const  offset,
                       double * const norm, t_split_buf3D * const part)
{

  __m256d vx0, vx1, vy0, vy1,  vz0, vz1, vq;

  __m256d vqnx, vqny, vqnz;
  __m256d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP];

  __m256d vwl1[NP], vwl2[NP], vwl3[NP];
  __m256d vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];

  DECLARE_ALIGNED_32( int idx[VEC_WIDTH] );

  unsigned const Dy = 3 * size[0];
  unsigned const Dz = Dy * size[1];

  __m128i const vDy = _mm_set1_epi32( Dy );
  __m128i const vDz = _mm_set1_epi32( Dz );

  double* const pj = current + ( offset[0] - OFFSET  ) * 3 +
  ( offset[1] - OFFSET  ) * Dy +
  ( offset[2] - OFFSET  ) * Dz;

  __m256d const vnorm1 = _mm256_set1_pd(norm[0]);
  __m256d const vnorm2 = _mm256_set1_pd(norm[1]);
  __m256d const vnorm3 = _mm256_set1_pd(norm[2]);



  unsigned i;
  unsigned np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end
  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
      unsigned k = np + i;

      part -> x0[k] = part -> x1[k] = 0.;
      part -> y0[k] = part -> y1[k] = 0.;
      part -> z0[k] = part -> z1[k] = 0.;
      part -> q[k] = 0.;
      part -> ix[k] = part -> iy[k] = part -> iz[k] = 1;
    }
  }

  for( i = 0; i < np; i+=VEC_WIDTH ) {

    unsigned k1, k2, k3;
    __m128i vix, viy, viz;

    // load VEC_WIDTH particles
    LOAD4P3D( part, i, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz )

    // Store cell index
    vix = _mm_add_epi32( vix, _mm_add_epi32( vix, vix ) );
    viy = _mm_mullo_epi32( viy, vDy );
    viz = _mm_mullo_epi32( viz, vDz );
    _mm_store_si128( (__m128i *) idx, _mm_add_epi32( _mm_add_epi32( vix, viy ), viz ) );

    // Get spline weights
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );
    SPLINE( vz0, vs0z );
    SPLINE( vz1, vs1z );

    vqnx = _mm256_mul_pd( vq, vnorm1);
    vqny = _mm256_mul_pd( vq, vnorm2);
    vqnz = _mm256_mul_pd( vq, vnorm3);

    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 ); vwl1[NP-1] = _mm256_setzero_pd();
    WL( vqny, vy0, vy1, vwl2 ); vwl2[NP-1] = _mm256_setzero_pd();
    WL( vqnz, vz0, vz1, vwl3 ); vwl3[NP-1] = _mm256_setzero_pd();

    // get perpendicular weights
    for( k2 = 0; k2 < NP; k2++ ) {
      for( k1 = 0; k1 < NP; k1++ ) {

        const __m256d c1_2 = _mm256_set1_pd(0.5f);
        __m256d tmp1, tmp2;

        // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] +
        //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )

        tmp1 = _mm256_mul_pd( vs0y[k1], vs0z[k2] );
        tmp2 = _mm256_mul_pd( vs0y[k1], vs1z[k2] );
        tmp1 = _mm256_fmadd_pd( vs1y[k1], vs1z[k2], tmp1 );
        tmp2 = _mm256_fmadd_pd( vs1y[k1], vs0z[k2], tmp2 );

        vwp1[k2][k1] = _mm256_fmadd_pd( c1_2, tmp2, tmp1 );

        // wp2[k2][k1] = s0x[k1]*s0z[k2] + s1x[k1]*s1z[k2] +
        //               0.5*( s0x[k1]*s1z[k2] + s1x[k1]*s0z[k2] )

        tmp1 = _mm256_mul_pd( vs0x[k1], vs0z[k2] );
        tmp2 = _mm256_mul_pd( vs0x[k1], vs1z[k2] );
        tmp1 = _mm256_fmadd_pd( vs1x[k1], vs1z[k2], tmp1 );
        tmp2 = _mm256_fmadd_pd( vs1x[k1], vs0z[k2], tmp2 );

        vwp2[k2][k1] = _mm256_fmadd_pd( c1_2, tmp2, tmp1 );

        // wp3[k2][k1] = s0x[k1]*s0y[k2] + s1x[k1]*s1y[k2] +
        //               0.5*( s0x[k1]*s1y[k2] + s1x[k1]*s0y[k2] )

        tmp1 = _mm256_mul_pd( vs0x[k1], vs0y[k2] );
        tmp2 = _mm256_mul_pd( vs0x[k1], vs1y[k2] );
        tmp1 = _mm256_fmadd_pd( vs1x[k1], vs1y[k2], tmp1 );
        tmp2 = _mm256_fmadd_pd( vs1x[k1], vs0y[k2], tmp2 );

        vwp3[k2][k1] = _mm256_fmadd_pd( c1_2, tmp2, tmp1 );

      }
    }

    // Accumulate current 1 plane at a time using the 2D algorithm
    for( k3 = 0; k3 < NP; k3++ ) {
      __m256d a[NP][NP], b[NP][NP], c[NP][NP], d[NP][NP];

      for( k2 = 0; k2 < NP; k2++ ) {

        for ( k1 = 0; k1 < NP; k1++ ) {
          __m256d j1, j2, j3;

          // Calculate current components
          j1 = _mm256_mul_pd( vwl1[k1] , vwp1[k3][k2] );
          j2 = _mm256_mul_pd( vwl2[k2] , vwp2[k3][k1] );
          j3 = _mm256_mul_pd( vwl3[k3] , vwp3[k2][k1] );

          // Do a 4x4 transpose, setting the last term to 0
          __m256d const zero = _mm256_setzero_pd();

          __m256d t0 = _mm256_shuffle_pd( j1, j2, 0x0 );
          __m256d t2 = _mm256_shuffle_pd( j1, j2, 0xF );
          __m256d t1 = _mm256_shuffle_pd( j3, zero, 0x0 );
          __m256d t3 = _mm256_shuffle_pd( j3, zero, 0xEE );

          a[k2][k1] = _mm256_permute2f128_pd( t0, t2, 0x20 );
          b[k2][k1] = _mm256_permute2f128_pd( t0, t2, 0x31 );
          c[k2][k1] = _mm256_permute2f128_pd( t1, t3, 0x20 );
          d[k2][k1] = _mm256_permute2f128_pd( t1, t3, 0x31 );
        }
      }

      // Accumulate electric current in k3 plane
      ACC_CURRENT3D( 0, a )
      ACC_CURRENT3D( 1, b )
      ACC_CURRENT3D( 2, c )
      ACC_CURRENT3D( 3, d )

    }

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








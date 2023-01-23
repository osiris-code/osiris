/*****************************************************************************************

Charge conserving current deposition, AVX-512 optimized version (double precision)

*****************************************************************************************/


/* if __TEMPLATE__ is defined then read the template definition at the end of the file */
#ifndef __TEMPLATE__

#include "os-spec-current-avx512d.h"
#include "splines-avx512.h"
/*****************************************************************************************
vwl_s1
*****************************************************************************************/

inline void vwl_s1( __m512d const vqn, __m512d const vx0, __m512d const vx1, __m512d vwl[] )
{

  vwl[0] = _mm512_mul_pd(vqn, _mm512_sub_pd(vx1, vx0));
}


/*****************************************************************************************
vwl_s2
*****************************************************************************************/

inline void vwl_s2( __m512d const vqn, __m512d const vx0, __m512d const vx1, __m512d vwl[] )
{

  const __m512d c1_2 = _mm512_set1_pd( 0.5 );

  __m512d d    = _mm512_sub_pd( vx1, vx0 );
  __m512d s1_2 = _mm512_sub_pd( c1_2, vx0 );
  __m512d p1_2 = _mm512_add_pd( c1_2, vx0 );

  __m512d n    = _mm512_mul_pd( vqn, d );

  vwl[0] = _mm512_mul_pd( n, _mm512_fnmadd_pd( c1_2, d, s1_2 ));
  vwl[1] = _mm512_mul_pd( n, _mm512_fmadd_pd(  c1_2, d, p1_2 ));

}


/*****************************************************************************************
vwl_s3
*****************************************************************************************/

inline void vwl_s3( __m512d const vqn, __m512d const vx0, __m512d const vx1, __m512d vwl[] )
{
  __m512d const c1_2 = _mm512_set1_pd( 0.5 );
  __m512d const c3_4 = _mm512_set1_pd( 0.75 );
  __m512d const c1_3 = _mm512_set1_pd( 1.0/3.0 );

  __m512d   d = _mm512_sub_pd(  vx1, vx0 );
  __m512d   s = _mm512_sub_pd( c1_2, vx0 );
  __m512d   p = _mm512_add_pd( c1_2, vx0 );
  __m512d d_3 = _mm512_mul_pd( c1_3,   d );

  vwl[0] = _mm512_mul_pd( c1_2, _mm512_fmsub_pd( s, s, _mm512_mul_pd( d, _mm512_sub_pd( s, d_3 )) ) );
  vwl[1] = _mm512_fnmadd_pd( d, _mm512_add_pd( vx0, d_3 ), _mm512_fnmadd_pd( vx0, vx0, c3_4 ));
  vwl[2] = _mm512_mul_pd( c1_2, _mm512_fmadd_pd( p, p, _mm512_mul_pd( d, _mm512_add_pd( p, d_3 )) ) );

  __m512d n = _mm512_mul_pd( vqn, d );
  vwl[0] = _mm512_mul_pd( n, vwl[0] );
  vwl[1] = _mm512_mul_pd( n, vwl[1] );
  vwl[2] = _mm512_mul_pd( n, vwl[2] );
}

/*****************************************************************************************
vwl_s4
*****************************************************************************************/

inline void vwl_s4( __m512d const vqn, __m512d const vx0, __m512d const vx1, __m512d vwl[] ) {

  __m512d const c1_2 = _mm512_set1_pd( 0.5 );
  __m512d const c1_4 = _mm512_set1_pd( 0.25 );
  __m512d const c1_6 = _mm512_set1_pd( 1.0/6.0 );
  __m512d const c3_2 = _mm512_set1_pd( 1.5 );
  __m512d const c3_4 = _mm512_set1_pd( 0.75 );
  __m512d const c1_3 = _mm512_set1_pd( 1.0/3.0 );
  __m512d const c2_3 = _mm512_set1_pd( 2.0/3.0 );

  __m512d d = _mm512_sub_pd( vx1, vx0 );
  __m512d s = _mm512_sub_pd( c1_2, vx0 );
  __m512d p = _mm512_add_pd( c1_2, vx0 );

  __m512d t = _mm512_sub_pd( vx0, c1_6 );
  __m512d u = _mm512_add_pd( vx0, c1_6 );

  __m512d s2 = _mm512_mul_pd( s, s );
  __m512d s3 = _mm512_mul_pd( s2, s );

  __m512d p2 = _mm512_mul_pd( p, p );
  __m512d p3 = _mm512_mul_pd( p2, p );

  __m512d d_2 = _mm512_mul_pd( d, c1_2 );
  __m512d d_4 = _mm512_mul_pd( d, c1_4 );

  vwl[0] = _mm512_mul_pd( c1_6, _mm512_fmadd_pd( d, _mm512_fmsub_pd( d, _mm512_sub_pd( s, d_4 ), _mm512_mul_pd( c3_2, s2 )  ), s3 ) );
  vwl[1] = _mm512_fmadd_pd( d, _mm512_fmadd_pd(d_2, _mm512_add_pd(t, d_4), _mm512_fmsub_pd(c3_4, _mm512_mul_pd(t, t), c1_3)),
                          _mm512_fmadd_pd(c1_2, p3, _mm512_sub_pd(c2_3, p2)));
  vwl[2] = _mm512_fnmadd_pd(d, _mm512_fmadd_pd(d_2, _mm512_add_pd(u, d_4), _mm512_fmsub_pd(c3_4, _mm512_mul_pd(u, u), c1_3)),
                          _mm512_fmadd_pd(c1_2, s3, _mm512_sub_pd(c2_3, s2)));
  vwl[3] = _mm512_mul_pd( c1_6, _mm512_fmadd_pd( d, _mm512_fmadd_pd( d, _mm512_add_pd( p, d_4 ), _mm512_mul_pd( c3_2, p2 )  ), p3 ) );

  __m512d n = _mm512_mul_pd( vqn, d );
  vwl[0] = _mm512_mul_pd( n, vwl[0] );
  vwl[1] = _mm512_mul_pd( n, vwl[1] );
  vwl[2] = _mm512_mul_pd( n, vwl[2] );
  vwl[3] = _mm512_mul_pd( n, vwl[3] );

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

#if 0

#warning Using reference 1D algorithm

void DEP_CURRENT_1D
(double * const current, int const * const size, int const * const offset,
 double * const norm, t_split_buf1D * const part)
{

  typedef struct Current { double j1, j2, j3; } t_current;
  t_current *p0;

  dvec j2[NP],j3[NP];
  dvec vwl1[ORDER];

  __m512d vx0, vx1, vq, vvy, vvz;

  DECLARE_ALIGNED_64( int idx[VEC_WIDTH] );

  __m512d vqnx, vqvy, vqvz;
  __m512d vs0x[NP], vs1x[NP];

  int i, np;

  __m512d const c1_2  = _mm512_set1_pd( 0.5 );

  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 );

  __m512d const vnorm1 = _mm512_set1_pd(norm[0]);

  np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       int k = np + i;

       part->x0[k] = part->x1[k] = 0.0;
       part-> q[k] = part->vy[k] = part->vz[k] = 0.0;

       part->ix[k] = 1;
    }
  }

  for( i = 0; i < np; i += VEC_WIDTH ) {

    int k, k1;
    __m256i vix;

    // load 8 particles
    LOAD8P1D( part, i, vx0, vx1, vq, vvy, vvz, vix )

    // Store cell index
    _mm256_store_si256( (__m256i *) idx, vix );

    // Get splines
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );

    vqnx = _mm512_mul_pd( vq, vnorm1 );
    vqvy = _mm512_mul_pd( vq, vvy );
    vqvz = _mm512_mul_pd( vq, vvz );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m512d *) vwl1 );

    // get j2,j3 current
    for ( k1 = 0; k1 < NP; k1++ ) {
      __m512d vwp1 = _mm512_mul_pd( c1_2, _mm512_add_pd(vs0x[k1], vs1x[k1]) );
      j2[k1].v8  = _mm512_mul_pd( vqvy, vwp1 );
      j3[k1].v8  = _mm512_mul_pd( vqvz, vwp1 );
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
        Small ~10% improvement from reference implementation
*/

#if 1

#warning Using mk. I 1D - algorithm

#if ( ORDER == 1 )

#define ACC_CURRENT1D( k, v, p ) { \
  double *p0 = pj + idx[k]; \
  __m512d va, ja; \
  va = _mm512_mask_expandloadu_pd( va, 0x77, (void *) (p0) ); \
  ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[0]), p )); \
  ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[1]), p )); \
  va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
  _mm512_mask_compressstoreu_pd( (void *) (p0), 0x77, va); \
}

#elif ( ORDER == 2 )

#define ACC_CURRENT1D( k, v, p ) { \
  double *p0 = pj + idx[k]; \
  __m512d va, vb, ja, jb; \
  va = _mm512_mask_expandloadu_pd( va, 0x77, (void *) (p0) ); \
  vb = _mm512_loadu_pd( (void *) (p0 + 6) ); \
  ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[0]), p )); \
  ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[1]), p )); \
  jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[2]), p )); \
  va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
  vb = _mm512_mask_add_pd( vb, 0x07 , jb, vb ); \
  _mm512_mask_compressstoreu_pd( (void *) (p0), 0x77, va); \
  _mm512_mask_storeu_pd( (void *) (p0 + 6), 0x07, vb) ; \
}

#elif ( ORDER == 3 )

#define ACC_CURRENT1D( k, v, p ) {                                      \
  double *p0 = pj + idx[k]; \
  __m512d va, vb, ja, jb; \
  va = _mm512_mask_expandloadu_pd( va, 0x77, (void *) (p0) ); \
  vb = _mm512_mask_expandloadu_pd( vb, 0x77, (void *) (p0 + 6) ); \
  ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[0]), p )); \
  ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[1]), p )); \
  jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[2]), p )); \
  jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x3F00, _mm512_castpd_ps(v[3]), p )); \
  va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
  vb = _mm512_mask_add_pd( vb, 0x77 , jb, vb ); \
  _mm512_mask_compressstoreu_pd( (void *) (p0), 0x77, va); \
  _mm512_mask_compressstoreu_pd( (void *) (p0 + 6), 0x77, vb); \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT1D( k, v, p ) {                                   \
  double *p0 = pj + idx[k]; \
  __m512d va, vb, vc, ja, jb, jc; \
  va = _mm512_mask_expandloadu_pd( va, 0x77, (void *) (p0) ); \
  vb = _mm512_mask_expandloadu_pd( vb, 0x77, (void *) (p0 + 6) ); \
  vc = _mm512_loadu_pd( (void *) (p0 + 12) ); \
  ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[0]), p )); \
  ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[1]), p )); \
  jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[2]), p )); \
  jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x3F00, _mm512_castpd_ps(v[3]), p )); \
  jc = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jc), 0x003F, _mm512_castpd_ps(v[4]), p )); \
  va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
  vb = _mm512_mask_add_pd( vb, 0x77 , jb, vb ); \
  vc = _mm512_mask_add_pd( vc, 0x07 , jc, vc ); \
  _mm512_mask_compressstoreu_pd( (void *) (p0),      0x77, va); \
  _mm512_mask_compressstoreu_pd( (void *) (p0 + 6),  0x77, vb); \
  _mm512_mask_storeu_pd(         (void *) (p0 + 12), 0x07, vc); \
}

#else
#error Unsupported interpolation order
#endif


void DEP_CURRENT_1D
(double * const current, int const * const size, int const * const offset,
                       double * const norm, t_split_buf1D * const part)
{


  __m512d vx0, vx1, vq, vvy, vvz;

  __m512d vqnx, vqvy, vqvz;
  __m512d vs0x[NP], vs1x[NP];

  int i, np;

  __m512d const c1_2  = _mm512_set1_pd( 0.5 );

  double* const pj =  current + ( offset[0] - OFFSET ) * 3 ;

  __m512d const vnorm1 = _mm512_set1_pd(norm[0]);

  np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       int k = np + i;

       part->x0[k] = part->x1[k] = 0.0;
       part-> q[k] = part->vy[k] = part->vz[k] = 0.0;

       part->ix[k] = 1;
    }
  }

  for( i = 0; i < np; i += VEC_WIDTH ) {

    int k1;

    __m512d vwl1[NP];

    __m256i vix;
    DECLARE_ALIGNED_64( int idx[VEC_WIDTH] );
    __m512d a[NP], b[NP], c[NP], d[NP];

    // load VEC_WIDTH particles
    LOAD8P1D( part, i, vx0, vx1, vq, vvy, vvz, vix )

    // Store cell index
    vix = _mm256_add_epi32( vix, _mm256_add_epi32( vix, vix ) );
    _mm256_store_si256( (__m256i *) idx, vix );

    // Get splines
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );

    vqnx = _mm512_mul_pd( vq, vnorm1 );
    vqvy = _mm512_mul_pd( vq, vvy );
    vqvz = _mm512_mul_pd( vq, vvz );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 ); vwl1[NP-1] = _mm512_setzero_pd();

    // get currents and transpose vectors
    for ( k1 = 0; k1 < NP; k1++ ) {
      __m512d vwp1, j1, j2, j3;

      vwp1 = _mm512_mul_pd( c1_2, _mm512_add_pd(vs0x[k1], vs1x[k1]) );

      j1  = vwl1[k1];
      j2  = _mm512_mul_pd( vqvy, vwp1 );
      j3  = _mm512_mul_pd( vqvz, vwp1 );

      // Do a 4x4 transpose within 256 bit lanes ignoring 4th current component
      __m512d t0 = _mm512_mask_swizzle_pd( j1, 0xAA, j2, _MM_SWIZ_REG_CDAB );
      __m512d t1 = _mm512_mask_swizzle_pd( j3, 0xAA, j3, _MM_SWIZ_REG_CDAB );
      __m512d t2 = _mm512_mask_swizzle_pd( j2, 0x55, j1, _MM_SWIZ_REG_CDAB );
      __m512d t3 = _mm512_mask_swizzle_pd( j3, 0x55, j3, _MM_SWIZ_REG_CDAB );

      a[k1] = _mm512_mask_swizzle_pd( t0, 0xCC, t1, _MM_SWIZ_REG_BADC );
      b[k1] = _mm512_mask_swizzle_pd( t2, 0xCC, t3, _MM_SWIZ_REG_BADC );
      c[k1] = _mm512_mask_swizzle_pd( t1, 0x33, t0, _MM_SWIZ_REG_BADC );
      d[k1] = _mm512_mask_swizzle_pd( t3, 0x33, t2, _MM_SWIZ_REG_BADC );

    }

    ACC_CURRENT1D(  0, a, _MM_PERM_BABA )
    ACC_CURRENT1D(  1, b, _MM_PERM_BABA )
    ACC_CURRENT1D(  2, c, _MM_PERM_BABA )
    ACC_CURRENT1D(  3, d, _MM_PERM_BABA )
    ACC_CURRENT1D(  4, a, _MM_PERM_DCDC )
    ACC_CURRENT1D(  5, b, _MM_PERM_DCDC )
    ACC_CURRENT1D(  6, c, _MM_PERM_DCDC )
    ACC_CURRENT1D(  7, d, _MM_PERM_DCDC )

  }

}

#undef ACC_CURRENT1D

#endif

/****************************************************************************************
  2D current deposition
*****************************************************************************************/

#if 0

/**
 * Reference implementation
 *
 * For each particle, deposit one current component at a time
 *
 * This is the slowest
 */

#warning Using reference 2D algorithm

void DEP_CURRENT_2D
(double * const current, int const * const size, int const * const offset,
                       double * const norm, t_split_buf2D * const part)
{

  typedef struct Current { double j1, j2, j3; } t_current;
  t_current *pjpart, *p0;


  dvec j3[NP][NP];
  dvec vwp1[NP], vwp2[NP];
  dvec vwl1[ORDER], vwl2[ORDER];

  __m512d vx0, vx1, vy0, vy1, vq, vvz;

  DECLARE_ALIGNED_64( int idx[VEC_WIDTH] );

  __m512d vqnx, vqny, vqvz;
  __m512d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP];

  int np;
  int i, k, k1, k2;

  __m512d const c1_3 = _mm512_set1_pd( 1.0/3.0 );
  __m512d const c1_2  = _mm512_set1_pd( 0.5 );

  int const Dy = size[0];
  t_current* const pj = (t_current *) (current + ( offset[0] - OFFSET ) * 3 +
                                                 ( offset[1] - OFFSET ) * 3 * Dy );

  // Setting this to 1/4 removes 1 multiplication from wl2
  __m512d const vnorm1 = _mm512_set1_pd(norm[0]);
  __m512d const vnorm2 = _mm512_set1_pd(norm[1]);

  np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH != 0 ) {
    for( i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       int k = np + i;

       part->x0[k] = part->x1[k] = 0.0;
       part->y0[k] = part->y1[k] = 0.0;
       part-> q[k] = part->vz[k] = 0.0;

       part->ix[k] = part->iy[k] = 1;
    }
  }

  // Each virtual particle uses 8 doubles, so 4 vp = 32 doubles
  for( i = 0; i < np; i += VEC_WIDTH ) {

    __m256i vix, viy;

    // load 8 particles
    LOAD8P2D( part, i, vx0, vx1, vy0, vy1, vq, vvz, vix, viy )

    // Store cell index
    _mm256_store_si256( (__m256i *) idx,
            _mm256_add_epi32( vix, _mm256_mullo_epi32( viy, _mm256_set1_epi32( Dy ) ) ) );

    // Get splines
    SPLINE( vx0, vs0x );
    SPLINE( vx1, vs1x );
    SPLINE( vy0, vs0y );
    SPLINE( vy1, vs1y );

    vqnx = _mm512_mul_pd( vq, vnorm1 );
    vqny = _mm512_mul_pd( vq, vnorm2 );
    vqvz = _mm512_mul_pd( vq, vvz);
    vqvz = _mm512_mul_pd( vqvz, c1_3 );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m512d *) vwl1 );
    WL( vqny, vy0, vy1, (__m512d *) vwl2 );

    // get perpendicular weights
    for( k = 0; k < NP; k++ ) {
      vwp1[k].v8 = _mm512_add_pd(vs0y[k], vs1y[k]);
      vwp2[k].v8 = _mm512_add_pd(vs0x[k], vs1x[k]);
    }

    // get j3 current
    for( k2 = 0; k2 < NP; k2++ ) {
      for ( k1 = 0; k1 < NP; k1++ ) {
        __m512d s00, s10, tmp1, tmp2;

    s00  = _mm512_mul_pd( vs0x[k1], vs0y[k2] );
    tmp1 = _mm512_fmadd_pd( vs1y[k2], vs1x[k1], s00 );  // tmp1 = s0x*s0y + s1x*s1y

    s10  = _mm512_mul_pd( vs0x[k1], vs1y[k2] );
    tmp2 = _mm512_fmadd_pd( vs0y[k2], vs1x[k1], s10 );  // tmp2 = s0x*s1y + s1x*s0y

    tmp1 = _mm512_fmadd_pd( tmp2, c1_2, tmp1 );       // tmp1 = tmp1 + 0.5*tmp2

    j3[k1][k2].v8 = _mm512_mul_pd( vqvz, tmp1 );          // j3 = vqvz * tmp1
      }
    }

    // Loop by particle on the outside loop
    for ( k = 0; k < VEC_WIDTH; k ++ ) {

      // This is not vectorially because it's a 64bit op
      pjpart = pj + idx[k];

    // accumulate j1
    for( k2 = 0; k2 < NP; k2++ ) {
    p0 = pjpart + k2*Dy;
    for ( k1 = 0; k1 < ORDER; k1++ ) {
      p0[k1].j1 += vwl1[k1].v[k] * vwp1[k2].v[k];
    }
    }

    // accumulate j2 - making k2 the outside loop gives marginal perf. gain
    for( k2 = 0; k2 < ORDER; k2++ ) {
    p0 = pjpart + k2*Dy;
    for ( k1 = 0; k1 < NP; k1++ ) {
       p0[k1].j2 += vwl2[k2].v[k] * vwp2[k1].v[k];
    }
    }

    // accumulate j3
    for( k2 = 0; k2 < NP; k2++ ) {
    p0 = pjpart + k2*Dy;
    for ( k1 = 0; k1 < NP; k1++ ) {
       p0[k1].j3 += (j3[k1][k2]).v[k];
    }
    }

    }

  }

}
#endif


#if 1

/**
 * mk. I
 *
 * Transpose currents, then for each particle deposit one current line at a time
 *
 * This is the recommended version
 */


#warning Using mk. I 2D algorithm

#if ( ORDER == 1 )

#define ACC_CURRENT2D( k, v, p ) { \
  double * const line = pj + idx[k]; \
  for( int k2 = 0; k2 < 2; k2++) { \
    __m512d va, ja; \
    va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy ] ); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k2][0]), p )); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k2][1]), p )); \
    va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy ], 0x77, va); \
  } \
}



#elif ( ORDER == 2 )

#define ACC_CURRENT2D( k, v, p ) { \
  double * const line = pj + idx[k]; \
  for( int k2 = 0; k2 < 3; k2++) { \
    __m512d va, vb, ja, jb; \
    va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy     ] ); \
    vb = _mm512_loadu_pd(                      &line[ k2*Dy + 6 ] ); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k2][0]), p )); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k2][1]), p )); \
    jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[k2][2]), p )); \
    va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
    vb = _mm512_mask_add_pd( vb, 0x07 , jb, vb ); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy     ], 0x77, va); \
    _mm512_mask_storeu_pd(         &line[ k2*Dy + 6 ], 0x07, vb) ; \
  } \
}

#elif ( ORDER == 3 )

#define ACC_CURRENT2D( k, v, p ) { \
  double * const line = pj + idx[k]; \
  for( int k2 = 0; k2 < 4; k2++) { \
    __m512d va, vb, ja, jb; \
    va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy     ] ); \
    vb = _mm512_mask_expandloadu_pd( vb, 0x77, &line[ k2*Dy + 6 ] ); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k2][0]), p )); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k2][1]), p )); \
    jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[k2][2]), p )); \
    jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x3F00, _mm512_castpd_ps(v[k2][3]), p )); \
    va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
    vb = _mm512_mask_add_pd( vb, 0x77 , jb, vb ); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy     ], 0x77, va); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy + 6 ], 0x77, vb); \
  } \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT2D( k, v, p ) { \
  double * const line = pj + idx[k]; \
  for( int k2 = 0; k2 < 5; k2++) { \
    __m512d va, vb, vc, ja, jb, jc; \
    va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy      ] ); \
    vb = _mm512_mask_expandloadu_pd( vb, 0x77, &line[ k2*Dy +  6 ] ); \
    vc = _mm512_loadu_pd(                      &line[ k2*Dy + 12 ] ); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k2][0]), p )); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k2][1]), p )); \
    jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[k2][2]), p )); \
    jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x3F00, _mm512_castpd_ps(v[k2][3]), p )); \
    jc = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jc), 0x003F, _mm512_castpd_ps(v[k2][4]), p )); \
    va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
    vb = _mm512_mask_add_pd( vb, 0x77 , jb, vb ); \
    vc = _mm512_mask_add_pd( vc, 0x07 , jc, vc ); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy      ],  0x77, va); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy +  6 ],  0x77, vb); \
    _mm512_mask_storeu_pd(         &line[ k2*Dy + 12 ],  0x07, vc); \
  } \
}

#else
#error Unsupported interpolation order
#endif


void DEP_CURRENT_2D
(double * restrict const current, int const * restrict const size, int const * restrict const offset,
                       double * restrict const norm, t_split_buf2D * restrict const part)
{

  __m512d vwp1[NP], vwp2[NP];
  __m512d vwl1[NP], vwl2[NP];

  __m512d vx0, vx1, vy0, vy1, vq, vvz;

  DECLARE_ALIGNED_64( int idx[VEC_WIDTH] );

  __m512d vqnx, vqny, vqvz;
  __m512d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP];

  __m512d const c1_3  = _mm512_set1_pd( 1.0/3.0 );
  __m512d const c1_2  = _mm512_set1_pd( 0.5 );

  // Current array has 3 field components
  double* const pj = current + 3 * ( ( offset[0] - OFFSET ) +
                                     ( offset[1] - OFFSET ) * size[0] );

  // Normalization for current
  __m512d const vnorm1 = _mm512_set1_pd(norm[0]);
  __m512d const vnorm2 = _mm512_set1_pd(norm[1]);

  int     const Dy = 3*size[0];
  __m256i const vDy = _mm256_set1_epi32( Dy );

  int np = part -> np;

  // If number of particles in buffer is not multiple of VEC_WIDTH add dummy particles to the end

  if ( np % VEC_WIDTH != 0 ) {
    for( int i = 0; i < VEC_WIDTH - np%VEC_WIDTH; i ++ ) {
       int k = np + i;

       part->x0[k] = part->x1[k] = 0.0;
       part->y0[k] = part->y1[k] = 0.0;
       part-> q[k] = part->vz[k] = 0.0;

       part->ix[k] = part->iy[k] = 1;
    }
  }

  for( int i = 0; i < np; i += VEC_WIDTH ) {

    __m256i vix, viy;

    // load 8 particles
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

    vqnx = _mm512_mul_pd( vq, vnorm1 );
    vqny = _mm512_mul_pd( vq, vnorm2 );
    vqvz = _mm512_mul_pd( vq, vvz);
    vqvz = _mm512_mul_pd( vqvz, c1_3 );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 ); vwl1[NP-1] = _mm512_setzero_pd();
    WL( vqny, vy0, vy1, vwl2 ); vwl2[NP-1] = _mm512_setzero_pd();

    // get perpendicular weights
    for( int k = 0; k < NP; k++ ) {
      vwp1[k] = _mm512_add_pd(vs0y[k], vs1y[k]);
      vwp2[k] = _mm512_add_pd(vs0x[k], vs1x[k]);
    }


    __m512d a[NP][NP], b[NP][NP], c[NP][NP], d[NP][NP];

    for( int k2 = 0; k2 < NP; k2++ ) {
      for ( int k1 = 0; k1 < NP; k1++ ) {
         __m512d j1, j2, j3;

        j1 = _mm512_mul_pd( vwl1[k1] , vwp1[k2] );
        j2 = _mm512_mul_pd( vwl2[k2] , vwp2[k1] );

        __m512d s00, s10, tmp1, tmp2;

        s00  = _mm512_mul_pd( vs0x[k1], vs0y[k2] );
        tmp1 = _mm512_fmadd_pd( vs1y[k2], vs1x[k1], s00 );

        s10  = _mm512_mul_pd( vs0x[k1], vs1y[k2] );
        tmp2 = _mm512_fmadd_pd( vs0y[k2], vs1x[k1], s10 );

        tmp1 = _mm512_fmadd_pd( tmp2, c1_2, tmp1 );

        j3 = _mm512_mul_pd( vqvz, tmp1 );

        // Do a 4x4 transpose within 256 bit lanes ignoring 4th current component
        __m512d t0 = _mm512_mask_swizzle_pd( j1, 0xAA, j2, _MM_SWIZ_REG_CDAB );
        __m512d t1 = _mm512_mask_swizzle_pd( j3, 0xAA, j3, _MM_SWIZ_REG_CDAB );
        __m512d t2 = _mm512_mask_swizzle_pd( j2, 0x55, j1, _MM_SWIZ_REG_CDAB );
        __m512d t3 = _mm512_mask_swizzle_pd( j3, 0x55, j3, _MM_SWIZ_REG_CDAB );

        a[k2][k1] = _mm512_mask_swizzle_pd( t0, 0xCC, t1, _MM_SWIZ_REG_BADC );
        b[k2][k1] = _mm512_mask_swizzle_pd( t2, 0xCC, t3, _MM_SWIZ_REG_BADC );
        c[k2][k1] = _mm512_mask_swizzle_pd( t1, 0x33, t0, _MM_SWIZ_REG_BADC );
        d[k2][k1] = _mm512_mask_swizzle_pd( t3, 0x33, t2, _MM_SWIZ_REG_BADC );
      }
    }


    ACC_CURRENT2D(  0, a, _MM_PERM_BABA )
    ACC_CURRENT2D(  1, b, _MM_PERM_BABA )
    ACC_CURRENT2D(  2, c, _MM_PERM_BABA )
    ACC_CURRENT2D(  3, d, _MM_PERM_BABA )
    ACC_CURRENT2D(  4, a, _MM_PERM_DCDC )
    ACC_CURRENT2D(  5, b, _MM_PERM_DCDC )
    ACC_CURRENT2D(  6, c, _MM_PERM_DCDC )
    ACC_CURRENT2D(  7, d, _MM_PERM_DCDC )

  }

}

#undef ACC_CURRENT2D

#endif

/****************************************************************************************
  3D current deposition
****************************************************************************************/

/**
 * Reference implementation
 *
 * For each particle, accumulate 1 field component at a time.
 *
 * This is the fastest version that has the same roundoff as the Fortran version.
 */

#if 0

#warning Using reference 3D algorithm

void DEP_CURRENT_3D
(double * const  current, int const * const  size, int const * const  offset,
                       double * const norm, t_split_buf3D * const part)
{

  typedef struct Current { double j1, j2, j3; } t_current;

  __m512d vx0, vx1, vy0, vy1,  vz0, vz1, vq;

  __m512d vqnx, vqny, vqnz;
  __m512d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP];

  dvec vwl1[ORDER], vwl2[ORDER], vwl3[ORDER];
  dvec vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];

  DECLARE_ALIGNED_64( int idx[VEC_WIDTH] );

  int const Dy = size[0];
  int const Dz = Dy * size[1];

  t_current*  pj = (t_current *) ( current + ( offset[0] - OFFSET  ) * 3 +
                                             ( offset[1] - OFFSET  ) * 3 * Dy +
                                             ( offset[2] - OFFSET  ) * 3 * Dz );

  __m512d const vnorm1 = _mm512_set1_pd(norm[0]);
  __m512d const vnorm2 = _mm512_set1_pd(norm[1]);
  __m512d const vnorm3 = _mm512_set1_pd(norm[2]);

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

    // load 8 particles
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

    vqnx = _mm512_mul_pd( vq, vnorm1 );
    vqny = _mm512_mul_pd( vq, vnorm2 );
    vqnz = _mm512_mul_pd( vq, vnorm3 );

    // get longitudinal weights
    WL( vqnx, vx0, vx1, (__m512d *) vwl1 );
    WL( vqny, vy0, vy1, (__m512d *) vwl2 );
    WL( vqnz, vz0, vz1, (__m512d *) vwl3 );

    // get perpendicular weights
    for( int k2 = 0; k2 < NP; k2++ ) {
      for( int k1 = 0; k1 < NP; k1++ ) {

         const __m512d c1_2 = _mm512_set1_pd(0.5);
         __m512d tmp1, tmp2;

         // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] +
         //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )

         tmp1 = _mm512_mul_pd( vs0y[k1], vs0z[k2] );
         tmp2 = _mm512_mul_pd( vs0y[k1], vs1z[k2] );
         tmp1 = _mm512_fmadd_pd( vs1y[k1], vs1z[k2], tmp1 );
         tmp2 = _mm512_fmadd_pd( vs1y[k1], vs0z[k2], tmp2 );

         vwp1[k2][k1].v8 = _mm512_fmadd_pd( c1_2, tmp2, tmp1 );

         // wp2[k2][k1] = s0x[k1]*s0z[k2] + s1x[k1]*s1z[k2] +
         //               0.5*( s0x[k1]*s1z[k2] + s1x[k1]*s0z[k2] )

         tmp1 = _mm512_mul_pd( vs0x[k1], vs0z[k2] );
         tmp2 = _mm512_mul_pd( vs0x[k1], vs1z[k2] );
         tmp1 = _mm512_fmadd_pd( vs1x[k1], vs1z[k2], tmp1 );
         tmp2 = _mm512_fmadd_pd( vs1x[k1], vs0z[k2], tmp2 );

         vwp2[k2][k1].v8 = _mm512_fmadd_pd( c1_2, tmp2, tmp1 );

         // wp3[k2][k1] = s0x[k1]*s0y[k2] + s1x[k1]*s1y[k2] +
         //               0.5*( s0x[k1]*s1y[k2] + s1x[k1]*s0y[k2] )

         tmp1 = _mm512_mul_pd( vs0x[k1], vs0y[k2] );
         tmp2 = _mm512_mul_pd( vs0x[k1], vs1y[k2] );
         tmp1 = _mm512_fmadd_pd( vs1x[k1], vs1y[k2], tmp1 );
         tmp2 = _mm512_fmadd_pd( vs1x[k1], vs0y[k2], tmp2 );

         vwp3[k2][k1].v8 = _mm512_fmadd_pd( c1_2, tmp2, tmp1 );

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

/**
 * mk. I
 *
 * For each vector of particles loop over z coordinate and deposit one current plane
 * at a time using the same algorithm as 2D mk. I
 *
 * Note: This will give different roundoff from the other versions because we are not
 * adding 1 particle at a time.
 *
 * This is the fastest version
 */

#warning Using mk. I 3D algorithm

#if ( ORDER == 1 )

#define ACC_CURRENT3D( k, v, p ) { \
  double * const line = pj + idx[k] + k3*Dz; \
  for( int k2 = 0; k2 < 2; k2++) { \
    __m512d va, ja; \
    va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy ] ); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k2][0]), p )); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k2][1]), p )); \
    va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy ], 0x77, va); \
  } \
}

#elif ( ORDER == 2 )

#define ACC_CURRENT3D( k, v, p ) { \
  double * const line = pj + idx[k] + k3*Dz; \
  for( int k2 = 0; k2 < 3; k2++) { \
    __m512d va, vb, ja, jb; \
    va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy     ] ); \
    vb = _mm512_loadu_pd(                      &line[ k2*Dy + 6 ] ); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k2][0]), p )); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k2][1]), p )); \
    jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[k2][2]), p )); \
    va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
    vb = _mm512_mask_add_pd( vb, 0x07 , jb, vb ); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy     ], 0x77, va); \
    _mm512_mask_storeu_pd(         &line[ k2*Dy + 6 ], 0x07, vb) ; \
  } \
}

#elif ( ORDER == 3 )

#define ACC_CURRENT3D( k, v, p ) { \
  double * const line = pj + idx[k] + k3*Dz; \
  for( int k2 = 0; k2 < 4; k2++) { \
    __m512d va, vb, ja, jb; \
    va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy     ] ); \
    vb = _mm512_mask_expandloadu_pd( vb, 0x77, &line[ k2*Dy + 6 ] ); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k2][0]), p )); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k2][1]), p )); \
    jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[k2][2]), p )); \
    jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x3F00, _mm512_castpd_ps(v[k2][3]), p )); \
    va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
    vb = _mm512_mask_add_pd( vb, 0x77 , jb, vb ); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy     ], 0x77, va); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy + 6 ], 0x77, vb); \
  } \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT3D( k, v, p ) { \
  double * const line = pj + idx[k] + k3 * Dz; \
  for( int k2 = 0; k2 < 5; k2++) { \
    __m512d va, vb, vc, ja, jb, jc; \
    va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy      ] ); \
    vb = _mm512_mask_expandloadu_pd( vb, 0x77, &line[ k2*Dy +  6 ] ); \
    vc = _mm512_loadu_pd(                      &line[ k2*Dy + 12 ] ); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k2][0]), p )); \
    ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k2][1]), p )); \
    jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[k2][2]), p )); \
    jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x3F00, _mm512_castpd_ps(v[k2][3]), p )); \
    jc = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jc), 0x003F, _mm512_castpd_ps(v[k2][4]), p )); \
    va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
    vb = _mm512_mask_add_pd( vb, 0x77 , jb, vb ); \
    vc = _mm512_mask_add_pd( vc, 0x07 , jc, vc ); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy      ],  0x77, va); \
    _mm512_mask_compressstoreu_pd( &line[ k2*Dy +  6 ],  0x77, vb); \
    _mm512_mask_storeu_pd(         &line[ k2*Dy + 12 ],  0x07, vc); \
  } \
}

#else
#error Unsupported interpolation order
#endif



void DEP_CURRENT_3D
(double * const  current, int const * const  size, int const * const  offset,
                       double * const norm, t_split_buf3D * const part)
{

  __m512d vx0, vx1, vy0, vy1,  vz0, vz1, vq;

  __m512d vqnx, vqny, vqnz;
  __m512d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP];

  __m512d vwl1[NP], vwl2[NP], vwl3[NP];
  __m512d vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];

  DECLARE_ALIGNED_64( int idx[VEC_WIDTH] );

  int const Dy = 3 * size[0];
  int const Dz = Dy * size[1];

  __m256i const vDy = _mm256_set1_epi32( Dy );
  __m256i const vDz = _mm256_set1_epi32( Dz );

  double* const pj = current + ( offset[0] - OFFSET  ) * 3 +
                               ( offset[1] - OFFSET  ) * Dy +
                               ( offset[2] - OFFSET  ) * Dz;

  __m512d const vnorm1 = _mm512_set1_pd(norm[0]);
  __m512d const vnorm2 = _mm512_set1_pd(norm[1]);
  __m512d const vnorm3 = _mm512_set1_pd(norm[2]);



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

    vqnx = _mm512_mul_pd( vq, vnorm1);
    vqny = _mm512_mul_pd( vq, vnorm2);
    vqnz = _mm512_mul_pd( vq, vnorm3);

    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 ); vwl1[NP-1] = _mm512_setzero_pd();
    WL( vqny, vy0, vy1, vwl2 ); vwl2[NP-1] = _mm512_setzero_pd();
    WL( vqnz, vz0, vz1, vwl3 ); vwl3[NP-1] = _mm512_setzero_pd();

    // get perpendicular weights
    for( int k2 = 0; k2 < NP; k2++ ) {
      for( int k1 = 0; k1 < NP; k1++ ) {

        const __m512d c1_2 = _mm512_set1_pd(0.5);
        __m512d tmp1, tmp2;

        // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] +
        //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )

        tmp1 = _mm512_mul_pd( vs0y[k1], vs0z[k2] );
        tmp2 = _mm512_mul_pd( vs0y[k1], vs1z[k2] );
        tmp1 = _mm512_fmadd_pd( vs1y[k1], vs1z[k2], tmp1 );
        tmp2 = _mm512_fmadd_pd( vs1y[k1], vs0z[k2], tmp2 );

        vwp1[k2][k1] = _mm512_fmadd_pd( c1_2, tmp2, tmp1 );

        // wp2[k2][k1] = s0x[k1]*s0z[k2] + s1x[k1]*s1z[k2] +
        //               0.5*( s0x[k1]*s1z[k2] + s1x[k1]*s0z[k2] )

        tmp1 = _mm512_mul_pd( vs0x[k1], vs0z[k2] );
        tmp2 = _mm512_mul_pd( vs0x[k1], vs1z[k2] );
        tmp1 = _mm512_fmadd_pd( vs1x[k1], vs1z[k2], tmp1 );
        tmp2 = _mm512_fmadd_pd( vs1x[k1], vs0z[k2], tmp2 );

        vwp2[k2][k1] = _mm512_fmadd_pd( c1_2, tmp2, tmp1 );

        // wp3[k2][k1] = s0x[k1]*s0y[k2] + s1x[k1]*s1y[k2] +
        //               0.5*( s0x[k1]*s1y[k2] + s1x[k1]*s0y[k2] )

        tmp1 = _mm512_mul_pd( vs0x[k1], vs0y[k2] );
        tmp2 = _mm512_mul_pd( vs0x[k1], vs1y[k2] );
        tmp1 = _mm512_fmadd_pd( vs1x[k1], vs1y[k2], tmp1 );
        tmp2 = _mm512_fmadd_pd( vs1x[k1], vs0y[k2], tmp2 );

        vwp3[k2][k1] = _mm512_fmadd_pd( c1_2, tmp2, tmp1 );

      }
    }

    // Accumulate current 1 plane at a time using the 2D algorithm
    for( int k3 = 0; k3 < NP; k3++ ) {
      __m512d a[NP][NP], b[NP][NP], c[NP][NP], d[NP][NP];

      for( int k2 = 0; k2 < NP; k2++ ) {

        for ( int k1 = 0; k1 < NP; k1++ ) {
          __m512d j1, j2, j3;

          // Calculate current components
          j1 = _mm512_mul_pd( vwl1[k1] , vwp1[k3][k2] );
          j2 = _mm512_mul_pd( vwl2[k2] , vwp2[k3][k1] );
          j3 = _mm512_mul_pd( vwl3[k3] , vwp3[k2][k1] );

          // Do a 4x4 transpose within 256 bit lanes ignoring 4th current component
          __m512d t0 = _mm512_mask_swizzle_pd( j1, 0xAA, j2, _MM_SWIZ_REG_CDAB );
          __m512d t1 = _mm512_mask_swizzle_pd( j3, 0xAA, j3, _MM_SWIZ_REG_CDAB );
          __m512d t2 = _mm512_mask_swizzle_pd( j2, 0x55, j1, _MM_SWIZ_REG_CDAB );
          __m512d t3 = _mm512_mask_swizzle_pd( j3, 0x55, j3, _MM_SWIZ_REG_CDAB );

          a[k2][k1] = _mm512_mask_swizzle_pd( t0, 0xCC, t1, _MM_SWIZ_REG_BADC );
          b[k2][k1] = _mm512_mask_swizzle_pd( t2, 0xCC, t3, _MM_SWIZ_REG_BADC );
          c[k2][k1] = _mm512_mask_swizzle_pd( t1, 0x33, t0, _MM_SWIZ_REG_BADC );
          d[k2][k1] = _mm512_mask_swizzle_pd( t3, 0x33, t2, _MM_SWIZ_REG_BADC );
        }
      }

      // Accumulate electric current in k3 plane
      ACC_CURRENT3D(  0, a, _MM_PERM_BABA )
      ACC_CURRENT3D(  1, b, _MM_PERM_BABA )
      ACC_CURRENT3D(  2, c, _MM_PERM_BABA )
      ACC_CURRENT3D(  3, d, _MM_PERM_BABA )
      ACC_CURRENT3D(  4, a, _MM_PERM_DCDC )
      ACC_CURRENT3D(  5, b, _MM_PERM_DCDC )
      ACC_CURRENT3D(  6, c, _MM_PERM_DCDC )
      ACC_CURRENT3D(  7, d, _MM_PERM_DCDC )

    }

  }
}


#undef ACC_CURRENT3D

#endif

#if 0

/**
 * mk. II
 *
 * Pre calculate all currents, deposit 1 particle at a time
 *
 * This is actually slower than the reference version and should not be used
 *
 */

#warning Using mk. II 3D algorithm

#if ( ORDER == 1 )

#define ACC_CURRENT3D( k, v, p ) { \
  double * const line = pj + idx[k]; \
  for( int k3 = 0; k3 < 2; k3++) { \
    for( int k2 = 0; k2 < 2; k2++) { \
      __m512d va, ja; \
      va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy + k3*Dz ] ); \
      ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k3][k2][0]), p )); \
      ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k3][k2][1]), p )); \
      va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
      _mm512_mask_compressstoreu_pd( &line[ k2*Dy + k3*Dz ], 0x77, va); \
    } \
  } \
}

#elif ( ORDER == 2 )

#define ACC_CURRENT3D( k, v, p ) { \
  double * const line = pj + idx[k]; \
  for( int k3 = 0; k3 < 3; k3++) { \
    for( int k2 = 0; k2 < 3; k2++) { \
      __m512d va, vb, ja, jb; \
      va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy + k3*Dz     ] ); \
      vb = _mm512_loadu_pd(                      &line[ k2*Dy + k3*Dz + 6 ] ); \
      ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k3][k2][0]), p )); \
      ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k3][k2][1]), p )); \
      jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[k3][k2][2]), p )); \
      va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
      vb = _mm512_mask_add_pd( vb, 0x07 , jb, vb ); \
      _mm512_mask_compressstoreu_pd( &line[ k2*Dy + k3*Dz     ], 0x77, va); \
      _mm512_mask_storeu_pd(         &line[ k2*Dy + k3*Dz + 6 ], 0x07, vb) ; \
    } \
  } \
}

#elif ( ORDER == 3 )

#define ACC_CURRENT3D( k, v, p ) { \
  double * const line = pj + idx[k]; \
  for( int k3 = 0; k3 < 4; k3++) { \
    for( int k2 = 0; k2 < 4; k2++) { \
      __m512d va, vb, ja, jb; \
      va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy + k3*Dz     ] ); \
      vb = _mm512_mask_expandloadu_pd( vb, 0x77, &line[ k2*Dy + k3*Dz + 6 ] ); \
      ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k3][k2][0]), p )); \
      ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k3][k2][1]), p )); \
      jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[k3][k2][2]), p )); \
      jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x3F00, _mm512_castpd_ps(v[k3][k2][3]), p )); \
      va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
      vb = _mm512_mask_add_pd( vb, 0x77 , jb, vb ); \
      _mm512_mask_compressstoreu_pd( &line[ k2*Dy + k3*Dz     ], 0x77, va); \
      _mm512_mask_compressstoreu_pd( &line[ k2*Dy + k3*Dz + 6 ], 0x77, vb); \
    } \
  } \
}

#elif ( ORDER == 4 )

#define ACC_CURRENT3D( k, v, p ) { \
  double * const line = pj + idx[k]; \
  for( int k3 = 0; k3 < 5; k3++) { \
    for( int k2 = 0; k2 < 5; k2++) { \
      __m512d va, vb, vc, ja, jb, jc; \
      va = _mm512_mask_expandloadu_pd( va, 0x77, &line[ k2*Dy + k3*Dz      ] ); \
      vb = _mm512_mask_expandloadu_pd( vb, 0x77, &line[ k2*Dy + k3*Dz +  6 ] ); \
      vc = _mm512_loadu_pd(                      &line[ k2*Dy + k3*Dz + 12 ] ); \
      ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x003F, _mm512_castpd_ps(v[k3][k2][0]), p )); \
      ja = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(ja), 0x3F00, _mm512_castpd_ps(v[k3][k2][1]), p )); \
      jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x003F, _mm512_castpd_ps(v[k3][k2][2]), p )); \
      jb = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jb), 0x3F00, _mm512_castpd_ps(v[k3][k2][3]), p )); \
      jc = _mm512_castps_pd(_mm512_mask_permute4f128_ps( _mm512_castpd_ps(jc), 0x003F, _mm512_castpd_ps(v[k3][k2][4]), p )); \
      va = _mm512_mask_add_pd( va, 0x77 , ja, va ); \
      vb = _mm512_mask_add_pd( vb, 0x77 , jb, vb ); \
      vc = _mm512_mask_add_pd( vc, 0x07 , jc, vc ); \
      _mm512_mask_compressstoreu_pd( &line[ k2*Dy + k3*Dz      ],  0x77, va); \
      _mm512_mask_compressstoreu_pd( &line[ k2*Dy + k3*Dz +  6 ],  0x77, vb); \
      _mm512_mask_storeu_pd(         &line[ k2*Dy + k3*Dz + 12 ],  0x07, vc); \
    } \
  } \
}

#else
#error Unsupported interpolation order
#endif


void DEP_CURRENT_3D
(double * const  current, int const * const  size, int const * const  offset,
                       double * const norm, t_split_buf3D * const part)
{

  __m512d vx0, vx1, vy0, vy1,  vz0, vz1, vq;

  __m512d vqnx, vqny, vqnz;
  __m512d vs0x[NP], vs1x[NP], vs0y[NP], vs1y[NP], vs0z[NP], vs1z[NP];

  __m512d vwl1[NP], vwl2[NP], vwl3[NP];
  __m512d vwp1[NP][NP], vwp2[NP][NP], vwp3[NP][NP];

  DECLARE_ALIGNED_64( int idx[VEC_WIDTH] );

  int const Dy =  3 * size[0];
  int const Dz = Dy * size[1];

  __m256i const vDy = _mm256_set1_epi32( Dy );
  __m256i const vDz = _mm256_set1_epi32( Dz );

  double* const pj = current + ( offset[0] - OFFSET  ) * 3 +
                               ( offset[1] - OFFSET  ) * Dy +
                               ( offset[2] - OFFSET  ) * Dz;

  __m512d const vnorm1 = _mm512_set1_pd(norm[0]);
  __m512d const vnorm2 = _mm512_set1_pd(norm[1]);
  __m512d const vnorm3 = _mm512_set1_pd(norm[2]);



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
    // idx = 3 * ix + Dy * iy + Dz * iz
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

    vqnx = _mm512_mul_pd( vq, vnorm1);
    vqny = _mm512_mul_pd( vq, vnorm2);
    vqnz = _mm512_mul_pd( vq, vnorm3);

    // get longitudinal weights
    WL( vqnx, vx0, vx1, vwl1 ); vwl1[NP-1] = _mm512_setzero_pd();
    WL( vqny, vy0, vy1, vwl2 ); vwl2[NP-1] = _mm512_setzero_pd();
    WL( vqnz, vz0, vz1, vwl3 ); vwl3[NP-1] = _mm512_setzero_pd();

    // get perpendicular weights
    for( int k2 = 0; k2 < NP; k2++ ) {
      for( int k1 = 0; k1 < NP; k1++ ) {

        const __m512d c1_2 = _mm512_set1_pd(0.5);
        __m512d tmp1, tmp2;

        // wp1[k2][k1] = s0y[k1]*s0z[k2] + s1y[k1]*s1z[k2] +
        //               0.5*( s0y[k1]*s1z[k2] + s1y[k1]*s0z[k2] )

        tmp1 = _mm512_mul_pd( vs0y[k1], vs0z[k2] );
        tmp2 = _mm512_mul_pd( vs0y[k1], vs1z[k2] );
        tmp1 = _mm512_fmadd_pd( vs1y[k1], vs1z[k2], tmp1 );
        tmp2 = _mm512_fmadd_pd( vs1y[k1], vs0z[k2], tmp2 );

        vwp1[k2][k1] = _mm512_fmadd_pd( c1_2, tmp2, tmp1 );

        // wp2[k2][k1] = s0x[k1]*s0z[k2] + s1x[k1]*s1z[k2] +
        //               0.5*( s0x[k1]*s1z[k2] + s1x[k1]*s0z[k2] )

        tmp1 = _mm512_mul_pd( vs0x[k1], vs0z[k2] );
        tmp2 = _mm512_mul_pd( vs0x[k1], vs1z[k2] );
        tmp1 = _mm512_fmadd_pd( vs1x[k1], vs1z[k2], tmp1 );
        tmp2 = _mm512_fmadd_pd( vs1x[k1], vs0z[k2], tmp2 );

        vwp2[k2][k1] = _mm512_fmadd_pd( c1_2, tmp2, tmp1 );

        // wp3[k2][k1] = s0x[k1]*s0y[k2] + s1x[k1]*s1y[k2] +
        //               0.5*( s0x[k1]*s1y[k2] + s1x[k1]*s0y[k2] )

        tmp1 = _mm512_mul_pd( vs0x[k1], vs0y[k2] );
        tmp2 = _mm512_mul_pd( vs0x[k1], vs1y[k2] );
        tmp1 = _mm512_fmadd_pd( vs1x[k1], vs1y[k2], tmp1 );
        tmp2 = _mm512_fmadd_pd( vs1x[k1], vs0y[k2], tmp2 );

        vwp3[k2][k1] = _mm512_fmadd_pd( c1_2, tmp2, tmp1 );

      }
    }

    __m512d a[NP][NP][NP], b[NP][NP][NP], c[NP][NP][NP], d[NP][NP][NP];
    for( int k3 = 0; k3 < NP; k3++ ) {
      for( int k2 = 0; k2 < NP; k2++ ) {
        for ( int k1 = 0; k1 < NP; k1++ ) {
          __m512d j1, j2, j3;

          // Calculate current components
          j1 = _mm512_mul_pd( vwl1[k1] , vwp1[k3][k2] );
          j2 = _mm512_mul_pd( vwl2[k2] , vwp2[k3][k1] );
          j3 = _mm512_mul_pd( vwl3[k3] , vwp3[k2][k1] );

          // Do a 4x4 transpose within 256 bit lanes ignoring 4th current component
          __m512d t0 = _mm512_mask_swizzle_pd( j1, 0xAA, j2, _MM_SWIZ_REG_CDAB );
          __m512d t1 = _mm512_mask_swizzle_pd( j3, 0xAA, j3, _MM_SWIZ_REG_CDAB );
          __m512d t2 = _mm512_mask_swizzle_pd( j2, 0x55, j1, _MM_SWIZ_REG_CDAB );
          __m512d t3 = _mm512_mask_swizzle_pd( j3, 0x55, j3, _MM_SWIZ_REG_CDAB );

          a[k3][k2][k1] = _mm512_mask_swizzle_pd( t0, 0xCC, t1, _MM_SWIZ_REG_BADC );
          b[k3][k2][k1] = _mm512_mask_swizzle_pd( t2, 0xCC, t3, _MM_SWIZ_REG_BADC );
          c[k3][k2][k1] = _mm512_mask_swizzle_pd( t1, 0x33, t0, _MM_SWIZ_REG_BADC );
          d[k3][k2][k1] = _mm512_mask_swizzle_pd( t3, 0x33, t2, _MM_SWIZ_REG_BADC );
        }
      }
    }

    // Accumulate electric current in k3 plane
    ACC_CURRENT3D(  0, a, _MM_PERM_BABA )
    ACC_CURRENT3D(  1, b, _MM_PERM_BABA )
    ACC_CURRENT3D(  2, c, _MM_PERM_BABA )
    ACC_CURRENT3D(  3, d, _MM_PERM_BABA )
    ACC_CURRENT3D(  4, a, _MM_PERM_DCDC )
    ACC_CURRENT3D(  5, b, _MM_PERM_DCDC )
    ACC_CURRENT3D(  6, c, _MM_PERM_DCDC )
    ACC_CURRENT3D(  7, d, _MM_PERM_DCDC )


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








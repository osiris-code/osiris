/*

 Spline function definitions for x86 AVX extensions

*/

#ifndef _SPLINES_H
#define _SPLINES_H


#include "vector-avx.h"

#if defined( PRECISION_SINGLE )

/*
 Single precision variants
 */

/********************************* Linear Interpolation ********************************** 
vspline_s1

1st order (linear) splines.

w[0] = 1/2 - x
w[1] = 1/2 + x
*****************************************************************************************/

static inline void vspline_s1( __m256 const dx, __m256 w[] ) {
  
  __m256 const c1_2 = _mm256_set1_ps( 0.5f );
  
  w[0] = _mm256_sub_ps( c1_2, dx );
  w[1] = _mm256_add_ps( c1_2, dx );
}

static inline void vsplineh_s1( __m256 const dx, __m256 const h, __m256 w[] ) {
  
  __m256 const c1   = _mm256_set1_ps( 1.0f );
  
  w[0] = _mm256_sub_ps( _mm256_sub_ps( c1, h ), dx );
  w[1] = _mm256_add_ps(                    h  , dx );
}

/******************************** Quadratic Interpolation *******************************/ 

static inline void vspline_s2( __m256 const dx, __m256 w[] ) {
  
  __m256 const c1_2 = _mm256_set1_ps( 0.5f );

  __m256 t0 = _mm256_sub_ps( c1_2, dx );
  __m256 t1 = _mm256_add_ps( c1_2, dx );

  w[0] = _mm256_mul_ps( c1_2, _mm256_mul_ps( t0, t0 ) );
  w[1] = _mm256_add_ps( c1_2, _mm256_mul_ps( t0, t1 ) );
  w[2] = _mm256_mul_ps( c1_2, _mm256_mul_ps( t1, t1 ) );  
}

static inline void vsplineh_s2( __m256 const dx, __m256 const h, __m256 w[] ) {
  
  __m256 const c1   = _mm256_set1_ps( 1.0f );
  __m256 const c1_2 = _mm256_set1_ps( 0.5f );

  __m256 t0 = _mm256_sub_ps( _mm256_sub_ps( c1, h ), dx );
  __m256 t1 = _mm256_add_ps(                    h  , dx );

  w[0] = _mm256_mul_ps( c1_2, _mm256_mul_ps( t0, t0 ) );
  w[1] = _mm256_add_ps( c1_2, _mm256_mul_ps( t0, t1 ) );
  w[2] = _mm256_mul_ps( c1_2, _mm256_mul_ps( t1, t1 ) );  
}


/********************************** Cubic Interpolation *********************************/ 

static inline void vspline_s3( __m256 const dx, __m256 w[] ) {
  
  __m256 const c1_2 = _mm256_set1_ps( 1.0f/2.0f );
  __m256 const c1_6 = _mm256_set1_ps( 1.0f/6.0f );
  __m256 const c2_3 = _mm256_set1_ps( 2.0f/3.0f );

  __m256 t0 = _mm256_sub_ps( c1_2, dx );
  __m256 t1 = _mm256_add_ps( c1_2, dx );
  
  __m256 t2 = _mm256_mul_ps( t0, t0 );  // t2 = (1/2 - dx)^2
  __m256 t3 = _mm256_mul_ps( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm256_mul_ps( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm256_mul_ps( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm256_mul_ps( c1_6, t0 );
  w[1] = _mm256_add_ps( _mm256_sub_ps( c2_3, t3 ), _mm256_mul_ps( c1_2, t1 ) );
  w[2] = _mm256_add_ps( _mm256_sub_ps( c2_3, t2 ), _mm256_mul_ps( c1_2, t0 ) );
  w[3] = _mm256_mul_ps( c1_6, t1 ); 

}

static inline void vsplineh_s3( __m256 const dx, __m256 const h, __m256 w[] ) {
  
  __m256 const c1   = _mm256_set1_ps( 1.0f );
  __m256 const c1_2 = _mm256_set1_ps( 1.0f/2.0f );
  __m256 const c1_6 = _mm256_set1_ps( 1.0f/6.0f );
  __m256 const c2_3 = _mm256_set1_ps( 2.0f/3.0f );

  __m256 t0 = _mm256_sub_ps( _mm256_sub_ps( c1, h ), dx );
  __m256 t1 = _mm256_add_ps(                    h  , dx );
  
  __m256 t2 = _mm256_mul_ps( t0, t0 );  // t2 = (1/2 - dx)^2
  __m256 t3 = _mm256_mul_ps( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm256_mul_ps( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm256_mul_ps( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm256_mul_ps( c1_6, t0 );
  w[1] = _mm256_add_ps( _mm256_sub_ps( c2_3, t3 ), _mm256_mul_ps( c1_2, t1 ) );
  w[2] = _mm256_add_ps( _mm256_sub_ps( c2_3, t2 ), _mm256_mul_ps( c1_2, t0 ) );
  w[3] = _mm256_mul_ps( c1_6, t1 ); 

}


/********************************* Quartic Interpolation ********************************/ 

static inline void vspline_s4( __m256 const dx, __m256 w[] ) {
  
  __m256 const c1_2  = _mm256_set1_ps( 1.0f/2.0f );
  __m256 const c1_24 = _mm256_set1_ps( 1.0f/24.0f );
  __m256 const c1_4  = _mm256_set1_ps( 1.0f/4.0f );
  __m256 const c1    = _mm256_set1_ps( 1.0f );
  __m256 const c3_2  = _mm256_set1_ps( 3.0f/2.0f );
  __m256 const c1_6  = _mm256_set1_ps( 1.0f/6.0f );
  __m256 const c11_24 = _mm256_set1_ps( 11.0f/24.0f );
  
  __m256 t0 = _mm256_sub_ps( c1_2, dx );
  __m256 t1 = _mm256_add_ps( c1_2, dx );
  __m256 t2 = _mm256_mul_ps( t0, t1 );

  __m256 t02 = _mm256_mul_ps( t0, t0 );
  __m256 t12 = _mm256_mul_ps( t1, t1 );

  
  w[0] = _mm256_mul_ps( c1_24, _mm256_mul_ps( t02, t02 ) );
  w[1] = _mm256_mul_ps( c1_6, _mm256_add_ps( c1_4, _mm256_mul_ps( t0, _mm256_add_ps( c1, _mm256_mul_ps( t0, _mm256_add_ps( c3_2, t2 ) ) ) ) ) );
  w[2] = _mm256_add_ps( c11_24, _mm256_mul_ps( t2, _mm256_add_ps( c1_2, _mm256_mul_ps( c1_4, t2 ) ) ) );
  w[3] = _mm256_mul_ps( c1_6, _mm256_add_ps( c1_4, _mm256_mul_ps( t1, _mm256_add_ps( c1, _mm256_mul_ps( t1, _mm256_add_ps( c3_2, t2 ) ) ) ) ) );
  w[4] = _mm256_mul_ps( c1_24, _mm256_mul_ps( t12, t12 ) );
  
}

static inline void vsplineh_s4( __m256 const dx, __m256 const h, __m256 w[] ) {
  
  __m256 const c1_2  = _mm256_set1_ps( 1.0f/2.0f );
  __m256 const c1_24 = _mm256_set1_ps( 1.0f/24.0f );
  __m256 const c1_4  = _mm256_set1_ps( 1.0f/4.0f );
  __m256 const c1    = _mm256_set1_ps( 1.0f );
  __m256 const c3_2  = _mm256_set1_ps( 3.0f/2.0f );
  __m256 const c1_6  = _mm256_set1_ps( 1.0f/6.0f );
  __m256 const c11_24 = _mm256_set1_ps( 11.0f/24.0f );

  __m256 t0 = _mm256_sub_ps( _mm256_sub_ps( c1, h ), dx );
  __m256 t1 = _mm256_add_ps(                    h  , dx );
  __m256 t2 = _mm256_mul_ps( t0, t1 );

  __m256 t02 = _mm256_mul_ps( t0, t0 );
  __m256 t12 = _mm256_mul_ps( t1, t1 );
  
  w[0] = _mm256_mul_ps( c1_24, _mm256_mul_ps( t02, t02 ) );
  w[1] = _mm256_mul_ps( c1_6, _mm256_add_ps( c1_4, _mm256_mul_ps( t0, _mm256_add_ps( c1, _mm256_mul_ps( t0, _mm256_add_ps( c3_2, t2 ) ) ) ) ) );
  w[2] = _mm256_add_ps( c11_24, _mm256_mul_ps( t2, _mm256_add_ps( c1_2, _mm256_mul_ps( c1_4, t2 ) ) ) );
  w[3] = _mm256_mul_ps( c1_6, _mm256_add_ps( c1_4, _mm256_mul_ps( t1, _mm256_add_ps( c1, _mm256_mul_ps( t1, _mm256_add_ps( c3_2, t2 ) ) ) ) ) );
  w[4] = _mm256_mul_ps( c1_24, _mm256_mul_ps( t12, t12 ) );
  
}

#elif defined( PRECISION_DOUBLE )

/*
 Double precision variants
 */

/********************************* Linear Interpolation ********************************** 
vspline_s1

1st order (linear) splines.

w[0] = 1/2 - x
w[1] = 1/2 + x
*****************************************************************************************/

static inline void vsplined_s1( __m256d const dx, __m256d w[] ) {
  
  __m256d const c1_2 = _mm256_set1_pd( 0.5 );
  
  w[0] = _mm256_sub_pd( c1_2, dx );
  w[1] = _mm256_add_pd( c1_2, dx );
}

static inline void vsplinehd_s1( __m256d const dx, __m256d const h, __m256d w[] ) {

  __m256d const c1   = _mm256_set1_pd( 1.0 );
  
  w[0] = _mm256_sub_pd( _mm256_sub_pd( c1, h ), dx );
  w[1] = _mm256_add_pd(                    h  , dx );
}

/******************************** Quadratic Interpolation ********************************/

static inline void vsplined_s2( __m256d const dx, __m256d w[] ) {
  
  __m256d const c1_2 = _mm256_set1_pd( 0.5 );

  __m256d t0 = _mm256_sub_pd( c1_2, dx );
  __m256d t1 = _mm256_add_pd( c1_2, dx );

  w[0] = _mm256_mul_pd( c1_2, _mm256_mul_pd( t0, t0 ) );
  w[1] = _mm256_add_pd( c1_2, _mm256_mul_pd( t0, t1 ) );
  w[2] = _mm256_mul_pd( c1_2, _mm256_mul_pd( t1, t1 ) );  
}

static inline void vsplinehd_s2( __m256d const dx, __m256d const h, __m256d w[] ) {
  
  __m256d const c1   = _mm256_set1_pd( 1.0 );
  __m256d const c1_2 = _mm256_set1_pd( 0.5 );

  __m256d t0 = _mm256_sub_pd( _mm256_sub_pd( c1, h ), dx );
  __m256d t1 = _mm256_add_pd(                    h  , dx );

  w[0] = _mm256_mul_pd( c1_2, _mm256_mul_pd( t0, t0 ) );
  w[1] = _mm256_add_pd( c1_2, _mm256_mul_pd( t0, t1 ) );
  w[2] = _mm256_mul_pd( c1_2, _mm256_mul_pd( t1, t1 ) );  
}


/********************************** Cubic Interpolation *********************************/


static inline void vsplined_s3( __m256d const dx, __m256d w[] ) {
  
  __m256d const c1_2 = _mm256_set1_pd( 1.0/2.0 );
  __m256d const c1_6 = _mm256_set1_pd( 1.0/6.0 );
  __m256d const c2_3 = _mm256_set1_pd( 2.0/3.0 );

  __m256d t0 = _mm256_sub_pd( c1_2, dx );
  __m256d t1 = _mm256_add_pd( c1_2, dx );
  
  __m256d t2 = _mm256_mul_pd( t0, t0 );  // t2 = (1/2 - dx)^2
  __m256d t3 = _mm256_mul_pd( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm256_mul_pd( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm256_mul_pd( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm256_mul_pd( c1_6, t0 );
  w[1] = _mm256_add_pd( _mm256_sub_pd( c2_3, t3 ), _mm256_mul_pd( c1_2, t1 ) );
  w[2] = _mm256_add_pd( _mm256_sub_pd( c2_3, t2 ), _mm256_mul_pd( c1_2, t0 ) );
  w[3] = _mm256_mul_pd( c1_6, t1 ); 

}

static inline void vsplinehd_s3( __m256d const dx, __m256d const h, __m256d w[] ) {
  
  __m256d const c1   = _mm256_set1_pd( 1.0 );
  __m256d const c1_2 = _mm256_set1_pd( 1.0/2.0 );
  __m256d const c1_6 = _mm256_set1_pd( 1.0/6.0 );
  __m256d const c2_3 = _mm256_set1_pd( 2.0/3.0 );

  __m256d t0 = _mm256_sub_pd( _mm256_sub_pd( c1, h ), dx );
  __m256d t1 = _mm256_add_pd(                    h  , dx );
  
  __m256d t2 = _mm256_mul_pd( t0, t0 );  // t2 = (1/2 - dx)^2
  __m256d t3 = _mm256_mul_pd( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm256_mul_pd( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm256_mul_pd( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm256_mul_pd( c1_6, t0 );
  w[1] = _mm256_add_pd( _mm256_sub_pd( c2_3, t3 ), _mm256_mul_pd( c1_2, t1 ) );
  w[2] = _mm256_add_pd( _mm256_sub_pd( c2_3, t2 ), _mm256_mul_pd( c1_2, t0 ) );
  w[3] = _mm256_mul_pd( c1_6, t1 ); 

}


/********************************* Quartic Interpolation ********************************/

static inline void vsplined_s4( __m256d const dx, __m256d w[] ) {
  
  __m256d const c1_2  = _mm256_set1_pd( 1.0/2.0 );
  __m256d const c1_24 = _mm256_set1_pd( 1.0/24.0 );
  __m256d const c1_4  = _mm256_set1_pd( 1.0/4.0 );
  __m256d const c1    = _mm256_set1_pd( 1.0 );
  __m256d const c3_2  = _mm256_set1_pd( 3.0/2.0 );
  __m256d const c1_6  = _mm256_set1_pd( 1.0/6.0 );
  __m256d const c11_24 = _mm256_set1_pd( 11.0/24.0 );
  
  __m256d t0 = _mm256_sub_pd( c1_2, dx );
  __m256d t1 = _mm256_add_pd( c1_2, dx );
  __m256d t2 = _mm256_mul_pd( t0, t1 );

  __m256d t02 = _mm256_mul_pd( t0, t0 );
  __m256d t12 = _mm256_mul_pd( t1, t1 );

  
  w[0] = _mm256_mul_pd( c1_24, _mm256_mul_pd( t02, t02 ) );
  w[1] = _mm256_mul_pd( c1_6, _mm256_add_pd( c1_4, _mm256_mul_pd( t0, _mm256_add_pd( c1, _mm256_mul_pd( t0, _mm256_add_pd( c3_2, t2 ) ) ) ) ) );
  w[2] = _mm256_add_pd( c11_24, _mm256_mul_pd( t2, _mm256_add_pd( c1_2, _mm256_mul_pd( c1_4, t2 ) ) ) );
  w[3] = _mm256_mul_pd( c1_6, _mm256_add_pd( c1_4, _mm256_mul_pd( t1, _mm256_add_pd( c1, _mm256_mul_pd( t1, _mm256_add_pd( c3_2, t2 ) ) ) ) ) );
  w[4] = _mm256_mul_pd( c1_24, _mm256_mul_pd( t12, t12 ) );
  
}

static inline void vsplinehd_s4( __m256d const dx, __m256d const h, __m256d w[] ) {
  
  __m256d const c1_2  = _mm256_set1_pd( 1.0/2.0 );
  __m256d const c1_24 = _mm256_set1_pd( 1.0/24.0 );
  __m256d const c1_4  = _mm256_set1_pd( 1.0/4.0 );
  __m256d const c1    = _mm256_set1_pd( 1.0 );
  __m256d const c3_2  = _mm256_set1_pd( 3.0/2.0 );
  __m256d const c1_6  = _mm256_set1_pd( 1.0/6.0 );
  __m256d const c11_24 = _mm256_set1_pd( 11.0/24.0 );

  __m256d t0 = _mm256_sub_pd( _mm256_sub_pd( c1, h ), dx );
  __m256d t1 = _mm256_add_pd(                    h  , dx );
  __m256d t2 = _mm256_mul_pd( t0, t1 );

  __m256d t02 = _mm256_mul_pd( t0, t0 );
  __m256d t12 = _mm256_mul_pd( t1, t1 );
  
  w[0] = _mm256_mul_pd( c1_24, _mm256_mul_pd( t02, t02 ) );
  w[1] = _mm256_mul_pd( c1_6, _mm256_add_pd( c1_4, _mm256_mul_pd( t0, _mm256_add_pd( c1, _mm256_mul_pd( t0, _mm256_add_pd( c3_2, t2 ) ) ) ) ) );
  w[2] = _mm256_add_pd( c11_24, _mm256_mul_pd( t2, _mm256_add_pd( c1_2, _mm256_mul_pd( c1_4, t2 ) ) ) );
  w[3] = _mm256_mul_pd( c1_6, _mm256_add_pd( c1_4, _mm256_mul_pd( t1, _mm256_add_pd( c1, _mm256_mul_pd( t1, _mm256_add_pd( c3_2, t2 ) ) ) ) ) );
  w[4] = _mm256_mul_pd( c1_24, _mm256_mul_pd( t12, t12 ) );
  
}

#else
#error PRECISION_SINGLE or PRECISION_DOUBLE must be defined
#endif 


#endif /* _SPLINES_H */

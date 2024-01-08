/*

 Spline function definitions for Intel MIC 

*/

#ifndef _SPLINES_H
#define _SPLINES_H

#include "vector-mic.h"

/********************************* Linear Interpolation ********************************** 
vspline_s1

1st order (linear) splines.

w[0] = 1/2 - x
w[1] = 1/2 + x
*****************************************************************************************/

static inline void vspline_s1( __m512 const dx, __m512 w[] ) {
  
  __m512 const c1_2 = _mm512_set1_ps( 0.5f );
  
  w[0] = _mm512_sub_ps( c1_2, dx );
  w[1] = _mm512_add_ps( c1_2, dx );
}

static inline void vsplineh_s1( __m512 const dx, __m512 const h, __m512 w[] ) {
  
  __m512 const c1 = _mm512_set1_ps( 1.0f );
  
  w[0] = _mm512_sub_ps( _mm512_sub_ps( c1, h ), dx );
  w[1] = _mm512_add_ps(                    h  , dx );
}

/******************************** Quadratic Interpolation ********************************/

static inline void vspline_s2( __m512 const dx, __m512 w[] ) {
  
  __m512 const c1_2 = _mm512_set1_ps( 0.5f );

  __m512 t0 = _mm512_sub_ps( c1_2, dx );
  __m512 t1 = _mm512_add_ps( c1_2, dx );

  w[0] = _mm512_mul_ps( _mm512_mul_ps( t0, t0 ), c1_2 );
  w[1] = _mm512_fmadd_ps(              t0, t1,   c1_2 );
  w[2] = _mm512_mul_ps( _mm512_mul_ps( t1, t1 ), c1_2 );  
}

static inline void vsplineh_s2( __m512 const dx, __m512 const h, __m512 w[] ) {
  
  __m512 const c1   = _mm512_set1_ps( 1.0f );
  __m512 const c1_2 = _mm512_set1_ps( 0.5f );

  __m512 t0 = _mm512_sub_ps( _mm512_sub_ps( c1, h ), dx );
  __m512 t1 = _mm512_add_ps(                    h  , dx );

  w[0] = _mm512_mul_ps( _mm512_mul_ps( t0, t0 ), c1_2 );
  w[1] = _mm512_fmadd_ps(              t0, t1,   c1_2 );
  w[2] = _mm512_mul_ps( _mm512_mul_ps( t1, t1 ), c1_2 );  
}


/********************************** Cubic Interpolation *********************************/


static inline void vspline_s3( __m512 const dx, __m512 w[] ) {
  
  __m512 const c1_2 = _mm512_set1_ps( 1.0f/2.0f );
  __m512 const c1_6 = _mm512_set1_ps( 1.0f/6.0f );
  __m512 const c2_3 = _mm512_set1_ps( 2.0f/3.0f );

  __m512 t0 = _mm512_sub_ps( c1_2, dx );
  __m512 t1 = _mm512_add_ps( c1_2, dx );
  
  __m512 t2 = _mm512_mul_ps( t0, t0 );  // t2 = (1/2 - dx)^2
  __m512 t3 = _mm512_mul_ps( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm512_mul_ps( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm512_mul_ps( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm512_mul_ps( c1_6, t0 );
  w[1] = _mm512_fmadd_ps( c1_2, t1, _mm512_sub_ps( c2_3, t3 ) );
  w[2] = _mm512_fmadd_ps( c1_2, t0, _mm512_sub_ps( c2_3, t2 ) );
  w[3] = _mm512_mul_ps( c1_6, t1 ); 


}

static inline void vsplineh_s3( __m512 const dx, __m512 const h, __m512 w[] ) {
  
  __m512 const c1   = _mm512_set1_ps( 1.0f );
  __m512 const c1_2 = _mm512_set1_ps( 1.0f/2.0f );
  __m512 const c1_6 = _mm512_set1_ps( 1.0f/6.0f );
  __m512 const c2_3 = _mm512_set1_ps( 2.0f/3.0f );

  __m512 t0 = _mm512_sub_ps( _mm512_sub_ps( c1, h ), dx );
  __m512 t1 = _mm512_add_ps(                    h  , dx );
  
  __m512 t2 = _mm512_mul_ps( t0, t0 );  // t2 = (1/2 - dx)^2
  __m512 t3 = _mm512_mul_ps( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm512_mul_ps( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm512_mul_ps( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm512_mul_ps( c1_6, t0 );
  w[1] = _mm512_fmadd_ps( c1_2, t1, _mm512_sub_ps( c2_3, t3 ) );
  w[2] = _mm512_fmadd_ps( c1_2, t0, _mm512_sub_ps( c2_3, t2 ) );
  w[3] = _mm512_mul_ps( c1_6, t1 ); 

}


/********************************* Quartic Interpolation ********************************/

static inline void vspline_s4( __m512 const dx, __m512 w[] ) {
/*
  t0 = 0.5 - x
  t1 = 0.5 + x
  t2 = t0 * t1

  s(-2) = t0**4/24.
  s(-1) = ( 0.25 + t0 * ( 1 + t0 * ( 1.5 + t2 ) ) ) / 6.
  s( 0) = 11./24. + t2 * ( 0.5 + 0.25 * t2 )
  s( 1) = ( 0.25 + t1 * ( 1 + t1 * ( 1.5 + t2 ) ) ) / 6.
  s( 2) = t1**4/24.
*/
  
  __m512 const c1_2  = _mm512_set1_ps( 1.0f/2.0f );
  __m512 const c1_24 = _mm512_set1_ps( 1.0f/24.0f );
  __m512 const c1_4  = _mm512_set1_ps( 1.0f/4.0f );
  __m512 const c1    = _mm512_set1_ps( 1.0f );
  __m512 const c3_2  = _mm512_set1_ps( 3.0f/2.0f );
  __m512 const c1_6  = _mm512_set1_ps( 1.0f/6.0f );
  __m512 const c11_24 = _mm512_set1_ps( 11.0f/24.0f );
  
  __m512 t0 = _mm512_sub_ps( c1_2, dx );
  __m512 t1 = _mm512_add_ps( c1_2, dx );
  __m512 t2 = _mm512_mul_ps( t0, t1 );

  __m512 t02 = _mm512_mul_ps( t0, t0 );
  __m512 t12 = _mm512_mul_ps( t1, t1 );
  
  w[0] = _mm512_mul_ps( c1_24, _mm512_mul_ps( t02, t02 ) );
  w[1] = _mm512_mul_ps( c1_6, _mm512_fmadd_ps( t0, _mm512_fmadd_ps( t0, _mm512_add_ps( c3_2, t2 ), c1), c1_4 ) );
  w[2] = _mm512_fmadd_ps( t2, _mm512_fmadd_ps( c1_4, t2, c1_2 ), c11_24 );
  w[3] = _mm512_mul_ps( c1_6, _mm512_fmadd_ps( t1, _mm512_fmadd_ps( t1, _mm512_add_ps( c3_2, t2 ), c1), c1_4 ) );
  w[4] = _mm512_mul_ps( c1_24, _mm512_mul_ps( t12, t12 ) );
  
}

static inline void vsplineh_s4( __m512 const dx, __m512 const h, __m512 w[] ) {
/*
  t0 = 0.5 - x
  t1 = 0.5 + x
  t2 = t0 * t1

  s(-2) = t0**4/24.
  s(-1) = ( 0.25 + t0 * ( 1 + t0 * ( 1.5 + t2 ) ) ) / 6.
  s( 0) = 11./24. + t2 * ( 0.5 + 0.25 * t2 )
  s( 1) = ( 0.25 + t1 * ( 1 + t1 * ( 1.5 + t2 ) ) ) / 6.
  s( 2) = t1**4/24.
*/
  
  __m512 const c1_2  = _mm512_set1_ps( 1.0f/2.0f );
  __m512 const c1_24 = _mm512_set1_ps( 1.0f/24.0f );
  __m512 const c1_4  = _mm512_set1_ps( 1.0f/4.0f );
  __m512 const c1    = _mm512_set1_ps( 1.0f );
  __m512 const c3_2  = _mm512_set1_ps( 3.0f/2.0f );
  __m512 const c1_6  = _mm512_set1_ps( 1.0f/6.0f );
  __m512 const c11_24 = _mm512_set1_ps( 11.0f/24.0f );
  
  __m512 t0 = _mm512_sub_ps( _mm512_sub_ps( c1, h ), dx );
  __m512 t1 = _mm512_add_ps(                    h  , dx );
  __m512 t2 = _mm512_mul_ps( t0, t1 );

  __m512 t02 = _mm512_mul_ps( t0, t0 );
  __m512 t12 = _mm512_mul_ps( t1, t1 );
  
  w[0] = _mm512_mul_ps( c1_24, _mm512_mul_ps( t02, t02 ) );
  w[1] = _mm512_mul_ps( c1_6, _mm512_fmadd_ps( t0, _mm512_fmadd_ps( t0, _mm512_add_ps( c3_2, t2 ), c1), c1_4 ) );
  w[2] = _mm512_fmadd_ps( t2, _mm512_fmadd_ps( c1_4, t2, c1_2 ), c11_24 );
  w[3] = _mm512_mul_ps( c1_6, _mm512_fmadd_ps( t1, _mm512_fmadd_ps( t1, _mm512_add_ps( c3_2, t2 ), c1), c1_4 ) );
  w[4] = _mm512_mul_ps( c1_24, _mm512_mul_ps( t12, t12 ) );
  
}


#endif /* _SPLINES_H */

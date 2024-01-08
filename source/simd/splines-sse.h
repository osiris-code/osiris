/*

 Spline function definitions.

*/

#ifndef _SPLINES_H
#define _SPLINES_H


#include "vector-sse.h"

/********************************* Linear Interpolation *********************************/ 

static inline void vspline_s1( __m128 const dx, __m128 w[] ) {
  
  __m128 const c1_2 = _mm_set1_ps( 0.5f );
  
  w[0] = _mm_sub_ps( c1_2, dx );
  w[1] = _mm_add_ps( c1_2, dx );
}

static inline void vsplineh_s1( __m128 const dx, __m128 const h, __m128 w[] ) {
  
  __m128 const c1   = _mm_set1_ps( 1.0f );
  
  w[0] = _mm_sub_ps( _mm_sub_ps( c1, h ), dx );
  w[1] = _mm_add_ps(                 h  , dx );
}


static inline void vsplined_s1( __m128d const dx, __m128d w[] ) {

  __m128d const c1_2 = _mm_set1_pd( 0.5 );
  
  w[0] = _mm_sub_pd( c1_2, dx );
  w[1] = _mm_add_pd( c1_2, dx );
}

static inline void vsplinehd_s1( __m128d const dx, __m128d const h, __m128d w[] ) {

  __m128d const c1   = _mm_set1_pd( 1.0 );
  
  w[0] = _mm_sub_pd( _mm_sub_pd( c1, h ), dx );
  w[1] = _mm_add_pd(                 h  , dx );
}


/******************************** Quadratic Interpolation *******************************/ 

static inline void vspline_s2( __m128 const dx, __m128 w[] ) {
  
  __m128 const c1_2 = _mm_set1_ps( 0.5f );

  __m128 t0 = _mm_sub_ps( c1_2, dx );
  __m128 t1 = _mm_add_ps( c1_2, dx );

  w[0] = _mm_mul_ps( c1_2, _mm_mul_ps( t0, t0 ) );
  w[1] = _mm_add_ps( c1_2, _mm_mul_ps( t0, t1 ) );
  w[2] = _mm_mul_ps( c1_2, _mm_mul_ps( t1, t1 ) );  
}

static inline void vsplineh_s2( __m128 const dx, __m128 const h, __m128 w[] ) {
  
  __m128 const c1   = _mm_set1_ps( 1.0f );
  __m128 const c1_2 = _mm_set1_ps( 0.5f );

  __m128 t0 = _mm_sub_ps( _mm_sub_ps( c1, h ), dx );
  __m128 t1 = _mm_add_ps(                 h  , dx );

  w[0] = _mm_mul_ps( c1_2, _mm_mul_ps( t0, t0 ) );
  w[1] = _mm_add_ps( c1_2, _mm_mul_ps( t0, t1 ) );
  w[2] = _mm_mul_ps( c1_2, _mm_mul_ps( t1, t1 ) );  
}


static inline void vsplined_s2( __m128d const dx, __m128d w[] ) {
  
  __m128d const c1_2 = _mm_set1_pd( 0.5 );

  __m128d t0 = _mm_sub_pd( c1_2, dx );
  __m128d t1 = _mm_add_pd( c1_2, dx );

  w[0] = _mm_mul_pd( c1_2, _mm_mul_pd( t0, t0 ) );
  w[1] = _mm_add_pd( c1_2, _mm_mul_pd( t0, t1 ) );
  w[2] = _mm_mul_pd( c1_2, _mm_mul_pd( t1, t1 ) );  
}

static inline void vsplinehd_s2( __m128d const dx, __m128d const h, __m128d w[] ) {
  
  __m128d const c1   = _mm_set1_pd( 1.0 );
  __m128d const c1_2 = _mm_set1_pd( 0.5 );

  __m128d t0 = _mm_sub_pd( _mm_sub_pd( c1, h ), dx );
  __m128d t1 = _mm_add_pd(                 h  , dx );

  w[0] = _mm_mul_pd( c1_2, _mm_mul_pd( t0, t0 ) );
  w[1] = _mm_add_pd( c1_2, _mm_mul_pd( t0, t1 ) );
  w[2] = _mm_mul_pd( c1_2, _mm_mul_pd( t1, t1 ) );  
}


/********************************** Cubic Interpolation *********************************/ 

static inline void vspline_s3( __m128 const dx, __m128 w[] ) {
  
  __m128 const c1_2 = _mm_set1_ps( 1.0f/2.0f );
  __m128 const c1_6 = _mm_set1_ps( 1.0f/6.0f );
  __m128 const c2_3 = _mm_set1_ps( 2.0f/3.0f );

  __m128 t0 = _mm_sub_ps( c1_2, dx );
  __m128 t1 = _mm_add_ps( c1_2, dx );
  
  __m128 t2 = _mm_mul_ps( t0, t0 );  // t2 = (1/2 - dx)^2
  __m128 t3 = _mm_mul_ps( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm_mul_ps( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm_mul_ps( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm_mul_ps( c1_6, t0 );
  w[1] = _mm_add_ps( _mm_sub_ps( c2_3, t3 ), _mm_mul_ps( c1_2, t1 ) );
  w[2] = _mm_add_ps( _mm_sub_ps( c2_3, t2 ), _mm_mul_ps( c1_2, t0 ) );
  w[3] = _mm_mul_ps( c1_6, t1 ); 

}

static inline void vsplineh_s3( __m128 const dx, __m128 const h, __m128 w[] ) {
  
  __m128 const c1   = _mm_set1_ps( 1.0f );
  __m128 const c1_2 = _mm_set1_ps( 1.0f/2.0f );
  __m128 const c1_6 = _mm_set1_ps( 1.0f/6.0f );
  __m128 const c2_3 = _mm_set1_ps( 2.0f/3.0f );

  __m128 t0 = _mm_sub_ps( _mm_sub_ps( c1, h ), dx );
  __m128 t1 = _mm_add_ps(                 h  , dx );
  
  __m128 t2 = _mm_mul_ps( t0, t0 );  // t2 = (1/2 - dx)^2
  __m128 t3 = _mm_mul_ps( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm_mul_ps( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm_mul_ps( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm_mul_ps( c1_6, t0 );
  w[1] = _mm_add_ps( _mm_sub_ps( c2_3, t3 ), _mm_mul_ps( c1_2, t1 ) );
  w[2] = _mm_add_ps( _mm_sub_ps( c2_3, t2 ), _mm_mul_ps( c1_2, t0 ) );
  w[3] = _mm_mul_ps( c1_6, t1 ); 

}


static inline void vsplined_s3( __m128d const dx, __m128d w[] ) {
  
  __m128d const c1_2 = _mm_set1_pd( 1.0/2.0 );
  __m128d const c1_6 = _mm_set1_pd( 1.0/6.0 );
  __m128d const c2_3 = _mm_set1_pd( 2.0/3.0 );

  __m128d t0 = _mm_sub_pd( c1_2, dx );
  __m128d t1 = _mm_add_pd( c1_2, dx );
  
  __m128d t2 = _mm_mul_pd( t0, t0 );  // t2 = (1/2 - dx)^2
  __m128d t3 = _mm_mul_pd( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm_mul_pd( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm_mul_pd( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm_mul_pd( c1_6, t0 );
  w[1] = _mm_add_pd( _mm_sub_pd( c2_3, t3 ), _mm_mul_pd( c1_2, t1 ) );
  w[2] = _mm_add_pd( _mm_sub_pd( c2_3, t2 ), _mm_mul_pd( c1_2, t0 ) );
  w[3] = _mm_mul_pd( c1_6, t1 ); 

}

static inline void vsplinehd_s3( __m128d const dx, __m128d const h, __m128d w[] ) {
  
  __m128d const c1   = _mm_set1_pd( 1.0 );
  __m128d const c1_2 = _mm_set1_pd( 1.0/2.0 );
  __m128d const c1_6 = _mm_set1_pd( 1.0/6.0 );
  __m128d const c2_3 = _mm_set1_pd( 2.0/3.0 );

  __m128d t0 = _mm_sub_pd( _mm_sub_pd( c1, h ), dx );
  __m128d t1 = _mm_add_pd(                 h  , dx );
  
  __m128d t2 = _mm_mul_pd( t0, t0 );  // t2 = (1/2 - dx)^2
  __m128d t3 = _mm_mul_pd( t1, t1 );  // t3 = (1/2 + dx)^2
  
  t0 = _mm_mul_pd( t0, t2 );  // t0 = (1/2 - dx)^3
  t1 = _mm_mul_pd( t1, t3 );  // t1 = (1/2 + dx)^3 
  
  w[0] = _mm_mul_pd( c1_6, t0 );
  w[1] = _mm_add_pd( _mm_sub_pd( c2_3, t3 ), _mm_mul_pd( c1_2, t1 ) );
  w[2] = _mm_add_pd( _mm_sub_pd( c2_3, t2 ), _mm_mul_pd( c1_2, t0 ) );
  w[3] = _mm_mul_pd( c1_6, t1 ); 

}


/********************************* Quartic Interpolation ********************************/ 


static inline void vspline_s4( __m128 const dx, __m128 w[] ) {
  
  __m128 const c1_2  = _mm_set1_ps( 1.0f/2.0f );
  __m128 const c1_24 = _mm_set1_ps( 1.0f/24.0f );
  __m128 const c1_4  = _mm_set1_ps( 1.0f/4.0f );
  __m128 const c1    = _mm_set1_ps( 1.0f );
  __m128 const c3_2  = _mm_set1_ps( 3.0f/2.0f );
  __m128 const c1_6  = _mm_set1_ps( 1.0f/6.0f );
  __m128 const c11_24 = _mm_set1_ps( 11.0f/24.0f );
  
  __m128 t0 = _mm_sub_ps( c1_2, dx );
  __m128 t1 = _mm_add_ps( c1_2, dx );
  __m128 t2 = _mm_mul_ps( t0, t1 );

  __m128 t02 = _mm_mul_ps( t0, t0 );
  __m128 t12 = _mm_mul_ps( t1, t1 );

  
  w[0] = _mm_mul_ps( c1_24, _mm_mul_ps( t02, t02 ) );
  w[1] = _mm_mul_ps( c1_6, _mm_add_ps( c1_4, _mm_mul_ps( t0, _mm_add_ps( c1, _mm_mul_ps( t0, _mm_add_ps( c3_2, t2 ) ) ) ) ) );
  w[2] = _mm_add_ps( c11_24, _mm_mul_ps( t2, _mm_add_ps( c1_2, _mm_mul_ps( c1_4, t2 ) ) ) );
  w[3] = _mm_mul_ps( c1_6, _mm_add_ps( c1_4, _mm_mul_ps( t1, _mm_add_ps( c1, _mm_mul_ps( t1, _mm_add_ps( c3_2, t2 ) ) ) ) ) );
  w[4] = _mm_mul_ps( c1_24, _mm_mul_ps( t12, t12 ) );
  
}

static inline void vsplineh_s4( __m128 const dx, __m128 const h, __m128 w[] ) {
  
  __m128 const c1_2  = _mm_set1_ps( 1.0f/2.0f );
  __m128 const c1_24 = _mm_set1_ps( 1.0f/24.0f );
  __m128 const c1_4  = _mm_set1_ps( 1.0f/4.0f );
  __m128 const c1    = _mm_set1_ps( 1.0f );
  __m128 const c3_2  = _mm_set1_ps( 3.0f/2.0f );
  __m128 const c1_6  = _mm_set1_ps( 1.0f/6.0f );
  __m128 const c11_24 = _mm_set1_ps( 11.0f/24.0f );

  __m128 t0 = _mm_sub_ps( _mm_sub_ps( c1, h ), dx );
  __m128 t1 = _mm_add_ps(                 h  , dx );
  __m128 t2 = _mm_mul_ps( t0, t1 );

  __m128 t02 = _mm_mul_ps( t0, t0 );
  __m128 t12 = _mm_mul_ps( t1, t1 );
  
  w[0] = _mm_mul_ps( c1_24, _mm_mul_ps( t02, t02 ) );
  w[1] = _mm_mul_ps( c1_6, _mm_add_ps( c1_4, _mm_mul_ps( t0, _mm_add_ps( c1, _mm_mul_ps( t0, _mm_add_ps( c3_2, t2 ) ) ) ) ) );
  w[2] = _mm_add_ps( c11_24, _mm_mul_ps( t2, _mm_add_ps( c1_2, _mm_mul_ps( c1_4, t2 ) ) ) );
  w[3] = _mm_mul_ps( c1_6, _mm_add_ps( c1_4, _mm_mul_ps( t1, _mm_add_ps( c1, _mm_mul_ps( t1, _mm_add_ps( c3_2, t2 ) ) ) ) ) );
  w[4] = _mm_mul_ps( c1_24, _mm_mul_ps( t12, t12 ) );
  
}


static inline void vsplined_s4( __m128d const dx, __m128d w[] ) {
  
  __m128d const c1_2  = _mm_set1_pd( 1.0/2.0 );
  __m128d const c1_24 = _mm_set1_pd( 1.0/24.0 );
  __m128d const c1_4  = _mm_set1_pd( 1.0/4.0 );
  __m128d const c1    = _mm_set1_pd( 1.0 );
  __m128d const c3_2  = _mm_set1_pd( 3.0/2.0 );
  __m128d const c1_6  = _mm_set1_pd( 1.0/6.0 );
  __m128d const c11_24 = _mm_set1_pd( 11.0/24.0 );
  
  __m128d t0 = _mm_sub_pd( c1_2, dx );
  __m128d t1 = _mm_add_pd( c1_2, dx );
  __m128d t2 = _mm_mul_pd( t0, t1 );

  __m128d t02 = _mm_mul_pd( t0, t0 );
  __m128d t12 = _mm_mul_pd( t1, t1 );

  
  w[0] = _mm_mul_pd( c1_24, _mm_mul_pd( t02, t02 ) );
  w[1] = _mm_mul_pd( c1_6, _mm_add_pd( c1_4, _mm_mul_pd( t0, _mm_add_pd( c1, _mm_mul_pd( t0, _mm_add_pd( c3_2, t2 ) ) ) ) ) );
  w[2] = _mm_add_pd( c11_24, _mm_mul_pd( t2, _mm_add_pd( c1_2, _mm_mul_pd( c1_4, t2 ) ) ) );
  w[3] = _mm_mul_pd( c1_6, _mm_add_pd( c1_4, _mm_mul_pd( t1, _mm_add_pd( c1, _mm_mul_pd( t1, _mm_add_pd( c3_2, t2 ) ) ) ) ) );
  w[4] = _mm_mul_pd( c1_24, _mm_mul_pd( t12, t12 ) );
  
}

static inline void vsplinehd_s4( __m128d const dx, __m128d const h, __m128d w[] ) {
  
  __m128d const c1_2  = _mm_set1_pd( 1.0/2.0 );
  __m128d const c1_24 = _mm_set1_pd( 1.0/24.0 );
  __m128d const c1_4  = _mm_set1_pd( 1.0/4.0 );
  __m128d const c1    = _mm_set1_pd( 1.0 );
  __m128d const c3_2  = _mm_set1_pd( 3.0/2.0 );
  __m128d const c1_6  = _mm_set1_pd( 1.0/6.0 );
  __m128d const c11_24 = _mm_set1_pd( 11.0/24.0 );

  __m128d t0 = _mm_sub_pd( _mm_sub_pd( c1, h ), dx );
  __m128d t1 = _mm_add_pd(                 h  , dx );
  __m128d t2 = _mm_mul_pd( t0, t1 );

  __m128d t02 = _mm_mul_pd( t0, t0 );
  __m128d t12 = _mm_mul_pd( t1, t1 );
  
  w[0] = _mm_mul_pd( c1_24, _mm_mul_pd( t02, t02 ) );
  w[1] = _mm_mul_pd( c1_6, _mm_add_pd( c1_4, _mm_mul_pd( t0, _mm_add_pd( c1, _mm_mul_pd( t0, _mm_add_pd( c3_2, t2 ) ) ) ) ) );
  w[2] = _mm_add_pd( c11_24, _mm_mul_pd( t2, _mm_add_pd( c1_2, _mm_mul_pd( c1_4, t2 ) ) ) );
  w[3] = _mm_mul_pd( c1_6, _mm_add_pd( c1_4, _mm_mul_pd( t1, _mm_add_pd( c1, _mm_mul_pd( t1, _mm_add_pd( c3_2, t2 ) ) ) ) ) );
  w[4] = _mm_mul_pd( c1_24, _mm_mul_pd( t12, t12 ) );
  
}


#endif /* _SPLINES_H */

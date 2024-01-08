/*****************************************************************************************

 Spline function definitions, BG/Q optimized version

*****************************************************************************************/

#ifndef _SPLINES_H
#define _SPLINES_H

// #include "vector-bgq.h"


/*****************************************************************************************
vspline_n1

1st order (linear) splines.

w[0] = 1/2 - x
w[1] = 1/2 + x
*****************************************************************************************/

inline void vspline_s1( vector const dx, vector w[] ) {
  
  vector const c1_2 = vec_splats(0.5);
  
  w[0] = vec_sub( c1_2, dx );
  w[1] = vec_add( c1_2, dx );
}

inline void vsplineh_s1( vector const dx, vector const h, vector w[] ) {
  
  vector const c1 = vec_splats(1.0);
  
  w[0] = vec_sub( vec_sub( c1, h ), dx );
  w[1] = vec_add(              h  , dx );
}



/*****************************************************************************************
vspline_n2

2nd order (quadratic) splines.


w[0] = (1 - 2*dx)^2/8.
w[1] = 0.75 - dx^2
w[2] = (1 + 2*dx)^2/8.

*****************************************************************************************/


 inline void vspline_s2( vector const dx, vector w[] ) {
  
  vector const c1_2      = vec_splats(0.5);
  
  vector t0 = vec_sub( c1_2, dx );
  vector t1 = vec_add( c1_2, dx );
  
  w[0] = vec_mul( c1_2, vec_mul( t0, t0 ) );
  w[1] = vec_madd( t0, t1, c1_2 );
  w[2] = vec_mul( c1_2, vec_mul( t1, t1 ) );
}

inline void vsplineh_s2( vector const dx, vector const h, vector w[] ) {

  vector const c1        = vec_splats(1.0);
  vector const c1_2      = vec_splats(0.5);
  
  vector t0 = vec_sub( vec_sub( c1, h ), dx );
  vector t1 = vec_add(              h  , dx );
  
  w[0] = vec_mul( c1_2, vec_mul( t0, t0 ) );
  w[1] = vec_madd( t0, t1, c1_2 );
  w[2] = vec_mul( c1_2, vec_mul( t1, t1 ) );

}

/*****************************************************************************************
vspline_n3

3nd order (cubic) splines.

w[0] =  1/6 (1/2 - dx)^3
w[1] = 1/6 (4 - 6 (1/2 + dx)^2 + 3 (1/2 + dx)^3)
w[2] = 1/6 (4 - 6 (1/2 - dx)^2 + 3 (1/2 - dx)^3)
w[3] = 1/6 (1/2 + dx)^3
*****************************************************************************************/


inline void vspline_s3( vector const dx, vector w[] ) {
  
  vector const c1_2 = vec_splats( 1.0/2.0 );
  vector const c2_3 = vec_splats( 2.0/3.0 );
  vector const c1_6 = vec_splats( 1.0/6.0 );

  vector t0 = vec_sub( c1_2, dx );
  vector t1 = vec_add( c1_2, dx );

  vector t2 = vec_mul( t0, t0 );
  vector t3 = vec_mul( t1, t1 );
    
  t0 = vec_mul( t0, t2 );
  t1 = vec_mul( t1, t3 );
  
  w[0] = vec_mul( t0, c1_6 );
  w[1] = vec_madd( t1, c1_2, vec_sub( c2_3, t3 ) );
  w[2] = vec_madd( t0, c1_2, vec_sub( c2_3, t2 ) );
  w[3] = vec_mul( t1, c1_6 );
}

inline void vsplineh_s3( vector const dx, vector const h, vector w[] ) {

  vector const c1   = vec_splats( 1.0 );
  vector const c1_2 = vec_splats( 1.0/2.0 );
  vector const c2_3 = vec_splats( 2.0/3.0 );
  vector const c1_6 = vec_splats( 1.0/6.0 );

  vector t0 = vec_sub( vec_sub( c1, h ), dx );
  vector t1 = vec_add(              h  , dx );

  vector t2 = vec_mul( t0, t0 );
  vector t3 = vec_mul( t1, t1 );
    
  t0 = vec_mul( t0, t2 );
  t1 = vec_mul( t1, t3 );
  
  w[0] = vec_mul( t0, c1_6 );
  w[1] = vec_madd( t1, c1_2, vec_sub( c2_3, t3 ) );
  w[2] = vec_madd( t0, c1_2, vec_sub( c2_3, t2 ) );
  w[3] = vec_mul( t1, c1_6 );

}

/*****************************************************************************************
vspline_n4

4th order (quartic) splines.

w[0] = (1/384) (1 - 2*x)^4
w[1] = (1/96) (19 - 44*x + 24*x^2 + 16*x^3 - 16*x^4)
w[2] = (115/192) - (5*x^2)/8 + x^4/4
w[3] = (1/96) (19 + 44*x + 24*x^2 - 16*x^3 - 16*x^4)
w[4] = (1/384) (1 + 2*x)^4
*****************************************************************************************/


inline void vspline_s4( vector const dx, vector w[] ) {
  
  vector const c1_2  = vec_splats( 1.0f/2.0f );
  vector const c1_24 = vec_splats( 1.0f/24.0f );
  vector const c1_4  = vec_splats( 1.0f/4.0f );
  vector const c1    = vec_splats( 1.0f );
  vector const c3_2  = vec_splats( 3.0f/2.0f );
  vector const c1_6  = vec_splats( 1.0f/6.0f );
  vector const c11_24 = vec_splats( 11.0f/24.0f );
  
  vector t0 = vec_sub( c1_2, dx );
  vector t1 = vec_add( c1_2, dx );
  vector t2 = vec_mul( t0, t1 );

  vector t02 = vec_mul( t0, t0 );
  vector t12 = vec_mul( t1, t1 );
  
  w[0] = vec_mul( c1_24, vec_mul( t02, t02 ) );
  w[1] = vec_mul( c1_6, vec_madd( vec_madd( vec_add( c3_2, t2 ), t0, c1 ), t0, c1_4 ) );
  w[2] = vec_madd( vec_madd( c1_4, t2, c1_2 ), t2, c11_24 );
  w[3] = vec_mul( c1_6, vec_madd( vec_madd( vec_add( c3_2, t2 ), t1, c1 ), t1, c1_4 ) );
  w[4] = vec_mul( c1_24, vec_mul( t12, t12 ) );
  
}

inline void vsplineh_s4( vector const dx, vector const h, vector w[] ) {

  vector const c1_2  = vec_splats( 1.0f/2.0f );
  vector const c1_24 = vec_splats( 1.0f/24.0f );
  vector const c1_4  = vec_splats( 1.0f/4.0f );
  vector const c1    = vec_splats( 1.0f );
  vector const c3_2  = vec_splats( 3.0f/2.0f );
  vector const c1_6  = vec_splats( 1.0f/6.0f );
  vector const c11_24 = vec_splats( 11.0f/24.0f );
  
  vector t0  = vec_sub( vec_sub( c1, h ), dx );
  vector t1  = vec_add(              h  , dx );

  vector t2  = vec_mul( t0, t1 );
  vector t02 = vec_mul( t0, t0 );
  vector t12 = vec_mul( t1, t1 );

  w[0] = vec_mul( c1_24, vec_mul( t02, t02 ) );
  w[1] = vec_mul( c1_6, vec_madd( vec_madd( vec_add( c3_2, t2 ), t0, c1 ), t0, c1_4 ) );
  w[2] = vec_madd( vec_madd( c1_4, t2, c1_2 ), t2, c11_24 );
  w[3] = vec_mul( c1_6, vec_madd( vec_madd( vec_add( c3_2, t2 ), t1, c1 ), t1, c1_4 ) );
  w[4] = vec_mul( c1_24, vec_mul( t12, t12 ) );

}

#endif /* _SPLINES_H */

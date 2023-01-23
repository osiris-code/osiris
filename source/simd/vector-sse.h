/****************************************************************************************/
/* SIMD utilites.                                                                       */
/* Definitions for Intel SSE                                                            */
/****************************************************************************************/

#ifndef _VECTOR_SSE_H
#define _VECTOR_SSE_H

#ifdef __SSE__

#include <stdio.h>

/* SSE definitions */
#include <xmmintrin.h>
#include <pmmintrin.h>


#if defined(__GNUC__) || defined(__PGI)
  /* gcc or Portland Group pgcc */
  #define DECLARE_ALIGNED_16(v)       v __attribute__ ((aligned (16)))
#else
  /* Intel icc */
  #define DECLARE_ALIGNED_16(v)       __declspec(align(16)) v
#endif

#include <stdint.h>
#define IS_ALIGNED_16(addr) (((uintptr_t)addr & 0xF) == 0)

#if defined( PRECISION_SINGLE )
#define VEC_WIDTH 4
#elif defined( PRECISION_DOUBLE )
#define VEC_WIDTH 2
#else
#error PRECISION_SINGLE or PRECISION_DOUBLE must be defined
#endif 

typedef union vecDouble {double v[2]; __m128d v2;} dvec;
typedef union vecFloat {float v[4]; __m128 v4;} fvec;
typedef union vecInt {int v[4]; __m128i v4;} ivec;



/* Multiply an integer vector by a scalar */

#ifdef __SSE4_1__

#include <smmintrin.h>
  
  // Use SSE4 DWORD multiplication 
  static inline __m128i vsmul( const __m128i a, const unsigned b) {
	return  _mm_mullo_epi32( a, _mm_set1_epi32(b) );
  }
#else
  /* Use SSE3 instructions */
  static inline __m128i vsmul( const __m128i a, const unsigned b) {
	register __m128i vb, tmp1, tmp2;
	vb = _mm_set1_epi32( b );  
	tmp1 = _mm_mul_epu32( a, vb );
	tmp2 = _mm_mul_epu32( _mm_shuffle_epi32( a,  _MM_SHUFFLE( 2, 3, 0, 1 )) , vb );
	return  _mm_unpacklo_epi32( _mm_shuffle_epi32( tmp1, _MM_SHUFFLE( 3, 1, 2, 0 )), 
								_mm_shuffle_epi32( tmp2, _MM_SHUFFLE( 3, 1, 2, 0 )) );
  }

#endif

/* Blendv Operation 

The SSE4 blendv operation does:

__m128 _mm_blendv_ps(__m128 a, __m128 b, __m128 m)

Does:

r0 = ( m0 == 0 )? a0 : b0
r1 = ( m1 == 0 )? a1 : b1
r2 = ( m2 == 0 )? a2 : b2
r3 = ( m3 == 0 )? a3 : b3

and m must be set to either 0x0 or 0xffffffff. This can be done using the result of a
logical operation, e.g.:

m = _mm_cmplt_ps( c, d );  // Compares for less than

Or a custom mask can be created using an integer cast, e.g.: 

// Does c0 = a0, c1 = b1, c2 = b2, c3 = a3

m = _mm_castsi128_ps( _mm_set_epi32( 0, -1, -1, 0 ) );
c = _mm_blendv_ps( a, b, m ); 

This mimics the behavior of the select operation available in other SIMD architectures
(e.g. Power). 

If SSE4 is not available then this could be done doing (b & m) | (a & ~m). A more 
efficient form can be found in:

http://markplusplus.wordpress.com/2007/03/14/fast-sse-select-operation/

*/

#ifndef __SSE4_1__

static inline __m128 _mm_blendv_ps(const __m128 a, const __m128 b, const __m128 mask)
{
    // (((b ^ a) & mask)^a)
    return _mm_xor_ps( a, _mm_and_ps( mask, _mm_xor_ps( b, a ) ) );
}

static inline __m128d _mm_blendv_pd(const __m128d a, const __m128d b, const __m128d mask)
{
    // (((b ^ a) & mask)^a)
    return _mm_xor_pd( a, _mm_and_pd( mask, _mm_xor_pd( b, a ) ) );
}

#endif

/*****************************************************************************************
  Add all vector elements and return result
*****************************************************************************************/

static inline float _mm_reduce_add_ps( const __m128 a )
{
  __m128 sum;
  sum = _mm_hadd_ps( a, a );
  sum = _mm_hadd_ps( sum, sum );
  return _mm_cvtss_f32( sum );
} 

static inline double _mm_reduce_add_pd( const __m128d a )
{
  __m128d sum;
  sum = _mm_hadd_pd( a, a );
  return _mm_cvtsd_f64( sum );
} 

/*****************************************************************************************

IEEE accurate single precision reciprocal and reciprocal square root using fast estimates
followed by a Newton-Raphson iteration.

*****************************************************************************************/

static inline __m128 rsqrt_ps( __m128 a ) {
  
  __m128 r0 = a;
  __m128 r1 = _mm_rsqrt_ps( r0 );  
  r0 = _mm_mul_ps( _mm_mul_ps( r0, r1 ), r1 );
  r0 = _mm_sub_ps( r0, _mm_set1_ps( 3.0f ) );
  r1 = _mm_mul_ps( r0, r1 );
  r1 = _mm_mul_ps( r1, _mm_set1_ps( -0.5f ) );
  
  return r1;
}


static inline __m128 rcp_ps( __m128 a ) {
  
  __m128 r0, r1;
  r0 = a;
  r1 = _mm_rcp_ps( r0 );  
  r0 = _mm_mul_ps( r0, r1 );
  r0 = _mm_mul_ps( r0, r1 );
  r1 = _mm_add_ps( r1, r1 );
  r1 = _mm_sub_ps( r1, r0 );
  
  return r1;
}


/*****************************************************************************************
_MM_LOAD4v3_PS, _MM_LOAD4v3_EPI32

Loads 4 x 3 elements vectors stored sequentially in memory into 3 vector
variables:

b = [x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4] -> _MM_LOAD4v3_PS(v1,v2,v3,b)

-> v1 = [x1,x2,x3,x4]
-> v2 = [y1,y2,y3,y4]
-> v3 = [z1,z2,z3,z4]

********************************************************************************************/

#define _MM_LOAD4v3_PS(v1, v2, v3, b) {                                 \
	   register __m128 tmpA, tmpB, tmpC, tmp0, tmp1;                    \
	                                                                    \
	   tmpA = _mm_load_ps( (b) );                                       \
	   tmpB = _mm_load_ps( (b) + 4 );                                   \
	   tmpC = _mm_load_ps( (b) + 8 );                                   \
                                                                        \
	   tmp0 = _mm_shuffle_ps( tmpA, tmpB, _MM_SHUFFLE( 1, 0, 2, 1 ) );  \
       tmp1 = _mm_shuffle_ps( tmpB, tmpC, _MM_SHUFFLE( 2, 1, 3, 2 ) );  \
	                                                                    \
	   (v1) = _mm_shuffle_ps( tmpA, tmp1, _MM_SHUFFLE( 2, 0, 3, 0 ) );  \
	   (v2) = _mm_shuffle_ps( tmp0, tmp1, _MM_SHUFFLE( 3, 1, 2, 0 ) );  \
	   (v3) = _mm_shuffle_ps( tmp0, tmpC, _MM_SHUFFLE( 3, 0, 3, 1) );   \
}

// Using PS shuffles saves some cycles
#define _MM_LOAD4v3_EPI32(v1, v2, v3, b) {                              \
	   register __m128 tmpA, tmpB, tmpC, tmp0, tmp1;                    \
	                                                                    \
	   tmpA = _mm_load_ps( (float *)((b)) );                            \
	   tmpB = _mm_load_ps( (float *)((b) + 4) );                        \
	   tmpC = _mm_load_ps( (float *)((b) + 8) );                        \
                                                                        \
	   tmp0 = _mm_shuffle_ps( tmpA, tmpB, _MM_SHUFFLE( 1, 0, 2, 1 ) );  \
       tmp1 = _mm_shuffle_ps( tmpB, tmpC, _MM_SHUFFLE( 2, 1, 3, 2 ) );  \
	                                                                    \
	   (v1) = _mm_castps_si128( _mm_shuffle_ps( tmpA, tmp1, _MM_SHUFFLE( 2, 0, 3, 0 ) ) );  \
	   (v2) = _mm_castps_si128( _mm_shuffle_ps( tmp0, tmp1, _MM_SHUFFLE( 3, 1, 2, 0 ) ) );  \
	   (v3) = _mm_castps_si128( _mm_shuffle_ps( tmp0, tmpC, _MM_SHUFFLE( 3, 0, 3, 1 ) ) );   \
}


/*****************************************************************************************
_MM_STORE4v3_PS, _MM_STORE4v3_EPI32

Stores 3 vector variables into memory putting each vector element 
sequentially:

v1 = [x1,x2,x3,x4], v2 = [y1,y2,y3,y4], v3 = [z1,z2,z3,z4]

-> _MM_STORE4v3_PS(b, v1, v2, v3, b)

-> b = [x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4]

*****************************************************************************************/

#define _MM_STORE4v3_PS(b, v1, v2, v3)      {							\
       register __m128 tmp0, tmp1, tmp2;								\
       register __m128 tmpA, tmpB, tmpC;								\
     																	\
       tmp0 = _mm_shuffle_ps( (v1), (v2), _MM_SHUFFLE( 2, 0, 2, 0 ) );	\
       tmp1 = _mm_shuffle_ps( (v3), (v1), _MM_SHUFFLE( 3, 1, 2, 0 ) );	\
       tmp2 = _mm_shuffle_ps( (v2), (v3), _MM_SHUFFLE( 3, 1, 3, 1 ) );	\
																		\
       tmpA = _mm_shuffle_ps( tmp0, tmp1, _MM_SHUFFLE( 2, 0, 2, 0 ) );	\
       tmpB = _mm_shuffle_ps( tmp2, tmp0, _MM_SHUFFLE( 3, 1, 2, 0 ) );	\
       tmpC = _mm_shuffle_ps( tmp1, tmp2, _MM_SHUFFLE( 3, 1, 3, 1 ) );	\
																		\
       _mm_store_ps((b),   tmpA );										\
       _mm_store_ps((b)+4, tmpB );										\
	   _mm_store_ps((b)+8, tmpC );    									\
}

// Using PS shuffles saves some cycles
#define _MM_STORE4v3_EPI32(b, v1, v2, v3)      {												  \
  register __m128 tmp0, tmp1, tmp2;																  \
  register __m128 tmpA, tmpB, tmpC;																  \
																   								  \
  tmp0 = _mm_shuffle_ps( _mm_castsi128_ps(v1), _mm_castsi128_ps(v2), _MM_SHUFFLE( 2, 0, 2, 0 ) ); \
  tmp1 = _mm_shuffle_ps( _mm_castsi128_ps(v3), _mm_castsi128_ps(v1), _MM_SHUFFLE( 3, 1, 2, 0 ) ); \
  tmp2 = _mm_shuffle_ps( _mm_castsi128_ps(v2), _mm_castsi128_ps(v3), _MM_SHUFFLE( 3, 1, 3, 1 ) ); \
																   								  \
  tmpA = _mm_shuffle_ps( tmp0, tmp1, _MM_SHUFFLE( 2, 0, 2, 0 ) );								  \
  tmpB = _mm_shuffle_ps( tmp2, tmp0, _MM_SHUFFLE( 3, 1, 2, 0 ) );								  \
  tmpC = _mm_shuffle_ps( tmp1, tmp2, _MM_SHUFFLE( 3, 1, 3, 1 ) );								  \
																   								  \
  _mm_store_ps((float *)((b)),   tmpA );														  \
  _mm_store_ps((float *)((b)+4), tmpB );														  \
  _mm_store_ps((float *)((b)+8), tmpC );    													  \
}

/*****************************************************************************************
_MM_LOAD4v2_PS, _MM_LOAD4v2_EPI32

Loads 4 x 2 elements vectors stored sequentially in memory into 2 vector
variables:

b = [x1,y1,x2,y2,x3,y3,x4,y4] -> _MM_LOAD4v2_PS(v1,v2,b)

-> v1 = [x1,x2,x3,x4]
-> v2 = [y1,y2,y3,y4]

*****************************************************************************************/

#define _MM_LOAD4v2_PS(v1, v2, b) {                                     \
	   register __m128 tmp0, tmp1;                     		 			\
	                                                                    \
	   tmp0 = _mm_load_ps( (b) );                                       \
	   tmp1 = _mm_load_ps( (b) + 4 );                                   \
	                                                                    \
	   (v1) = _mm_shuffle_ps( tmp0, tmp1 , _MM_SHUFFLE( 2, 0, 2, 0 ) ); \
	   (v2) = _mm_shuffle_ps( tmp0, tmp1 , _MM_SHUFFLE( 3, 1, 3, 1 ) ); \
}

#define _MM_LOAD4v2_EPI32(v1, v2, b) {                                 	\
	   register __m128i tmp0, tmp1, tmpA, tmpB;							\
																		\
	   tmp0 = _mm_load_si128( (__m128i *)(b) );							\
	   tmp1 = _mm_load_si128( (__m128i *)((b)+4));   					\
																		\
	   tmpA = _mm_unpacklo_epi32( tmp0, tmp1 );							\
	   tmpB = _mm_unpackhi_epi32( tmp0, tmp1 );							\
																		\
	   (v1) = _mm_unpacklo_epi32( tmpA, tmpB );							\
	   (v2) = _mm_unpackhi_epi32( tmpA, tmpB );							\
}


/*****************************************************************************************
_MM_STORE4v2_PS, _MM_STORE4v2_EPI32

Stores 2 vector variables into memory putting each vector element 
sequentially:

v1 = [x1,x2,x3,x4], v2 = [y1,y2,y3,y4]

-> _MM_STORE4v2_PS(b, v1, v2)

-> b = [x1,y1,x2,y2,x3,y3,x4,y4]

*****************************************************************************************/

#define _MM_STORE4v2_PS(b, v1, v2) {                                    \
       _mm_store_ps((b),   _mm_unpacklo_ps( (v1), (v2) ) );				\
       _mm_store_ps((b)+4, _mm_unpackhi_ps( (v1), (v2) ) );				\
}

#define _MM_STORE4v2_EPI32(b, v1, v2) {                                        \
       _mm_store_si128((__m128i *)(b),     _mm_unpacklo_epi32( (v1), (v2) ) ); \
       _mm_store_si128((__m128i *)((b)+4), _mm_unpackhi_epi32( (v1), (v2) ) ); \
}

/*****************************************************************************************
_MM_LOAD2v3_PD

Loads 2 x 3 elements vectors stored sequentially in memory into 3 vector
variables:

b = [x1,y1,z1,x2,y2,z2] -> _MM_LOAD2v3_PS(v1,v2,v3,b)

-> v1 = [x1,x2]
-> v2 = [y1,y2]
-> v3 = [z1,z2]

*****************************************************************************************/

#define _MM_LOAD2v3_PD(v1, v2, v3, buf) {                           \
	   register __m128d tmp0, tmp1, tmp2;               	      \
	   tmp0 = _mm_load_pd( &(buf)[0] );                           \
	   tmp1 = _mm_load_pd( &(buf)[2] );                           \
	   tmp2 = _mm_load_pd( &(buf)[4] );                           \
	   (v1) = _mm_shuffle_pd( tmp0, tmp1, _MM_SHUFFLE2( 1, 0 ) ); \
	   (v2) = _mm_shuffle_pd( tmp0, tmp2, _MM_SHUFFLE2( 0, 1 ) ); \
	   (v3) = _mm_shuffle_pd( tmp1, tmp2, _MM_SHUFFLE2( 1, 0 ) ); \
}

/*****************************************************************************************
_MM_STORE2v3_PD

Stores 3 vector variables into memory putting each vector element 
sequentially:

v1 = [x1,x2], v2 = [y1,y2], v3 = [z1,z2]

-> _MM_STORE2v3_PD(b, v1, v2, v3)

-> b = [x1,y1,z1,x2,y2,z2]

*****************************************************************************************/

#define _MM_STORE2v3_PD(b, v1, v2, v3) {                                             \
       _mm_store_pd( &(b)[0], _mm_unpacklo_pd( (v1), (v2) ) );				         \
       _mm_store_pd( &(b)[2], _mm_shuffle_pd( (v3), (v1), _MM_SHUFFLE2( 1, 0 ) ) );	 \
       _mm_store_pd( &(b)[4], _mm_unpackhi_pd( (v2), (v3) ) );				         \
}

/*****************************************************************************************
_MM_LOAD2v2_PD

Loads 2 x 2 elements vectors stored sequentially in memory into 2 vector
variables:

b = [x1,y1,x2,y2] -> _MM_LOAD2v2_PD(v1,v2,b)

-> v1 = [x1,x2]
-> v2 = [y1,y2]

*****************************************************************************************/

#define _MM_LOAD2v2_PD(v1, v2, buf) {                             \
	   register __m128d tmp0, tmp1;                    		      \
	   tmp0 = _mm_load_pd( &(buf)[0] );                           \
	   tmp1 = _mm_load_pd( &(buf)[2] );                           \
	   (v1) = _mm_unpacklo_pd( tmp0, tmp1 );                      \
	   (v2) = _mm_unpackhi_pd( tmp0, tmp1 );                      \
}

/*****************************************************************************************
_MM_STORE2v2_PD

Stores 2 vector variables into memory putting each vector element 
sequentially:

v1 = [x1,x2], v2 = [y1,y2]

-> _MM_STORE2v2_PD(b, v1, v2)

-> b = [x1,y1,x2,y2]

*****************************************************************************************/


#define _MM_STORE2v2_PD(b, v1, v2) {                                        \
       _mm_store_pd( &(b)[0],  _mm_unpacklo_pd( (v1), (v2) ) );				\
       _mm_store_pd( &(b)[2],  _mm_unpackhi_pd( (v1), (v2) ) );				\
}


/*****************************************************************************************
_MM_LOAD2v2_EPI32

Loads 2 x 2 elements integer vectors stored sequentially in memory into 2 vector
variables:

b = [x1,y1,x2,y2] -> _MM_LOAD2v2_EPI32(v1,v2,b)

-> v1 = [x1,x2,X,X]
-> v2 = [y1,y2,X,X]


*****************************************************************************************/

#define _MM_LOAD2v2_EPI32( vix, viy, buf ) {          \
  __m128i t0  = _mm_load_si128( (__m128i *) buf );  \
  __m128i t1 = _mm_unpacklo_epi32( t0, t0 );        \
  t0 = _mm_unpackhi_epi32( t0, t0 );                \
  vix = _mm_unpacklo_epi32( t1, t0 );               \
  viy = _mm_unpackhi_epi32( t1, t0 );               \
}

/*****************************************************************************************
_MM_STORE2v2_EPI32

Stores 2x 2 element vector variables into memory putting each vector element 
sequentially:

v1 = [x1,x2,X,X], v2 = [y1,y2,X,X]

-> _MM_STORE2v2_EPI32(b, v1, v2)

-> b = [x1,y1,x2,y2]

*****************************************************************************************/

#define _MM_STORE2v2_EPI32( buf, va, vb ) { \
  _mm_store_si128( (__m128i *) &buf[0] , _mm_unpacklo_epi32( va, vb ) ); \
}

/*****************************************************************************************
_MM_LOAD2v3_EPI32

Loads 3 x 2 elements integer vectors stored sequentially in memory into 2 vector
variables:

b = [x1,y1,z1,x2,y2,z2] -> _MM_LOAD2v3_EPI32(v1,v2,v3,b)

-> v1 = [x1,x2,X,X]
-> v2 = [y1,y2,X,X]
-> v3 = [z1,z2,X,X]


*****************************************************************************************/

#define _MM_LOAD2v3_EPI32( vix, viy, viz, buf ) {                                     \
  __m128i t0, t1;                                                                     \
  t0 = _mm_loadu_si128( (__m128i *) &buf[0] );                                        \
  t1 = _mm_loadl_epi64( (__m128i *) &buf[4] );                                        \
  vix = _mm_shuffle_epi32( t0, _MM_SHUFFLE( 0, 2, 3, 0 ) );                           \
  viy = _mm_shuffle_epi32( _mm_unpacklo_epi32( t1, t0 ), _MM_SHUFFLE( 0, 2, 0, 3 ) ); \
  viz = _mm_unpackhi_epi32( vix, viy );                                               \
}


/*****************************************************************************************
_MM_STORE2v3_EPI32

Stores 3x 2 element vector variables into memory putting each vector element 
sequentially:

v1 = [x1,x2,X,X], v2 = [y1,y2,X,X], v3 = [z1,z2,X,X]

-> _MM_STORE2v3_EPI32(b, v1, v2, v3)

-> b = [x1,y1,z1,x2,y2,z3]

*****************************************************************************************/

#define _MM_STORE2v3_EPI32( buf, vix, viy, viz ) {                                       \
  __m128i t0, t1;                                                                        \
  t0 = _mm_unpacklo_epi32( vix, viy );                                                   \
  t1 = _mm_unpacklo_epi64( viy, viz );                                                   \
  _mm_storeu_si128( (__m128i *) &buf[0],                                                 \
                     _mm_unpacklo_epi64( t0, _mm_unpackhi_epi32( t1, t0 ) ) );           \
  _mm_storel_epi64( (__m128i *) &buf[4],                                                 \
                     _mm_shuffle_epi32( t1, _MM_SHUFFLE( 0, 2, 3, 1 ) ) );               \
}


/*****************************************************************************************

Convert 2 x 2 double precision vectors into / from 1 x 4 single precision vector

*****************************************************************************************/

#define _MM_CVTPD2_PS( a, b ) \
_mm_movelh_ps( _mm_cvtpd_ps( a ), _mm_cvtpd_ps( b ))

#define _MM_CVTPS2_PD( a, b, f ) {             \
  a = _mm_cvtps_pd( f );                       \
  b = _mm_cvtps_pd( _mm_movehl_ps( f, f ) );   \
}


/*****************************************************************************************

Utilities for printing vector values

*****************************************************************************************/

static inline void printfv(char *tag, __m128 v)
{
  fvec p;
  p.v4 = v;
  printf("%s = %g, %g, %g, %g\n", tag, p.v[0], p.v[1], p.v[2], p.v[3] );
  
}

static inline  void printiv(char *tag, __m128i v)
{
  ivec p;
  p.v4 = v;
  printf("%s = %i, %i, %i, %i\n", tag, p.v[0], p.v[1], p.v[2], p.v[3] );
  
}

static inline void printdv(char *tag, __m128d v)
{
  dvec p;
  p.v2 = v;
  printf("%s = %f, %f\n", tag, p.v[0], p.v[1] );
  
}

#include <math.h>
static inline int testfv( __m128 vec)
{
  DECLARE_ALIGNED_16(float v[4]);
  
  _mm_store_ps( v, vec );

  int i;
  for(i = 0; i < 4; i++) {
    if ( isnan( v[i] ) || isinf( v[i] ) ) return -1;
  } 
  return 0;
}

#define TEST_VECTOR( a ) {	                                                             \
  if ( testfv(a) ) {                                                                     \
     DECLARE_ALIGNED_16(float v[4]);                                                     \
     _mm_store_ps( v, a );                                                               \
 	 printf(" %s:%d %s = %.3f", __FILE__, __LINE__, #a, v[0] );                          \
	 for(i=1; i<4; i++) printf(", %.3f", v[i]);                                          \
	 printf("\n");                                                                       \
  }                                                                                      \
}

#endif

/* End of definitions for Intel SSE */

#endif

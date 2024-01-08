/****************************************************************************************/
/* SIMD utilites.                                                                       */
/* Definitions for Intel AVX                                                            */
/****************************************************************************************/

#ifndef _VECTOR_AVX_H
#define _VECTOR_AVX_H


#include <stdio.h>
#include <immintrin.h>

#if defined(__GNUC__) || defined(__PGI)
  /* gcc or Portland Group pgcc */
  #define DECLARE_ALIGNED_16(v)       v __attribute__ ((aligned (16)))
  #define DECLARE_ALIGNED_32(v)       v __attribute__ ((aligned (32)))
#else
  /* Intel icc */
  #define DECLARE_ALIGNED_16(v)       __declspec(align(16)) v
  #define DECLARE_ALIGNED_32(v)       __declspec(align(32)) v
#endif

#include <stdint.h>
#define IS_ALIGNED_16(addr) (((uintptr_t)addr & 0x0F) == 0)
#define IS_ALIGNED_32(addr) (((uintptr_t)addr & 0x1F) == 0)

#if defined( PRECISION_SINGLE )
#define VEC_WIDTH 8
#elif defined( PRECISION_DOUBLE )
#define VEC_WIDTH 4
#else
#error PRECISION_SINGLE or PRECISION_DOUBLE must be defined
#endif 

typedef union vecDouble {double v[4]; __m256d v4;} dvec;
typedef union vecFloat  {float v[8];  __m256  v8;} fvec;
typedef union vecInt    {int v[8];    __m256i v8;} ivec;

/*****************************************************************************************
  Integer arithmetic support instructions.
*****************************************************************************************/

// Multiply a vector by a scalar
static inline __m256i vimul_s( const __m256i a, const unsigned b) {
  
  register __m256i m;
  register __m128i c0;
  
  c0 = _mm_set1_epi32(b);
  m = _mm256_castsi128_si256( _mm_mullo_epi32( _mm256_castsi256_si128(a), c0 ) );
  m = _mm256_insertf128_si256( m, _mm_mullo_epi32( _mm256_extractf128_si256(a,1), c0 ), 1 );
  
  return  m;
}

// vector addition
static inline __m256i viadd( const __m256i a, const __m256i b) {
  
  register __m256i v; 
  
  v = _mm256_castsi128_si256( _mm_add_epi32( _mm256_castsi256_si128(a), 
                                             _mm256_castsi256_si128(b) ) );
  v = _mm256_insertf128_si256( v, _mm_add_epi32( _mm256_extractf128_si256(a,1), 
                                                 _mm256_extractf128_si256(b,1) ), 1 );
  
  return  v;
}

// 3 vector addition
static inline __m256i viadd3( const __m256i a, const __m256i b, const __m256i c) {
  
  register __m256i v; 
  
  v = _mm256_castsi128_si256( _mm_add_epi32( 
                                 _mm_add_epi32( _mm256_castsi256_si128(a), 
                                                _mm256_castsi256_si128(b) ),
                                 _mm256_castsi256_si128(c) ) );
  v = _mm256_insertf128_si256( v, _mm_add_epi32( 
                                     _mm_add_epi32( _mm256_extractf128_si256(a,1), 
                                                    _mm256_extractf128_si256(b,1) ),
                                     _mm256_extractf128_si256(c,1) ), 1 );
  return  v;
}

static inline void viadd_store( const int *buf, const __m256i a, const __m256i b) {
  
  _mm_store_si128( (__m128i *) &buf[0], _mm_add_epi32( _mm256_castsi256_si128(a), 
                                         _mm256_castsi256_si128(b) ) );
  _mm_store_si128( (__m128i *) &buf[4], _mm_add_epi32( _mm256_extractf128_si256(a,1), 
                                           _mm256_extractf128_si256(b,1) ) );
}

static inline void viadd3_store( const int buf[], const __m256i a, const __m256i b, const __m256i c) {
  
  _mm_store_si128((__m128i *)  &buf[0], _mm_add_epi32( 
                                 _mm_add_epi32( _mm256_castsi256_si128(a), 
                                                _mm256_castsi256_si128(b) ),
                                                _mm256_castsi256_si128(c) ) );
  _mm_store_si128((__m128i *) &buf[4], _mm_add_epi32( 
                                 _mm_add_epi32( _mm256_extractf128_si256(a,1), 
                                                _mm256_extractf128_si256(b,1) ),
                                                _mm256_extractf128_si256(c,1) ) );
}

// vector subtraction
static inline __m256i visub( const __m256i a, const __m256i b) {
  
  register __m256i v; 
  
  v = _mm256_castsi128_si256( _mm_sub_epi32( _mm256_castsi256_si128(a), 
                                             _mm256_castsi256_si128(b) ) );
  v = _mm256_insertf128_si256( v, _mm_sub_epi32( _mm256_extractf128_si256(a,1), 
                                                 _mm256_extractf128_si256(b,1) ), 1 );
  
  return  v;
}


/*****************************************************************************************

IEEE accurate single precision reciprocal and reciprocal square root using fast estimates
followed by a Newton-Raphson iteration.

*****************************************************************************************/

static inline __m256 rsqrt_ps( const __m256 a ) {
  
  __m256 r0 = a;
  __m256 r1 = _mm256_rsqrt_ps( r0 );  
  r0 = _mm256_mul_ps( _mm256_mul_ps( r0, r1 ), r1 );
  r0 = _mm256_sub_ps( r0, _mm256_set1_ps( 3.0f ) );
  r1 = _mm256_mul_ps( r0, r1 );
  r1 = _mm256_mul_ps( r1, _mm256_set1_ps( -0.5f ) );
  
  return r1;
}


static inline __m256 rcp_ps( const __m256 a ) {
  
  __m256 r0, r1;
  r0 = a;
  r1 = _mm256_rcp_ps( r0 );  
  r0 = _mm256_mul_ps( r0, r1 );
  r0 = _mm256_mul_ps( r0, r1 );
  r1 = _mm256_add_ps( r1, r1 );
  r1 = _mm256_sub_ps( r1, r0 );
  
  return r1;
}

/*****************************************************************************************

Reduction operations

*****************************************************************************************/

static inline double _mm256_reduce_add_pd( const __m256d a ) {
  
  __m256d r; 
 
  r =  _mm256_hadd_pd( a, _mm256_permute2f128_pd( a, a, 1 ) );
  r =  _mm256_hadd_pd( r, r );

  return _mm_cvtsd_f64( _mm256_castpd256_pd128( r ) );
}


static inline float _mm256_reduce_add_ps( const __m256 a ) {

   // Do it with SSE instructions
   __m128 r = _mm_add_ps( _mm256_extractf128_ps(a, 1), _mm256_castps256_ps128(a));
   r = _mm_hadd_ps( r, r );
   r = _mm_hadd_ps( r, r );
   return _mm_cvtss_f32( r );
}



/*****************************************************************************************
_MM256_LOAD8v3_PS, _MM256_LOAD4v3_EPI32

Loads 8 x 3 elements vectors stored sequentially in memory into 3 vector
variables:

b = [x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7] 
-> _MM256_LOAD8v3_PS(v1,v2,v3,b)

-> v1 = [x0,x1,x2,x3,x4,x5,x6,x7]
-> v2 = [y0,y1,y2,y3,y4,y5,y6,y7]
-> v3 = [z0,z1,z2,z3,z4,z5,z6,z7]

********************************************************************************************/

#define _MM256_LOAD8v3_PS(v1, v2, v3, b) {                                               \
   register __m256 m03, m14, m25;                                                        \
   m03 = _mm256_castps128_ps256(_mm_load_ps(&b[0]));                                     \
   m14 = _mm256_castps128_ps256(_mm_load_ps(&b[4]));                                     \
   m25 = _mm256_castps128_ps256(_mm_load_ps(&b[8]));                                     \
   m03 = _mm256_insertf128_ps(m03 ,_mm_load_ps(&b[12]),1);                               \
   m14 = _mm256_insertf128_ps(m14 ,_mm_load_ps(&b[16]),1);                               \
   m25 = _mm256_insertf128_ps(m25 ,_mm_load_ps(&b[20]),1);                               \
   __m256 xy = _mm256_shuffle_ps(m14, m25, _MM_SHUFFLE( 2,1,3,2));                       \
   __m256 yz = _mm256_shuffle_ps(m03, m14, _MM_SHUFFLE( 1,0,2,1));                       \
   v1 = _mm256_shuffle_ps(m03, xy , _MM_SHUFFLE( 2,0,3,0));                              \
   v2 = _mm256_shuffle_ps(yz , xy , _MM_SHUFFLE( 3,1,2,0));                              \
   v3 = _mm256_shuffle_ps(yz , m25, _MM_SHUFFLE( 3,0,3,1));                              \
}

#define _MM256_LOAD8v3_EPI32(v1, v2, v3, b) {                                            \
   register __m256 m03, m14, m25;                                                        \
   const float* buf = (float *) b;                                                       \
   m03 = _mm256_castps128_ps256(_mm_load_ps(&buf[0]));                                   \
   m14 = _mm256_castps128_ps256(_mm_load_ps(&buf[4]));                                   \
   m25 = _mm256_castps128_ps256(_mm_load_ps(&buf[8]));                                   \
   m03 = _mm256_insertf128_ps(m03 ,_mm_load_ps(&buf[12]),1);                             \
   m14 = _mm256_insertf128_ps(m14 ,_mm_load_ps(&buf[16]),1);                             \
   m25 = _mm256_insertf128_ps(m25 ,_mm_load_ps(&buf[20]),1);                             \
   __m256 xy = _mm256_shuffle_ps(m14, m25, _MM_SHUFFLE( 2,1,3,2));                       \
   __m256 yz = _mm256_shuffle_ps(m03, m14, _MM_SHUFFLE( 1,0,2,1));                       \
   v1 = _mm256_castps_si256( _mm256_shuffle_ps(m03, xy , _MM_SHUFFLE( 2,0,3,0)));        \
   v2 = _mm256_castps_si256( _mm256_shuffle_ps(yz , xy , _MM_SHUFFLE( 3,1,2,0)));        \
   v3 = _mm256_castps_si256( _mm256_shuffle_ps(yz , m25, _MM_SHUFFLE( 3,0,3,1)));        \
}


/*****************************************************************************************
_MM256_STORE8v3_PS, _MM256_STORE8v3_EPI32

Stores 3 vector variables into memory putting each vector element 
sequentially:

v1 = [x0,x1,x2,x3,x4,x5,x6,x7]
v2 = [y0,y1,y2,y3,y4,y5,y6,y7]
v3 = [z0,z1,z2,z3,z4,z5,z6,z7]

-> _MM_STORE8v3_PS(b, v1, v2, v3, b)

-> b = [x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7] 

*****************************************************************************************/

#define _MM256_STORE8v3_PS(b, x, y, z) {       							                 \
  __m256 rxy = _mm256_shuffle_ps(x,y, _MM_SHUFFLE(2,0,2,0));                             \
  __m256 ryz = _mm256_shuffle_ps(y,z, _MM_SHUFFLE(3,1,3,1));                             \
  __m256 rzx = _mm256_shuffle_ps(z,x, _MM_SHUFFLE(3,1,2,0));                             \
  __m256 r03 = _mm256_shuffle_ps(rxy, rzx, _MM_SHUFFLE(2,0,2,0));                        \
  __m256 r14 = _mm256_shuffle_ps(ryz, rxy, _MM_SHUFFLE(3,1,2,0));                        \
  __m256 r25 = _mm256_shuffle_ps(rzx, ryz, _MM_SHUFFLE(3,1,3,1));                        \
                                                                                         \
  _mm_store_ps( &b[ 0], _mm256_castps256_ps128( r03 ) );                                 \
  _mm_store_ps( &b[ 4], _mm256_castps256_ps128( r14 ) );                                 \
  _mm_store_ps( &b[ 8], _mm256_castps256_ps128( r25 ) );                                 \
  _mm_store_ps( &b[12], _mm256_extractf128_ps( r03 ,1) );                                \
  _mm_store_ps( &b[16], _mm256_extractf128_ps( r14 ,1) );                                \
  _mm_store_ps( &b[20], _mm256_extractf128_ps( r25 ,1) );                                \
}

#define _MM256_STORE8v3_EPI32(b, x, y, z) {	         	  				                 \
  __m256 rxy = _mm256_shuffle_ps(_mm256_castsi256_ps(x),                                 \
                                 _mm256_castsi256_ps(y), _MM_SHUFFLE(2,0,2,0));          \
  __m256 ryz = _mm256_shuffle_ps(_mm256_castsi256_ps(y),                                 \
                                 _mm256_castsi256_ps(z), _MM_SHUFFLE(3,1,3,1));          \
  __m256 rzx = _mm256_shuffle_ps(_mm256_castsi256_ps(z),                                 \
                                 _mm256_castsi256_ps(x), _MM_SHUFFLE(3,1,2,0));          \
  __m256 r03 = _mm256_shuffle_ps(rxy, rzx, _MM_SHUFFLE(2,0,2,0));                        \
  __m256 r14 = _mm256_shuffle_ps(ryz, rxy, _MM_SHUFFLE(3,1,2,0));                        \
  __m256 r25 = _mm256_shuffle_ps(rzx, ryz, _MM_SHUFFLE(3,1,3,1));                        \
                                                                                         \
  _mm_store_ps( (float*)&b[ 0], _mm256_castps256_ps128( r03 ) );                         \
  _mm_store_ps( (float*)&b[ 4], _mm256_castps256_ps128( r14 ) );                         \
  _mm_store_ps( (float*)&b[ 8], _mm256_castps256_ps128( r25 ) );                         \
  _mm_store_ps( (float*)&b[12], _mm256_extractf128_ps( r03 ,1) );                        \
  _mm_store_ps( (float*)&b[16], _mm256_extractf128_ps( r14 ,1) );                        \
  _mm_store_ps( (float*)&b[20], _mm256_extractf128_ps( r25 ,1) );                        \
}

/*****************************************************************************************
_MM256_LOAD8v2_PS, _MM256_LOAD8v2_EPI32

Loads 8 x 2 elements vectors stored sequentially in memory into 2 vector
variables:

b = [x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7]
-> _MM_LOAD8v2_PS( v1, v2, b )

-> v1 = [x0,x1,x2,x3,x4,x5,x6,x7]
-> v2 = [y0,y1,y2,y3,y4,y5,y6,y7]

*****************************************************************************************/

#define _MM256_LOAD8v2_PS(v1, v2, b) {                                                      \
  register __m256 m02, m13;                                                              \
  m02 = _mm256_castps128_ps256(_mm_load_ps(&b[0]));                                      \
  m13 = _mm256_castps128_ps256(_mm_load_ps(&b[4]));                                      \
  m02 = _mm256_insertf128_ps(m02 ,_mm_load_ps(&b[8]),1);                                 \
  m13 = _mm256_insertf128_ps(m13 ,_mm_load_ps(&b[12]),1);                                \
  (v1) = _mm256_shuffle_ps( m02, m13, _MM_SHUFFLE( 2, 0, 2, 0 ) );                       \
  (v2) = _mm256_shuffle_ps( m02, m13, _MM_SHUFFLE( 3, 1, 3, 1 ) );                       \
}


#define _MM256_LOAD8v2_EPI32(v1, v2, b) {                                 	                 \
  register __m256 m02, m13;                                                              \
  m02 = _mm256_castps128_ps256(_mm_load_ps((float*) &b[0]));                             \
  m13 = _mm256_castps128_ps256(_mm_load_ps((float*) &b[4]));                             \
  m02 = _mm256_insertf128_ps(m02 ,_mm_load_ps((float*) &b[8]),1);                        \
  m13 = _mm256_insertf128_ps(m13 ,_mm_load_ps((float*) &b[12]),1);                       \
  (v1) = _mm256_castps_si256( _mm256_shuffle_ps( m02, m13, _MM_SHUFFLE( 2, 0, 2, 0 ) )); \
  (v2) = _mm256_castps_si256( _mm256_shuffle_ps( m02, m13, _MM_SHUFFLE( 3, 1, 3, 1 ) )); \
}


/*****************************************************************************************
_MM256_STORE8v2_PS, _MM256_STORE8v2_EPI32

Stores 2 vector variables into memory putting each vector element 
sequentially:

v1 = [x0,x1,x2,x3,x4,x5,x6,x7], v2 = [y0,y1,y2,y3,y4,y5,y6,y7]

-> _MM256_STORE8v2_PS(b, v1, v2)

-> b = [x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7]



*****************************************************************************************/

#define _MM256_STORE8v2_PS(buf, x, y) {                                                     \
  register __m256 r02, r13;                                                              \
  r02 = _mm256_unpacklo_ps( x, y );                                                      \
  r13 = _mm256_unpackhi_ps( x, y );                                                      \
  _mm_store_ps( &buf[ 0], _mm256_castps256_ps128( r02 ) );                               \
  _mm_store_ps( &buf[ 4], _mm256_castps256_ps128( r13 ) );                               \
  _mm_store_ps( &buf[ 8], _mm256_extractf128_ps( r02 ,1) );                              \
  _mm_store_ps( &buf[12], _mm256_extractf128_ps( r13 ,1) );                              \
}

#define _MM256_STORE8v2_EPI32(buf, x, y) {                                                  \
  register __m256 r02, r13;                                                              \
  r02 = _mm256_unpacklo_ps( _mm256_castsi256_ps(x), _mm256_castsi256_ps(y) );            \
  r13 = _mm256_unpackhi_ps( _mm256_castsi256_ps(x), _mm256_castsi256_ps(y) );            \
  _mm_store_ps( (float *) &buf[ 0], _mm256_castps256_ps128( r02 ) );                     \
  _mm_store_ps( (float *) &buf[ 4], _mm256_castps256_ps128( r13 ) );                     \
  _mm_store_ps( (float *) &buf[ 8], _mm256_extractf128_ps( r02 ,1) );                    \
  _mm_store_ps( (float *) &buf[12], _mm256_extractf128_ps( r13 ,1) );                    \
}

/*****************************************************************************************

Convert 2 x 4 double precision vectors to / from 1 x 8 single precision vector

*****************************************************************************************/

#define _MM256_CVTPS_PD( a, b, f ) {                              \
  (a) = _mm256_cvtps_pd( _mm256_castps256_ps128( f ) );           \
  (b) = _mm256_cvtps_pd( _mm256_extractf128_ps( f , 1 ) );        \
}

#define _MM256_CVTPD_PS( f, a, b ) {                              \
  (f) = _mm256_castps128_ps256( _mm256_cvtpd_ps( a ) );           \
  (f) = _mm256_insertf128_ps( (f), _mm256_cvtpd_ps( b ), 1 );     \
}


/*****************************************************************************************
_MM256_LOAD4v3_PD

Loads 4 x 3 elements vectors stored sequentially in memory into 3 vector
variables:

b = [x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4] -> _MM_LOA42v3_PD(v1,v2,v3,b)

-> v1 = [x1,x2,x3,x4]
-> v2 = [y1,y2,y3,y4]
-> v3 = [z1,z2,z3,z4]

*****************************************************************************************/

#define _MM256_LOAD4v3_PD(v1, v2, v3, buf) { \
   register __m256d m03, m14, m25; \
   m03 = _mm256_castpd128_pd256(_mm_load_pd(&buf[0])); \
   m14 = _mm256_castpd128_pd256(_mm_load_pd(&buf[2])); \
   m25 = _mm256_castpd128_pd256(_mm_load_pd(&buf[4])); \
   m03 = _mm256_insertf128_pd(m03 ,_mm_load_pd(&buf[6]),1); \
   m14 = _mm256_insertf128_pd(m14 ,_mm_load_pd(&buf[8]),1); \
   m25 = _mm256_insertf128_pd(m25 ,_mm_load_pd(&buf[10]),1); \
   (v1) = _mm256_shuffle_pd(m03, m14, 0x0A); \
   (v2) = _mm256_shuffle_pd(m03, m25, 0x05); \
   (v3) = _mm256_shuffle_pd(m14, m25, 0x0A); \
}

#define _MM256_STORE4v3_PD(buf, v1, v2, v3) { \
   register __m256d r03, r14, r25; \
   r03 = _mm256_unpacklo_pd( (v1), (v2) ); \
   r14 = _mm256_shuffle_pd( (v3), (v1), 0x0A ); \
   r25 = _mm256_unpackhi_pd( (v2), (v3) ); \
   _mm_store_pd( &buf[ 0], _mm256_castpd256_pd128( r03 ) ); \
   _mm_store_pd( &buf[ 2], _mm256_castpd256_pd128( r14 ) ); \
   _mm_store_pd( &buf[ 4], _mm256_castpd256_pd128( r25 ) ); \
   _mm_store_pd( &buf[ 6], _mm256_extractf128_pd( r03, 1 ) ); \
   _mm_store_pd( &buf[ 8], _mm256_extractf128_pd( r14, 1 ) ); \
   _mm_store_pd( &buf[10], _mm256_extractf128_pd( r25, 1 ) ); \
}

#define _MM256_LOAD4v2_PD(v1, v2, buf) { \
   register __m256d m02, m13; \
   m02 = _mm256_castpd128_pd256(_mm_load_pd(&buf[0])); \
   m13 = _mm256_castpd128_pd256(_mm_load_pd(&buf[2])); \
   m02 = _mm256_insertf128_pd(m02, _mm_load_pd(&buf[4]),1); \
   m13 = _mm256_insertf128_pd(m13, _mm_load_pd(&buf[6]),1); \
   (v1) = _mm256_unpacklo_pd( m02, m13 ); \
   (v2) = _mm256_unpackhi_pd( m02, m13 ); \
}

#define _MM256_STORE4v2_PD(buf, v1, v2) { \
   register __m256d r02, r13; \
   r02 = _mm256_unpacklo_pd( (v1), (v2) ); \
   r13 = _mm256_unpackhi_pd( (v1), (v2) ); \
   _mm_store_pd( &buf[ 0], _mm256_castpd256_pd128( r02 ) ); \
   _mm_store_pd( &buf[ 2], _mm256_castpd256_pd128( r13 ) ); \
   _mm_store_pd( &buf[ 4], _mm256_extractf128_pd( r02, 1 ) ); \
   _mm_store_pd( &buf[ 6], _mm256_extractf128_pd( r13, 1 ) ); \
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

Utilities for printing / testing vector values

*****************************************************************************************/

static inline void printfv( char* s, __m256 vec )
{
  int i;
  DECLARE_ALIGNED_32(float v[8]);
  
  _mm256_store_ps( v, vec );
  
  printf(" %s = %g", s, v[0] );
  for(i=1; i<8; i++) printf(", %g", v[i]);
  printf("\n");  
}

static inline void printiv( char* s, __m256i vec )
{
  int i;
  DECLARE_ALIGNED_32(int v[8]);
  
  _mm256_store_si256( ( __m256i *) v, vec );
  
  printf(" %s = %i", s, v[0] );
  for(i=1; i<8; i++) printf(", %i", v[i]);
  printf("\n");  
}

static inline void printiv128( char* s, __m128i vec )
{
  int i;
  DECLARE_ALIGNED_16(int v[4]);
  
  _mm_store_si128( ( __m128i *) v, vec );
  
  printf(" %s = %i", s, v[0] );
  for(i=1; i<4; i++) printf(", %i", v[i]);
  printf("\n");  
}

static inline void printdv( char* s, __m256d vec )
{
  int i;
  DECLARE_ALIGNED_32(double v[4]);
  
  _mm256_store_pd( v, vec );
  
  printf(" %s = %.6f", s, v[0] );
  for(i=1; i<4; i++) printf(", %.6f", v[i]);
  printf("\n");  
}

#include <math.h>
static inline int testfv( __m256 vec)
{
  DECLARE_ALIGNED_32(float v[8]);
  
  _mm256_store_ps( v, vec );
  int i;
  for(i = 0; i < 8; i++) {
    if ( isnan( v[i] ) || isinf( v[i] ) ) return -1;
  } 
  return 0;
}

#define TEST_VECTOR( a ) {	                                                             \
  if ( testfv(a) ) {                                                                     \
     DECLARE_ALIGNED_32(float v[8]);                                                     \
     _mm256_store_ps( v, a );                                                            \
 	 printf(" %s:%d %s = %.3f", __FILE__, __LINE__, #a, v[0] );                          \
	 for(i=1; i<8; i++) printf(", %.3f", v[i]);                                          \
	 printf("\n");                                                                       \
  }                                                                                      \
}


/* End of definitions for Intel AVX */

#endif

/****************************************************************************************/
/* SIMD utilites.                                                                       */
/* Definitions for Intel MIC                                                            */
/****************************************************************************************/

#ifndef _VECTOR_MIC_H
#define _VECTOR_MIC_H


#include <stdio.h>
#include <immintrin.h>

#if defined(__GNUC__)
  /* gcc */
  #define DECLARE_ALIGNED_64(v)       v __attribute__ ((aligned (64)))
#else
  /* Intel icc */
  #define DECLARE_ALIGNED_64(v)       __declspec(align(64)) v
#endif

#include <stdint.h>

#define IS_ALIGNED_64(addr) (((uintptr_t)addr & 0x3F) == 0)
#define CHECK_ALIGN(addr) { if ((uintptr_t)addr & 0x3F) fprintf( stderr, "%s is not 64 aligned.\n", #addr ); } 

#if defined( PRECISION_SINGLE )
#define VEC_WIDTH 16
#elif defined( PRECISION_DOUBLE )
#define VEC_WIDTH 8
#else
#error PRECISION_SINGLE or PRECISION_DOUBLE must be defined
#endif 

typedef union vecDouble {double v[8]; __m512d v8;} dvec;
typedef union vecFloat  {float v[16]; __m512  v16;} fvec;
typedef union vecInt    {int v[16];   __m512i v16;} ivec;


/*****************************************************************************************
  Unaligned store / load.
*****************************************************************************************/

#define  _MM512_STOREU_PS(buf, v1) {  \
  _mm512_packstorelo_ps( (void*)(buf),    (v1)); \
  _mm512_packstorehi_ps( (void*)(buf+16), (v1)); \
}

#define  _MM512_LOADU_PS(v1, buf) {  \
  (v1) = _mm512_loadunpacklo_ps( (v1), (void*)(buf) ); \
  (v1) = _mm512_loadunpackhi_ps( (v1), (void*)(buf+16) ); \
}

#define  _MM512_STOREU_EPI32(buf, v1) {  \
  _mm512_packstorelo_epi32( (void*)(buf),    (v1)); \
  _mm512_packstorehi_epi32( (void*)(buf+16), (v1)); \
}

#define  _MM512_LOADU_EPI32(v1, buf) {  \
  (v1) = _mm512_loadunpacklo_epi32( (v1), (void*)(buf) ); \
  (v1) = _mm512_loadunpackhi_epi32( (v1), (void*)(buf+16) ); \
}

/*****************************************************************************************
Transpose ops: 16x3, 3x16


********************************************************************************************/

#define _MM512_TRANSPOSE16v3_PS(a, b, c) { \
  const __m512i perm0 = _mm512_set_epi32(14,11,8,5,2,13,10,7,4,1,15,12,9,6,3,0); \
  __m512i t0 = _mm512_permutevar_epi32( perm0, _mm512_castps_si512(a) ); \
  __m512i t1 = _mm512_permutevar_epi32( perm0, _mm512_castps_si512(b) ); \
  __m512i t2 = _mm512_permutevar_epi32( perm0, _mm512_castps_si512(c) ); \
  __m512i vx = _mm512_mask_alignr_epi32( t0, 0x07C0, t1, t1, 5 ); \
  __m512i vy = _mm512_mask_alignr_epi32( t2, 0x001F, t0, t0, 6 ); \
  __m512i vz = _mm512_mask_alignr_epi32( _mm512_alignr_epi32( t0, t0, 11 ), 0x03E0,t1, t1, 1 ); \
  a = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( vx, 0xF800, t2, t2, 11 ) ); \
  b = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( vy, 0x07E0, t1, t1, 11 ) ); \
  c = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( vz, 0xFC00, t2, t2,  6 ) ); \
}

#define _MM512_TRANSPOSE16v3_EPI32(a, b, c) { \
  const __m512i perm0 = _mm512_set_epi32(14,11,8,5,2,13,10,7,4,1,15,12,9,6,3,0); \
  __m512i t0 = _mm512_permutevar_epi32( perm0, a ); \
  __m512i t1 = _mm512_permutevar_epi32( perm0, b ); \
  __m512i t2 = _mm512_permutevar_epi32( perm0, c ); \
  a = _mm512_mask_alignr_epi32( t0, 0x07C0, t1, t1, 5 ); \
  b = _mm512_mask_alignr_epi32( t2, 0x001F, t0, t0, 6 ); \
  c = _mm512_mask_alignr_epi32( _mm512_alignr_epi32( t0, t0, 11 ), 0x03E0,t1, t1, 1 ); \
  a = _mm512_mask_alignr_epi32( a, 0xF800, t2, t2, 11 ); \
  b = _mm512_mask_alignr_epi32( b, 0x07E0, t1, t1, 11 ); \
  c = _mm512_mask_alignr_epi32( c, 0xFC00, t2, t2,  6 ); \
}

#define _MM512_TRANSPOSE3v16_PS( a, b, c) { \
  const __m512i perm0 =  _mm512_set_epi32( 5,10,15, 4, 9,14, 3, 8,13, 2, 7,12, 1, 6,11, 0); \
  __m512i vx = _mm512_permutevar_epi32( perm0, _mm512_castps_si512( a ) ); \
  __m512i vy = _mm512_permutevar_epi32( perm0, _mm512_castps_si512( b ) ); \
  __m512i vz = _mm512_permutevar_epi32( perm0, _mm512_castps_si512( c ) ); \
  vy = _mm512_alignr_epi32( vy, vy, 15 ); \
  vz = _mm512_alignr_epi32( vz, vz, 14 ); \
  a = _mm512_castsi512_ps(_mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x2492, vy ), 0x4924, vz )); \
  b = _mm512_castsi512_ps(_mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x9249, vy ), 0x2492, vz )); \
  c = _mm512_castsi512_ps(_mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x4924, vy ), 0x9249, vz )); \
}

#define _MM512_TRANSPOSE3v16_EPI32( a, b, c) { \
  const __m512i perm0 =  _mm512_set_epi32( 5,10,15, 4, 9,14, 3, 8,13, 2, 7,12, 1, 6,11, 0); \
  __m512i vx = _mm512_permutevar_epi32( perm0, a ); \
  __m512i vy = _mm512_permutevar_epi32( perm0, b ); \
  __m512i vz = _mm512_permutevar_epi32( perm0, c ); \
  vy = _mm512_alignr_epi32( vy, vy, 15 ); \
  vz = _mm512_alignr_epi32( vz, vz, 14 ); \
  a = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x2492, vy ), 0x4924, vz ); \
  b = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x9249, vy ), 0x2492, vz ); \
  c = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x4924, vy ), 0x9249, vz ); \
}

#define _MM512_TRANSPOSE16v2_PS( a, b ) { \
  const __m512i perm =  _mm512_set_epi32(15,13,11, 9, 7, 5, 3, 1,14,12,10, 8, 6, 4, 2, 0); \
  __m512i t0 = _mm512_permutevar_epi32( perm, _mm512_castps_si512( a ) ); \
  __m512i t1 = _mm512_permutevar_epi32( perm, _mm512_castps_si512( b ) ); \
  a = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( t0, 0xFF00, t1, t1, 8 ) ); \
  b = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( t1, 0x00FF, t0, t0, 8 ) ); \
}

#define _MM512_TRANSPOSE16v2_EPI32( a, b ) { \
  const __m512i perm =  _mm512_set_epi32(15,13,11, 9, 7, 5, 3, 1,14,12,10, 8, 6, 4, 2, 0); \
  __m512i t0 = _mm512_permutevar_epi32( perm, a ); \
  __m512i t1 = _mm512_permutevar_epi32( perm, b ); \
  a = _mm512_mask_alignr_epi32( t0, 0xFF00, t1, t1, 8 ); \
  b = _mm512_mask_alignr_epi32( t1, 0x00FF, t0, t0, 8 ); \
}

#define _MM512_TRANSPOSE2v16_PS( a, b ) {  \
  const __m512i perm0 =  _mm512_set_epi32(15, 7,14, 6,13, 5,12, 4,11, 3,10, 2, 9, 1, 8, 0); \
  __m512i t0 = _mm512_permutevar_epi32( perm0, _mm512_castps_si512(a) ); \
  __m512i t1 = _mm512_permutevar_epi32( perm0, _mm512_castps_si512(b) ); \
  a = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( t0, 0xAAAA, t1, t1, 15 ) ); \
  b = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( t1, 0x5555, t0, t0, 1 ) ); \
}

#define _MM512_TRANSPOSE2v16_EPI32( a, b ) {  \
  const __m512i perm0 =  _mm512_set_epi32(15, 7,14, 6,13, 5,12, 4,11, 3,10, 2, 9, 1, 8, 0); \
  __m512i t0 = _mm512_permutevar_epi32( perm0, a ); \
  __m512i t1 = _mm512_permutevar_epi32( perm0, b ); \
  a = ( _mm512_mask_alignr_epi32( t0, 0xAAAA, t1, t1, 15 ) ); \
  b = ( _mm512_mask_alignr_epi32( t1, 0x5555, t0, t0, 1 ) ); \
}


/*****************************************************************************************
_MM512_LOAD16v3_PS, _MM512_LOAD16v3_EPI32

Loads 16 x 3 elements vectors stored sequentially in memory into 3 vector
variables:

b = [x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,
     x8,y8,z8,x9,y9,z9,xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD,xE,yE,zE,xF,yF,zF] 
-> _MM512_LOAD16v3_PS(v1,v2,v3,b)

-> v1 = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xA,xB,xC,xD,xE,xF]
-> v2 = [y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,yA,yB,yC,yD,yE,yF]
-> v3 = [z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,zA,zB,zC,zD,zE,zF]

********************************************************************************************/

#define _MM512_LOAD16v3_PS(a, b, c, buf) { \
  const __m512i perm0 = _mm512_set_epi32(14,11,8,5,2,13,10,7,4,1,15,12,9,6,3,0); \
  __m512i t0 = _mm512_permutevar_epi32( perm0, _mm512_load_epi32(&buf[ 0]) ); \
  __m512i t1 = _mm512_permutevar_epi32( perm0, _mm512_load_epi32(&buf[16]) ); \
  __m512i t2 = _mm512_permutevar_epi32( perm0, _mm512_load_epi32(&buf[32]) ); \
  __m512i vx = _mm512_mask_alignr_epi32( t0, 0x07C0, t1, t1, 5 ); \
  __m512i vy = _mm512_mask_alignr_epi32( t2, 0x001F, t0, t0, 6 ); \
  __m512i vz = _mm512_mask_alignr_epi32( _mm512_alignr_epi32( t0, t0, 11 ), 0x03E0,t1, t1, 1 ); \
  a = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( vx, 0xF800, t2, t2, 11 ) ); \
  b = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( vy, 0x07E0, t1, t1, 11 ) ); \
  c = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( vz, 0xFC00, t2, t2,  6 ) ); \
}


#define _MM512_LOAD16v3_EPI32(a, b, c, buf) { \
  const __m512i perm0 = _mm512_set_epi32(14,11,8,5,2,13,10,7,4,1,15,12,9,6,3,0); \
  __m512i t0 = _mm512_permutevar_epi32( perm0, _mm512_load_epi32(&buf[ 0]) ); \
  __m512i t1 = _mm512_permutevar_epi32( perm0, _mm512_load_epi32(&buf[16]) ); \
  __m512i t2 = _mm512_permutevar_epi32( perm0, _mm512_load_epi32(&buf[32]) ); \
  a = _mm512_mask_alignr_epi32( t0, 0x07C0, t1, t1, 5 ); \
  b = _mm512_mask_alignr_epi32( t2, 0x001F, t0, t0, 6 ); \
  c = _mm512_mask_alignr_epi32( _mm512_alignr_epi32( t0, t0, 11 ), 0x03E0,t1, t1, 1 ); \
  a = _mm512_mask_alignr_epi32( a, 0xF800, t2, t2, 11 ); \
  b = _mm512_mask_alignr_epi32( b, 0x07E0, t1, t1, 11 ); \
  c = _mm512_mask_alignr_epi32( c, 0xFC00, t2, t2,  6 ); \
}

/*
// Use gather instructions instead - slower
#define _MM512_LOAD16v3_PS(a, b, c, buf) { \                                              
   __m512i const stride3 = _mm512_set_epi32(45,42,39,36,33,30,27,24,21,18,15,12,9,6,3,0); \
   a = _mm512_i32gather_ps( stride3, &buf[0], _MM_SCALE_4 ); \
   b = _mm512_i32gather_ps( stride3, &buf[1], _MM_SCALE_4 ); \
   c = _mm512_i32gather_ps( stride3, &buf[2], _MM_SCALE_4 ); \
}

#define _MM512_LOAD16v3_EPI32(a, b, c, buf) { \                                              
   __m512i const stride3 = _mm512_set_epi32(45,42,39,36,33,30,27,24,21,18,15,12,9,6,3,0); \
   a = _mm512_i32gather_epi32( stride3, &buf[0], _MM_SCALE_4 ); \
   b = _mm512_i32gather_epi32( stride3, &buf[1], _MM_SCALE_4 ); \
   c = _mm512_i32gather_epi32( stride3, &buf[2], _MM_SCALE_4 ); \
}
*/

/*****************************************************************************************
_MM512_STORE16v3_PS, _MM512_STORE16v3_EPI32

Stores 3 vector variables into memory putting each vector element 
sequentially:

v1 = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xA,xB,xC,xD,xE,xF]
v2 = [y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,yA,yB,yC,yD,yE,yF]
v3 = [z0,z1,z2,z3,z4,z5,z6,z7,z8,z9,zA,zB,zC,zD,zE,zF]

-> _MM512_STORE16v3_PS(b, v1, v2, v3)

-> b = [x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,
        x8,y8,z8,x9,y9,z9,xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD,xE,yE,zE,xF,yF,zF] 

*****************************************************************************************/


#define _MM512_STORE16v3_PS(buf, a, b, c) { \
  const __m512i perm0 =  _mm512_set_epi32( 5,10,15, 4, 9,14, 3, 8,13, 2, 7,12, 1, 6,11, 0); \
  __m512i vx = _mm512_permutevar_epi32( perm0, _mm512_castps_si512( a ) ); \
  __m512i vy = _mm512_permutevar_epi32( perm0, _mm512_castps_si512( b ) ); \
  __m512i vz = _mm512_permutevar_epi32( perm0, _mm512_castps_si512( c ) ); \
  vy = _mm512_alignr_epi32( vy, vy, 15 ); \
  vz = _mm512_alignr_epi32( vz, vz, 14 ); \
  __m512i t0 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x2492, vy ), 0x4924, vz ); \
  __m512i t1 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x9249, vy ), 0x2492, vz ); \
  __m512i t2 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x4924, vy ), 0x9249, vz ); \
  _mm512_store_epi32( &buf[ 0], t0 ); \
  _mm512_store_epi32( &buf[16], t1 ); \
  _mm512_store_epi32( &buf[32], t2 ); \
}




#define _MM512_STORE16v3_EPI32(buf, a, b, c) { \
  const __m512i perm0 =  _mm512_set_epi32( 5,10,15, 4, 9,14, 3, 8,13, 2, 7,12, 1, 6,11, 0); \
  __m512i vx = _mm512_permutevar_epi32( perm0, a ); \
  __m512i vy = _mm512_permutevar_epi32( perm0, b ); \
  __m512i vz = _mm512_permutevar_epi32( perm0, c ); \
  vy = _mm512_alignr_epi32( vy, vy, 15 ); \
  vz = _mm512_alignr_epi32( vz, vz, 14 ); \
  __m512i t0 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x2492, vy ), 0x4924, vz ); \
  __m512i t1 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x9249, vy ), 0x2492, vz ); \
  __m512i t2 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x4924, vy ), 0x9249, vz ); \
  _mm512_store_epi32( &buf[ 0], t0 ); \
  _mm512_store_epi32( &buf[16], t1 ); \
  _mm512_store_epi32( &buf[32], t2 ); \
}

/*
#define _MM512_STORE16v3_PS(buf, a, b, c) { \
  const __m512i perm0 =  _mm512_set_epi32( 5,10,15, 4, 9,14, 3, 8,13, 2, 7,12, 1, 6,11, 0); \
  _mm_prefetch( (const char *) &buf[ 0], _MM_HINT_T0 ); \
  _mm_prefetch( (const char *) &buf[16], _MM_HINT_T0 ); \
  _mm_prefetch( (const char *) &buf[32], _MM_HINT_T0 ); \
  __m512i vx = _mm512_permutevar_epi32( perm0, _mm512_castps_si512( a ) ); \
  __m512i vy = _mm512_permutevar_epi32( perm0, _mm512_castps_si512( b ) ); \
  __m512i vz = _mm512_permutevar_epi32( perm0, _mm512_castps_si512( c ) ); \
  vy = _mm512_alignr_epi32( vy, vy, 15 ); \
  vz = _mm512_alignr_epi32( vz, vz, 14 ); \
  __m512i t0 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x2492, vy ), 0x4924, vz ); \   
  __m512i t1 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x9249, vy ), 0x2492, vz ); \   
  __m512i t2 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x4924, vy ), 0x9249, vz ); \   
  _mm512_store_epi32( &buf[ 0], t0 ); \
  _mm512_store_epi32( &buf[16], t1 ); \
  _mm512_store_epi32( &buf[32], t2 ); \
}

#define _MM512_STORE16v3_EPI32(buf, a, b, c) { \
  const __m512i perm0 =  _mm512_set_epi32( 5,10,15, 4, 9,14, 3, 8,13, 2, 7,12, 1, 6,11, 0); \
  _mm_prefetch( (const char *) &buf[ 0], _MM_HINT_T0 ); \
  _mm_prefetch( (const char *) &buf[16], _MM_HINT_T0 ); \
  _mm_prefetch( (const char *) &buf[32], _MM_HINT_T0 ); \
  __m512i vx = _mm512_permutevar_epi32( perm0, a ); \
  __m512i vy = _mm512_permutevar_epi32( perm0, b ); \
  __m512i vz = _mm512_permutevar_epi32( perm0, c ); \
  vy = _mm512_alignr_epi32( vy, vy, 15 ); \
  vz = _mm512_alignr_epi32( vz, vz, 14 ); \
  __m512i t0 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x2492, vy ), 0x4924, vz ); \   
  __m512i t1 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x9249, vy ), 0x2492, vz ); \   
  __m512i t2 = _mm512_mask_mov_epi32( _mm512_mask_mov_epi32( vx, 0x4924, vy ), 0x9249, vz ); \   
  _mm512_store_epi32( &buf[ 0], t0 ); \
  _mm512_store_epi32( &buf[16], t1 ); \
  _mm512_store_epi32( &buf[32], t2 ); \
}
*/

/*
// Use scatter instructions instead
#define _MM512_STORE16v3_PS(buf, a, b, c) { \
  __m512i const stride3 = _mm512_set_epi32(45,42,39,36,33,30,27,24,21,18,15,12,9,6,3,0); \
  _mm512_i32scatter_ps( &buf[0], stride3, a, _MM_SCALE_4 ); \
  _mm512_i32scatter_ps( &buf[1], stride3, b, _MM_SCALE_4 ); \
  _mm512_i32scatter_ps( &buf[2], stride3, c, _MM_SCALE_4 ); \
}

#define _MM512_STORE16v3_EPI32(buf, a, b, c) { \
  __m512i const stride3 = _mm512_set_epi32(45,42,39,36,33,30,27,24,21,18,15,12,9,6,3,0); \
  _mm512_i32scatter_epi32( &buf[0], stride3, a, _MM_SCALE_4 ); \
  _mm512_i32scatter_epi32( &buf[1], stride3, b, _MM_SCALE_4 ); \
  _mm512_i32scatter_epi32( &buf[2], stride3, c, _MM_SCALE_4 ); \
}
*/

/*****************************************************************************************
_MM512_LOAD16v2_PS, _MM512_LOAD16v2_EPI32

Loads 16 x 2 elements vectors stored sequentially in memory into 2 vector
variables:

b = [x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7
     x8,y8,x9,y9,xA,yA,xB,yB,xC,yC,xD,yD,xE,yE,xF,yF]
-> _MM512_LOAD16v2_PS( v1, v2, b )

-> v1 = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xA,xB,xC,xD,xE,xF]
-> v2 = [y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,yA,yB,yC,yD,yE,yF]

*****************************************************************************************/

#define _MM512_LOAD16v2_PS( a, b, buf ) { \
  const __m512i perm =  _mm512_set_epi32(15,13,11, 9, 7, 5, 3, 1,14,12,10, 8, 6, 4, 2, 0); \
  __m512i t0 = _mm512_permutevar_epi32( perm, _mm512_load_epi32(&buf[ 0]) ); \
  __m512i t1 = _mm512_permutevar_epi32( perm, _mm512_load_epi32(&buf[16]) ); \
  a = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( t0, 0xFF00, t1, t1, 8 ) ); \
  b = _mm512_castsi512_ps( _mm512_mask_alignr_epi32( t1, 0x00FF, t0, t0, 8 ) ); \
}

#define _MM512_LOAD16v2_EPI32( a, b, buf ) { \
  const __m512i perm =  _mm512_set_epi32(15,13,11, 9, 7, 5, 3, 1,14,12,10, 8, 6, 4, 2, 0); \
  __m512i t0 = _mm512_permutevar_epi32( perm, _mm512_load_epi32(&buf[ 0]) ); \
  __m512i t1 = _mm512_permutevar_epi32( perm, _mm512_load_epi32(&buf[16]) ); \
  a = _mm512_mask_alignr_epi32( t0, 0xFF00, t1, t1, 8 ); \
  b = _mm512_mask_alignr_epi32( t1, 0x00FF, t0, t0, 8 ); \
}

/*
// Use gather instructions instead
#define _MM512_LOAD16v2_PS(a, b, buf) { \                                              
   __m512i const stride2 = _mm512_set_epi32(30,28,26,24,22,20,18,16,14,12,10,8,6,4,2,0); \
   a = _mm512_i32gather_ps( stride2, &buf[0], _MM_SCALE_4 ); \
   b = _mm512_i32gather_ps( stride2, &buf[1], _MM_SCALE_4 ); \
}

#define _MM512_LOAD16v3_EPI32(a, b, buf) { \                                              
   __m512i const stride2 = _mm512_set_epi32(30,28,26,24,22,20,18,16,14,12,10,8,6,4,2,0); \
   a = _mm512_i32gather_epi32( stride2, &buf[0], _MM_SCALE_4 ); \
   b = _mm512_i32gather_epi32( stride2, &buf[1], _MM_SCALE_4 ); \
}
*/

/*****************************************************************************************
_MM512_STORE16v2_PS, _MM256_STORE16v2_EPI32

Stores 2 vector variables into memory putting each vector element 
sequentially:

v1 = [x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,xA,xB,xC,xD,xE,xF]
v2 = [y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,yA,yB,yC,yD,yE,yF]

-> _MM512_STORE8v2_PS(b, v1, v2)

-> b = [x0,y0,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7
        x8,y8,x9,y9,xA,yA,xB,yB,xC,yC,xD,yD,xE,yE,xF,yF]



*****************************************************************************************/

#define _MM512_STORE16v2_PS(buf, a, b) {  \
  const __m512i perm0 =  _mm512_set_epi32(15, 7,14, 6,13, 5,12, 4,11, 3,10, 2, 9, 1, 8, 0); \
  __m512i t0 = _mm512_permutevar_epi32( perm0, _mm512_castps_si512(a) ); \
  __m512i t1 = _mm512_permutevar_epi32( perm0, _mm512_castps_si512(b) ); \
  _mm512_store_epi32( &buf[ 0], _mm512_mask_alignr_epi32( t0, 0xAAAA, t1, t1, 15 ) ); \
  _mm512_store_epi32( &buf[16], _mm512_mask_alignr_epi32( t1, 0x5555, t0, t0, 1 ) ); \
}

#define _MM512_STORE16v2_EPI32(buf, a, b) {  \
  const __m512i perm0 =  _mm512_set_epi32(15, 7,14, 6,13, 5,12, 4,11, 3,10, 2, 9, 1, 8, 0); \
  __m512i t0 = _mm512_permutevar_epi32( perm0, a ); \
  __m512i t1 = _mm512_permutevar_epi32( perm0, b ); \
  _mm512_store_epi32( &buf[ 0], _mm512_mask_alignr_epi32( t0, 0xAAAA, t1, t1, 15 ) ); \
  _mm512_store_epi32( &buf[16], _mm512_mask_alignr_epi32( t1, 0x5555, t0, t0, 1 ) ); \
}

/*
// Use scatter instructions instead
#define _MM512_STORE16v2_PS(buf, a, b, c) { \
  __m512i const stride2 = _mm512_set_epi32(30,28,26,24,22,20,18,16,14,12,10,8,6,4,2,0); \
  _mm512_i32scatter_ps( &buf[0], stride2, v0, _MM_SCALE_4 ); \
  _mm512_i32scatter_ps( &buf[1], stride2, v1, _MM_SCALE_4 ); \
}

#define _MM512_STORE16v2_EPI32(buf, a, b, c) { \
  __m512i const stride2 = _mm512_set_epi32(30,28,26,24,22,20,18,16,14,12,10,8,6,4,2,0); \
  _mm512_i32scatter_epi32( &buf[0], stride2, v0, _MM_SCALE_4 ); \
  _mm512_i32scatter_epi32( &buf[1], stride2, v1, _MM_SCALE_4 ); \
}
*/


/*****************************************************************************************

Utilities for printing / testing vector values

*****************************************************************************************/

static inline void printfv( char* s, __m512 vec )
{
  int i;
  DECLARE_ALIGNED_64(float v[16]);
  
  _mm512_store_ps( v, vec );
  
  printf(" %s = %g", s, v[0] );
  for(i=1; i<16; i++) printf(", %g", v[i]);
  printf("\n");  
}

static inline void printiv( char* s, __m512i vec )
{
  int i;
  DECLARE_ALIGNED_64(int v[16]);
  
  _mm512_store_epi32( v, vec );
  
  printf(" %s = %i", s, v[0] );
  for(i=1; i<16; i++) printf(", %i", v[i]);
  printf("\n");  
}

#include <math.h>
static inline int testfv( __m512 vec)
{
  DECLARE_ALIGNED_64(float v[16]);
  
  _mm512_store_ps( v, vec );
  int i;
  for(i = 0; i < 8; i++) {
    if ( isnan( v[i] ) || isinf( v[i] ) ) return -1;
  } 
  return 0;
}

#define TEST_VECTOR( a ) {	                                                             \
  if ( testfv(a) ) {                                                                     \
     DECLARE_ALIGNED_64(float v[16]);                                                    \
     _mm512_store_ps( v, a );                                                            \
 	 printf(" %s:%d %s = %.3f", __FILE__, __LINE__, #a, v[0] );                          \
	 for(i=1; i<16; i++) printf(", %.3f", v[i]);                                         \
	 printf("\n");                                                                       \
  }                                                                                      \
}


/* End of definitions for Intel MIC */

#endif

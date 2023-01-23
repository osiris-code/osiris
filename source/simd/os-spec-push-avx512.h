#ifndef _OS_SPEC_PUSH_H_AVX512
#define _OS_SPEC_PUSH_H_AVX512

/* Size of particle buffer */
/* These must be a multiple of 16 */
#define p_cache_size_1D 64
#define p_cache_size_2D 464
#define p_cache_size_3D 160

#include "vector-avx512.h"

/*****************************************************************************************
t_split_buf1D

Buffer to hold virtual particles for current deposition in 1D.
*****************************************************************************************/

typedef struct {
  // 2 splits maximum
  DECLARE_ALIGNED_64( float x0[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_64( float x1[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_64( float  q[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_64( float vy[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_64( float vz[ 2 * p_cache_size_1D ] );

  DECLARE_ALIGNED_64( int   ix[ 2 * p_cache_size_1D ] );

  unsigned int np;
} t_split_buf1D;


// Because of particle splitting the stores are not always aligned
#define STOREU16P1D( buf, idx, vx0, vx1, vq, vvy, vvz, vix ) { \
  _mm512_storeu_ps( &buf.x0[ idx ], vx0); \
  _mm512_storeu_ps( &buf.x1[ idx ], vx1); \
  _mm512_storeu_ps(  &buf.q[ idx ],  vq); \
  _mm512_storeu_ps( &buf.vy[ idx ], vvy); \
  _mm512_storeu_ps( &buf.vz[ idx ], vvz); \
\
   _mm512_storeu_si512( &buf.ix[ idx ], vix); \
}

// For the current deposition routines the loads will always be aligned

#define LOAD16P1D( pbuf, idx, vx0, vx1, vq, vvy, vvz, vix ) { \
   vx0 = _mm512_load_ps( &(pbuf->x0[ idx ]) ); \
   vx1 = _mm512_load_ps( &(pbuf->x1[ idx ]) ); \
   vq  = _mm512_load_ps( &(pbuf->q[ idx ]) ); \
   vvy = _mm512_load_ps( &(pbuf->vy[ idx ]) ); \
   vvz = _mm512_load_ps( &(pbuf->vz[ idx ]) ); \
\
   vix = _mm512_load_epi32( (const __m512i *) &(pbuf->ix[ idx ]) ); \
}


/*****************************************************************************************
t_split_buf2D

Buffer to hold virtual particles for current deposition in 2D.
*****************************************************************************************/

typedef struct {
  // 3 splits maximum
  DECLARE_ALIGNED_64( float x0[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_64( float x1[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_64( float y0[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_64( float y1[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_64( float  q[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_64( float vz[ 3 * p_cache_size_2D ] );

  DECLARE_ALIGNED_64( int   ix[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_64( int   iy[ 3 * p_cache_size_2D ] );

  unsigned int np;
} t_split_buf2D;


// Because of particle splitting the stores are not always aligned
#define STOREU16P2D( buf, idx, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
   _mm512_storeu_ps( &buf.x0[ idx ], vx0); \
   _mm512_storeu_ps( &buf.x1[ idx ], vx1); \
   _mm512_storeu_ps( &buf.y0[ idx ], vy0); \
   _mm512_storeu_ps( &buf.y1[ idx ], vy1); \
   _mm512_storeu_ps(  &buf.q[ idx ],  vq); \
   _mm512_storeu_ps( &buf.vz[ idx ], vvz); \
\
   _mm512_storeu_si512( &buf.ix[ idx ], vix); \
   _mm512_storeu_si512( &buf.iy[ idx ], viy); \
}

// For the current deposition routines the loads will always be aligned

#define LOAD16P2D( pbuf, idx, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
   vx0 = _mm512_load_ps( &(pbuf->x0[ idx ]) ); \
   vx1 = _mm512_load_ps( &(pbuf->x1[ idx ]) ); \
   vy0 = _mm512_load_ps( &(pbuf->y0[ idx ]) ); \
   vy1 = _mm512_load_ps( &(pbuf->y1[ idx ]) ); \
   vq  = _mm512_load_ps( &(pbuf->q[ idx ]) ); \
   vvz = _mm512_load_ps( &(pbuf->vz[ idx ]) ); \
\
   vix = _mm512_load_epi32( (const __m512i *) &(pbuf->ix[ idx ]) ); \
   viy = _mm512_load_epi32( (const __m512i *) &(pbuf->iy[ idx ]) ); \
}

/*****************************************************************************************
t_split_buf3D

Buffer to hold virtual particles for current deposition in 3D.
*****************************************************************************************/

typedef struct {
  // 4 splits maximum
  DECLARE_ALIGNED_64( float x0[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_64( float x1[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_64( float y0[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_64( float y1[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_64( float z0[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_64( float z1[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_64( float  q[ 4 * p_cache_size_3D ] );

  DECLARE_ALIGNED_64( int   ix[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_64( int   iy[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_64( int   iz[ 4 * p_cache_size_3D ] );

  unsigned int np;
} t_split_buf3D;



// Because of particle splitting the stores are not always aligned



// Reference

#define STOREU16P3D( buf, idx, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { \
   _mm512_storeu_ps( &buf.x0[ idx ], vx0); \
   _mm512_storeu_ps( &buf.x1[ idx ], vx1); \
   _mm512_storeu_ps( &buf.y0[ idx ], vy0); \
   _mm512_storeu_ps( &buf.y1[ idx ], vy1); \
   _mm512_storeu_ps( &buf.z0[ idx ], vz0); \
   _mm512_storeu_ps( &buf.z1[ idx ], vz1); \
   _mm512_storeu_ps(  &buf.q[ idx ],  vq); \
   _mm512_storeu_si512( &buf.ix[ idx ], vix); \
   _mm512_storeu_si512( &buf.iy[ idx ], viy); \
   _mm512_storeu_si512( &buf.iz[ idx ], viz); \
}

/*

// Software prefetch

#define STOREU16P3D( buf, idx, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { \
_mm_prefetch( (void*)(&buf.y0[ idx ]), _MM_HINT_T0 ); \
_mm_prefetch( (void*)(&buf.y0[ idx ]) + 64, _MM_HINT_T0 ); \
   _mm512_storeu_ps( &buf.x0[ idx ], vx0); \
_mm_prefetch( (void*)(&buf.y1[ idx ]), _MM_HINT_T0 ); \
_mm_prefetch( (void*)(&buf.y1[ idx ]) + 64, _MM_HINT_T0 ); \
   _mm512_storeu_ps( &buf.x1[ idx ], vx1); \
_mm_prefetch( (void*)(&buf.z0[ idx ]), _MM_HINT_T0 ); \
_mm_prefetch( (void*)(&buf.z0[ idx ]) + 64, _MM_HINT_T0 ); \
   _mm512_storeu_ps( &buf.y0[ idx ], vy0); \
_mm_prefetch( (void*)(&buf.z1[ idx ]), _MM_HINT_T0 ); \
_mm_prefetch( (void*)(&buf.z1[ idx ]) + 64, _MM_HINT_T0 ); \
   _mm512_storeu_ps( &buf.y1[ idx ], vy1); \
_mm_prefetch( (void*)(&buf.q[ idx ]), _MM_HINT_T0 ); \
_mm_prefetch( (void*)(&buf.q[ idx ]) + 64, _MM_HINT_T0 ); \
   _mm512_storeu_ps( &buf.z0[ idx ], vz0); \
_mm_prefetch( (void*)(&buf.ix[ idx ]), _MM_HINT_T0 ); \
_mm_prefetch( (void*)(&buf.ix[ idx ]) + 64, _MM_HINT_T0 ); \
   _mm512_storeu_ps( &buf.z1[ idx ], vz1); \
_mm_prefetch( (void*)(&buf.iy[ idx ]), _MM_HINT_T0 ); \
_mm_prefetch( (void*)(&buf.iy[ idx ]) + 64, _MM_HINT_T0 ); \
   _mm512_storeu_ps(  &buf.q[ idx ],  vq); \
_mm_prefetch( (void*)(&buf.iz[ idx ]), _MM_HINT_T0 ); \
_mm_prefetch( (void*)(&buf.iz[ idx ]) + 64, _MM_HINT_T0 ); \
   _mm512_storeu_si512( &buf.ix[ idx ], vix); \
   _mm512_storeu_si512( &buf.iy[ idx ], viy); \
   _mm512_storeu_si512( &buf.iz[ idx ], viz); \
}
*/



// For the current deposition routines the loads will always be aligned

#define LOAD16P3D( pbuf, idx, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { \
   vx0 = _mm512_load_ps( &(pbuf->x0[ idx ]) ); \
   vx1 = _mm512_load_ps( &(pbuf->x1[ idx ]) ); \
   vy0 = _mm512_load_ps( &(pbuf->y0[ idx ]) ); \
   vy1 = _mm512_load_ps( &(pbuf->y1[ idx ]) ); \
   vz0 = _mm512_load_ps( &(pbuf->z0[ idx ]) ); \
   vz1 = _mm512_load_ps( &(pbuf->z1[ idx ]) ); \
   vq  = _mm512_load_ps( &(pbuf->q[ idx ]) ); \
\
   vix = _mm512_load_epi32( (const __m512i *) &(pbuf->ix[ idx ]) ); \
   viy = _mm512_load_epi32( (const __m512i *) &(pbuf->iy[ idx ]) ); \
   viz = _mm512_load_epi32( (const __m512i *) &(pbuf->iz[ idx ]) ); \
}

#endif

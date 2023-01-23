#ifndef _OS_SPEC_PUSH_H_SSE
#define _OS_SPEC_PUSH_H_SSE

/* Size of particle buffer */
/* for values between 32 and 8192 this has little impact on performance */
#define p_cache_size_1D 64
#define p_cache_size_2D 64
#define p_cache_size_3D 1024

#include "vector-sse.h"

/*****************************************************************************************
t_split_buf1D

Buffer to hold virtual particles for current deposition in 1D. 
*****************************************************************************************/

typedef struct { 
  // 2 splits maximum
  DECLARE_ALIGNED_16( float x0[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_16( float x1[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_16( float  q[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_16( float vy[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_16( float vz[ 2 * p_cache_size_1D ] );
  
  DECLARE_ALIGNED_16( int ix[ 2 * p_cache_size_1D ] );

  unsigned int np;
} t_split_buf1D;

#define STOREU4P1D( buf, idx, vx0, vx1, vq, vvy, vvz, vix ) { \
   _mm_storeu_ps( &buf.x0[ idx ], vx0); \
   _mm_storeu_ps( &buf.x1[ idx ], vx1); \
   _mm_storeu_ps(  &buf.q[ idx ],  vq); \
   _mm_storeu_ps( &buf.vy[ idx ], vvy); \
   _mm_storeu_ps( &buf.vz[ idx ], vvz); \
\
   _mm_storeu_si128( (__m128i *) &buf.ix[ idx ], vix); \
}

// For the current deposition routines the loads will always be aligned

#define LOAD4P1D( pbuf, idx, vx0, vx1, vq, vvy, vvz, vix ) { \
   vx0 = _mm_load_ps( &(pbuf->x0[ idx ]) ); \
   vx1 = _mm_load_ps( &(pbuf->x1[ idx ]) ); \
   vq  = _mm_load_ps( &(pbuf->q[ idx ]) ); \
   vvy = _mm_load_ps( &(pbuf->vy[ idx ]) ); \
   vvz = _mm_load_ps( &(pbuf->vz[ idx ]) ); \
\
   vix = _mm_load_si128( (const __m128i *) &(pbuf->ix[ idx ]) ); \
}


/*****************************************************************************************
t_split_buf2D

Buffer to hold virtual particles for current deposition in 2D. 
*****************************************************************************************/

typedef struct { 
  // 3 splits maximum
  DECLARE_ALIGNED_16( float x0[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_16( float x1[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_16( float y0[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_16( float y1[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_16( float  q[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_16( float vz[ 3 * p_cache_size_2D ] );
  
  DECLARE_ALIGNED_16( int ix[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_16( int iy[ 3 * p_cache_size_2D ] );

  unsigned int np;
} t_split_buf2D;

#define STOREU4P2D( buf, idx, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
   _mm_storeu_ps( &buf.x0[ idx ], vx0); \
   _mm_storeu_ps( &buf.x1[ idx ], vx1); \
   _mm_storeu_ps( &buf.y0[ idx ], vy0); \
   _mm_storeu_ps( &buf.y1[ idx ], vy1); \
   _mm_storeu_ps(  &buf.q[ idx ],  vq); \
   _mm_storeu_ps( &buf.vz[ idx ], vvz); \
\
   _mm_storeu_si128( (__m128i *) &buf.ix[ idx ], vix); \
   _mm_storeu_si128( (__m128i *) &buf.iy[ idx ], viy); \
}

// For the current deposition routines the loads will always be aligned

#define LOAD4P2D( pbuf, idx, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
   vx0 = _mm_load_ps( &(pbuf->x0[ idx ]) ); \
   vx1 = _mm_load_ps( &(pbuf->x1[ idx ]) ); \
   vy0 = _mm_load_ps( &(pbuf->y0[ idx ]) ); \
   vy1 = _mm_load_ps( &(pbuf->y1[ idx ]) ); \
   vq  = _mm_load_ps( &(pbuf->q[ idx ]) ); \
   vvz = _mm_load_ps( &(pbuf->vz[ idx ]) ); \
\
   vix = _mm_load_si128( (const __m128i *) &(pbuf->ix[ idx ]) ); \
   viy = _mm_load_si128( (const __m128i *) &(pbuf->iy[ idx ]) ); \
}

/*****************************************************************************************
t_split_buf3D

Buffer to hold virtual particles for current deposition in 3D. 
*****************************************************************************************/

typedef struct { 
  // 4 splits maximum
  DECLARE_ALIGNED_16( float x0[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_16( float x1[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_16( float y0[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_16( float y1[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_16( float z0[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_16( float z1[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_16( float  q[ 4 * p_cache_size_3D ] );
  
  DECLARE_ALIGNED_16( int ix[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_16( int iy[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_16( int iz[ 4 * p_cache_size_3D ] );

  unsigned int np;
} t_split_buf3D;

#define STOREU4P3D( buf, idx, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { \
   _mm_storeu_ps( &buf.x0[ idx ], vx0); \
   _mm_storeu_ps( &buf.x1[ idx ], vx1); \
   _mm_storeu_ps( &buf.y0[ idx ], vy0); \
   _mm_storeu_ps( &buf.y1[ idx ], vy1); \
   _mm_storeu_ps( &buf.z0[ idx ], vz0); \
   _mm_storeu_ps( &buf.z1[ idx ], vz1); \
   _mm_storeu_ps(  &buf.q[ idx ],  vq); \
   _mm_storeu_si128( (__m128i *) &buf.ix[ idx ], vix); \
   _mm_storeu_si128( (__m128i *) &buf.iy[ idx ], viy); \
   _mm_storeu_si128( (__m128i *) &buf.iz[ idx ], viz); \
}

#define LOAD4P3D( pbuf, idx, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { \
   vx0 = _mm_load_ps( &(pbuf->x0[ idx ]) ); \
   vx1 = _mm_load_ps( &(pbuf->x1[ idx ]) ); \
   vy0 = _mm_load_ps( &(pbuf->y0[ idx ]) ); \
   vy1 = _mm_load_ps( &(pbuf->y1[ idx ]) ); \
   vz0 = _mm_load_ps( &(pbuf->z0[ idx ]) ); \
   vz1 = _mm_load_ps( &(pbuf->z1[ idx ]) ); \
   vq  = _mm_load_ps( &(pbuf->q[ idx ]) ); \
\
   vix = _mm_load_si128( (const __m128i *) &(pbuf->ix[ idx ]) ); \
   viy = _mm_load_si128( (const __m128i *) &(pbuf->iy[ idx ]) ); \
   viz = _mm_load_si128( (const __m128i *) &(pbuf->iz[ idx ]) ); \
}

#endif

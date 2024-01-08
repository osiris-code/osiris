#ifndef _OS_SPEC_PUSH_H_BGQ
#define _OS_SPEC_PUSH_H_BGQ



/* Size of particle buffer */

#define p_cache_size_1D 32
#define p_cache_size_2D 32
#define p_cache_size_3D 512

#include "vector-bgq.h"

/*****************************************************************************************
LOADFLD4

Loads 4 field values corresponding to 4 particle positions into a vector
variable. 

*****************************************************************************************/


#define LOADFLD4( fp, shift ) \
 vec_perm( vec_perm( vec_lds( 0, &(fp[0])[shift] ), vec_lds( 0, &(fp[1])[shift] ), vec_gpci( 00415 ) ), \
           vec_perm( vec_lds( 0, &(fp[2])[shift] ), vec_lds( 0, &(fp[3])[shift] ), vec_gpci( 00415 ) ), \
           vec_gpci( 00145 ) )

/*****************************************************************************************
t_split_buf1D

Buffer to hold virtual particles for current deposition in 1D. 
*****************************************************************************************/

typedef struct { 
  // 2 splits maximum
  DECLARE_ALIGNED_32( double x0[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_32( double x1[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_32( double  q[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_32( double vy[ 2 * p_cache_size_1D ] );
  DECLARE_ALIGNED_32( double vz[ 2 * p_cache_size_1D ] );
  
  DECLARE_ALIGNED_32( int ix[ 2 * p_cache_size_1D ] );

  unsigned int np;
} t_split_buf1D;


#define STOREU4P1D( buf, idx, vx0, vx1, vq, vvy, vvz, vix ) { \
  buf.x0[idx  ] = vec_extract( vx0, 0 ); \
  buf.x0[idx+1] = vec_extract( vx0, 1 ); \
  buf.x0[idx+2] = vec_extract( vx0, 2 ); \
  buf.x0[idx+3] = vec_extract( vx0, 3 ); \
\
  buf.x1[idx  ] = vec_extract( vx1, 0 ); \
  buf.x1[idx+1] = vec_extract( vx1, 1 ); \
  buf.x1[idx+2] = vec_extract( vx1, 2 ); \
  buf.x1[idx+3] = vec_extract( vx1, 3 ); \
\
  buf.q[idx  ] = vec_extract( vq, 0 ); \
  buf.q[idx+1] = vec_extract( vq, 1 ); \
  buf.q[idx+2] = vec_extract( vq, 2 ); \
  buf.q[idx+3] = vec_extract( vq, 3 ); \
\
  buf.vy[idx  ] = vec_extract( vvy, 0 ); \
  buf.vy[idx+1] = vec_extract( vvy, 1 ); \
  buf.vy[idx+2] = vec_extract( vvy, 2 ); \
  buf.vy[idx+3] = vec_extract( vvy, 3 ); \
\
  buf.vz[idx  ] = vec_extract( vvz, 0 ); \
  buf.vz[idx+1] = vec_extract( vvz, 1 ); \
  buf.vz[idx+2] = vec_extract( vvz, 2 ); \
  buf.vz[idx+3] = vec_extract( vvz, 3 ); \
\
  buf.ix[idx]=vix[0]; buf.ix[idx+1]=vix[1]; buf.ix[idx+2]=vix[2]; buf.ix[idx+3]=vix[3]; \
\
}

// For the current deposition routines the loads will always be aligned

#define LOAD4P1D( pbuf, idx, vx0, vx1, vq, vvy, vvz, vix ) { \
  register vector _vix; \
\
  _vix = vec_ldiaa( 0x00, &(pbuf->ix[ idx ]) ); \
  vx0  =   vec_lda( 0x00, &(pbuf->x0[ idx ]) ); \
  vx1  =   vec_lda( 0x00, &(pbuf->x1[ idx ]) ); \
  vq   =   vec_lda( 0x00, &(pbuf->q[ idx ]) ); \
  vvy  =   vec_lda( 0x00, &(pbuf->vy[ idx ]) ); \
  vvz  =   vec_lda( 0x00, &(pbuf->vz[ idx ]) ); \
\
  vec_sta( _vix, 0x00, vix ); \
}


/*****************************************************************************************
t_split_buf2D

Buffer to hold virtual particles for current deposition in 2D. 
*****************************************************************************************/

typedef struct { 
  // 3 splits maximum
  DECLARE_ALIGNED_32( double x0[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_32( double x1[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_32( double y0[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_32( double y1[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_32( double  q[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_32( double vz[ 3 * p_cache_size_2D ] );
  
  DECLARE_ALIGNED_32( int ix[ 3 * p_cache_size_2D ] );
  DECLARE_ALIGNED_32( int iy[ 3 * p_cache_size_2D ] );

  unsigned int np;
} t_split_buf2D;

/* 
 BG/Q has no (single) instructions for unaligned memory access, so this is a bit more
 complicated than usual. I have found that scalar access is slightly faster than
 vector unaligned access (and much simpler) so we use that instead. Manual unrolling
 the loops also helps.
*/

#define STOREU4P2D( buf, idx, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
  buf.x0[idx  ] = vec_extract( vx0, 0 ); \
  buf.x0[idx+1] = vec_extract( vx0, 1 ); \
  buf.x0[idx+2] = vec_extract( vx0, 2 ); \
  buf.x0[idx+3] = vec_extract( vx0, 3 ); \
\
  buf.x1[idx  ] = vec_extract( vx1, 0 ); \
  buf.x1[idx+1] = vec_extract( vx1, 1 ); \
  buf.x1[idx+2] = vec_extract( vx1, 2 ); \
  buf.x1[idx+3] = vec_extract( vx1, 3 ); \
\
  buf.y0[idx  ] = vec_extract( vy0, 0 ); \
  buf.y0[idx+1] = vec_extract( vy0, 1 ); \
  buf.y0[idx+2] = vec_extract( vy0, 2 ); \
  buf.y0[idx+3] = vec_extract( vy0, 3 ); \
\
  buf.y1[idx  ] = vec_extract( vy1, 0 ); \
  buf.y1[idx+1] = vec_extract( vy1, 1 ); \
  buf.y1[idx+2] = vec_extract( vy1, 2 ); \
  buf.y1[idx+3] = vec_extract( vy1, 3 ); \
\
  buf.q[idx  ] = vec_extract( vq, 0 ); \
  buf.q[idx+1] = vec_extract( vq, 1 ); \
  buf.q[idx+2] = vec_extract( vq, 2 ); \
  buf.q[idx+3] = vec_extract( vq, 3 ); \
\
  buf.vz[idx  ] = vec_extract( vvz, 0 ); \
  buf.vz[idx+1] = vec_extract( vvz, 1 ); \
  buf.vz[idx+2] = vec_extract( vvz, 2 ); \
  buf.vz[idx+3] = vec_extract( vvz, 3 ); \
\
  buf.ix[idx]=vix[0]; buf.ix[idx+1]=vix[1]; buf.ix[idx+2]=vix[2]; buf.ix[idx+3]=vix[3]; \
  buf.iy[idx]=viy[0]; buf.iy[idx+1]=viy[1]; buf.iy[idx+2]=viy[2]; buf.iy[idx+3]=viy[3]; \
\
}

#if 0

// This is slower than the serial memory access version above

#define STOREU4P2D( buf, idx, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
  register vector _perm, _ins, _va, _vb, _vc, _vix, _viy; \
  \
  unsigned int _idxa = idx & 0xFFFFFFFCL; \
  \
  _vix = vec_ldiza( 0x00, vix ); \
  _viy = vec_ldiza( 0x00, viy ); \
  \
  _perm = vec_lvsr( 0x00, &( buf.x0[ idx ] ) ); \
  _ins  = vec_perm( vec_cmplt(vx0,vx0), vec_cmpeq(vx0,vx0), _perm ); \
  \
  /* vx0 */ \
  _va = vec_lda( 0x00, &( buf.x0[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.x0[ _idxa ] ) ); \
  _vc = vec_perm( vx0, vx0, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.x0[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.x0[ _idxa ] ) ); \
  \
  /* vx1 */ \
  _va = vec_lda( 0x00, &( buf.x1[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.x1[ _idxa ] ) ); \
  _vc = vec_perm( vx1, vx1, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.x1[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.x1[ _idxa ] ) ); \
  \
  /* vy0 */ \
  _va = vec_lda( 0x00, &( buf.y0[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.y0[ _idxa ] ) ); \
  _vc = vec_perm( vy0, vy0, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.y0[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.y0[ _idxa ] ) ); \
  \
  /* vy1 */ \
  _va = vec_lda( 0x00, &( buf.y1[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.y1[ _idxa ] ) ); \
  _vc = vec_perm( vy1, vy1, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.y1[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.y1[ _idxa ] ) ); \
  \
  /* q */ \
  _va = vec_lda( 0x00, &( buf.q[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.q[ _idxa ] ) ); \
  _vc = vec_perm( vq, vq, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.q[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.q[ _idxa ] ) ); \
  \ 
  /* vz */ \
  _va = vec_lda( 0x00, &( buf.vz[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.vz[ _idxa ] ) ); \
  _vc = vec_perm( vvz, vvz, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.vz[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.vz[ _idxa ] ) ); \
  \
  /* vix */ \
  _va = vec_ldiza( 0x00, &( buf.ix[ _idxa ] ) ); _vb = vec_ldiza( 0x10, &( buf.ix[ _idxa ] ) ); \
  _vc = vec_perm( _vix, _vix, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.ix[ _idxa ] ) ); vec_sta( _vb, 0x10, &( buf.ix[ _idxa ] ) ); \
  \
  /* viy */ \
  _va = vec_ldiza( 0x00, &( buf.iy[ _idxa ] ) ); _vb = vec_ldiza( 0x10, &( buf.iy[ _idxa ] ) ); \
  _vc = vec_perm( _viy, _viy, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.iy[ _idxa ] ) ); vec_sta( _vb, 0x10, &( buf.iy[ _idxa ] ) ); \
}

#endif

// For the current deposition routines the loads will always be aligned

#define LOAD4P2D( pbuf, idx, vx0, vx1, vy0, vy1, vq, vvz, vix, viy ) { \
   register vector _vix; \
   register vector _viy; \
\
   _vix = vec_ldiaa( 0x00, &(pbuf->ix[ idx ]) ); \
   _viy = vec_ldiaa( 0x00, &(pbuf->iy[ idx ]) ); \
   vx0  =   vec_lda( 0x00, &(pbuf->x0[ idx ]) ); \
   vx1  =   vec_lda( 0x00, &(pbuf->x1[ idx ]) ); \
   vy0  =   vec_lda( 0x00, &(pbuf->y0[ idx ]) ); \
   vy1  =   vec_lda( 0x00, &(pbuf->y1[ idx ]) ); \
   vq   =   vec_lda( 0x00, &(pbuf->q[ idx ]) ); \
   vvz  =   vec_lda( 0x00, &(pbuf->vz[ idx ]) ); \
\
   vec_sta( _vix, 0x00, vix ); \
   vec_sta( _viy, 0x00, viy ); \
}


/*****************************************************************************************
t_split_buf3D

Buffer to hold virtual particles for current deposition in 3D. 
*****************************************************************************************/

typedef struct { 
  // 4 splits maximum
  DECLARE_ALIGNED_32( double x0[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_32( double x1[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_32( double y0[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_32( double y1[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_32( double z0[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_32( double z1[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_32( double  q[ 4 * p_cache_size_3D ] );
  
  DECLARE_ALIGNED_32( int ix[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_32( int iy[ 4 * p_cache_size_3D ] );
  DECLARE_ALIGNED_32( int iz[ 4 * p_cache_size_3D ] );

  unsigned int np;
} t_split_buf3D;

#define STOREU4P3D( buf, idx, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { \
  buf.x0[idx  ] = vec_extract( vx0, 0 ); \
  buf.x0[idx+1] = vec_extract( vx0, 1 ); \
  buf.x0[idx+2] = vec_extract( vx0, 2 ); \
  buf.x0[idx+3] = vec_extract( vx0, 3 ); \
\
  buf.x1[idx  ] = vec_extract( vx1, 0 ); \
  buf.x1[idx+1] = vec_extract( vx1, 1 ); \
  buf.x1[idx+2] = vec_extract( vx1, 2 ); \
  buf.x1[idx+3] = vec_extract( vx1, 3 ); \
\
  buf.y0[idx  ] = vec_extract( vy0, 0 ); \
  buf.y0[idx+1] = vec_extract( vy0, 1 ); \
  buf.y0[idx+2] = vec_extract( vy0, 2 ); \
  buf.y0[idx+3] = vec_extract( vy0, 3 ); \
\
  buf.y1[idx  ] = vec_extract( vy1, 0 ); \
  buf.y1[idx+1] = vec_extract( vy1, 1 ); \
  buf.y1[idx+2] = vec_extract( vy1, 2 ); \
  buf.y1[idx+3] = vec_extract( vy1, 3 ); \
\
  buf.z0[idx  ] = vec_extract( vz0, 0 ); \
  buf.z0[idx+1] = vec_extract( vz0, 1 ); \
  buf.z0[idx+2] = vec_extract( vz0, 2 ); \
  buf.z0[idx+3] = vec_extract( vz0, 3 ); \
\
  buf.z1[idx  ] = vec_extract( vz1, 0 ); \
  buf.z1[idx+1] = vec_extract( vz1, 1 ); \
  buf.z1[idx+2] = vec_extract( vz1, 2 ); \
  buf.z1[idx+3] = vec_extract( vz1, 3 ); \
\
  buf.q[idx  ] = vec_extract( vq, 0 ); \
  buf.q[idx+1] = vec_extract( vq, 1 ); \
  buf.q[idx+2] = vec_extract( vq, 2 ); \
  buf.q[idx+3] = vec_extract( vq, 3 ); \
\
  buf.ix[idx]=vix[0]; buf.ix[idx+1]=vix[1]; buf.ix[idx+2]=vix[2]; buf.ix[idx+3]=vix[3]; \
  buf.iy[idx]=viy[0]; buf.iy[idx+1]=viy[1]; buf.iy[idx+2]=viy[2]; buf.iy[idx+3]=viy[3]; \
  buf.iz[idx]=viz[0]; buf.iz[idx+1]=viz[1]; buf.iz[idx+2]=viz[2]; buf.iz[idx+3]=viz[3]; \
\
}

#if 0

#define STOREU4P3D( buf, idx, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { \
  register vector _perm, _ins, _va, _vb, _vc, _vix, _viy, _viz; \
  \
  unsigned _idxa = idx & 0xFFFFFFFCL; \
  \
  _vix = vec_ldiza( 0x00, vix ); \
  _viy = vec_ldiza( 0x00, viy ); \
  _viz = vec_ldiza( 0x00, viz ); \
  \
  _perm = vec_lvsr( 0x00, &( buf.x0[ idx ] ) ); \
  _ins  = vec_perm( vec_cmplt(vx0,vx0), vec_cmpeq(vx0,vx0), _perm ); \
  \
  /* vx0 */ \
  _va = vec_lda( 0x00, &( buf.x0[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.x0[ _idxa ] ) ); \
  _vc = vec_perm( vx0, vx0, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.x0[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.x0[ _idxa ] ) ); \
  \
  /* vx1 */ \
  _va = vec_lda( 0x00, &( buf.x1[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.x1[ _idxa ] ) ); \
  _vc = vec_perm( vx1, vx1, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.x1[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.x1[ _idxa ] ) ); \
  \
  /* vy0 */ \
  _va = vec_lda( 0x00, &( buf.y0[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.y0[ _idxa ] ) ); \
  _vc = vec_perm( vy0, vy0, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.y0[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.y0[ _idxa ] ) ); \
  \
  /* vy1 */ \
  _va = vec_lda( 0x00, &( buf.y1[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.y1[ _idxa ] ) ); \
  _vc = vec_perm( vy1, vy1, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.y1[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.y1[ _idxa ] ) ); \
  \
  /* vz0 */ \
  _va = vec_lda( 0x00, &( buf.z0[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.z0[ _idxa ] ) ); \
  _vc = vec_perm( vz0, vz0, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.z0[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.z0[ _idxa ] ) ); \
  \
  /* vz1 */ \
  _va = vec_lda( 0x00, &( buf.z1[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.z1[ _idxa ] ) ); \
  _vc = vec_perm( vz1, vz1, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.z1[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.z1[ _idxa ] ) ); \
  \
  /* q */ \
  _va = vec_lda( 0x00, &( buf.q[ _idxa ] ) ); _vb = vec_lda( 0x20, &( buf.q[ _idxa ] ) ); \
  _vc = vec_perm( vq, vq, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.q[ _idxa ] ) ); vec_sta( _vb, 0x20, &( buf.q[ _idxa ] ) ); \
  \ 
  /* vix */ \
  _va = vec_ldiza( 0x00, &( buf.ix[ _idxa ] ) ); _vb = vec_ldiza( 0x10, &( buf.ix[ _idxa ] ) ); \
  _vc = vec_perm( _vix, _vix, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.ix[ _idxa ] ) ); vec_sta( _vb, 0x10, &( buf.ix[ _idxa ] ) ); \
  \
  /* viy */ \
  _va = vec_ldiza( 0x00, &( buf.iy[ _idxa ] ) ); _vb = vec_ldiza( 0x10, &( buf.iy[ _idxa ] ) ); \
  _vc = vec_perm( _viy, _viy, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.iy[ _idxa ] ) ); vec_sta( _vb, 0x10, &( buf.iy[ _idxa ] ) ); \
  \
  /* viz */ \
  _va = vec_ldiza( 0x00, &( buf.iz[ _idxa ] ) ); _vb = vec_ldiza( 0x10, &( buf.iz[ _idxa ] ) ); \
  _vc = vec_perm( _viz, _viz, _perm ); \
  _va = vec_sel( _va, _vc, _ins ); _vb = vec_sel( _vc, _vb, _ins ); \
  vec_sta( _va, 0x00, &( buf.iz[ _idxa ] ) ); vec_sta( _vb, 0x10, &( buf.iz[ _idxa ] ) ); \
}

#endif

#define LOAD4P3D( pbuf, idx, vx0, vx1, vy0, vy1, vz0, vz1, vq, vix, viy, viz ) { \
   register vector _vix; \
   register vector _viy; \
   register vector _viz; \
\
   _vix = vec_ldiaa( 0x00, &(pbuf->ix[ idx ]) ); \
   _viy = vec_ldiaa( 0x00, &(pbuf->iy[ idx ]) ); \
   _viz = vec_ldiaa( 0x00, &(pbuf->iz[ idx ]) ); \
   vx0 = vec_lda( 0x00, &(pbuf->x0[ idx ]) ); \
   vx1 = vec_lda( 0x00, &(pbuf->x1[ idx ]) ); \
   vy0 = vec_lda( 0x00, &(pbuf->y0[ idx ]) ); \
   vy1 = vec_lda( 0x00, &(pbuf->y1[ idx ]) ); \
   vq  = vec_lda( 0x00, &(pbuf->q[ idx ]) ); \
\
   vec_sta( _vix, 0x00, vix ); \
   vec_sta( _viy, 0x00, viy ); \
   vec_sta( _viz, 0x00, viz ); \
}

/*****************************************************************************************
VRGAMMA

Get 1 / Lorentz gamma. 
*****************************************************************************************/

#define VRGAMMA( vrg, vu1, vu2, vu3 ) {                     \
   const vector one = vec_splats( 1.0 );                    \
   (vrg) = vec_rsqrt( vec_madd( (vu3), (vu3),               \
                      vec_madd( (vu2), (vu2),               \
                      vec_madd( (vu1), (vu1), one ) ) ) );  \
}


inline vector vntrim( const vector vx );


#endif

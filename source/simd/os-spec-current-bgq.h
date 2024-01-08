#ifndef _OS_SPEC_CURRENT_H
#define _OS_SPEC_CURRENT_H

#include "vector-bgq.h"
#include "os-spec-push-bgq.h"

inline void vwl_s1( vector const vqn, vector const vx0, vector const vx1, vector vwl[] );
inline void vwl_s2( vector const vqn, vector const vx0, vector const vx1, vector vwl[] );
inline void vwl_s3( vector const vqn, vector const vx0, vector const vx1, vector vwl[] );
inline void vwl_s4( vector const vqn, vector const vx0, vector const vx1, vector vwl[] );


/* restrict Current deposition routines for BG/Q */

inline void vdepcurrent_1d_s1(t_real * restrict const current, int const * restrict const size, int const * restrict const offset, 
                       double * restrict const norm, t_split_buf1D * restrict const part);
inline void vdepcurrent_1d_s2(t_real * restrict const current, int const * restrict const size, int const * restrict const offset, 
                       double * restrict const norm, t_split_buf1D * restrict const part);
inline void vdepcurrent_1d_s3(t_real * restrict const current, int const * restrict const size, int const * restrict const offset, 
                       double * restrict const norm, t_split_buf1D * restrict const part);
inline void vdepcurrent_1d_s4(t_real * restrict const current, int const * restrict const size, int const * restrict const offset, 
                       double * restrict const norm, t_split_buf1D * restrict const part);


inline void vdepcurrent_2d_s1(t_real * restrict const current, int const * restrict const size, int const * restrict const offset, 
                       double * restrict const norm, t_split_buf2D * restrict const part);
inline void vdepcurrent_2d_s2(t_real * restrict const current, int const * restrict const size, int const * restrict const offset, 
                       double * restrict const norm, t_split_buf2D * restrict const part);
inline void vdepcurrent_2d_s3(t_real * restrict const current, int const * restrict const size, int const * restrict const offset, 
                       double * restrict const norm, t_split_buf2D * restrict const part);
inline void vdepcurrent_2d_s4(t_real * restrict const current, int const * restrict const size, int const * restrict const offset, 
                       double * restrict const norm, t_split_buf2D * restrict const part);

inline void vdepcurrent_3d_s1(t_real * restrict const  current, int const * restrict const  size, int const * restrict const  offset,
                       double * restrict const norm, t_split_buf3D * restrict const part);
inline void vdepcurrent_3d_s2(t_real * restrict const  current, int const * restrict const  size, int const * restrict const  offset,
                       double * restrict const norm, t_split_buf3D * restrict const part);
inline void vdepcurrent_3d_s3(t_real * restrict const  current, int const * restrict const  size, int const * restrict const  offset,
                       double * restrict const norm, t_split_buf3D * restrict const part);
inline void vdepcurrent_3d_s4(t_real * restrict const  current, int const * restrict const  size, int const * restrict const  offset,
                       double * restrict const norm, t_split_buf3D * restrict const part);

#endif

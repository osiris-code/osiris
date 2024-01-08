#ifndef _OS_SPEC_CURRENT_H
#define _OS_SPEC_CURRENT_H

#include "vector-sse.h"
#include "os-spec-push-ssed.h"

void vwl_s1( __m128d const vqn, __m128d const vx0, __m128d const vx1, __m128d vwl[] );
void vwl_s2( __m128d const vqn, __m128d const vx0, __m128d const vx1, __m128d vwl[] );
void vwl_s3( __m128d const vqn, __m128d const vx0, __m128d const vx1, __m128d vwl[] );
void vwl_s4( __m128d const vqn, __m128d const vx0, __m128d const vx1, __m128d vwl[] );

/* Current deposition routines for x86 SSE */

void vdepcurrent_1d_s1(double * const current, int const * const size, int const * const offset, 
                       double * const norm, t_split_buf1D * const part);
void vdepcurrent_1d_s2(double * const current, int const * const size, int const * const offset, 
                       double * const norm, t_split_buf1D * const part);
void vdepcurrent_1d_s3(double * const current, int const * const size, int const * const offset, 
                       double * const norm, t_split_buf1D * const part);
void vdepcurrent_1d_s4(double * const current, int const * const size, int const * const offset, 
                       double * const norm, t_split_buf1D * const part);

void vdepcurrent_2d_s1(double * const current, int const * const size, int const * const offset, 
                       double * const norm, t_split_buf2D * const part);
void vdepcurrent_2d_s2(double * const current, int const * const size, int const * const offset, 
                       double * const norm, t_split_buf2D * const part);
void vdepcurrent_2d_s3(double * const current, int const * const size, int const * const offset, 
                       double * const norm, t_split_buf2D * const part);
void vdepcurrent_2d_s4(double * const current, int const * const size, int const * const offset, 
                       double * const norm, t_split_buf2D * const part);

void vdepcurrent_3d_s1(double * const current, int const * const  size, int const * const  offset,
                       double * const norm, t_split_buf3D * const part);
void vdepcurrent_3d_s2(double * const current, int const * const  size, int const * const  offset,
                       double * const norm, t_split_buf3D * const part);
void vdepcurrent_3d_s3(double * const current, int const * const  size, int const * const  offset,
                       double * const norm, t_split_buf3D * const part);
void vdepcurrent_3d_s4(double * const current, int const * const  size, int const * const  offset,
                       double * const norm, t_split_buf3D * const part);


#endif
